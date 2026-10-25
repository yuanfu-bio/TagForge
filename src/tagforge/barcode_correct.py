from __future__ import annotations

import csv
import json
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timedelta
from functools import lru_cache
from itertools import chain
from pathlib import Path
from typing import Dict, Iterable

from . import __version__
from .columnar import (ParquetFragmentWriter, _pa, parquet_files, parquet_size, record_batches,
                       resolve_table, write_parquet_columns)
from .config import CorrectionConfig, TagForgeConfig, ConfigError
from .extract import decode_method_payload, decode_segment_payload, hamming
from .io_utils import atomic_text, sample_dirs, write_progress, write_json, write_tsv
from .logging_utils import sample_logger


_WORKER_CONFIG = None
_WORKER_CORRECTORS = None
_WORKER_ANNOTATION = None


@dataclass(frozen=True)
class CorrectionResult:
    raw_sequence: str
    shifted_sequence: str
    corrected_sequence: str
    correction_status: str
    shift_distance: int
    mismatch_distance: int
    correction_type: str
    ambiguous: bool = False

    @property
    def success(self):
        return self.correction_status in {"valid", "disabled"}


class WhitelistCorrector:
    def __init__(self, values: Iterable[str], config: CorrectionConfig, length: int):
        self.values = tuple(sorted({x.strip().upper() for x in values if x.strip()}))
        self.value_set = set(self.values)
        if any(len(value) != length for value in self.values):
            raise ConfigError(f"Whitelist sequences must all have configured length {length}")
        self.config = config
        self.length = length
        # Substitution neighbours can be enumerated directly when every
        # whitelist sequence is plain ACGT; otherwise fall back to a full scan.
        self._acgt_whitelist = all(set(value) <= set("ACGT") for value in self.values)

    def _neighbours(self, shifted: str):
        """Whitelist sequences at Hamming distance exactly 1, in sorted order."""
        found = []
        value_set = self.value_set
        for index, current in enumerate(shifted):
            prefix, suffix = shifted[:index], shifted[index + 1:]
            for base in "ACGT":
                if base != current:
                    candidate = prefix + base + suffix
                    if candidate in value_set:
                        found.append(candidate)
        found.sort()
        return found

    @lru_cache(maxsize=200000)
    def correct(self, raw: str, anchor: int = 0) -> CorrectionResult:
        raw = raw.upper()
        if not self.config.enabled:
            sequence = raw[anchor:anchor + self.length]
            if len(sequence) != self.length:
                return CorrectionResult(raw, sequence, "", "failed", 0, -1, "failed")
            return CorrectionResult(raw, sequence, sequence, "disabled", 0, 0, "disabled")
        original = raw[anchor:anchor + self.length]
        if original in self.value_set:
            # The original exact hit has the unique minimum candidate score.
            return CorrectionResult(raw, original, original, "valid", 0, 0, "exact")
        candidates = []
        max_shift = self.config.max_shift if self.config.allow_shift else 0
        max_mismatch = self.config.max_mismatch if self.config.allow_mismatch else 0
        shifts = [0] + [signed for distance in range(1, max_shift + 1) for signed in (-distance, distance)]
        for shift in shifts:
            start = anchor + shift
            if start < 0:
                continue
            shifted = raw[start:start + self.length]
            if len(shifted) != self.length:
                continue
            if shifted in self.value_set:
                candidates.append((shift, 0, shifted, shifted))
                continue
            if max_mismatch == 1 and self._acgt_whitelist:
                # shifted is not a member, so the nearest entries are exactly the
                # distance-1 neighbours; anything farther exceeds the limit.
                for known in self._neighbours(shifted):
                    candidates.append((shift, 1, known, shifted))
            elif max_mismatch:
                best_distance = max_mismatch + 1
                best = []
                for known in self.values:
                    distance = hamming(shifted, known)
                    if distance < best_distance:
                        best_distance, best = distance, [known]
                    elif distance == best_distance:
                        best.append(known)
                if best_distance <= max_mismatch:
                    for known in best:
                        candidates.append((shift, best_distance, known, shifted))
        if not candidates:
            return CorrectionResult(raw, raw[anchor:anchor + self.length], "", "failed", 0, -1, "failed")
        score = min(
            ((x[0] != 0) + (x[1] > 0), abs(x[0]) + x[1], x[1], abs(x[0]))
            for x in candidates
        )
        best = [x for x in candidates if ((x[0] != 0) + (x[1] > 0), abs(x[0]) + x[1], x[1], abs(x[0])) == score]
        sequences = {x[2] for x in best}
        if len(sequences) != 1:
            return CorrectionResult(raw, best[0][3], "", "ambiguous", best[0][0], best[0][1], "failed", True)
        shift, mismatch, corrected, shifted = best[0]
        kind = "exact" if not shift and not mismatch else (
            "shift_and_mismatch" if shift and mismatch else "shift_only" if shift else "mismatch_only")
        return CorrectionResult(raw, shifted, corrected, "valid", shift, mismatch, kind)


def load_whitelist(path: Path):
    return list(load_whitelist_names(path))


def load_whitelist_names(path: Path) -> Dict[str, str]:
    """Read a sequence whitelist, optionally with ``name<TAB>sequence`` rows."""
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        sequence = fields[-1].upper()
        name = fields[0] if len(fields) > 1 else sequence
        if sequence in values:
            raise ConfigError(f"Whitelist contains duplicate entry: {sequence}")
        values[sequence] = name
    if not values:
        raise ConfigError(f"Whitelist is empty: {path}")
    return values


def load_fb_annotation(config: TagForgeConfig) -> Dict[str, str]:
    if getattr(config, "barcode2_sequence_only", False):
        return {}
    if config.fb_info is None:
        raise ConfigError("barcode2 annotation is required unless barcode2.sequence_only is true")
    with open(config.fb_info, encoding="utf-8", newline="") as handle:
        required = {config.fb_id_column, config.fb_sequence_column, config.fb_name_column}
        reader = csv.reader(handle, delimiter="\t")
        first_row = next(reader, None)
        if first_row is None:
            raise ConfigError(f"FB annotation is empty: {config.fb_info}")

        has_header = required.issubset(set(first_row))
        if has_header:
            rows = csv.DictReader(handle, fieldnames=first_row, delimiter="\t")
        else:
            if len(first_row) < 2:
                raise ConfigError("Headerless FB annotation requires at least ID and sequence columns")
            # Historical 10x-style files have no header: ID, sequence, name.
            # Use the sequence as the name when the optional third column is absent.
            rows = chain((first_row,), reader)
        mapping, names = {}, set()
        for row in rows:
            if has_header:
                sequence = row[config.fb_sequence_column].strip().upper()
                name = row[config.fb_name_column].strip()
            else:
                sequence = row[1].strip().upper()
                name = row[2].strip() if len(row) > 2 else sequence
            if sequence in mapping:
                raise ConfigError(f"Duplicate FB sequence in annotation: {sequence}")
            if name in names and not config.allow_duplicate_names:
                raise ConfigError(f"Duplicate antibody name in annotation: {name}")
            mapping[sequence] = name; names.add(name)
    return mapping


def _build_correctors(config: TagForgeConfig, whitelist_values_by_segment=None):
    correctors = {}
    whitelist_values_by_segment = whitelist_values_by_segment or {}
    for segment in config.segments:
        if segment.target == "umi":
            continue
        if segment.whitelist:
            values = whitelist_values_by_segment.get(segment.name)
            if values is None:
                values = load_whitelist(segment.whitelist)
            correctors[segment.name] = WhitelistCorrector(
                values, segment.correction, segment.length
            )
        elif segment.correction.enabled:
            raise ConfigError(f"Segment {segment.name}: correction enabled but whitelist is missing")
        else:
            correctors[segment.name] = WhitelistCorrector([], segment.correction, segment.length)
    return correctors


def _correction_worker_init(config, whitelist_values_by_segment, annotation):
    global _WORKER_CONFIG, _WORKER_CORRECTORS, _WORKER_ANNOTATION
    _WORKER_CONFIG = config
    _WORKER_CORRECTORS = _build_correctors(config, whitelist_values_by_segment)
    _WORKER_ANNOTATION = annotation


def _correction_worker_ready():
    return True


def _merge_counter_maps(target, source):
    for key, value in source.items():
        target[key].update(value)


def _target_display(config, role: str) -> str:
    for segment in config.segments:
        if segment.target == role:
            return segment.target_name or role
    return role


def _correct_batch(batch, config, correctors, annotation):
    counters = {
        segment.name: Counter() for segment in config.segments if segment.target != "umi"
    }
    summary = Counter()
    valid_rows = []
    trace_rows = []
    for row in batch:
        valid_row, rows = _correct_row(
            row, config, correctors, annotation, counters, summary, config.trace_enabled
        )
        if valid_row is not None:
            valid_rows.append(valid_row)
        if config.trace_enabled:
            trace_rows.extend(rows)
    return (
        valid_rows,
        trace_rows,
        {name: dict(counter) for name, counter in counters.items()},
        dict(summary),
    )


def _correct_batch_worker(batch):
    if _WORKER_CONFIG is None or _WORKER_CORRECTORS is None or _WORKER_ANNOTATION is None:
        raise RuntimeError("Correction worker was not initialized")
    return _correct_batch(batch, _WORKER_CONFIG, _WORKER_CORRECTORS, _WORKER_ANNOTATION)


def _read_part_rows(path):
    _, pq = _pa()
    table = pq.read_table(path)
    names = table.column_names
    columns = [table.column(name).to_pylist() for name in names]
    return [
        {name: "" if value is None else str(value) for name, value in zip(names, values)}
        for values in zip(*columns)
    ]


def _correct_part(task, config, correctors, annotation):
    """Correct one extracted fragment and commit its output fragments."""
    index, path, valid_dir, trace_dir, compression, compression_level = task
    valid_rows, trace_rows, local_counters, local_summary = _correct_batch(
        _read_part_rows(path), config, correctors, annotation)
    name = f"part-{index:08d}.parquet"
    fields = valid_fields(config)
    write_parquet_columns(
        Path(valid_dir) / name, {field: [row[field] for row in valid_rows] for field in fields},
        fields, compression_level=compression_level, compression=compression)
    if trace_dir is not None:
        write_parquet_columns(
            Path(trace_dir) / name,
            {field: [row[field] for row in trace_rows] for field in TRACE_FIELDS},
            TRACE_FIELDS, compression_level=compression_level, compression=compression)
    return len(valid_rows), len(trace_rows), local_counters, local_summary


def _correct_part_worker(task):
    if _WORKER_CONFIG is None or _WORKER_CORRECTORS is None or _WORKER_ANNOTATION is None:
        raise RuntimeError("Correction worker was not initialized")
    return _correct_part(task, _WORKER_CONFIG, _WORKER_CORRECTORS, _WORKER_ANNOTATION)


def _fragment_dataset(source: Path):
    """Row-counted fragment list for a Parquet dataset directory, else None."""
    if not source.is_dir():
        return None
    _, pq = _pa()
    files = parquet_files(source)
    return [(part, pq.ParquetFile(part).metadata.num_rows) for part in files]


TRACE_FIELDS = ["read_id", "segment_name", "target_type", "raw_sequence", "shifted_sequence",
                "corrected_sequence", "whitelist_hit", "correction_status", "shift_distance",
                "mismatch_distance", "correction_type"]

REJECTION_REASONS = ("extraction_failed", "barcode_ambiguous", "barcode_invalid",
                     "annotation_missing", "umi_invalid")


def valid_fields(config: TagForgeConfig):
    return [
        "read_id", config.target_name("barcode1"),
        f"{config.target_name('barcode2')}_sequence",
        f"{config.target_name('barcode2')}_name",
        config.target_name("umi"), "correction_summary",
    ]


def _correction_fingerprint(config: TagForgeConfig, source: Path, sample_name: str) -> str:
    from .checkpoint import resume_fingerprint
    return resume_fingerprint(config, sample_name, "correct")


def evaluate_extracted_read(row, config, correctors, annotation, counters, summary, trace=True):
    """Shared formal/QC read eligibility and exclusive rejection counters."""
    summary["total_reads"] += 1
    segment_columns = {
        "barcode1": config.segment_column("barcode1"),
        "barcode2": config.segment_column("barcode2"),
        "umi": config.segment_column("umi"),
    }
    required = {
        "read_id", *segment_columns.values(),
        "methods", "status",
    }
    missing = required - set(row)
    if missing:
        raise ConfigError(
            "Incompatible extracted table schema; missing column(s): "
            + ", ".join(sorted(missing))
            + f". Re-run 'tagforge extract --overwrite' with TagForge {__version__}."
        )
    if row["status"] != "success":
        summary["failure_extraction_failed"] += 1
        return None, []
    summary["extracted_reads"] += 1
    extraction_methods = decode_method_payload(row["methods"], config.segments)
    values = {}
    trace_rows = []
    ambiguous = False
    for target, col in (("barcode1", segment_columns["barcode1"]), ("barcode2", segment_columns["barcode2"])):
        target_segments = [segment for segment in config.segments if segment.target == target]
        raw_values = decode_segment_payload(row[col], target_segments)
        corrected_values = []
        for segment in target_segments:
            source_method = extraction_methods.get(segment.name, "unknown")
            extra = (
                segment.correction.max_shift
                if source_method == "fixed" and segment.correction.enabled
                and segment.correction.allow_shift else 0
            )
            anchor = min(extra, segment.start or 0) if source_method == "fixed" else 0
            result = correctors[segment.name].correct(raw_values.get(segment.name, ""), anchor)
            ambiguous |= result.ambiguous
            corrected_values.append(result.corrected_sequence if result.success else "")
            count = counters[segment.name]
            count["extracted_reads"] += 1
            count[result.correction_type] += 1
            count[f"{source_method}_extracted"] += 1
            count[f"mismatch_{max(result.mismatch_distance, 0)}"] += int(result.mismatch_distance >= 0)
            count[f"shift_{abs(result.shift_distance)}"] += int(result.success)
            if result.success and result.shift_distance:
                direction = "left" if result.shift_distance < 0 else "right"
                count[f"shift_{direction}_{abs(result.shift_distance)}"] += 1
            if result.success:
                count["valid_reads"] += 1
                count[f"{source_method}_valid"] += 1
            else:
                count["invalid_reads"] += 1
            if trace:
                trace_rows.append({
                "read_id": row["read_id"], "segment_name": segment.name,
                "target_type": segment.target_name or target, "raw_sequence": result.raw_sequence,
                "shifted_sequence": result.shifted_sequence,
                "corrected_sequence": result.corrected_sequence,
                "whitelist_hit": str(result.success).lower(),
                "correction_status": result.correction_status,
                "shift_distance": result.shift_distance,
                "mismatch_distance": result.mismatch_distance,
                "correction_type": result.correction_type,
            })
        values[target] = "".join(corrected_values) if all(corrected_values) else ""
    umi_segments = [segment for segment in config.segments if segment.target == "umi"]
    umi_values = decode_segment_payload(row[segment_columns["umi"]], umi_segments)
    umi = "".join(umi_values.get(segment.name, "") for segment in umi_segments)
    fb_name = values["barcode2"] if getattr(config, "barcode2_sequence_only", False) else annotation.get(values["barcode2"], "")
    if values["barcode2"]:
        summary["barcode2_whitelist_reads"] += 1
        if not fb_name:
            summary["missing_annotation:" + values["barcode2"]] += 1
    if values["barcode1"]:
        summary["barcode1_valid"] += 1
    if values["barcode2"] and fb_name:
        summary["barcode2_valid"] += 1
    umi_valid = all(len(umi_values.get(segment.name, "")) == segment.length
                    for segment in umi_segments) and bool(umi) and set(umi) <= set("ACGT")
    if not values["barcode1"] or not values["barcode2"]:
        summary["failure_barcode_ambiguous" if ambiguous else "failure_barcode_invalid"] += 1
    elif not fb_name:
        summary["failure_annotation_missing"] += 1
    elif not umi_valid:
        summary["failure_umi_invalid"] += 1
    else:
        summary["combined_valid"] += 1
        return {
            "read_id": row["read_id"], config.target_name("barcode1"): values["barcode1"],
            f"{config.target_name('barcode2')}_sequence": values["barcode2"],
            f"{config.target_name('barcode2')}_name": fb_name,
            config.target_name("umi"): umi, "correction_summary": "valid",
        }, trace_rows
    return None, trace_rows


# Keep the existing batch entry point while QC uses the same evaluator.
_correct_row = evaluate_extracted_read


def _correction_stats(config, counters, summary):
    stat_rows = []
    total = summary["total_reads"]
    for segment in (s for s in config.segments if s.target != "umi"):
        c = counters[segment.name]
        row = {"scope": segment.name, "target_type": segment.target_name or segment.target, "total_reads": total,
               "extracted_reads": c["extracted_reads"], "valid_reads": c["valid_reads"], "invalid_reads": c["invalid_reads"],
               "valid_rate": c["valid_reads"] / total if total else 0,
               "exact_count": c["exact"], "mismatch_only_count": c["mismatch_only"],
               "shift_only_count": c["shift_only"], "shift_and_mismatch_count": c["shift_and_mismatch"],
               "failed_count": c["failed"], "disabled_count": c["disabled"],
               "linker_extracted_count": c["linker_extracted"],
               "linker_barcode_valid_count": c["linker_valid"],
               "linker_barcode_valid_rate": c["linker_valid"] / c["linker_extracted"] if c["linker_extracted"] else "",
               "fixed_extracted_count": c["fixed_extracted"],
               "fixed_barcode_valid_count": c["fixed_valid"],
               "fixed_barcode_valid_rate": c["fixed_valid"] / c["fixed_extracted"] if c["fixed_extracted"] else ""}
        for i in range(max(segment.correction.max_mismatch, 2) + 1): row[f"mismatch_{i}_count"] = c[f"mismatch_{i}"]
        for i in range(max(segment.correction.max_shift, 2) + 1):
            row[f"shift_{i}_count"] = c[f"shift_{i}"]
            if i:
                row[f"shift_left_{i}_count"] = c[f"shift_left_{i}"]
                row[f"shift_right_{i}_count"] = c[f"shift_right_{i}"]
        stat_rows.append(row)
    stat_rows.extend([
        {"scope": f"final_{config.target_name('barcode1')}", "target_type": _target_display(config, "barcode1"), "total_reads": total, "valid_reads": summary["barcode1_valid"], "valid_rate": summary["barcode1_valid"] / total if total else 0},
        {"scope": f"final_{config.target_name('barcode2')}", "target_type": _target_display(config, "barcode2"), "total_reads": total, "valid_reads": summary["barcode2_valid"], "valid_rate": summary["barcode2_valid"] / total if total else 0},
        {"scope": "combined", "target_type": "all", "total_reads": total, "valid_reads": summary["combined_valid"], "valid_rate": summary["combined_valid"] / total if total else 0},
    ])
    stat_rows.extend({"scope": f"failure:{reason}", "target_type": "all", "total_reads": total,
                      "invalid_reads": summary[f"failure_{reason}"]}
                     for reason in REJECTION_REASONS)
    all_fields = ["scope", "target_type", "total_reads", "extracted_reads", "valid_reads", "invalid_reads", "valid_rate",
                  "exact_count", "mismatch_only_count", "shift_only_count", "shift_and_mismatch_count", "failed_count", "disabled_count"]
    dynamic = sorted({k for row in stat_rows for k in row if k not in all_fields})
    result_summary = dict(summary)
    result_summary["segment_method_qc"] = [
        {key: row.get(key, "") for key in (
            "scope", "target_type", "linker_extracted_count", "linker_barcode_valid_count",
            "linker_barcode_valid_rate", "fixed_extracted_count", "fixed_barcode_valid_count",
            "fixed_barcode_valid_rate",
        )}
        for row in stat_rows if not row.get("scope", "").startswith("failure:") and row.get("scope") not in {
            f"final_{config.target_name('barcode1')}",
            f"final_{config.target_name('barcode2')}",
            "combined",
        }
    ]
    return all_fields + dynamic, stat_rows, result_summary


def correct_sample(config: TagForgeConfig, sample_name: str, resume: bool = True, *, keep_resume: bool = False):
    dirs = sample_dirs(config.output_dir, sample_name)
    source = resolve_table(dirs["extracted"] / f"{sample_name}.extracted.parquet",
                           dirs["extracted"] / f"{sample_name}.extracted.tsv.gz")
    trace_parquet = dirs["corrected"] / f"{sample_name}.barcode_correction_trace.parquet"
    stats_path = dirs["corrected"] / f"{sample_name}.correction_qc.json"
    valid_parquet = dirs["corrected"] / f"{sample_name}.valid_reads.parquet"
    progress_path = dirs["logs"] / f"{sample_name}.progress.json"
    resume_path = dirs["checkpoint"] / f"{sample_name}.correct.resume.json"
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    annotation = load_fb_annotation(config)
    whitelist_values_by_segment = {
        segment.name: load_whitelist(segment.whitelist)
        for segment in config.segments
        if segment.target != "umi" and segment.whitelist
    }
    correctors = _build_correctors(config, whitelist_values_by_segment)
    counters = {
        segment.name: Counter() for segment in config.segments if segment.target != "umi"
    }
    summary = Counter()
    started = time.monotonic()
    previous_elapsed = 0.0
    last_input_fraction = 0.0
    fingerprint = _correction_fingerprint(config, source, sample_name)
    requested_workers = config.barcode_workers or config.threads
    workers = requested_workers
    executor = None
    fields = valid_fields(config)

    def save_resume():
        state = {
            "schema": 4, "tagforge_version": __version__, "fingerprint": fingerprint,
            "trace_enabled": config.trace_enabled,
            "reads_completed": summary["total_reads"],
            "safe_valid_parts": valid_writer.part,
            "safe_trace_parts": trace_writer.part if config.trace_enabled else 0,
            "elapsed_seconds": previous_elapsed + time.monotonic() - started,
            "input_fraction": last_input_fraction,
            "summary": dict(summary),
            "segment_counters": {
                name: dict(counter) for name, counter in counters.items()
            },
        }
        with atomic_text(resume_path) as state_handle:
            json.dump(state, state_handle, separators=(",", ":"))


    resumed = False
    if resume and resume_path.exists():
        state = json.loads(resume_path.read_text(encoding="utf-8"))
        if (
            state.get("schema") != 4
            or state.get("fingerprint") != fingerprint
            or bool(state.get("trace_enabled")) != config.trace_enabled
        ):
            raise ConfigError(
                f"Correction resume state does not match current inputs/config: {resume_path}. "
                "Use --overwrite to restart safely."
            )
        valid_writer = ParquetFragmentWriter(
            valid_parquet, fields, config.chunk_size, resume=True,
            parts=int(state["safe_valid_parts"]),
         compression_level=getattr(config, "compression_level", 3), compression=getattr(config, "parquet_compression", "snappy"))
        trace_writer = (
            ParquetFragmentWriter(trace_parquet, TRACE_FIELDS, config.chunk_size, resume=True,
                                  parts=int(state["safe_trace_parts"]), compression_level=getattr(config, "compression_level", 3), compression=getattr(config, "parquet_compression", "snappy"))
            if config.trace_enabled else None
        )
        summary.update(state.get("summary", {}))
        for name, values in state.get("segment_counters", {}).items():
            if name in counters:
                counters[name].update(values)
        if valid_writer.rows != summary["combined_valid"] or (trace_writer is not None
            and trace_writer.rows != sum(counter["extracted_reads"] for counter in counters.values())):
            raise ConfigError("Committed correction rows disagree with resume counters; use --overwrite")
        previous_elapsed = float(state.get("elapsed_seconds", 0.0))
        last_input_fraction = float(state.get("input_fraction", 0.0))
        resumed = True
        logger.info(
            "correction_resume\treads=%s\tsafe_valid_parts=%s\tsafe_trace_parts=%s",
            summary["total_reads"], valid_writer.part, trace_writer.part if trace_writer else 0,
        )
    else:
        if resume and not resume_path.exists() and (valid_parquet.exists() or trace_parquet.exists()):
            raise ConfigError(
                "Untracked correction Parquet dataset cannot be resumed safely. Use --overwrite to restart."
            )
        valid_writer = ParquetFragmentWriter(valid_parquet, fields, config.chunk_size, compression_level=getattr(config, "compression_level", 3), compression=getattr(config, "parquet_compression", "snappy"))
        trace_writer = (ParquetFragmentWriter(trace_parquet, TRACE_FIELDS, config.chunk_size, compression_level=getattr(config, "compression_level", 3), compression=getattr(config, "parquet_compression", "snappy"))
                        if config.trace_enabled else None)
        resume_path.unlink(missing_ok=True)
        save_resume()


    def update_progress(status, input_fraction, resume_skip=None):
        elapsed = previous_elapsed + time.monotonic() - started
        rate = summary["total_reads"] / elapsed if elapsed else 0.0
        eta = (
            elapsed * (1.0 - input_fraction) / input_fraction
            if input_fraction is not None and 0 < input_fraction < 1 else 0.0
        )
        finish = (
            (datetime.now().astimezone() + timedelta(seconds=eta)).isoformat(timespec="seconds")
            if eta else ""
        )
        progress = {
            "status": status, "reads_completed": summary["total_reads"],
            "extracted_reads": summary["extracted_reads"],
            "combined_valid_reads": summary["combined_valid"],
            "input_percent": f"{input_fraction * 100:.2f}" if input_fraction is not None else "NA",
            "elapsed_seconds": f"{elapsed:.2f}", "reads_per_second": f"{rate:.2f}",
            "eta_seconds": f"{eta:.2f}" if eta else "", "estimated_finish": finish,
            "resume_skip_percent": f"{resume_skip * 100:.2f}" if resume_skip is not None else "",
            "valid_tmp_bytes": parquet_size(valid_parquet),
            "trace_tmp_bytes": parquet_size(trace_parquet) if config.trace_enabled else 0,
        }
        write_progress(progress_path, sample_name, "correct", progress)
        logger.info(
            "correction_progress\treads=%s\textracted_reads=%s\tcombined_valid_reads=%s\t"
            "input_percent=%s\treads_per_second=%s\teta_seconds=%s\t"
            "estimated_finish=%s\tresume_skip_percent=%s\tstatus=%s",
            progress["reads_completed"], progress["extracted_reads"],
            progress["combined_valid_reads"], progress["input_percent"],
            progress["reads_per_second"], progress["eta_seconds"] or "NA",
            progress["estimated_finish"] or "NA", progress["resume_skip_percent"] or "NA",
            status,
        )
        if status == "running":
            print(
                f"{sample_name} correct: {summary['total_reads']:,} reads, "
                f"input≈{progress['input_percent']}%, speed={float(progress['reads_per_second']):,.0f}/s, "
                f"ETA={progress['eta_seconds'] + 's' if progress['eta_seconds'] else 'NA'}",
                flush=True,
            )

    def append_correction_member(result, committed=False):
        valid_rows, trace_rows, local_counters, local_summary = result
        if committed:
            # Workers already wrote their fragments; only advance the counters.
            valid_writer.part += 1
            valid_writer.rows += valid_rows
            if trace_writer is not None:
                trace_writer.part += 1
                trace_writer.rows += trace_rows
        else:
            valid_writer.writerows(valid_rows)
            valid_writer.flush()
            if trace_writer is not None:
                trace_writer.writerows(trace_rows)
                trace_writer.flush()
        summary.update(local_summary)
        _merge_counter_maps(counters, local_counters)

    try:
        if requested_workers > 1:
            try:
                executor = ProcessPoolExecutor(
                    max_workers=requested_workers,
                    initializer=_correction_worker_init,
                    initargs=(config, whitelist_values_by_segment, annotation),
                )
                executor.submit(_correction_worker_ready).result()
            except Exception as exc:
                if executor is not None:
                    executor.shutdown(wait=False, cancel_futures=True)
                executor = None
                workers = 1
                logger.info(
                    "correction_parallel_fallback\trequested_workers=%s\tworkers=1\treason=%s",
                    requested_workers, type(exc).__name__,
                )
        logger.info(
            "correction_parallel_start\tbackend=tagforge-process-pool\trequested_workers=%s\t"
            "workers=%s\tchunk_size=%s\tcompression_level=%s\ttrace_enabled=%s",
            requested_workers, workers, config.chunk_size,
            config.compression_level, str(config.trace_enabled).lower(),
        )
        update_progress("resuming" if resumed else "running", last_input_fraction, 0.0 if resumed else None)
        reads_to_skip = summary["total_reads"]
        reads_seen = 0
        pending = {}
        completed = {}
        next_sequence = 0
        next_commit = 0

        def collect_finished(futures):
            for future in futures:
                sequence, input_fraction = pending.pop(future)
                completed[sequence] = (future.result(), input_fraction)

        def commit_ready():
            nonlocal last_input_fraction, next_commit
            while next_commit in completed:
                result, input_fraction = completed.pop(next_commit)
                append_correction_member(result, isinstance(result[0], int))
                last_input_fraction = max(last_input_fraction, input_fraction)
                save_resume()
                update_progress("running", last_input_fraction)
                next_commit += 1

        parts_checked = False
        fragments = _fragment_dataset(source)
        if fragments is not None:
            total_rows = sum(count for _, count in fragments) or 1
            rows_before = 0
            work_items = []
            for index, (part, count) in enumerate(fragments):
                rows_before += count
                work_items.append((count, min(1.0, rows_before / total_rows), index, part))
        else:
            work_items = None

        def source_items():
            if work_items is not None:
                yield from work_items
            else:
                for rows, fraction in record_batches(source, config.chunk_size):
                    yield len(rows), fraction, None, rows

        for item_rows, input_fraction, part_index, payload in source_items():
            if reads_seen < reads_to_skip:
                remaining = reads_to_skip - reads_seen
                if remaining >= item_rows:
                    reads_seen += item_rows
                    update_progress(
                        "resuming", last_input_fraction,
                        min(1.0, reads_seen / reads_to_skip) if reads_to_skip else 1.0,
                    )
                    continue
                if part_index is not None:
                    raise ConfigError(
                        "Correction resume state is not aligned to extracted fragments; "
                        "use --overwrite to restart safely."
                    )
                payload = payload[remaining:]
                reads_seen += remaining
            if part_index is not None:
                if not parts_checked and part_index != valid_writer.part:
                    raise ConfigError(
                        "Committed correction fragments do not match the extracted fragments; "
                        "use --overwrite to restart safely."
                    )
                parts_checked = True
                task = (part_index, str(payload), str(valid_parquet),
                        str(trace_parquet) if config.trace_enabled else None,
                        getattr(config, "parquet_compression", "snappy"),
                        getattr(config, "compression_level", 3))
            if executor is None:
                if part_index is not None:
                    append_correction_member(_correct_part(task, config, correctors, annotation), True)
                else:
                    append_correction_member(_correct_batch(payload, config, correctors, annotation))
                last_input_fraction = max(last_input_fraction, input_fraction)
                save_resume()
                update_progress("running", last_input_fraction)
                continue
            while len(pending) >= workers:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                collect_finished(done)
                commit_ready()
            if part_index is not None:
                future = executor.submit(_correct_part_worker, task)
            else:
                future = executor.submit(_correct_batch_worker, payload)
            pending[future] = (next_sequence, input_fraction)
            next_sequence += 1
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            collect_finished(done)
            commit_ready()
        valid_writer.close()
        if trace_writer is not None:
            trace_writer.close()
        save_resume()
        fields, stat_rows, result_summary = _correction_stats(config, counters, summary)
        missing_annotation = [{"barcode2_sequence": key.split(":", 1)[1], "reads_count": count}
                              for key, count in summary.items() if key.startswith("missing_annotation:")]
        missing_annotation.sort(key=lambda row: (-row["reads_count"], row["barcode2_sequence"]))
        write_json(stats_path, {"schema": 1, "sample": sample_name, "fields": fields,
                              "barcode_stats": stat_rows, "annotation": {
                                  "whitelist_reads": summary["barcode2_whitelist_reads"],
                                  "not_in_annotation_reads": sum(row["reads_count"] for row in missing_annotation)}})
        write_tsv(dirs["corrected"] / f"{sample_name}.barcode2_not_in_annotation.tsv.gz",
                  ["barcode2_sequence", "reads_count"], missing_annotation, config.compression_level)
        if not keep_resume:
            resume_path.unlink(missing_ok=True)
        update_progress("completed", 1.0)
    except Exception:
        update_progress("failed", last_input_fraction)
        raise
    finally:
        if executor is not None:
            executor.shutdown()
    outputs = [valid_parquet, stats_path,
               dirs["corrected"] / f"{sample_name}.barcode2_not_in_annotation.tsv.gz"] + ([trace_parquet] if config.trace_enabled else [])
    result_summary["wall_seconds"] = previous_elapsed + time.monotonic() - started
    result_summary["reads_per_second"] = (
        summary["total_reads"] / result_summary["wall_seconds"]
        if result_summary["wall_seconds"] else 0.0
    )
    result_summary["requested_workers"] = requested_workers
    result_summary["workers"] = workers
    result_summary["chunk_size"] = config.chunk_size
    return outputs, result_summary


def correction_qc_path(directory: Path, sample: str):
    current = directory / f"{sample}.correction_qc.json"
    legacy = directory / f"{sample}.barcode_correction_stats.tsv"
    return legacy if not current.exists() and legacy.exists() else current


def load_correction_stats(path: Path):
    if path.suffix != ".json":
        from .io_utils import open_tsv
        return list(open_tsv(path))
    payload = json.loads(path.read_text())
    if payload.get("schema") != 1 or not isinstance(payload.get("barcode_stats"), list):
        raise ConfigError(f"Invalid correction QC: {path}")
    return [{key: str(row.get(key, "")) for key in payload["fields"]}
            for row in payload["barcode_stats"]]
