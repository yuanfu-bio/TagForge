from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timedelta
from functools import lru_cache
from itertools import chain
from pathlib import Path
from typing import Dict, Iterable, Optional

from . import __version__
from .config import CorrectionConfig, SegmentConfig, TagForgeConfig, ConfigError
from .extract import decode_method_payload, decode_segment_payload, hamming
from .io_utils import atomic_text, sample_dirs, tsv_batches, write_tsv
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
        self.config = config
        self.length = length

    @lru_cache(maxsize=200000)
    def correct(self, raw: str, anchor: int = 0) -> CorrectionResult:
        raw = raw.upper()
        if not self.config.enabled:
            sequence = raw[anchor:anchor + self.length]
            return CorrectionResult(raw, sequence, sequence, "disabled", 0, 0, "disabled")
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
            if max_mismatch:
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
            row, config, correctors, annotation, counters, summary
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


TRACE_FIELDS = ["read_id", "segment_name", "target_type", "raw_sequence", "shifted_sequence",
                "corrected_sequence", "whitelist_hit", "correction_status", "shift_distance",
                "mismatch_distance", "correction_type"]


def valid_fields(config: TagForgeConfig):
    return [
        "read_id", config.target_name("barcode1"),
        f"{config.target_name('barcode2')}_sequence",
        f"{config.target_name('barcode2')}_name",
        config.target_name("umi"), "correction_summary",
    ]


def _correction_fingerprint(config: TagForgeConfig, source: Path) -> str:
    digest = hashlib.sha256(config.path.read_bytes())
    dependencies = [source]
    if config.fb_info is not None:
        dependencies.append(config.fb_info)
    dependencies.extend(
        segment.whitelist for segment in config.segments if segment.whitelist is not None
    )
    for path in dependencies:
        stat = path.stat()
        digest.update(f"{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    return digest.hexdigest()


def _correct_row(row, config, correctors, annotation, counters, summary):
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
        return None, []
    summary["extracted_reads"] += 1
    extraction_methods = decode_method_payload(row["methods"], config.segments)
    values = {}
    trace_rows = []
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
            corrected_values.append(result.corrected_sequence)
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
    if values["barcode1"]:
        summary["barcode1_valid"] += 1
    if values["barcode2"] and fb_name:
        summary["barcode2_valid"] += 1
    if values["barcode1"] and fb_name and umi and "N" not in umi:
        summary["combined_valid"] += 1
        return {
            "read_id": row["read_id"], config.target_name("barcode1"): values["barcode1"],
            f"{config.target_name('barcode2')}_sequence": values["barcode2"],
            f"{config.target_name('barcode2')}_name": fb_name,
            config.target_name("umi"): umi, "correction_summary": "valid",
        }, trace_rows
    return None, trace_rows


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
        for row in stat_rows if row.get("scope") not in {
            f"final_{config.target_name('barcode1')}",
            f"final_{config.target_name('barcode2')}",
            "combined",
        }
    ]
    return all_fields + dynamic, stat_rows, result_summary


def correct_sample(config: TagForgeConfig, sample_name: str, resume: bool = True):
    dirs = sample_dirs(config.output_dir, sample_name)
    source = dirs["extracted"] / f"{sample_name}.extracted.tsv.gz"
    if not source.is_file():
        raise FileNotFoundError(f"Extraction output missing for {sample_name}: {source}")
    trace_path = dirs["corrected"] / f"{sample_name}.barcode_correction_trace.tsv.gz"
    stats_path = dirs["corrected"] / f"{sample_name}.barcode_correction_stats.tsv"
    valid_path = dirs["detail"] / f"{sample_name}.valid_reads.tsv.gz"
    progress_path = dirs["logs"] / f"{sample_name}.correction_progress.tsv"
    resume_path = dirs["checkpoint"] / f"{sample_name}.correct.resume.json"
    valid_tmp = valid_path.with_name(valid_path.name + ".tmp")
    trace_tmp = trace_path.with_name(trace_path.name + ".tmp")
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
    fingerprint = _correction_fingerprint(config, source)
    requested_workers = config.barcode_workers or config.threads
    workers = requested_workers
    executor = None
    fields = valid_fields(config)

    def save_resume():
        state = {
            "schema": 2, "tagforge_version": __version__, "fingerprint": fingerprint,
            "trace_enabled": config.trace_enabled,
            "reads_completed": summary["total_reads"],
            "safe_valid_bytes": valid_tmp.stat().st_size,
            "safe_trace_bytes": trace_tmp.stat().st_size if config.trace_enabled else 0,
            "elapsed_seconds": previous_elapsed + time.monotonic() - started,
            "input_fraction": last_input_fraction,
            "summary": dict(summary),
            "segment_counters": {
                name: dict(counter) for name, counter in counters.items()
            },
        }
        with atomic_text(resume_path) as state_handle:
            json.dump(state, state_handle, separators=(",", ":"))

    def restore_output(tmp_path, final_path, safe_bytes, label):
        if not tmp_path.exists() and final_path.exists() and final_path.stat().st_size == safe_bytes:
            os.replace(final_path, tmp_path)
        if not tmp_path.exists():
            raise ConfigError(
                f"Correction resume {label} output is missing: {tmp_path}. Use --overwrite to restart."
            )
        with open(tmp_path, "r+b") as raw_handle:
            raw_handle.truncate(safe_bytes)

    resumed = False
    if resume and resume_path.exists():
        state = json.loads(resume_path.read_text(encoding="utf-8"))
        if (
            state.get("schema") != 2
            or state.get("tagforge_version") != __version__
            or state.get("fingerprint") != fingerprint
            or bool(state.get("trace_enabled")) != config.trace_enabled
        ):
            raise ConfigError(
                f"Correction resume state does not match current inputs/config: {resume_path}. "
                "Use --overwrite to restart safely."
            )
        restore_output(valid_tmp, valid_path, int(state["safe_valid_bytes"]), "valid-read")
        if config.trace_enabled:
            restore_output(trace_tmp, trace_path, int(state["safe_trace_bytes"]), "trace")
        summary.update(state.get("summary", {}))
        for name, values in state.get("segment_counters", {}).items():
            if name in counters:
                counters[name].update(values)
        previous_elapsed = float(state.get("elapsed_seconds", 0.0))
        last_input_fraction = float(state.get("input_fraction", 0.0))
        resumed = True
        logger.info(
            "correction_resume\treads=%s\tsafe_valid_bytes=%s\tsafe_trace_bytes=%s",
            summary["total_reads"], state["safe_valid_bytes"], state.get("safe_trace_bytes", 0),
        )
    else:
        if resume and not resume_path.exists() and (valid_tmp.exists() or trace_tmp.exists()):
            raise ConfigError(
                "Untracked correction temporary output cannot be resumed safely. Use --overwrite to restart."
            )
        valid_tmp.unlink(missing_ok=True)
        trace_tmp.unlink(missing_ok=True)
        resume_path.unlink(missing_ok=True)
        with gzip.open(
            valid_tmp, "wt", encoding="utf-8", newline="", compresslevel=config.compression_level
        ) as handle:
            csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            ).writeheader()
        if config.trace_enabled:
            with gzip.open(
                trace_tmp, "wt", encoding="utf-8", newline="",
                compresslevel=config.compression_level,
            ) as handle:
                csv.DictWriter(
                    handle, fieldnames=TRACE_FIELDS, delimiter="\t", lineterminator="\n"
                ).writeheader()
        save_resume()

    progress_fields = [
        "status", "reads_completed", "extracted_reads", "combined_valid_reads",
        "input_percent", "elapsed_seconds", "reads_per_second", "eta_seconds",
        "estimated_finish", "resume_skip_percent", "valid_tmp_bytes", "trace_tmp_bytes",
    ]

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
            "valid_tmp_bytes": valid_tmp.stat().st_size if valid_tmp.exists() else 0,
            "trace_tmp_bytes": trace_tmp.stat().st_size if trace_tmp.exists() else 0,
        }
        write_tsv(progress_path, progress_fields, [progress])
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

    def append_correction_member(result):
        valid_rows, trace_rows, local_counters, local_summary = result
        with gzip.open(
            valid_tmp, "at", encoding="utf-8", newline="",
            compresslevel=config.compression_level,
        ) as valid_handle:
            valid_writer = csv.DictWriter(
                valid_handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            valid_writer.writerows(valid_rows)
        if config.trace_enabled:
            with gzip.open(
                trace_tmp, "at", encoding="utf-8", newline="",
                compresslevel=config.compression_level,
            ) as trace_handle:
                trace_writer = csv.DictWriter(
                    trace_handle, fieldnames=TRACE_FIELDS, delimiter="\t", lineterminator="\n"
                )
                trace_writer.writerows(trace_rows)
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
                append_correction_member(result)
                last_input_fraction = max(last_input_fraction, input_fraction)
                save_resume()
                update_progress("running", last_input_fraction)
                next_commit += 1

        for batch, input_fraction in tsv_batches(source, config.chunk_size):
            if reads_seen < reads_to_skip:
                remaining = reads_to_skip - reads_seen
                if remaining >= len(batch):
                    reads_seen += len(batch)
                    update_progress(
                        "resuming", last_input_fraction,
                        min(1.0, reads_seen / reads_to_skip) if reads_to_skip else 1.0,
                    )
                    continue
                batch = batch[remaining:]
                reads_seen += remaining
            if executor is None:
                append_correction_member(
                    _correct_batch(batch, config, correctors, annotation)
                )
                last_input_fraction = max(last_input_fraction, input_fraction)
                save_resume()
                update_progress("running", last_input_fraction)
                continue
            while len(pending) >= workers:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                collect_finished(done)
                commit_ready()
            future = executor.submit(_correct_batch_worker, batch)
            pending[future] = (next_sequence, input_fraction)
            next_sequence += 1
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            collect_finished(done)
            commit_ready()
        os.replace(valid_tmp, valid_path)
        if config.trace_enabled:
            os.replace(trace_tmp, trace_path)
        fields, stat_rows, result_summary = _correction_stats(config, counters, summary)
        write_tsv(stats_path, fields, stat_rows)
        resume_path.unlink(missing_ok=True)
        update_progress("completed", 1.0)
    except Exception:
        update_progress("failed", last_input_fraction)
        raise
    finally:
        if executor is not None:
            executor.shutdown()
    outputs = [valid_path, stats_path] + ([trace_path] if config.trace_enabled else [])
    result_summary["wall_seconds"] = previous_elapsed + time.monotonic() - started
    result_summary["reads_per_second"] = (
        summary["total_reads"] / result_summary["wall_seconds"]
        if result_summary["wall_seconds"] else 0.0
    )
    result_summary["requested_workers"] = requested_workers
    result_summary["workers"] = workers
    result_summary["chunk_size"] = config.chunk_size
    return outputs, result_summary
