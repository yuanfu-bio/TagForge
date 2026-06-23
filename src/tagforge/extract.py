from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from collections import Counter
from dataclasses import dataclass, replace
from datetime import datetime, timedelta
from functools import lru_cache
from itertools import repeat
import time
from typing import Optional

from . import __version__
from .config import ConfigError, SegmentConfig, TagForgeConfig
from .fastq import paired_fastq_batches
from .io_utils import atomic_text, sample_dirs, write_tsv
from .logging_utils import sample_logger


@dataclass(frozen=True)
class SegmentExtractionResult:
    segment_name: str
    raw_sequence: str
    success: bool
    failure_reason: str = ""
    extraction_method: str = ""
    linker_failure_reason: str = ""
    linker_elapsed_seconds: float = 0.0
    left_linker_matches: int = 0
    right_linker_matches: int = 0
    linker_candidate_pairs: int = 0
    selected_linker_gap: int = -1


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


@lru_cache(maxsize=512)
def _cutadapt_adapter(linker: str, max_mismatch: int):
    try:
        from cutadapt.adapters import AnywhereAdapter
    except ImportError as exc:
        raise RuntimeError(
            "cutadapt is required for linker extraction. Activate the TagForge Conda environment."
        ) from exc
    return AnywhereAdapter(
        sequence=linker,
        # A rate is portable across Cutadapt releases. Full-length overlap and
        # disabled indels make this equivalent to max_mismatch substitutions.
        max_errors=max_mismatch / len(linker),
        min_overlap=len(linker),
        indels=False,
    )


def _find_linker(sequence: str, linker: str, max_mismatch: int, start: int = 0) -> Optional[int]:
    """Locate a linker with cutadapt's matching engine."""
    match = _cutadapt_adapter(linker, max_mismatch).match_to(sequence[start:])
    return None if match is None else start + match.rstart


def _find_all_linkers(sequence: str, linker: str, max_mismatch: int) -> list[int]:
    """Enumerate all full-length substitution-only matches, verified by Cutadapt."""
    adapter = _cutadapt_adapter(linker, max_mismatch)
    positions = []
    width = len(linker)
    for position in range(0, len(sequence) - width + 1):
        window = sequence[position:position + width]
        if hamming(window, linker) > max_mismatch:
            continue
        match = adapter.match_to(window)
        if match is not None and match.rstart == 0 and match.rstop == width:
            positions.append(position)
    return positions


def _closest_linker_pair(sequence: str, left: str, right: str, max_mismatch: int):
    left_positions = _find_all_linkers(sequence, left, max_mismatch)
    right_positions = _find_all_linkers(sequence, right, max_mismatch)
    pairs = [
        (right_position - (left_position + len(left)), left_position, right_position)
        for left_position in left_positions
        for right_position in right_positions
        if right_position >= left_position + len(left)
    ]
    selected = min(pairs, key=lambda item: (item[0], item[1], item[2])) if pairs else None
    return selected, len(left_positions), len(right_positions), len(pairs)


def _extract_fixed(sequence: str, segment: SegmentConfig, extra: int) -> SegmentExtractionResult:
    assert segment.start is not None
    required_end = segment.start + segment.length
    if required_end > len(sequence):
        return SegmentExtractionResult(
            segment.name, "", False,
            f"coordinate_out_of_range:{segment.start}:{required_end}",
            "fixed",
        )
    window_start = max(0, segment.start - extra)
    window_end = min(len(sequence), required_end + extra)
    return SegmentExtractionResult(
        segment.name, sequence[window_start:window_end], True, extraction_method="fixed"
    )


def _extract_linker_impl(sequence: str, segment: SegmentConfig) -> SegmentExtractionResult:
    left_end = 0
    right_start = len(sequence)
    left_matches = right_matches = candidate_pairs = 0
    selected_gap = -1
    if segment.left_linker and segment.right_linker:
        selected, left_matches, right_matches, candidate_pairs = _closest_linker_pair(
            sequence, segment.left_linker, segment.right_linker, segment.linker_max_mismatch
        )
        if not left_matches:
            return SegmentExtractionResult(
                segment.name, "", False, "left_linker_not_found", "linker",
                left_linker_matches=left_matches, right_linker_matches=right_matches,
                linker_candidate_pairs=candidate_pairs,
            )
        if selected is None:
            return SegmentExtractionResult(
                segment.name, "", False, "right_linker_not_found_after_left", "linker",
                left_linker_matches=left_matches, right_linker_matches=right_matches,
                linker_candidate_pairs=candidate_pairs,
            )
        selected_gap, left_position, right_start = selected
        left_end = left_position + len(segment.left_linker)
    elif segment.left_linker:
        pos = _find_linker(sequence, segment.left_linker, segment.linker_max_mismatch)
        if pos is None:
            return SegmentExtractionResult(segment.name, "", False, "left_linker_not_found", "linker")
        left_end = pos + len(segment.left_linker)
        left_matches = 1
    elif segment.right_linker:
        pos = _find_linker(sequence, segment.right_linker, segment.linker_max_mismatch, left_end)
        if pos is None:
            return SegmentExtractionResult(segment.name, "", False, "right_linker_not_found", "linker")
        right_start = pos
        right_matches = 1
    if segment.right_linker:
        pos = right_start
        if segment.left_linker:
            raw = sequence[left_end:pos]
        else:
            raw = sequence[max(0, pos - segment.length):pos]
    else:
        raw = sequence[left_end:left_end + segment.length]
    if len(raw) != segment.length:
        return SegmentExtractionResult(
            segment.name, raw, False, f"unexpected_length:{len(raw)}", "linker",
            left_linker_matches=left_matches, right_linker_matches=right_matches,
            linker_candidate_pairs=candidate_pairs, selected_linker_gap=selected_gap,
        )
    return SegmentExtractionResult(
        segment.name, raw, True, extraction_method="linker",
        left_linker_matches=left_matches, right_linker_matches=right_matches,
        linker_candidate_pairs=candidate_pairs, selected_linker_gap=selected_gap,
    )


def _extract_linker(sequence: str, segment: SegmentConfig) -> SegmentExtractionResult:
    started = time.perf_counter()
    result = _extract_linker_impl(sequence, segment)
    elapsed = time.perf_counter() - started
    return replace(
        result,
        linker_failure_reason=result.failure_reason if not result.success else "",
        linker_elapsed_seconds=elapsed,
    )


def extract_segment(sequence: str, segment: SegmentConfig) -> SegmentExtractionResult:
    extra = (
        segment.correction.max_shift
        if segment.target != "umi" and segment.correction.enabled and segment.correction.allow_shift
        else 0
    )
    if segment.method == "fixed":
        return _extract_fixed(sequence, segment, extra)
    linker_result = _extract_linker(sequence, segment)
    if linker_result.success or segment.method == "linker":
        return linker_result
    # Combined mode is fallback, not coordinate composition: fixed extraction
    # is attempted on the original read only when linker extraction failed.
    fixed_result = _extract_fixed(sequence, segment, extra)
    if fixed_result.success:
        return replace(
            fixed_result,
            linker_failure_reason=linker_result.failure_reason,
            linker_elapsed_seconds=linker_result.linker_elapsed_seconds,
            left_linker_matches=linker_result.left_linker_matches,
            right_linker_matches=linker_result.right_linker_matches,
            linker_candidate_pairs=linker_result.linker_candidate_pairs,
            selected_linker_gap=linker_result.selected_linker_gap,
        )
    return SegmentExtractionResult(
        segment.name, fixed_result.raw_sequence, False,
        f"linker_failed:{linker_result.failure_reason};fixed_failed:{fixed_result.failure_reason}",
        "failed",
        linker_result.failure_reason,
        linker_result.linker_elapsed_seconds,
        linker_result.left_linker_matches,
        linker_result.right_linker_matches,
        linker_result.linker_candidate_pairs,
        linker_result.selected_linker_gap,
    )


EXTRACT_TRAILING_FIELDS = ["methods", "status", "failure_reason"]


def extract_fields(config: TagForgeConfig):
    return [
        "read_id", config.segment_column("barcode1"),
        config.segment_column("barcode2"), config.segment_column("umi"),
        *EXTRACT_TRAILING_FIELDS,
    ]

EXTRACTION_STATS_FIELDS = [
    "sample", "segment", "target", "read", "configured_mode", "total_reads",
    "linker_attempted", "linker_success", "linker_failed", "linker_success_rate",
    "left_linker_not_found", "right_linker_not_found", "unexpected_length",
    "reads_with_multiple_left_matches", "reads_with_multiple_right_matches",
    "reads_with_multiple_candidate_pairs", "linker_candidate_pairs_total",
    "selected_gap_min", "selected_gap_max",
    "fixed_attempted", "fixed_success", "fixed_failed", "fixed_rescue_rate",
    "final_success", "final_success_rate", "workers", "parallel_backend",
    "wall_seconds", "cutadapt_cpu_seconds", "reads_per_second",
]


def _extract_record(record, segments, segment_columns=None):
    segment_columns = segment_columns or {
        "barcode1": "barcode1_segments",
        "barcode2": "barcode2_segments",
        "umi": "umi_segments",
    }
    by_target = {"barcode1": {}, "barcode2": {}, "umi": {}}
    methods = {}
    runtime_stats = {}
    failures = []
    for segment in segments:
        sequence = record.r1_seq if segment.read == "R1" else record.r2_seq
        result = extract_segment(sequence, segment)
        by_target[segment.target][segment.name] = result.raw_sequence
        methods[segment.name] = result.extraction_method or "unknown"
        runtime_stats[segment.name] = {
            "success": result.success,
            "method": result.extraction_method,
            "linker_failure_reason": result.linker_failure_reason,
            "linker_elapsed_seconds": result.linker_elapsed_seconds,
            "left_linker_matches": result.left_linker_matches,
            "right_linker_matches": result.right_linker_matches,
            "linker_candidate_pairs": result.linker_candidate_pairs,
            "selected_linker_gap": result.selected_linker_gap,
        }
        if not result.success:
            failures.append(f"{segment.name}:{result.failure_reason}")
    valid = not failures
    row = {
        "read_id": record.read_id,
        # Comma-separated segment values and method codes follow config order.
        # DNA/IUPAC sequences cannot contain commas, so JSON adds no value here.
        segment_columns["barcode1"]: ",".join(by_target["barcode1"].values()),
        segment_columns["barcode2"]: ",".join(by_target["barcode2"].values()),
        segment_columns["umi"]: ",".join(by_target["umi"].values()),
        "methods": "".join(
            {"fixed": "F", "linker": "L", "failed": "X"}.get(method, "?")
            for method in methods.values()
        ),
        "status": "success" if valid else "failed",
        "failure_reason": ";".join(failures),
    }
    return row, runtime_stats


def decode_segment_payload(payload: str, segments) -> dict[str, str]:
    """Decode comma-separated segment sequences in configuration order."""
    if any(marker in payload for marker in '{}[]"'):
        raise ValueError(
            "Legacy JSON segment payload is unsupported; rerun 'tagforge extract --overwrite'"
        )
    values = payload.split(",")
    if len(values) != len(segments):
        raise ValueError(
            f"Segment payload length mismatch: expected {len(segments)}, found {len(values)}"
        )
    return {segment.name: value for segment, value in zip(segments, values)}


def decode_method_payload(payload: str, segments) -> dict[str, str]:
    """Decode compact F/L/X codes in configuration order."""
    stripped = (payload or "").strip()
    if len(stripped) != len(segments):
        raise ValueError(
            f"Extraction-method payload length mismatch: expected {len(segments)}, "
            f"found {len(stripped)}"
        )
    names = {"F": "fixed", "L": "linker", "X": "failed", "?": "unknown"}
    unknown = sorted(set(stripped) - set(names))
    if unknown:
        raise ValueError(f"Unknown extraction method code(s): {', '.join(unknown)}")
    return {segment.name: names[code] for segment, code in zip(segments, stripped)}


def _resume_fingerprint(config: TagForgeConfig, sample) -> str:
    digest = hashlib.sha256()
    digest.update(config.path.read_bytes())
    for path in (sample.r1, sample.r2):
        stat = path.stat()
        digest.update(f"{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    return digest.hexdigest()


def extract_sample(config: TagForgeConfig, sample_name: str, resume: bool = True):
    sample = config.sample(sample_name)
    dirs = sample_dirs(config.output_dir, sample_name)
    output = dirs["extracted"] / f"{sample_name}.extracted.tsv.gz"
    stats_output = dirs["extracted"] / f"{sample_name}.extraction_stats.tsv"
    preview_output = dirs["extracted"] / f"{sample_name}.extracted.preview.tsv"
    progress_output = dirs["logs"] / f"{sample_name}.extraction_progress.tsv"
    resume_path = dirs["checkpoint"] / f"{sample_name}.extract.resume.json"
    tmp_output = output.with_name(output.name + ".tmp")
    fields = extract_fields(config)
    segment_columns = {
        "barcode1": config.segment_column("barcode1"),
        "barcode2": config.segment_column("barcode2"),
        "umi": config.segment_column("umi"),
    }
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    started = time.monotonic()
    totals = {"total_reads": 0, "extracted_reads": 0}
    segment_counters = {segment.name: Counter() for segment in config.segments}
    previous_elapsed = 0.0
    last_input_fraction = 0.0
    resume_fraction_known = True
    fingerprint = _resume_fingerprint(config, sample)

    def save_resume():
        state = {
            "schema": 2, "tagforge_version": __version__, "fingerprint": fingerprint,
            "reads_completed": totals["total_reads"],
            "extracted_reads": totals["extracted_reads"],
            "safe_output_bytes": tmp_output.stat().st_size,
            "elapsed_seconds": previous_elapsed + time.monotonic() - started,
            "input_fraction": last_input_fraction,
            "parallel_backend": totals.get("parallel_backend", "unknown"),
            "segment_counters": {name: dict(counter) for name, counter in segment_counters.items()},
        }
        with atomic_text(resume_path) as state_handle:
            json.dump(state, state_handle, separators=(",", ":"))

    resumed = False
    if resume_path.exists() and resume:
        state = json.loads(resume_path.read_text(encoding="utf-8"))
        if state.get("schema") != 2 or state.get("tagforge_version") != __version__ or state.get("fingerprint") != fingerprint:
            raise ConfigError(
                f"Extraction resume state does not match current inputs/config: {resume_path}. "
                "Use --overwrite to restart safely."
            )
        if not tmp_output.exists():
            raise ConfigError(f"Extraction resume output is missing: {tmp_output}. Use --overwrite to restart.")
        safe_bytes = int(state["safe_output_bytes"])
        with open(tmp_output, "r+b") as raw_handle:
            raw_handle.truncate(safe_bytes)
        totals["total_reads"] = int(state["reads_completed"])
        totals["extracted_reads"] = int(state["extracted_reads"])
        previous_elapsed = float(state.get("elapsed_seconds", 0.0))
        resume_fraction_known = "input_fraction" in state
        last_input_fraction = float(state.get("input_fraction", 0.0))
        for name, values in state.get("segment_counters", {}).items():
            if name in segment_counters:
                segment_counters[name].update(values)
        resumed = True
        logger.info("extract_resume\treads=%s\tsafe_output_bytes=%s", totals["total_reads"], safe_bytes)
    else:
        if tmp_output.exists() and resume and not resume_path.exists():
            raise ConfigError(
                f"Untracked extraction temporary file cannot be resumed safely: {tmp_output}. "
                "Use --overwrite to restart."
            )
        tmp_output.unlink(missing_ok=True); resume_path.unlink(missing_ok=True)
        with gzip.open(tmp_output, "wt", encoding="utf-8", newline="", compresslevel=config.compression_level) as initial:
            csv.DictWriter(initial, fieldnames=fields, delimiter="\t", lineterminator="\n").writeheader()
        save_resume()
    progress_fields = [
        "status", "reads_completed", "input_percent", "elapsed_seconds",
        "reads_per_second", "eta_seconds", "estimated_finish", "resume_skip_percent",
        "output_tmp_bytes",
    ]

    def update_progress(status: str, input_fraction: float | None, resume_skip: float | None = None):
        elapsed = previous_elapsed + time.monotonic() - started
        rate = totals["total_reads"] / elapsed if elapsed else 0.0
        eta = (
            elapsed * (1.0 - input_fraction) / input_fraction
            if input_fraction is not None and 0 < input_fraction < 1 else 0.0
        )
        finish = (
            (datetime.now().astimezone() + timedelta(seconds=eta)).isoformat(timespec="seconds")
            if eta else ""
        )
        row = {
            "status": status, "reads_completed": totals["total_reads"],
            "input_percent": f"{input_fraction * 100:.2f}" if input_fraction is not None else "NA",
            "elapsed_seconds": f"{elapsed:.2f}", "reads_per_second": f"{rate:.2f}",
            "eta_seconds": f"{eta:.2f}" if eta else "", "estimated_finish": finish,
            "resume_skip_percent": f"{resume_skip * 100:.2f}" if resume_skip is not None else "",
            "output_tmp_bytes": tmp_output.stat().st_size if tmp_output.exists() else 0,
        }
        write_tsv(progress_output, progress_fields, [row])
        logger.info(
            "extract_progress\treads=%s\tinput_percent=%s\treads_per_second=%s\t"
            "eta_seconds=%s\testimated_finish=%s\tresume_skip_percent=%s\tstatus=%s",
            row["reads_completed"], row["input_percent"], row["reads_per_second"],
            row["eta_seconds"] or "NA", row["estimated_finish"] or "NA",
            row["resume_skip_percent"] or "NA", status,
        )
        if status == "running":
            eta_text = f"{row['eta_seconds']}s" if row["eta_seconds"] else "NA"
            print(
                f"{sample_name} extract: {row['reads_completed']:,} reads, "
                f"input≈{row['input_percent']}%, speed={float(row['reads_per_second']):,.0f}/s, "
                f"ETA={eta_text}",
                flush=True,
            )

    executor = None
    preview_mode = "a" if resumed and preview_output.exists() else "w"
    preview_handle = open(preview_output, preview_mode, encoding="utf-8", newline="")
    preview_writer = csv.DictWriter(
        preview_handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    if preview_mode == "w": preview_writer.writeheader()
    preview_handle.flush()
    try:
        if config.threads == 1:
            totals["parallel_backend"] = "serial"
        else:
            try:
                executor = ProcessPoolExecutor(max_workers=config.threads)
                totals["parallel_backend"] = "process"
            except (PermissionError, NotImplementedError):
                executor = ThreadPoolExecutor(max_workers=config.threads)
                totals["parallel_backend"] = "thread"
        update_progress(
            "resuming" if resumed else "running",
            last_input_fraction if (not resumed or resume_fraction_known) else None,
            0.0 if resumed else None,
        )
        reads_to_skip = totals["total_reads"]
        reads_seen = 0
        for batch, input_fraction in paired_fastq_batches(
            sample.r1, sample.r2, config.chunk_size
        ):
            if reads_seen < reads_to_skip:
                remaining = reads_to_skip - reads_seen
                if remaining >= len(batch):
                    reads_seen += len(batch)
                    skip_fraction = min(1.0, reads_seen / reads_to_skip) if reads_to_skip else 1.0
                    if reads_seen >= reads_to_skip:
                        last_input_fraction = max(last_input_fraction, input_fraction)
                        resume_fraction_known = True
                    update_progress(
                        "resuming",
                        last_input_fraction if resume_fraction_known else None,
                        skip_fraction,
                    )
                    continue
                batch = batch[remaining:]
                reads_seen += remaining
            if executor is None:
                extracted = (
                    _extract_record(record, config.segments, segment_columns)
                    for record in batch
                )
            else:
                map_chunksize = max(1, min(1000, len(batch) // config.threads))
                extracted = executor.map(
                    _extract_record, batch, repeat(tuple(config.segments), len(batch)),
                    repeat(segment_columns, len(batch)),
                    chunksize=map_chunksize,
                )
            with gzip.open(
                tmp_output, "at", encoding="utf-8", newline="",
                compresslevel=config.compression_level,
            ) as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                for row, runtime_stats in extracted:
                    totals["total_reads"] += 1
                    if row["status"] == "success":
                        totals["extracted_reads"] += 1
                    for segment in config.segments:
                        stat = runtime_stats[segment.name]
                        counter = segment_counters[segment.name]
                        counter["total_reads"] += 1
                        counter["final_success"] += int(stat["success"])
                        counter["cutadapt_microseconds"] += round(stat["linker_elapsed_seconds"] * 1_000_000)
                        if segment.method in {"linker", "linker_fixed"}:
                            counter["linker_attempted"] += 1
                            linker_success = stat["method"] == "linker" and stat["success"]
                            counter["linker_success"] += int(linker_success)
                            counter["linker_failed"] += int(not linker_success)
                            if stat["linker_failure_reason"]:
                                reason = stat["linker_failure_reason"].split(":", 1)[0]
                                if reason.startswith("right_linker_not_found"):
                                    reason = "right_linker_not_found"
                                counter[reason] += 1
                            counter["reads_with_multiple_left_matches"] += int(stat["left_linker_matches"] > 1)
                            counter["reads_with_multiple_right_matches"] += int(stat["right_linker_matches"] > 1)
                            counter["reads_with_multiple_candidate_pairs"] += int(stat["linker_candidate_pairs"] > 1)
                            counter["linker_candidate_pairs_total"] += stat["linker_candidate_pairs"]
                            if stat["selected_linker_gap"] >= 0:
                                counter["selected_gap_count"] += 1
                                counter["selected_gap_sum"] += stat["selected_linker_gap"]
                                gap_key = f"selected_gap_{stat['selected_linker_gap']}"
                                counter[gap_key] += 1
                        fixed_attempted = segment.method == "fixed" or (
                            segment.method == "linker_fixed" and stat["method"] != "linker"
                        )
                        if fixed_attempted:
                            counter["fixed_attempted"] += 1
                            fixed_success = stat["method"] == "fixed" and stat["success"]
                            counter["fixed_success"] += int(fixed_success)
                            counter["fixed_failed"] += int(not fixed_success)
                    writer.writerow(row)
                    if totals["total_reads"] <= config.extraction_preview_reads:
                        preview_writer.writerow(row)
            preview_handle.flush()
            last_input_fraction = max(last_input_fraction, input_fraction)
            save_resume()
            update_progress("running", last_input_fraction)
        os.replace(tmp_output, output)
        resume_path.unlink(missing_ok=True)
        update_progress("completed", 1.0)
    except Exception:
        update_progress("failed", 0.0)
        raise
    finally:
        preview_handle.close()
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    wall_seconds = previous_elapsed + time.monotonic() - started
    reads_per_second = totals["total_reads"] / wall_seconds if wall_seconds else 0.0
    stats_rows = []
    for segment in config.segments:
        counter = segment_counters[segment.name]
        linker_attempted = counter["linker_attempted"]
        fixed_attempted = counter["fixed_attempted"]
        total_reads = counter["total_reads"]
        stats_rows.append({
            "sample": sample_name, "segment": segment.name,
            "target": segment.target_name or segment.target,
            "read": segment.read, "configured_mode": segment.method,
            "total_reads": total_reads, "linker_attempted": linker_attempted,
            "linker_success": counter["linker_success"], "linker_failed": counter["linker_failed"],
            "linker_success_rate": counter["linker_success"] / linker_attempted if linker_attempted else "",
            "left_linker_not_found": counter["left_linker_not_found"],
            "right_linker_not_found": counter["right_linker_not_found"],
            "unexpected_length": counter["unexpected_length"],
            "reads_with_multiple_left_matches": counter["reads_with_multiple_left_matches"],
            "reads_with_multiple_right_matches": counter["reads_with_multiple_right_matches"],
            "reads_with_multiple_candidate_pairs": counter["reads_with_multiple_candidate_pairs"],
            "linker_candidate_pairs_total": counter["linker_candidate_pairs_total"],
            "selected_gap_min": min(
                (int(key.removeprefix("selected_gap_")) for key in counter if key.startswith("selected_gap_") and key != "selected_gap_count" and key != "selected_gap_sum"),
                default="",
            ),
            "selected_gap_max": max(
                (int(key.removeprefix("selected_gap_")) for key in counter if key.startswith("selected_gap_") and key != "selected_gap_count" and key != "selected_gap_sum"),
                default="",
            ),
            "fixed_attempted": fixed_attempted, "fixed_success": counter["fixed_success"],
            "fixed_failed": counter["fixed_failed"],
            "fixed_rescue_rate": (
                counter["fixed_success"] / fixed_attempted
                if segment.method == "linker_fixed" and fixed_attempted else ""
            ),
            "final_success": counter["final_success"],
            "final_success_rate": counter["final_success"] / total_reads if total_reads else 0,
            "workers": config.threads, "parallel_backend": totals["parallel_backend"],
            "wall_seconds": f"{wall_seconds:.6f}",
            "cutadapt_cpu_seconds": f"{counter['cutadapt_microseconds'] / 1_000_000:.6f}",
            "reads_per_second": f"{reads_per_second:.3f}",
        })
    write_tsv(stats_output, EXTRACTION_STATS_FIELDS, stats_rows)
    totals.update({
        "wall_seconds": wall_seconds, "reads_per_second": reads_per_second,
        "cutadapt_cpu_seconds": sum(c["cutadapt_microseconds"] for c in segment_counters.values()) / 1_000_000,
        "segment_stats": stats_rows,
    })
    return [output, stats_output], totals
