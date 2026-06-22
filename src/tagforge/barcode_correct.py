from __future__ import annotations

import csv
import json
from collections import Counter
from contextlib import nullcontext
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Optional

from .config import CorrectionConfig, SegmentConfig, TagForgeConfig, ConfigError
from .extract import hamming
from .io_utils import atomic_text, open_tsv, sample_dirs, write_tsv


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
    def correct(self, raw: str) -> CorrectionResult:
        raw = raw.upper()
        if not self.config.enabled:
            return CorrectionResult(raw, raw, raw, "disabled", 0, 0, "disabled")
        candidates = []
        max_shift = self.config.max_shift if self.config.allow_shift else 0
        max_mismatch = self.config.max_mismatch if self.config.allow_mismatch else 0
        for shift in range(max_shift + 1):
            shifted = raw[shift:shift + self.length]
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
            return CorrectionResult(raw, raw[:self.length], "", "failed", -1, -1, "failed")
        score = min(
            ((x[0] > 0) + (x[1] > 0), x[0] + x[1], x[1], x[0])
            for x in candidates
        )
        best = [x for x in candidates if ((x[0] > 0) + (x[1] > 0), x[0] + x[1], x[1], x[0]) == score]
        sequences = {x[2] for x in best}
        if len(sequences) != 1:
            return CorrectionResult(raw, best[0][3], "", "ambiguous", best[0][0], best[0][1], "failed", True)
        shift, mismatch, corrected, shifted = best[0]
        kind = "exact" if not shift and not mismatch else (
            "shift_and_mismatch" if shift and mismatch else "shift_only" if shift else "mismatch_only")
        return CorrectionResult(raw, shifted, corrected, "valid", shift, mismatch, kind)


def load_whitelist(path: Path):
    values = [line.strip().upper() for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not values:
        raise ConfigError(f"Whitelist is empty: {path}")
    duplicates = len(values) - len(set(values))
    if duplicates:
        raise ConfigError(f"Whitelist contains {duplicates} duplicate entries: {path}")
    return values


def load_fb_annotation(config: TagForgeConfig) -> Dict[str, str]:
    with open(config.fb_info, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {config.fb_id_column, config.fb_sequence_column, config.fb_name_column}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ConfigError(f"FB annotation missing column(s): {', '.join(sorted(missing))}")
        mapping, names = {}, set()
        for row in reader:
            sequence = row[config.fb_sequence_column].strip().upper()
            name = row[config.fb_name_column].strip()
            if sequence in mapping:
                raise ConfigError(f"Duplicate FB sequence in annotation: {sequence}")
            if name in names and not config.allow_duplicate_names:
                raise ConfigError(f"Duplicate antibody name in annotation: {name}")
            mapping[sequence] = name; names.add(name)
    return mapping


TRACE_FIELDS = ["read_id", "segment_name", "target_type", "raw_sequence", "shifted_sequence",
                "corrected_sequence", "whitelist_hit", "correction_status", "shift_distance",
                "mismatch_distance", "correction_type"]
VALID_FIELDS = ["read_id", "barcode1", "barcode2_sequence", "barcode2_name", "umi", "correction_summary"]


def correct_sample(config: TagForgeConfig, sample_name: str):
    dirs = sample_dirs(config.output_dir, sample_name)
    source = dirs["extracted"] / f"{sample_name}.extracted.tsv.gz"
    if not source.is_file():
        raise FileNotFoundError(f"Extraction output missing for {sample_name}: {source}")
    trace_path = dirs["corrected"] / f"{sample_name}.barcode_correction_trace.tsv.gz"
    stats_path = dirs["corrected"] / f"{sample_name}.barcode_correction_stats.tsv"
    valid_path = dirs["detail"] / f"{sample_name}.valid_reads.tsv.gz"
    annotation = load_fb_annotation(config)
    correctors = {}
    for segment in config.segments:
        if segment.target == "umi":
            continue
        if segment.whitelist:
            correctors[segment.name] = WhitelistCorrector(load_whitelist(segment.whitelist), segment.correction, segment.length)
        elif segment.correction.enabled:
            raise ConfigError(f"Segment {segment.name}: correction enabled but whitelist is missing")
        else:
            correctors[segment.name] = WhitelistCorrector([], segment.correction, segment.length)
    counters = {segment.name: Counter() for segment in config.segments if segment.target != "umi"}
    summary = Counter()
    with atomic_text(valid_path, config.compression_level) as valid_handle:
        valid_writer = csv.DictWriter(valid_handle, fieldnames=VALID_FIELDS, delimiter="\t", lineterminator="\n")
        valid_writer.writeheader()
        trace_context = atomic_text(trace_path, config.compression_level) if config.trace_enabled else nullcontext(None)
        with trace_context as trace_handle:
            trace_writer = csv.DictWriter(trace_handle, fieldnames=TRACE_FIELDS, delimiter="\t", lineterminator="\n") if trace_handle else None
            if trace_writer: trace_writer.writeheader()
            for row in open_tsv(source):
                summary["total_reads"] += 1
                if row["extraction_status"] != "success":
                    continue
                summary["extracted_reads"] += 1
                extraction_methods = json.loads(row.get("segment_extraction_methods") or "{}")
                values = {}
                for target, col in (("barcode1", "barcode1_segment_raw_values"), ("barcode2", "barcode2_segment_raw_values")):
                    raw_values = json.loads(row[col]); corrected_values = []
                    for segment in (s for s in config.segments if s.target == target):
                        result = correctors[segment.name].correct(raw_values.get(segment.name, ""))
                        corrected_values.append(result.corrected_sequence)
                        count = counters[segment.name]
                        count["extracted_reads"] += 1; count[result.correction_type] += 1
                        source_method = extraction_methods.get(segment.name, "unknown")
                        count[f"{source_method}_extracted"] += 1
                        count[f"mismatch_{max(result.mismatch_distance, 0)}"] += int(result.mismatch_distance >= 0)
                        count[f"shift_{max(result.shift_distance, 0)}"] += int(result.shift_distance >= 0)
                        if result.success:
                            count["valid_reads"] += 1
                            count[f"{source_method}_valid"] += 1
                        else: count["invalid_reads"] += 1
                        if trace_writer:
                            trace_writer.writerow({"read_id": row["read_id"], "segment_name": segment.name,
                                "target_type": target, "raw_sequence": result.raw_sequence,
                                "shifted_sequence": result.shifted_sequence, "corrected_sequence": result.corrected_sequence,
                                "whitelist_hit": str(result.success).lower(), "correction_status": result.correction_status,
                                "shift_distance": result.shift_distance, "mismatch_distance": result.mismatch_distance,
                                "correction_type": result.correction_type})
                    values[target] = "".join(corrected_values) if all(corrected_values) else ""
                umi_values = json.loads(row["umi_segment_raw_values"])
                umi = "".join(umi_values.get(s.name, "") for s in config.segments if s.target == "umi")
                fb_name = annotation.get(values["barcode2"], "")
                if values["barcode1"]: summary["barcode1_valid"] += 1
                if values["barcode2"] and fb_name: summary["barcode2_valid"] += 1
                if values["barcode1"] and fb_name and umi and "N" not in umi:
                    summary["combined_valid"] += 1
                    valid_writer.writerow({"read_id": row["read_id"], "barcode1": values["barcode1"],
                        "barcode2_sequence": values["barcode2"], "barcode2_name": fb_name, "umi": umi,
                        "correction_summary": "valid"})
    stat_rows = []
    total = summary["total_reads"]
    for segment in (s for s in config.segments if s.target != "umi"):
        c = counters[segment.name]
        row = {"scope": segment.name, "target_type": segment.target, "total_reads": total,
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
        for i in range(max(segment.correction.max_shift, 2) + 1): row[f"shift_{i}_count"] = c[f"shift_{i}"]
        stat_rows.append(row)
    stat_rows.extend([
        {"scope": "final_barcode1", "target_type": "barcode1", "total_reads": total, "valid_reads": summary["barcode1_valid"], "valid_rate": summary["barcode1_valid"] / total if total else 0},
        {"scope": "final_barcode2", "target_type": "barcode2", "total_reads": total, "valid_reads": summary["barcode2_valid"], "valid_rate": summary["barcode2_valid"] / total if total else 0},
        {"scope": "combined", "target_type": "all", "total_reads": total, "valid_reads": summary["combined_valid"], "valid_rate": summary["combined_valid"] / total if total else 0},
    ])
    all_fields = ["scope", "target_type", "total_reads", "extracted_reads", "valid_reads", "invalid_reads", "valid_rate",
                  "exact_count", "mismatch_only_count", "shift_only_count", "shift_and_mismatch_count", "failed_count", "disabled_count"]
    dynamic = sorted({k for row in stat_rows for k in row if k not in all_fields})
    write_tsv(stats_path, all_fields + dynamic, stat_rows)
    result_summary = dict(summary)
    result_summary["segment_method_qc"] = [
        {key: row.get(key, "") for key in (
            "scope", "target_type", "linker_extracted_count", "linker_barcode_valid_count",
            "linker_barcode_valid_rate", "fixed_extracted_count", "fixed_barcode_valid_count",
            "fixed_barcode_valid_rate",
        )}
        for row in stat_rows if row.get("scope") not in {"final_barcode1", "final_barcode2", "combined"}
    ]
    return [valid_path, stats_path] + ([trace_path] if config.trace_enabled else []), result_summary
