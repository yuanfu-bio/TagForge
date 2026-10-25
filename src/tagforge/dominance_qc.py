"""PB assignment diagnostics using read support after observed CB correction."""
from __future__ import annotations

import gzip
from functools import lru_cache
from pathlib import Path

from .config import ConfigError


class DominanceQC:
    def __init__(self, threshold):
        self.threshold = float(threshold)
        self.counts = dict.fromkeys((
            "total_pbs", "single_cb_pbs", "single_cb_retained_pbs", "multi_cb_pbs",
            "resolved_multi_cb_pbs", "unresolved_multi_cb_pbs", "unassigned_pbs", "retained_pbs",
            "total_relations", "removed_secondary_relations", "removed_unresolved_relations",
            "supporting_reads", "retained_supporting_reads"), 0)
        self.histogram = [{"lower": i / 20, "upper": (i + 1) / 20,
                           "passed_pbs": 0, "failed_pbs": 0} for i in range(20)]

    def add(self, cb_count, total_reads, dominant_reads):
        cb_count, total_reads, dominant_reads = int(cb_count), int(total_reads), int(dominant_reads)
        if cb_count < 1 or not 0 < dominant_reads <= total_reads:
            raise ConfigError("Invalid PB support while computing dominance diagnostics")
        ratio = dominant_reads / total_reads
        passed = ratio >= self.threshold
        c = self.counts
        c["total_pbs"] += 1
        c["total_relations"] += cb_count
        c["supporting_reads"] += total_reads
        if passed:
            c["retained_pbs"] += 1
            c["retained_supporting_reads"] += dominant_reads
            c["removed_secondary_relations"] += cb_count - 1
        else:
            c["unassigned_pbs"] += 1
            c["removed_unresolved_relations"] += cb_count
        if cb_count == 1:
            c["single_cb_pbs"] += 1
            c["single_cb_retained_pbs"] += int(passed)
        else:
            c["multi_cb_pbs"] += 1
            c["resolved_multi_cb_pbs" if passed else "unresolved_multi_cb_pbs"] += 1
            self.histogram[min(int(ratio * 20), 19)]["passed_pbs" if passed else "failed_pbs"] += 1

    def payload(self):
        return {"schema": 1, "basis": "post_cb_correction", "dominance_threshold": self.threshold,
                **self.counts, "dominant_ratio_histogram": self.histogram}


def mapping_sources(config, sample):
    from .io_utils import sample_paths
    directory = sample_paths(config.output_dir, sample)["matrix"]
    return (directory / f"{sample}.pb_cb_mapping.tsv.gz",
            directory / f"{sample}.cb_observed_correction.tsv.gz")


def source_signature(paths, threshold):
    files = []
    for path in paths:
        resolved = path.resolve()
        stat = resolved.stat() if resolved.exists() else None
        files.append({"path": str(resolved), "size": stat.st_size if stat else None,
                      "mtime_ns": stat.st_mtime_ns if stat else None})
    return {"schema": 1, "dominance_threshold": float(threshold), "files": files}


@lru_cache(maxsize=32)
def legacy_dominance_qc(mapping_path, mapping_size, mapping_mtime, alias_path, alias_size, alias_mtime, threshold):
    """Stream sorted legacy mappings, merging aliases once as in the mapping stage.

    Reads remain additive when aliases share a UMI, so no molecule reprocessing
    is needed to reproduce read-based dominance and corrected relation counts.
    File metadata forms the cache key; memory is bounded by aliases + one PB.
    """
    aliases = {}
    if alias_size is not None:
        with gzip.open(alias_path, "rb") as handle:
            next(handle, None)
            for line in handle:
                raw, corrected, *_ = line.rstrip(b"\r\n").split(b"\t")
                aliases[raw] = corrected
    qc = DominanceQC(threshold)
    previous, supports = None, {}
    with gzip.open(mapping_path, "rb") as handle:
        next(handle, None)
        for line in handle:
            pb, cb, _, reads = line.rstrip(b"\r\n").split(b"\t")
            if pb != previous:
                if previous is not None:
                    if pb < previous:
                        raise ConfigError(f"PB-CB mapping is not sorted by PB: {mapping_path}")
                    values = supports.values()
                    qc.add(len(supports), sum(values), max(values))
                    supports.clear()
                previous = pb
            cb = aliases.get(cb, cb)  # One-step mapping matches the SQL LEFT JOIN.
            supports[cb] = supports.get(cb, 0) + int(reads)
    if supports:
        values = supports.values()
        qc.add(len(supports), sum(values), max(values))
    return qc.payload()


def validate_dominance_qc(qc, summary):
    c = qc
    total, retained = c["total_pbs"], c["retained_pbs"]
    if (c.get("schema") != 1 or c.get("basis") != "post_cb_correction"
            or total != c["single_cb_pbs"] + c["multi_cb_pbs"]
            or retained != c["single_cb_retained_pbs"] + c["resolved_multi_cb_pbs"]
            or total != retained + c["unassigned_pbs"]
            or c["multi_cb_pbs"] != c["resolved_multi_cb_pbs"] + c["unresolved_multi_cb_pbs"]
            or c["total_relations"] != retained + c["removed_secondary_relations"] + c["removed_unresolved_relations"]
            or sum(b["passed_pbs"] + b["failed_pbs"] for b in c["dominant_ratio_histogram"]) != c["multi_cb_pbs"]):
        raise ConfigError("Inconsistent PB dominance diagnostics")
    for key in ("total_pbs", "retained_pbs", "supporting_reads", "retained_supporting_reads", "dominance_threshold"):
        if key in summary and summary[key] != c[key]:
            raise ConfigError(f"PB dominance diagnostics do not match mapping summary: {key}")
    return qc
