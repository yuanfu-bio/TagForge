from __future__ import annotations

import csv
import hashlib
import random
from pathlib import Path

from .config import TagForgeConfig
from .io_utils import atomic_text, open_tsv, sample_dirs, write_tsv
from .matrix import matrix_from_molecules


METRIC_FIELDS = ["sample", "downsample_ratio", "reads_sampled", "umi_types", "umi_detected_once",
                 "duplication_ratio", "sequencing_saturation", "repeat"]


def _seed(base: int, sample: str, ratio: float, repeat: int) -> int:
    digest = hashlib.sha256(f"{base}:{sample}:{ratio:.12g}:{repeat}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def _binomial(n: int, probability: float, rng: random.Random) -> int:
    if probability >= 1: return n
    return sum(rng.random() < probability for _ in range(n))


def calculate_metrics(supports):
    positive = [int(x) for x in supports if int(x) > 0]
    reads = sum(positive); molecules = len(positive); singles = sum(x == 1 for x in positive)
    duplication = ((reads - molecules) / reads * 100) if reads else 0.0
    saturation = ((1 - singles / molecules) * 100) if molecules else 0.0
    return reads, molecules, singles, duplication, saturation


def downsample_sample(config: TagForgeConfig, sample_name: str):
    dirs = sample_dirs(config.output_dir, sample_name)
    source = dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz"
    metrics_path = dirs["downsample"] / f"{sample_name}.downsample_metrics.tsv"
    point_path = dirs["downsample"] / f"{sample_name}.optimal_saturation_point.tsv"
    optimal_detail = dirs["detail"] / f"{sample_name}.optimal_saturation_molecule_detail.tsv.gz"
    optimal_matrix = dirs["matrix"] / f"{sample_name}.optimal_saturation_count_matrix.tsv.gz"
    barcode1_col = config.target_name("barcode1")
    barcode2_name_col = f"{config.target_name('barcode2')}_name"
    metric_rows = []
    for ratio in config.downsample_ratios:
        for repeat in range(1, config.downsample_repeats + 1):
            rng = random.Random(_seed(config.downsample_seed, sample_name, ratio, repeat))
            supports = (_binomial(int(row["reads_count"]), ratio, rng) for row in open_tsv(source))
            reads, molecules, singles, duplication, saturation = calculate_metrics(supports)
            metric_rows.append({"sample": sample_name, "downsample_ratio": ratio, "reads_sampled": reads,
                "umi_types": molecules, "umi_detected_once": singles, "duplication_ratio": f"{duplication:.6f}",
                "sequencing_saturation": f"{saturation:.6f}", "repeat": repeat})
    write_tsv(metrics_path, METRIC_FIELDS, metric_rows)
    optimal = max(metric_rows, key=lambda row: (float(row["sequencing_saturation"]), -float(row["downsample_ratio"]), -int(row["repeat"])))
    point_fields = [x for x in METRIC_FIELDS if x != "repeat"]
    point = dict(optimal); point["optimal_downsample_ratio"] = point.pop("downsample_ratio")
    point["max_sequencing_saturation"] = point.pop("sequencing_saturation")
    point_fields = ["sample", "optimal_downsample_ratio", "max_sequencing_saturation", "reads_sampled",
                    "umi_types", "umi_detected_once", "duplication_ratio"]
    write_tsv(point_path, point_fields, [point])
    ratio, repeat = float(optimal["downsample_ratio"]), int(optimal["repeat"])
    rng = random.Random(_seed(config.downsample_seed, sample_name, ratio, repeat))
    with atomic_text(optimal_detail, config.compression_level) as handle:
        fields = [barcode1_col, barcode2_name_col, "corrected_umi", "reads_count_at_optimal_downsample"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader()
        for row in open_tsv(source):
            count = _binomial(int(row["reads_count"]), ratio, rng)
            if count:
                writer.writerow({barcode1_col: row[barcode1_col], barcode2_name_col: row[barcode2_name_col],
                                 "corrected_umi": row["corrected_umi"], "reads_count_at_optimal_downsample": count})
    matrix_from_molecules(
        optimal_detail, optimal_matrix, config.compression_level,
        barcode_col=barcode1_col, feature_col=barcode2_name_col,
        row_header=barcode1_col,
    )
    return [metrics_path, point_path, optimal_detail, optimal_matrix], optimal
