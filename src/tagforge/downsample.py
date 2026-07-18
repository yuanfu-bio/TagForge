from __future__ import annotations

import csv
import hashlib
import html
import json
import multiprocessing
import random
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .config import TagForgeConfig
from .io_utils import atomic_text, sample_dirs, tsv_batches, write_tsv
from .logging_utils import sample_logger
from .matrix import matrix_from_molecules


METRIC_FIELDS = ["sample", "downsample_ratio", "reads_sampled", "umi_types", "umi_detected_once",
                 "duplication_ratio", "sequencing_saturation", "repeat"]
PROGRESS_FIELDS = [
    "status", "molecules_loaded", "total_reads", "ratios_completed", "total_ratios",
    "downsample_ratio", "repeat", "reads_sampled", "umi_types", "umi_detected_once",
    "duplication_ratio", "sequencing_saturation", "input_percent",
    "elapsed_seconds", "ratios_per_second", "eta_seconds", "estimated_finish",
]


# Forked workers inherit the molecule counts copy-on-write, avoiding repeated
# serialization of the complete vector for every ratio/repeat task.
_WORKER_COUNTS: list[int] | None = None


def _seed(base: int, sample: str, ratio: float, repeat: int) -> int:
    digest = hashlib.sha256(f"{base}:{sample}:{ratio:.12g}:{repeat}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def _binomial(n: int, probability: float, rng: random.Random) -> int:
    if probability >= 1: return n
    if probability <= 0: return 0
    if hasattr(rng, "binomialvariate"):
        return rng.binomialvariate(n, probability)
    return sum(rng.random() < probability for _ in range(n))


def calculate_metrics(supports):
    reads = 0
    molecules = 0
    singles = 0
    for value in supports:
        count = int(value)
        if count <= 0:
            continue
        reads += count
        molecules += 1
        if count == 1:
            singles += 1
    duplication = ((reads - molecules) / reads * 100) if reads else 0.0
    saturation = ((1 - singles / molecules) * 100) if molecules else 0.0
    return reads, molecules, singles, duplication, saturation


def _analysis_ratios(config: TagForgeConfig) -> list[float]:
    """Return the plotted/evaluated grid, including reference-style endpoints."""
    ratios = set(config.downsample_ratios)
    ratios.add(0.0)
    ratios.add(1.0)
    return sorted(ratios)


def _sample_counts(counts: list[int], ratio: float, sample_name: str, repeat: int, seed: int) -> list[int]:
    if ratio <= 0:
        return [0] * len(counts)
    if ratio >= 1:
        return list(counts)
    rng = random.Random(_seed(seed, sample_name, ratio, repeat))
    return [_binomial(count, ratio, rng) for count in counts]


def _sample_metrics(counts: list[int], ratio: float, sample_name: str, repeat: int, seed: int):
    """Sample and calculate metrics in one pass without retaining sampled counts."""
    if ratio <= 0:
        return 0, 0, 0, 0.0, 0.0
    if ratio >= 1:
        return calculate_metrics(counts)
    rng = random.Random(_seed(seed, sample_name, ratio, repeat))
    reads = molecules = singles = 0
    for count in counts:
        sampled = _binomial(count, ratio, rng)
        if sampled <= 0:
            continue
        reads += sampled
        molecules += 1
        singles += 1 if sampled == 1 else 0
    duplication = ((reads - molecules) / reads * 100) if reads else 0.0
    saturation = ((1 - singles / molecules) * 100) if molecules else 0.0
    return reads, molecules, singles, duplication, saturation


def _parallel_sample_metrics(job: tuple[float, int, str, int]):
    if _WORKER_COUNTS is None:
        raise RuntimeError("downsample worker counts are not initialized")
    ratio, repeat, sample_name, seed = job
    return _sample_metrics(_WORKER_COUNTS, ratio, sample_name, repeat, seed)


def _write_saturation_html(path: Path, sample_name: str, rows: list[dict]):
    ratios = [float(row["downsample_ratio"]) for row in rows]
    saturation = [float(row["sequencing_saturation"]) for row in rows]
    duplication = [float(row["duplication_ratio"]) for row in rows]
    umi_types = [int(row["umi_types"]) for row in rows]
    singletons = [int(row["umi_detected_once"]) for row in rows]
    max_saturation = max(saturation) if saturation else 0.0
    payload = {
        "ratios": ratios,
        "saturation": saturation,
        "duplication": duplication,
        "umi_types": umi_types,
        "singletons": singletons,
    }
    page = f'''<!doctype html><html><head><meta charset="utf-8">
<title>TagForge saturation · {html.escape(sample_name)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>body{{margin:0;background:#f8fafc;color:#0f172a;font:15px system-ui}}main{{max-width:980px;margin:auto;padding:32px}}.panel{{background:white;border:1px solid #e2e8f0;border-radius:16px;padding:20px;box-shadow:0 10px 30px #33415514}}.sub{{color:#64748b}}</style>
</head><body><main><p class="sub">TAGFORGE / DOWNSAMPLE SATURATION</p>
<h1>{html.escape(sample_name)} · max saturation {max_saturation:.2f}%</h1>
<section class="panel"><div id="saturation" style="height:540px"></div></section>
<script>
const data = {json.dumps(payload)};
Plotly.newPlot('saturation', [
  {{x:data.ratios,y:data.saturation,name:'Sequencing Saturation',mode:'lines+markers',yaxis:'y3',line:{{color:'#5963f5',width:3}}}},
  {{x:data.ratios,y:data.duplication,name:'Duplication Ratio',mode:'lines+markers',yaxis:'y3',line:{{color:'#e44c39',width:3}}}},
  {{x:data.ratios,y:data.umi_types,name:'UMI Types',mode:'lines+markers',yaxis:'y2',line:{{color:'#37c58d',width:3}}}},
  {{x:data.ratios,y:data.singletons,name:'UMI detected once',mode:'lines+markers',yaxis:'y',line:{{color:'#9c59f5',width:3}}}}
], {{
  template:'plotly_white',
  height:540,
  legend:{{orientation:'h',x:0.5,xanchor:'center',y:-0.12}},
  xaxis:{{title:'Downsample Ratio'}},
  yaxis:{{title:'UMI detected once',domain:[0.04,0.48],rangemode:'tozero',titlefont:{{color:'#9c59f5'}},tickfont:{{color:'#9c59f5'}}}},
  yaxis2:{{title:'UMI Types',overlaying:'y',side:'right',domain:[0.04,0.48],rangemode:'tozero',titlefont:{{color:'#37c58d'}},tickfont:{{color:'#37c58d'}}}},
  yaxis3:{{title:'Saturation / Duplication (%)',domain:[0.56,1],range:[0,100]}}
}}, {{responsive:true}});
</script></main></body></html>'''
    with atomic_text(path) as handle:
        handle.write(page)


def downsample_sample(config: TagForgeConfig, sample_name: str):
    dirs = sample_dirs(config.output_dir, sample_name)
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    source = dirs["detail"] / (f"{sample_name}.molecule_detail.rmMP.tsv.gz" if config.pi_seq_enabled else f"{sample_name}.molecule_detail.tsv.gz")
    metrics_path = dirs["downsample"] / f"{sample_name}.downsample_metrics.tsv"
    point_path = dirs["downsample"] / f"{sample_name}.optimal_saturation_point.tsv"
    html_path = dirs["downsample"] / f"{sample_name}.downsample.html"
    progress_path = dirs["logs"] / f"{sample_name}.downsample_progress.tsv"
    optimal_detail = dirs["detail"] / f"{sample_name}.optimal_saturation_molecule_detail.tsv.gz"
    optimal_matrix = dirs["matrix"] / f"{sample_name}.optimal_saturation_count_matrix.tsv.gz"
    barcode1_col = config.target_name("barcode1")
    barcode2_name_col = f"{config.target_name('barcode2')}_name"
    started = time.monotonic()
    progress_rows = []

    def record_progress(status: str, **values):
        now = time.monotonic()
        elapsed = now - started
        ratios_completed = int(values.get("ratios_completed", 0) or 0)
        total_ratios = int(values.get("total_ratios", 0) or 0)
        rate = ratios_completed / elapsed if elapsed and ratios_completed else 0.0
        remaining = total_ratios - ratios_completed
        eta = remaining / rate if rate and remaining > 0 else 0.0
        finish = (
            time.strftime(
                "%Y-%m-%dT%H:%M:%S%z",
                time.localtime(time.time() + eta),
            )
            if eta else ""
        )
        row = {
            "status": status,
            "molecules_loaded": values.get("molecules_loaded", ""),
            "total_reads": values.get("total_reads", ""),
            "ratios_completed": ratios_completed,
            "total_ratios": total_ratios,
            "downsample_ratio": values.get("downsample_ratio", ""),
            "repeat": values.get("repeat", ""),
            "reads_sampled": values.get("reads_sampled", ""),
            "umi_types": values.get("umi_types", ""),
            "umi_detected_once": values.get("umi_detected_once", ""),
            "duplication_ratio": values.get("duplication_ratio", ""),
            "sequencing_saturation": values.get("sequencing_saturation", ""),
            "input_percent": values.get("input_percent", ""),
            "elapsed_seconds": f"{elapsed:.2f}",
            "ratios_per_second": f"{rate:.3f}",
            "eta_seconds": f"{eta:.2f}" if eta else "",
            "estimated_finish": finish,
        }
        progress_rows.append(row)
        write_tsv(progress_path, PROGRESS_FIELDS, progress_rows)
        logger.info(
            "downsample_progress\tstatus=%s\tmolecules_loaded=%s\ttotal_reads=%s\t"
            "ratios_completed=%s\ttotal_ratios=%s\tdownsample_ratio=%s\trepeat=%s\t"
            "reads_sampled=%s\tumi_types=%s\tumi_detected_once=%s\t"
            "duplication_ratio=%s\tsequencing_saturation=%s\tinput_percent=%s\t"
            "elapsed_seconds=%s\tratios_per_second=%s\teta_seconds=%s\testimated_finish=%s",
            row["status"], row["molecules_loaded"], row["total_reads"],
            row["ratios_completed"], row["total_ratios"], row["downsample_ratio"],
            row["repeat"], row["reads_sampled"], row["umi_types"],
            row["umi_detected_once"], row["duplication_ratio"],
            row["sequencing_saturation"], row["input_percent"],
            row["elapsed_seconds"], row["ratios_per_second"],
            row["eta_seconds"] or "NA", row["estimated_finish"] or "NA",
        )

    records = []
    counts = []
    total_reads = 0
    for batch, fraction in tsv_batches(source, config.chunk_size):
        for row in batch:
            count = int(row["reads_count"])
            records.append((row[barcode1_col], row[barcode2_name_col], row["corrected_umi"], count))
            counts.append(count)
            total_reads += count
        record_progress(
            "loading", molecules_loaded=len(records), total_reads=total_reads,
            input_percent=f"{fraction * 100:.2f}",
        )

    ratios = _analysis_ratios(config)
    total_ratio_jobs = 0
    for ratio in ratios:
        total_ratio_jobs += 1 if ratio in {0.0, 1.0} else config.downsample_repeats
    logger.info(
        "downsample_start\tmolecules=%s\ttotal_reads=%s\tratios=%s\trepeats=%s\t"
        "sampling_backend=%s",
        len(records), total_reads, len(ratios), config.downsample_repeats,
        "random.binomialvariate" if hasattr(random.Random(), "binomialvariate") else "python-loop",
    )

    jobs = [(ratio, repeat, sample_name, config.downsample_seed) for ratio in ratios
            for repeat in ([0] if ratio in {0.0, 1.0} else range(1, config.downsample_repeats + 1))]
    requested_workers = config.downsample_workers or config.threads
    workers = min(requested_workers, len(jobs))
    executor = None
    if workers > 1 and "fork" in multiprocessing.get_all_start_methods():
        global _WORKER_COUNTS
        _WORKER_COUNTS = counts
        executor = ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("fork"))
        backend = "process-fork"
    else:
        workers = 1
        backend = "serial"
    logger.info("downsample_parallel_start\tbackend=%s\trequested_workers=%s\tworkers=%s\tratio_jobs=%s", backend, requested_workers, workers, len(jobs))
    metric_rows = []
    ratios_completed = 0
    try:
        results = executor.map(_parallel_sample_metrics, jobs) if executor else (_sample_metrics(counts, ratio, sample, repeat, seed) for ratio, repeat, sample, seed in jobs)
        for (ratio, repeat, _, _), (reads, molecules, singles, duplication, saturation) in zip(jobs, results):
            row = {"sample": sample_name, "downsample_ratio": ratio, "reads_sampled": reads, "umi_types": molecules, "umi_detected_once": singles, "duplication_ratio": f"{duplication:.6f}", "sequencing_saturation": f"{saturation:.6f}", "repeat": repeat}
            metric_rows.append(row)
            ratios_completed += 1
            record_progress("sampling", molecules_loaded=len(records), total_reads=total_reads, ratios_completed=ratios_completed, total_ratios=total_ratio_jobs, downsample_ratio=ratio, repeat=repeat, reads_sampled=reads, umi_types=molecules, umi_detected_once=singles, duplication_ratio=f"{duplication:.6f}", sequencing_saturation=f"{saturation:.6f}")
    finally:
        if executor:
            executor.shutdown()
        _WORKER_COUNTS = None
    write_tsv(metrics_path, METRIC_FIELDS, metric_rows)
    _write_saturation_html(html_path, sample_name, metric_rows)
    optimal = max(metric_rows, key=lambda row: (float(row["sequencing_saturation"]), -float(row["downsample_ratio"]), -int(row["repeat"])))
    point_fields = [x for x in METRIC_FIELDS if x != "repeat"]
    point = dict(optimal); point["optimal_downsample_ratio"] = point.pop("downsample_ratio")
    point["max_sequencing_saturation"] = point.pop("sequencing_saturation")
    point_fields = ["sample", "optimal_downsample_ratio", "max_sequencing_saturation", "reads_sampled",
                    "umi_types", "umi_detected_once", "duplication_ratio"]
    write_tsv(point_path, point_fields, [point])
    ratio, repeat = float(optimal["downsample_ratio"]), int(optimal["repeat"])
    optimal_counts = _sample_counts(counts, ratio, sample_name, repeat, config.downsample_seed)
    with atomic_text(optimal_detail, config.compression_level) as handle:
        fields = [barcode1_col, barcode2_name_col, "corrected_umi", "reads_count_at_optimal_downsample"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader()
        for record, count in zip(records, optimal_counts):
            if count:
                writer.writerow({barcode1_col: record[0], barcode2_name_col: record[1],
                                 "corrected_umi": record[2], "reads_count_at_optimal_downsample": count})
    matrix_from_molecules(
        optimal_detail, optimal_matrix, config.compression_level,
        barcode_col=barcode1_col, feature_col=barcode2_name_col,
        row_header=barcode1_col,
    )
    record_progress(
        "completed", molecules_loaded=len(records), total_reads=total_reads,
        ratios_completed=ratios_completed, total_ratios=total_ratio_jobs,
        downsample_ratio=ratio, repeat=repeat,
        reads_sampled=optimal["reads_sampled"], umi_types=optimal["umi_types"],
        umi_detected_once=optimal["umi_detected_once"],
        duplication_ratio=optimal["duplication_ratio"],
        sequencing_saturation=optimal["sequencing_saturation"],
    )
    logger.info(
        "downsample_summary\tmolecules=%s\ttotal_reads=%s\tratios=%s\t"
        "ratio_jobs=%s\toptimal_ratio=%s\tmax_sequencing_saturation=%s\t"
        "reads_sampled=%s\telapsed_time=%.3f",
        len(records), total_reads, len(ratios), total_ratio_jobs,
        point["optimal_downsample_ratio"], point["max_sequencing_saturation"],
        point["reads_sampled"], time.monotonic() - started,
    )
    return [metrics_path, point_path, html_path, optimal_detail, optimal_matrix], optimal
