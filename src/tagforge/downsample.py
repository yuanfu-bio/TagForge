from __future__ import annotations

from collections import Counter
import hashlib
import html
import json
import multiprocessing
import math
import os
import random
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .config import TagForgeConfig
from .columnar import ParquetFragmentWriter, open_records, resolve_table, record_batches
from .io_utils import atomic_text, sample_dirs, write_progress, write_tsv
from .logging_utils import sample_logger
from .selected import write_selection


METRIC_FIELDS = ["sample", "downsample_ratio", "reads_sampled", "umi_types", "umi_detected_once",
                 "sequencing_saturation", "non_singleton_molecule_fraction"]


# Forked workers inherit the molecule counts copy-on-write, avoiding repeated
# serialization of the complete vector for every ratio task.
_WORKER_COUNTS: list[int] | Counter[int] | None = None


def _seed(base: int, sample: str, ratio: float) -> int:
    # Keep the established single-draw seed stream stable across schema changes.
    digest = hashlib.sha256(f"{base}:{sample}:{ratio:.12g}:1".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def _binomial(n: int, probability: float, rng: random.Random) -> int:
    if probability >= 1: return n
    if probability <= 0: return 0
    if hasattr(rng, "binomialvariate"):
        return rng.binomialvariate(n, probability)
    return sum(rng.random() < probability for _ in range(n))


def calculate_metrics(supports):
    if isinstance(supports, Counter):
        reads = sum(int(value) * int(multiplicity) for value, multiplicity in supports.items())
        molecules = sum(int(multiplicity) for multiplicity in supports.values())
        singles = int(supports.get(1, 0))
        saturation = ((reads - molecules) / reads * 100) if reads else 0.0
        non_singleton_fraction = ((1 - singles / molecules) * 100) if molecules else 0.0
        return reads, molecules, singles, saturation, non_singleton_fraction
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
    saturation = ((reads - molecules) / reads * 100) if reads else 0.0
    non_singleton_fraction = ((1 - singles / molecules) * 100) if molecules else 0.0
    return reads, molecules, singles, saturation, non_singleton_fraction


CURVE_METRICS = ("reads_sampled", "umi_types", "umi_detected_once", "sequencing_saturation",
                 "non_singleton_molecule_fraction")
SAMPLING_NOTE = (
    "Each ratio is sampled once using a reproducible random seed. "
    "Reads are sampled from fixed full-depth corrected molecule groups; UMI correction and library filtering are not rerun. "
    "The selected ratio maximizes the observed non-singleton molecule fraction, with ties using the lower ratio. "
    "This diagnostic does not establish recommended sequencing depth. Analytical count expectations are stored separately."
)


def curve_rows(rows):
    """Return the single observed result at each ratio in ratio order."""
    result = []
    seen = set()
    for row in sorted(rows, key=lambda row: float(row["downsample_ratio"])):
        ratio = float(row["downsample_ratio"])
        if ratio in seen:
            raise ValueError("Expected one downsampling result per ratio")
        seen.add(ratio)
        result.append({"sample": row["sample"], "downsample_ratio": ratio,
                       **{metric: float(row[metric]) for metric in CURVE_METRICS}})
    return result


def choose_selection(rows):
    curve_rows(rows)
    selected = dict(max(rows, key=lambda row: (
        float(row["non_singleton_molecule_fraction"]), -float(row["downsample_ratio"]))))
    selected["selection_criterion"] = "max_non_singleton_molecule_fraction"
    return selected


def expected_counts(histogram, ratio):
    """Exact count expectations under thinning of fixed full-depth molecules.

    This does not report expectations of nonlinear fractions, re-cluster UMIs,
    or repeat library filtering at each depth.
    """
    if ratio == 0:
        return 0.0, 0.0, 0.0
    reads = sum(n * count for n, count in histogram.items()) * ratio
    if ratio == 1:
        return reads, float(sum(histogram.values())), float(histogram.get(1, 0))
    log_survival = math.log1p(-ratio)
    molecules = sum(-math.expm1(n * log_survival) * count for n, count in histogram.items())
    singletons = sum(n * ratio * math.exp((n - 1) * log_survival) * count for n, count in histogram.items())
    return reads, molecules, singletons


def _analysis_ratios(config: TagForgeConfig) -> list[float]:
    """Return the plotted/evaluated grid, including reference-style endpoints."""
    ratios = set(config.downsample_ratios)
    ratios.add(0.0)
    ratios.add(1.0)
    return sorted(ratios)


def _sample_counts(counts: list[int], ratio: float, sample_name: str, seed: int) -> list[int]:
    if ratio <= 0:
        return [0] * len(counts)
    if ratio >= 1:
        return list(counts)
    rng = random.Random(_seed(seed, sample_name, ratio))
    return [_binomial(count, ratio, rng) for count in counts]


def _sample_metrics(counts: list[int], ratio: float, sample_name: str, seed: int):
    """Sample and calculate metrics in one pass without retaining sampled counts."""
    if ratio <= 0:
        return 0, 0, 0, 0.0, 0.0
    if ratio >= 1:
        return calculate_metrics(counts)
    rng = random.Random(_seed(seed, sample_name, ratio))
    reads = molecules = singles = 0
    iterable = sorted(counts.items()) if isinstance(counts, Counter) else ((count, 1) for count in counts)
    for count, multiplicity in iterable:
        value_rng = random.Random(_seed(seed, f"{sample_name}:{count}", ratio)) if isinstance(counts, Counter) else rng
        for _ in range(multiplicity):
            sampled = _binomial(count, ratio, value_rng)
            if sampled <= 0:
                continue
            reads += sampled
            molecules += 1
            singles += 1 if sampled == 1 else 0
    saturation = ((reads - molecules) / reads * 100) if reads else 0.0
    non_singleton_fraction = ((1 - singles / molecules) * 100) if molecules else 0.0
    return reads, molecules, singles, saturation, non_singleton_fraction


def _parallel_sample_metrics(job: tuple[float, str, int]):
    if _WORKER_COUNTS is None:
        raise RuntimeError("downsample worker counts are not initialized")
    ratio, sample_name, seed = job
    return _sample_metrics(_WORKER_COUNTS, ratio, sample_name, seed)


def _write_saturation_html(path: Path, sample_name: str, rows: list[dict]):
    payload = {
        "sample": sample_name,
        "ratios": [float(row["downsample_ratio"]) for row in rows],
        "saturation": [float(row["sequencing_saturation"]) for row in rows],
        "non_singleton": [float(row["non_singleton_molecule_fraction"]) for row in rows],
        "umi_types": [float(row["umi_types"]) for row in rows],
        "singletons": [float(row["umi_detected_once"]) for row in rows],
    }
    full_depth = next(
        (row for row in rows if float(row["downsample_ratio"]) == 1.0),
        max(rows, key=lambda row: float(row["downsample_ratio"]), default=None),
    )
    title = json.dumps(f"<b>{html.escape(sample_name)}</b>")
    subtitle = json.dumps(
        f'<span style="color:#5963f5">Saturation {float(full_depth["sequencing_saturation"]):.2f}%</span>'
        f'<span style="color:#9ca3af">  |  </span>'
        f'<span style="color:#e44c39">Non-singleton Molecule Fraction {float(full_depth["non_singleton_molecule_fraction"]):.2f}%</span>'
    ) if full_depth is not None else '""'
    page = f'''<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(sample_name)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>html,body{{margin:0;background:#fff;font-family:Arial,sans-serif}}body{{min-height:100vh;display:grid;place-items:center}}#saturation{{width:100%;max-width:500px;height:540px}}@media print{{body{{min-height:0}}}}@page{{size:500px 540px;margin:0}}</style>
</head><body><div id="saturation"></div>
<script>
const data = {json.dumps(payload)};
Plotly.newPlot('saturation', [
  {{x:data.ratios,y:data.saturation,name:'Saturation',mode:'lines',yaxis:'y3',line:{{color:'#5963f5',width:3}}}},
  {{x:data.ratios,y:data.non_singleton,name:'Non-singleton Molecule Fraction',mode:'lines',yaxis:'y3',line:{{color:'#e44c39',width:3}}}},
  {{x:data.ratios,y:data.umi_types,name:'UMI Types',mode:'lines',yaxis:'y2',line:{{color:'#37c58d',width:3}}}},
  {{x:data.ratios,y:data.singletons,name:'UMI detected once',mode:'lines',yaxis:'y',line:{{color:'#9c59f5',width:3}}}}
], {{
  template:'plotly_white',
  height:540,
  font:{{family:'Arial'}},
  annotations:[
    {{text:{title},xref:'paper',yref:'paper',x:0,y:1.13,xanchor:'left',yanchor:'bottom',showarrow:false,font:{{family:'Arial',size:16,color:'#172033'}}}},
    {{text:{subtitle},xref:'paper',yref:'paper',x:0,y:1.075,xanchor:'left',yanchor:'bottom',showarrow:false,font:{{family:'Arial',size:11}}}}
  ],
  margin:{{l:58,r:58,t:75,b:90}},
  legend:{{orientation:'h',entrywidthmode:'fraction',entrywidth:0.5,font:{{size:10}},x:0.5,xanchor:'center',y:-0.12}},
  xaxis:{{title:'Downsample Ratio'}},
  yaxis:{{title:'UMI detected once',domain:[0.04,0.48],rangemode:'tozero',titlefont:{{color:'#9c59f5'}},tickfont:{{color:'#9c59f5'}},showline:true,linecolor:'#9c59f5',linewidth:1.5}},
  yaxis2:{{title:'UMI Types',overlaying:'y',side:'right',domain:[0.04,0.48],rangemode:'tozero',titlefont:{{color:'#37c58d'}},tickfont:{{color:'#37c58d'}},showline:true,linecolor:'#37c58d',linewidth:1.5}},
  yaxis3:{{title:'Saturation / non-singleton (%)',domain:[0.56,1],range:[0,100]}}
}}, {{responsive:true,displayModeBar:false,displaylogo:false}});
</script></body></html>'''
    with atomic_text(path) as handle:
        handle.write(page)


FONT_DIR = Path(os.environ.get("TAGFORGE_FONT_DIR", "/data/project/yanhuixu/yfu/01_source/fonts"))


def _write_saturation_pdf(path: Path, sample_name: str, rows: list[dict]):
    from matplotlib import font_manager, rc_context
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.ticker import FuncFormatter

    for font in ("Arial.ttf", "Arial Bold.ttf"):
        if (FONT_DIR / font).is_file():
            font_manager.fontManager.addfont(str(FONT_DIR / font))
    ratios = [float(row["downsample_ratio"]) for row in rows]
    saturation = [float(row["sequencing_saturation"]) for row in rows]
    non_singleton = [float(row["non_singleton_molecule_fraction"]) for row in rows]
    umi_types = [float(row["umi_types"]) for row in rows]
    singletons = [float(row["umi_detected_once"]) for row in rows]
    colors = {"sat": "#5963f5", "ns": "#e44c39", "types": "#37c58d", "once": "#9c59f5"}
    with rc_context({"font.family": ["Arial", "DejaVu Sans"], "pdf.fonttype": 42, "axes.unicode_minus": False}):
        fig = Figure(figsize=(5.2, 5.6))
        FigureCanvasAgg(fig)
        kilo = FuncFormatter(lambda v, _: f"{v / 1000:g}k" if v >= 1000 else f"{v:g}")
        top, bottom = fig.subplots(2, 1, gridspec_kw={"height_ratios": [1, 1], "hspace": 0.28})
        top.plot(ratios, saturation, color=colors["sat"], lw=2, label="Saturation")
        top.plot(ratios, non_singleton, color=colors["ns"], lw=2, label="Non-singleton Molecule Fraction")
        top.set_ylim(0, 100)
        top.set_ylabel("Saturation / non-singleton (%)")
        bottom.plot(ratios, singletons, color=colors["once"], lw=2, label="UMI detected once")
        bottom.set_ylabel("UMI detected once", color=colors["once"])
        bottom.tick_params(axis="y", colors=colors["once"])
        bottom.set_ylim(bottom=0)
        twin = bottom.twinx()
        twin.plot(ratios, umi_types, color=colors["types"], lw=2, label="UMI Types")
        twin.set_ylabel("UMI Types", color=colors["types"])
        twin.tick_params(axis="y", colors=colors["types"])
        twin.set_ylim(bottom=0)
        bottom.set_xlabel("Downsample Ratio")
        bottom.yaxis.set_major_formatter(kilo)
        twin.yaxis.set_major_formatter(kilo)
        top.spines["right"].set_visible(False)
        bottom.spines["left"].set_color(colors["once"])
        bottom.spines["right"].set_visible(False)
        twin.spines["right"].set_color(colors["types"])
        twin.spines["left"].set_visible(False)
        for axis in (top, bottom, twin):
            axis.spines["top"].set_visible(False)
            if axis is not twin:
                axis.grid(color="#e5e7eb", lw=0.6)
                axis.set_axisbelow(True)
            if axis is not top:
                axis.set_xlim(0, max(ratios, default=1))
        top.set_xlim(0, max(ratios, default=1))
        handles = top.get_legend_handles_labels()[0] + bottom.get_legend_handles_labels()[0] + twin.get_legend_handles_labels()[0]
        labels = top.get_legend_handles_labels()[1] + bottom.get_legend_handles_labels()[1] + twin.get_legend_handles_labels()[1]
        order = [0, 3, 1, 2]  # matplotlib fills columns first: rows read Sat | NS / Types | Once
        fig.legend([handles[i] for i in order], [labels[i] for i in order], loc="lower center", ncol=2,
                   frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.005))
        full_depth = next((r for r in rows if float(r["downsample_ratio"]) == 1.0),
                          max(rows, key=lambda r: float(r["downsample_ratio"]), default=None))
        fig.text(0.02, 0.975, sample_name, ha="left", va="top", fontsize=13, fontweight="bold")
        if full_depth is not None:
            renderer = fig.canvas.get_renderer()
            x = 0.02
            for text, color in ((f"Saturation {float(full_depth['sequencing_saturation']):.2f}%", colors["sat"]),
                                ("|", "#9ca3af"),
                                (f"Non-singleton Molecule Fraction {float(full_depth['non_singleton_molecule_fraction']):.2f}%", colors["ns"])):
                label = fig.text(x, 0.935, text, ha="left", va="top", fontsize=9, color=color)
                x = label.get_window_extent(renderer).transformed(fig.transFigure.inverted()).x1 + 0.015
        fig.subplots_adjust(left=0.15, right=0.85, top=0.88, bottom=0.17)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(str(path) + ".tmp")
        try:
            fig.savefig(tmp, format="pdf")
            os.replace(tmp, path)
        finally:
            tmp.unlink(missing_ok=True)


def downsample_sample(config: TagForgeConfig, sample_name: str):
    dirs = sample_dirs(config.output_dir, sample_name)
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    source = dirs["detail"] / (f"{sample_name}.molecule_detail.rmMP.parquet" if config.pi_seq_enabled else f"{sample_name}.molecule_detail.parquet")
    source = resolve_table(source, source.with_name(source.name.replace(".parquet", ".tsv.gz")))
    metrics_path = dirs["downsample"] / f"{sample_name}.downsample_metrics.tsv"
    html_path = dirs["downsample"] / f"{sample_name}.downsample.html"
    pdf_path = dirs["downsample"] / f"{sample_name}.downsample.pdf"
    progress_path = dirs["logs"] / f"{sample_name}.progress.json"
    optimal_parquet = dirs["detail"] / f"{sample_name}.optimal_saturation_molecule_detail.parquet"
    barcode1_col = config.target_name("barcode1")
    barcode2_name_col = f"{config.target_name('barcode2')}_name"
    started = time.monotonic()

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
            "reads_sampled": values.get("reads_sampled", ""),
            "umi_types": values.get("umi_types", ""),
            "umi_detected_once": values.get("umi_detected_once", ""),
            "sequencing_saturation": values.get("sequencing_saturation", ""),
            "non_singleton_molecule_fraction": values.get("non_singleton_molecule_fraction", ""),
            "input_percent": values.get("input_percent", ""),
            "elapsed_seconds": f"{elapsed:.2f}",
            "ratios_per_second": f"{rate:.3f}",
            "eta_seconds": f"{eta:.2f}" if eta else "",
            "estimated_finish": finish,
        }
        write_progress(progress_path, sample_name, "downsample", row)
        logger.info(
            "downsample_progress\tstatus=%s\tmolecules_loaded=%s\ttotal_reads=%s\t"
            "ratios_completed=%s\ttotal_ratios=%s\tdownsample_ratio=%s\t"
            "reads_sampled=%s\tumi_types=%s\tumi_detected_once=%s\t"
            "sequencing_saturation=%s\t"
            "non_singleton_molecule_fraction=%s\tinput_percent=%s\t"
            "elapsed_seconds=%s\tratios_per_second=%s\teta_seconds=%s\testimated_finish=%s",
            row["status"], row["molecules_loaded"], row["total_reads"],
            row["ratios_completed"], row["total_ratios"], row["downsample_ratio"],
            row["reads_sampled"], row["umi_types"],
            row["umi_detected_once"],
            row["sequencing_saturation"], row["non_singleton_molecule_fraction"], row["input_percent"],
            row["elapsed_seconds"], row["ratios_per_second"],
            row["eta_seconds"] or "NA", row["estimated_finish"] or "NA",
        )

    counts = Counter()
    molecules_loaded = 0
    total_reads = 0
    for batch, fraction in record_batches(source, config.chunk_size):
        for row in batch:
            count = int(row["reads_count"])
            counts[count] += 1
            molecules_loaded += 1
            total_reads += count
        record_progress(
            "loading", molecules_loaded=molecules_loaded, total_reads=total_reads,
            input_percent=f"{fraction * 100:.2f}",
        )

    ratios = _analysis_ratios(config)
    total_ratio_jobs = len(ratios)
    logger.info(
        "downsample_start\tmolecules=%s\ttotal_reads=%s\tratios=%s\t"
        "sampling_backend=%s",
        molecules_loaded, total_reads, len(ratios),
        "random.binomialvariate" if hasattr(random.Random(), "binomialvariate") else "python-loop",
    )

    jobs = [(ratio, sample_name, config.downsample_seed) for ratio in ratios]
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
        results = executor.map(_parallel_sample_metrics, jobs) if executor else (_sample_metrics(counts, ratio, sample, seed) for ratio, sample, seed in jobs)
        for (ratio, _, _), (reads, molecules, singles, saturation, non_singleton_fraction) in zip(jobs, results):
            row = {"sample": sample_name, "downsample_ratio": ratio, "reads_sampled": reads, "umi_types": molecules, "umi_detected_once": singles, "sequencing_saturation": f"{saturation:.6f}", "non_singleton_molecule_fraction": f"{non_singleton_fraction:.6f}"}
            metric_rows.append(row)
            ratios_completed += 1
            record_progress("sampling", molecules_loaded=molecules_loaded, total_reads=total_reads, ratios_completed=ratios_completed, total_ratios=total_ratio_jobs, downsample_ratio=ratio, reads_sampled=reads, umi_types=molecules, umi_detected_once=singles, sequencing_saturation=f"{saturation:.6f}", non_singleton_molecule_fraction=f"{non_singleton_fraction:.6f}")
    finally:
        if executor:
            executor.shutdown()
        _WORKER_COUNTS = None
    write_tsv(metrics_path, METRIC_FIELDS, metric_rows)
    curves = []
    for row in curve_rows(metric_rows):
        ratio = row["downsample_ratio"]
        curves.append({**{key if key in {"sample", "downsample_ratio"} else "observed_" + key: value
                         for key, value in row.items()},
                       **dict(zip(("expected_reads", "expected_molecules", "expected_singletons"), expected_counts(counts, ratio)))})
    write_tsv(dirs["downsample"] / f"{sample_name}.downsample_curve.tsv", list(curves[0]), curves)
    _write_saturation_html(html_path, sample_name, curve_rows(metric_rows))
    _write_saturation_pdf(pdf_path, sample_name, curve_rows(metric_rows))
    optimal = choose_selection(metric_rows)
    ratio = float(optimal["downsample_ratio"])
    owned = []
    selected_molecules = source
    matrix_paths = None
    if ratio < 1.0:
        rngs = {}
        fields = [barcode1_col, barcode2_name_col, "corrected_umi", "reads_count_at_optimal_downsample"]
        writer = ParquetFragmentWriter(optimal_parquet, fields, config.chunk_size, compression_level=getattr(config, "compression_level", 3), compression=getattr(config, "parquet_compression", "snappy"))
        selected_reads = selected_count = selected_singletons = 0
        for row in open_records(source):
            original_count = int(row["reads_count"])
            rng = rngs.setdefault(original_count, random.Random(_seed(config.downsample_seed, f"{sample_name}:{original_count}", ratio)))
            count = _binomial(original_count, ratio, rng)
            if count:
                selected_reads += count
                selected_count += 1
                selected_singletons += int(count == 1)
                writer.writerow({barcode1_col: row[barcode1_col], barcode2_name_col: row[barcode2_name_col],
                                 "corrected_umi": row["corrected_umi"], "reads_count_at_optimal_downsample": count})
        writer.close()
        if (selected_reads, selected_count, selected_singletons) != (
            int(optimal["reads_sampled"]), int(optimal["umi_types"]), int(optimal["umi_detected_once"])):
            raise RuntimeError("Selected molecule details do not match downsample metrics")
        selected_molecules = optimal_parquet
        owned.append(optimal_parquet)
    if not config.pb_cb_enabled:
        from .sparse_matrix import build_sparse_matrix
        if ratio >= 1.0:
            from .sparse_matrix import pair_path, summary_path
            matrix_paths = (pair_path(dirs["matrix"], sample_name), summary_path(dirs["matrix"], sample_name))
            if not all(path.exists() for path in matrix_paths):
                matrix_paths = build_sparse_matrix(source, dirs["matrix"], sample_name, barcode1_col, barcode2_name_col, prefix="optimal", threads=config.threads, memory_mb=config.matrix_memory_mb, compression_level=config.compression_level, compression=config.parquet_compression)
                owned.extend(matrix_paths)
        else:
            matrix_paths = build_sparse_matrix(selected_molecules, dirs["matrix"], sample_name, barcode1_col, barcode2_name_col, prefix="optimal", threads=config.threads, memory_mb=config.matrix_memory_mb, compression_level=config.compression_level, compression=config.parquet_compression)
            owned.extend(matrix_paths)
    write_selection(config, sample_name, optimal, selected_molecules, matrix_paths, owned)
    record_progress(
        "completed", molecules_loaded=molecules_loaded, total_reads=total_reads,
        ratios_completed=ratios_completed, total_ratios=total_ratio_jobs,
        downsample_ratio=ratio,
        reads_sampled=optimal["reads_sampled"], umi_types=optimal["umi_types"],
        umi_detected_once=optimal["umi_detected_once"],
        sequencing_saturation=optimal["sequencing_saturation"],
        non_singleton_molecule_fraction=optimal["non_singleton_molecule_fraction"],
    )
    logger.info(
        "downsample_summary\tmolecules=%s\ttotal_reads=%s\tratios=%s\t"
        "ratio_jobs=%s\toptimal_ratio=%s\tmax_non_singleton_molecule_fraction=%s\t"
        "reads_sampled=%s\telapsed_time=%.3f",
         molecules_loaded, total_reads, len(ratios), total_ratio_jobs,
        optimal["downsample_ratio"], optimal["non_singleton_molecule_fraction"],
        optimal["reads_sampled"], time.monotonic() - started,
    )
    return [metrics_path, dirs["downsample"] / f"{sample_name}.downsample_curve.tsv",
            dirs["downsample"] / f"{sample_name}.selection.json", html_path, pdf_path], optimal
