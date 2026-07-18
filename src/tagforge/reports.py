from __future__ import annotations

import csv
import fcntl
import gzip
import html
import json
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

from . import __version__
from .barcode_correct import load_fb_annotation, load_whitelist_names
from .config import TagForgeConfig
from .io_utils import atomic_text, open_tsv, sample_dirs, step_complete, write_tsv
from .xlsx import write_xlsx


def _rows(path: Path):
    records = list(open_tsv(path))
    return ([list(records[0])] + [[row[k] for k in records[0]] for row in records]) if records else [[]]


def _matrix_summary(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t"); header = next(reader, [])
        b1 = 0; totals = Counter()
        for row in reader:
            b1 += 1
            for name, value in zip(header[1:], row[1:]): totals[name] += int(value)
    return b1, totals


def report_sample(config: TagForgeConfig, sample_name: str):
    barcode1_name = config.target_name("barcode1")
    dirs = sample_dirs(config.output_dir, sample_name)
    extraction_stats = dirs["extracted"] / f"{sample_name}.extraction_stats.tsv"
    stats = dirs["corrected"] / f"{sample_name}.barcode_correction_stats.tsv"
    molecules = dirs["detail"] / (f"{sample_name}.molecule_detail.rmMP.tsv.gz" if config.pi_seq_enabled else f"{sample_name}.molecule_detail.tsv.gz")
    metrics = dirs["downsample"] / f"{sample_name}.downsample_metrics.tsv"
    point = dirs["downsample"] / f"{sample_name}.optimal_saturation_point.tsv"
    matrix = dirs["matrix"] / f"{sample_name}.raw_count_matrix.tsv.gz"
    xlsx = dirs["report"] / f"{sample_name}.report.xlsx"
    report_html = dirs["report"] / f"{sample_name}.report.html"
    stat_records = list(open_tsv(stats)); metric_records = list(open_tsv(metrics)); point_records = list(open_tsv(point))
    b1_count, features = _matrix_summary(matrix)
    molecule_count = sum(1 for _ in open_tsv(molecules)); read_count = sum(int(r["reads_count"]) for r in open_tsv(molecules))
    combined = next((r for r in stat_records if r["scope"] == "combined"), {})
    summary = [
        ["Metric", "Value"], ["Sample", sample_name], ["Total reads", combined.get("total_reads", 0)],
        ["Valid reads", combined.get("valid_reads", 0)], ["Combined valid rate", combined.get("valid_rate", 0)],
        ["Deduplicated molecules", molecule_count], ["Reads supporting molecules", read_count],
        [f"{barcode1_name} rows", b1_count], ["Features", len(features)],
    ]
    feature_rows = [["Feature", "Molecule count"]] + [[k, v] for k, v in features.most_common()]
    sheets = [("Summary", summary), ("Extraction Stats", _rows(extraction_stats)),
              ("Correction Stats", _rows(stats)), ("UMI Stats", [["Metric","Value"],["Molecules",molecule_count],["Supporting reads",read_count]]),
              ("Downsample Metrics", _rows(metrics)), ("Optimal Saturation", _rows(point)), ("Top Features", feature_rows[:1001])]
    write_xlsx(xlsx, sheets)
    ratios = [float(r["downsample_ratio"]) for r in metric_records]
    saturation = [float(r["sequencing_saturation"]) for r in metric_records]
    duplication = [float(r["duplication_ratio"]) for r in metric_records]
    umi_types = [int(r["umi_types"]) for r in metric_records]
    singletons = [int(r["umi_detected_once"]) for r in metric_records]
    valid_rate = float(combined.get("valid_rate", 0) or 0)
    flags = []
    if valid_rate < 0.5: flags.append("Combined barcode/UMI valid rate is below 50%.")
    if read_count and molecule_count / read_count < 0.1: flags.append("High duplication: fewer than 10% of supporting reads are unique molecules.")
    extraction_records = list(open_tsv(extraction_stats))
    extraction_rows = ''.join(
        f'<tr><td>{html.escape(row["segment"])}</td><td>{html.escape(row["configured_mode"])}</td>'
        f'<td>{float(row["linker_success_rate"]):.1%}</td><td>{row["fixed_success"]}</td>'
        f'<td>{float(row["final_success_rate"]):.1%}</td></tr>'
        for row in extraction_records if row["linker_attempted"] != "0"
    )
    page = f'''<!doctype html><html><head><meta charset="utf-8"><title>TagForge · {html.escape(sample_name)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script><style>
body{{margin:0;background:#07111f;color:#e5edf8;font:15px system-ui}}main{{max-width:1180px;margin:auto;padding:36px}}h1{{font-size:36px;margin:0}}.sub{{color:#94a3b8}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin:28px 0}}.card,.panel{{background:#0f1d30;border:1px solid #24344d;border-radius:14px;padding:18px}}.value{{font-size:28px;font-weight:750;color:#5eead4}}.panel{{margin:16px 0}}table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;border-bottom:1px solid #24344d;text-align:left}}th{{color:#93c5fd}}.flag{{color:#fbbf24}}</style></head><body><main>
<p class="sub">TAGFORGE / SAMPLE REPORT</p><h1>{html.escape(sample_name)}</h1><div class="cards">
<div class="card"><div class="sub">Total reads</div><div class="value">{int(combined.get('total_reads',0) or 0):,}</div></div>
<div class="card"><div class="sub">Valid rate</div><div class="value">{valid_rate:.1%}</div></div>
<div class="card"><div class="sub">Molecules</div><div class="value">{molecule_count:,}</div></div>
<div class="card"><div class="sub">{html.escape(barcode1_name)}</div><div class="value">{b1_count:,}</div></div></div>
<section class="panel"><h2>Cutadapt linker extraction</h2><table><tr><th>Segment</th><th>Mode</th><th>Linker success</th><th>Fixed rescues</th><th>Final success</th></tr>{extraction_rows or '<tr><td colspan="5">No linker segments configured.</td></tr>'}</table></section>
<section class="panel"><h2>Sequencing saturation</h2><div id="sat" style="height:390px"></div></section>
<section class="panel"><h2>Top features</h2><table><tr><th>Feature</th><th>Molecules</th></tr>{''.join(f'<tr><td>{html.escape(k)}</td><td>{v:,}</td></tr>' for k,v in features.most_common(25))}</table></section>
<section class="panel"><h2>QC flags</h2>{''.join(f'<p class="flag">⚠ {html.escape(f)}</p>' for f in flags) or '<p>No automatic QC warnings.</p>'}</section>
<script>Plotly.newPlot('sat',[
{{x:{json.dumps(ratios)},y:{json.dumps(saturation)},name:'Sequencing Saturation',mode:'lines+markers',yaxis:'y3',line:{{color:'#5eead4',width:3}}}},
{{x:{json.dumps(ratios)},y:{json.dumps(duplication)},name:'Duplication Ratio',mode:'lines+markers',yaxis:'y3',line:{{color:'#f97316',width:3}}}},
{{x:{json.dumps(ratios)},y:{json.dumps(umi_types)},name:'UMI Types',mode:'lines+markers',yaxis:'y2',line:{{color:'#22c55e',width:3}}}},
{{x:{json.dumps(ratios)},y:{json.dumps(singletons)},name:'UMI detected once',mode:'lines+markers',yaxis:'y',line:{{color:'#a78bfa',width:3}}}}
],{{paper_bgcolor:'#0f1d30',plot_bgcolor:'#0f1d30',font:{{color:'#cbd5e1'}},legend:{{orientation:'h',x:0.5,xanchor:'center',y:-0.15}},xaxis:{{title:'Downsample ratio',gridcolor:'#24344d'}},yaxis:{{title:'UMI detected once',domain:[0.05,0.48],gridcolor:'#24344d'}},yaxis2:{{title:'UMI Types',overlaying:'y',side:'right',domain:[0.05,0.48],gridcolor:'#24344d'}},yaxis3:{{title:'Saturation / Duplication (%)',domain:[0.56,1],range:[0,100],gridcolor:'#24344d'}}}},{{responsive:true}})</script></main></body></html>'''
    with atomic_text(report_html) as handle: handle.write(page)
    return [xlsx, report_html]


def batch_report(config: TagForgeConfig, sample_names):
    out = config.workdir / "00_report"; out.mkdir(parents=True, exist_ok=True)
    barcode1_name = config.target_name("barcode1")
    meta = [["sample", "total_reads", "valid_reads", "valid_rate", "molecules", f"{barcode1_name}_count", "feature_count"]]
    bulk = {}
    all_features = set()
    for sample in sample_names:
        dirs = sample_dirs(config.output_dir, sample)
        combined = next((r for r in open_tsv(dirs["corrected"] / f"{sample}.barcode_correction_stats.tsv") if r["scope"] == "combined"), {})
        b1_count, features = _matrix_summary(dirs["matrix"] / f"{sample}.raw_count_matrix.tsv.gz")
        molecules = sum(features.values()); bulk[sample] = features; all_features.update(features)
        meta.append([sample, combined.get("total_reads",0), combined.get("valid_reads",0), combined.get("valid_rate",0), molecules,b1_count,len(features)])
    feature_names = sorted(all_features)
    counts = [["sample"] + feature_names] + [[s] + [bulk[s].get(f,0) for f in feature_names] for s in sample_names]
    xlsx = out / "TagForge_batch_report.xlsx"; page_path = out / "TagForge_batch_report.html"
    write_xlsx(xlsx, [("meta",meta),("counts",counts)])
    cards = ''.join(f'<div class="card"><b>{html.escape(str(r[0]))}</b><span>{int(r[4]):,} molecules</span><span>{float(r[3]):.1%} valid</span></div>' for r in meta[1:])
    page = f'''<!doctype html><html><head><meta charset="utf-8"><title>TagForge batch</title><style>body{{background:#f6f1e8;color:#18202a;font:16px system-ui;margin:0}}main{{max-width:1150px;margin:auto;padding:50px}}h1{{font:700 48px Georgia}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:16px}}.card{{background:white;border-radius:16px;padding:22px;box-shadow:0 8px 28px #33415518;display:grid;gap:10px}}span{{color:#64748b}}</style></head><body><main><p>TAGFORGE · BATCH SUMMARY</p><h1>{len(sample_names)} samples, forged cleanly.</h1><div class="grid">{cards}</div></main></body></html>'''
    with atomic_text(page_path) as h: h.write(page)
    return [xlsx,page_path]


@contextmanager
def _summary_lock(path: Path):
    """Serialize rebuilds from independently finishing sample jobs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _summary_inputs(config: TagForgeConfig, sample: str):
    root = config.output_dir / sample
    return root / "01_checkpoint" / "downsample.done", [
        root / "03_corrected" / f"{sample}.barcode_correction_stats.tsv",
        root / "05_detail" / f"{sample}.optimal_saturation_molecule_detail.tsv.gz",
        root / "06_downsample" / f"{sample}.optimal_saturation_point.tsv",
    ]


def _barcode2_annotation_gap(config, sample, stats, barcode2_segments):
    """Return whitelist reads, unannotated reads, and their barcode counts.

    Existing correction traces are scanned at most once per trace/annotation
    version.  The compact cache keeps subsequent summary refreshes cheap.
    """
    root = config.output_dir / sample
    corrected = root / "03_corrected"
    trace = corrected / f"{sample}.barcode_correction_trace.tsv.gz"
    qc_path = corrected / f"{sample}.barcode2_annotation_qc.tsv"
    detail_path = corrected / f"{sample}.barcode2_not_in_annotation.tsv.gz"
    annotation_path = getattr(config, "fb_info", None)
    dependencies = [path for path in (trace, annotation_path) if path and path.is_file()]
    cache_current = (
        qc_path.is_file() and detail_path.is_file() and dependencies
        and min(qc_path.stat().st_mtime_ns, detail_path.stat().st_mtime_ns)
        >= max(path.stat().st_mtime_ns for path in dependencies)
    )
    if cache_current:
        qc = next(open_tsv(qc_path))
        return int(qc["whitelist_reads"]), int(qc["not_in_annotation_reads"]), list(open_tsv(detail_path))

    by_scope = {row["scope"]: row for row in stats}
    # A historical run without a trace can still report the rate for its usual
    # single barcode2 segment, but cannot recover the omitted sequences.
    if not trace.is_file():
        whitelist_reads = sum(int(by_scope.get(s.name, {}).get("valid_reads", 0) or 0) for s in barcode2_segments)
        annotated = int(by_scope.get(f"final_{config.target_name('barcode2')}", {}).get("valid_reads", 0) or 0)
        return whitelist_reads, max(0, whitelist_reads - annotated), []

    annotation = load_fb_annotation(config)
    segment_names = {segment.name for segment in barcode2_segments}
    segment_order = {segment.name: index for index, segment in enumerate(barcode2_segments)}
    counts = Counter(); whitelist_reads = 0
    current_id = None; current = {}

    def consume(values):
        nonlocal whitelist_reads
        if len(values) != len(segment_names):
            return
        sequence = "".join(values[name] for name in sorted(values, key=segment_order.__getitem__))
        whitelist_reads += 1
        if sequence not in annotation:
            counts[sequence] += 1

    with gzip.open(trace, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        columns = {name: index for index, name in enumerate(header)}
        read_id_index = columns["read_id"]
        segment_index = columns["segment_name"]
        hit_index = columns["whitelist_hit"]
        sequence_index = columns["corrected_sequence"]
        for row in reader:
            read_id = row[read_id_index]
            if current_id is not None and read_id != current_id:
                consume(current); current = {}
            current_id = read_id
            if row[segment_index] in segment_names and row[hit_index] == "true":
                current[row[segment_index]] = row[sequence_index]
    if current_id is not None:
        consume(current)
    details = [
        {"barcode2_sequence": sequence, "reads_count": count}
        for sequence, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    write_tsv(qc_path, ["whitelist_reads", "not_in_annotation_reads"], [{
        "whitelist_reads": whitelist_reads, "not_in_annotation_reads": sum(counts.values()),
    }])
    write_tsv(detail_path, ["barcode2_sequence", "reads_count"], details, config.compression_level)
    return whitelist_reads, sum(counts.values()), details


def _completed_summary_samples(config: TagForgeConfig):
    """Return samples in configuration order once their summary inputs exist."""
    completed = []
    for item in config.samples:
        checkpoint, outputs = _summary_inputs(config, item.sample)
        if step_complete(checkpoint, outputs, False, __version__):
            completed.append(item.sample)
    return completed


def write_summary(config: TagForgeConfig):
    """Atomically rebuild the dynamic multi-sample workbook without a checkpoint."""
    output = config.output_dir / "00_summary.xlsx"
    with _summary_lock(config.output_dir / ".summary.lock"):
        samples = _completed_summary_samples(config)
        barcode_segments = [s for s in config.segments if s.target in {"barcode1", "barcode2"}]
        barcode2_segments = [s for s in config.segments if s.target == "barcode2"]
        meta = [["Sample", "Total reads", *[f"{s.name} Valid rate" for s in barcode_segments],
                 "Annotated Rate",
                 "optimal_downsample_ratio", "sequencing_saturation", "reads_sampled",
                 "umi_types", "umi_detected_once", "duplication_ratio"] + (
                    ["MP Ratio", "Multi-PI ratio", "Dominant ratio", "Multi-PI FB_UMIs", "Retained molecules"]
                    if getattr(config, "pi_seq_enabled", False) else [])]
        feature_counts, all_features = {}, set()
        missing_annotation_sheets = []
        feature_column = f"{config.target_name('barcode2')}_name"
        barcode2_names = {}
        for segment in barcode2_segments:
            if getattr(segment, "whitelist", None):
                barcode2_names.update(load_whitelist_names(segment.whitelist))

        def publish():
            features = sorted(all_features)
            count = [["Sample"] + features] + [
                [sample] + [feature_counts[sample].get(feature, 0) for feature in features]
                for sample in feature_counts
            ]
            write_xlsx(output, [("meta", meta), ("count", count)])
            write_xlsx(config.output_dir / "01_FB_not_in_anno.xlsx", missing_annotation_sheets)

        for sample in samples:
            root = config.output_dir / sample
            stats = list(open_tsv(root / "03_corrected" / f"{sample}.barcode_correction_stats.tsv"))
            by_scope = {row["scope"]: row for row in stats}
            whitelist_reads, missing_reads, missing_details = _barcode2_annotation_gap(
                config, sample, stats, barcode2_segments)
            point = next(open_tsv(root / "06_downsample" / f"{sample}.optimal_saturation_point.tsv"))
            combined = by_scope.get("combined", {})
            pi_qc = next(open_tsv(root / "05_detail" / f"{sample}.pi_seq_qc.tsv"), {}) if getattr(config, "pi_seq_enabled", False) else {}
            meta.append([sample, combined.get("total_reads", ""),
                         *[by_scope.get(s.name, {}).get("valid_rate", "") for s in barcode_segments],
                         (whitelist_reads - missing_reads) / whitelist_reads if whitelist_reads else 0,
                         point.get("optimal_downsample_ratio", ""), point.get("max_sequencing_saturation", ""),
                         point.get("reads_sampled", ""), point.get("umi_types", ""),
                         point.get("umi_detected_once", ""), point.get("duplication_ratio", "")] + (
                            [pi_qc.get("mp_ratio", ""), pi_qc.get("multi_pi_ratio", ""),
                             pi_qc.get("dominant_ratio", ""), pi_qc.get("multi_pi_fb_umis", ""),
                             pi_qc.get("retained_molecules", "")]
                            if getattr(config, "pi_seq_enabled", False) else []))
            missing_annotation_sheets.append((sample, [["barcode2_name", "reads_count"]] + [
                [barcode2_names.get(row["barcode2_sequence"], row["barcode2_sequence"]), row["reads_count"]]
                for row in missing_details[:50]
            ]))
            # The optimal count matrix is the aggregation of the requested
            # molecule detail, and is substantially smaller to read.
            optimal_matrix = root / "04_matrix" / f"{sample}.optimal_saturation_count_matrix.tsv.gz"
            counts = (
                _matrix_summary(optimal_matrix)[1] if optimal_matrix.is_file() else
                Counter(row[feature_column] for row in open_tsv(
                    root / "05_detail" / f"{sample}.optimal_saturation_molecule_detail.tsv.gz"))
            )
            feature_counts[sample] = counts
            all_features.update(counts)
            publish()
    return output, samples
