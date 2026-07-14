from __future__ import annotations

import csv
import gzip
import html
import json
from collections import Counter
from pathlib import Path

from .config import TagForgeConfig
from .io_utils import atomic_text, open_tsv, sample_dirs
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
    molecules = dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz"
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
