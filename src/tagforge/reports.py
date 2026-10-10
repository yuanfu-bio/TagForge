from __future__ import annotations

import csv
import fcntl
import gzip
import html
import json
import math
import os
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from urllib.parse import quote

from . import __version__
from .barcode_correct import load_fb_annotation, load_whitelist_names, correction_qc_path, load_correction_stats
from .columnar import open_records, is_parquet_table, resolve_table
from .config import TagForgeConfig, ConfigError
from .io_utils import atomic_text, batch_report_dir, open_tsv, sample_dirs, sample_paths, write_tsv
from .xlsx import write_xlsx
from .checkpoint import stage_current
from .selected import selected_summary_path, selection_point
from .downsample import SAMPLING_NOTE


def primary_result(config, sample):
    paths = sample_paths(config.output_dir, sample)
    if getattr(config, "pb_cb_enabled", False):
        return "pb-cb", paths["matrix"] / f"{sample}.pb_cb_map.tsv.gz", "Final PB-to-CB mapping after dominance filtering; full-depth read support."
    if getattr(config, "pi_seq_enabled", False):
        return "pi-seq", paths["detail"] / f"{sample}.molecule_detail.rmMP.parquet", "Full-depth PI-filtered molecules for per-PB analysis."
    return "standard", paths["matrix"] / f"{sample}.raw_molecule_pairs.parquet", "Full-depth UMI-deduplicated barcode-by-feature counts in sparse long format."


def _result_rows(config, samples):
    rows = [["Sample", "Library", "Primary result (relative to output_dir)", "Description"]]
    for sample in samples:
        library, path, description = primary_result(config, sample)
        rows.append([sample, library, str(path.relative_to(config.output_dir)), description])
    return rows


def _result_html(config, samples, directory):
    rows = []
    for sample in samples:
        library, path, description = primary_result(config, sample)
        href = quote(os.path.relpath(path, directory))
        rows.append(f'<tr><td>{html.escape(sample)}</td><td>{html.escape(library)}</td>'
                    f'<td><a href="{html.escape(href, quote=True)}">{html.escape(path.name)}</a></td>'
                    f'<td>{html.escape(description)}</td></tr>')
    return '<section class="panel"><h2>Main results</h2><table><tr><th>Sample</th><th>Library</th><th>Primary result</th><th>Description</th></tr>' + ''.join(rows) + '</table></section>'


def _rows(path: Path):
    records = load_correction_stats(path) if path.suffix == ".json" else list(open_tsv(path))
    return ([list(records[0])] + [[row[k] for k in records[0]] for row in records]) if records else [[]]


def _matrix_summary(path: Path):
    if path.name.endswith(".json"):
        from .sparse_matrix import load_summary
        return load_summary(path)[:2]
    if not path.is_file():
        from .sparse_matrix import load_summary
        return load_summary(path.with_name(path.name.replace(".raw_count_matrix.tsv.gz", ".raw_matrix_summary.json")))[:2]
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t"); header = next(reader, [])
        b1 = 0; totals = Counter()
        for row in reader:
            b1 += 1
            for name, value in zip(header[1:], row[1:]): totals[name] += int(value)
    return b1, totals


def _pi_seq_pb_summary(path: Path):
    """Stream PB statistics, retaining only the two Top-20 views."""
    if path.name.endswith(".json") or not path.is_file():
        from .sparse_matrix import load_summary
        return load_summary(path if path.name.endswith(".json") else path.with_name(path.name.replace(".raw_count_matrix.tsv.gz", ".raw_matrix_summary.json")))[2]
    import heapq
    molecule_heap = []; protein_heap = []; molecule_bins = [0] * 8; protein_bins = [0] * 8; total = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t"); next(reader, [])
        for row in reader:
            counts = [int(value) for value in row[1:]]; molecules = sum(counts); proteins = sum(value > 0 for value in counts); total += 1
            molecule_bins[min(molecules, 8) - 1] += 1; protein_bins[min(proteins, 8) - 1] += 1
            for heap, value in ((molecule_heap, molecules), (protein_heap, proteins)):
                item = (value, row[0], molecules, proteins)
                if len(heap) < 20: heapq.heappush(heap, item)
                elif item > heap[0]: heapq.heapreplace(heap, item)
    top = lambda heap, key: [{"PB": item[1], "Molecules": item[2], "Proteins": item[3]} for item in sorted(heap, reverse=True)]
    return {"total": total, "molecule_bins": molecule_bins, "protein_bins": protein_bins, "top_molecules": top(molecule_heap, "Molecules"), "top_proteins": top(protein_heap, "Proteins")}


def _pi_seq_ratios(pi_qc):
    """Read both current and pre-Removed-Ratio PI-seq QC files."""
    mp_ratio = float(pi_qc.get("multi_pi_ratio", pi_qc.get("mp_ratio", 0)) or 0)
    dominant_ratio = float(pi_qc.get("dominant_ratio", 0) or 0)
    removed_ratio = float(pi_qc.get("removed_ratio", pi_qc.get("mp_ratio", 0)) or 0)
    return mp_ratio, dominant_ratio, removed_ratio


def _pi_seq_matrix_section(pb_records, pi_qc, dominance_threshold):
    top_molecule_pbs = pb_records["top_molecules"]
    top_protein_pbs = pb_records["top_proteins"]
    mp_ratio, dominant_ratio, removed_ratio = _pi_seq_ratios(pi_qc)
    total_pbs = pb_records["total"] or 1
    multi_molecule_pbs = total_pbs - pb_records["molecule_bins"][0]
    multi_protein_pbs = total_pbs - pb_records["protein_bins"][0]

    def distribution(field):
        counts = pb_records["molecule_bins" if field == "Molecules" else "protein_bins"]
        labels = [str(value) for value in range(1, 8)] + ["8+"]
        return [{"label": label, "count": count, "ratio": count / total_pbs * 100}
                for label, count in zip(labels, counts)]

    distributions = {"molecules": distribution("Molecules"), "proteins": distribution("Proteins")}
    return f"""<h2 class="section-title">PI-seq Matrix</h2><div class="report-grid"><div class="pi-column"><div class="pi-metrics"><div class="card"><div class="sub">MP Ratio</div><div class="value">{mp_ratio:.1%}</div></div><div class="card"><div class="sub">Dominant MP Ratio</div><div class="value">{dominant_ratio:.1%}</div></div><div class="card"><div class="sub">Removed Ratio</div><div class="value">{removed_ratio:.1%}</div></div></div><section class="panel"><h3>PB Distribution</h3><div class="mini-grid"><div class="card"><div class="sub">PBs with >=2 molecules</div><div class="value">{multi_molecule_pbs / total_pbs:.1%}</div><div class="sub">{multi_molecule_pbs:,} / {total_pbs:,} PBs</div></div><div class="card"><div class="sub">PBs with >=2 proteins</div><div class="value">{multi_protein_pbs / total_pbs:.1%}</div><div class="sub">{multi_protein_pbs:,} / {total_pbs:,} PBs</div></div></div><div id="pi-distribution" style="height:360px"></div></section></div><section class="panel pi-top-panel"><h3>Top 20 PBs</h3><h4>Molecules per PB</h4><div id="pi-molecules" style="height:260px"></div><h4>Proteins per PB</h4><div id="pi-proteins" style="height:260px"></div></section></div><script>const piMoleculePb={json.dumps(top_molecule_pbs)};const piProteinPb={json.dumps(top_protein_pbs)};const piDistribution={json.dumps(distributions)};const piMoleculeNames=piMoleculePb.map(row=>row.PB);const piMolecules=piMoleculePb.map(row=>row.Molecules);const piProteinNames=piProteinPb.map(row=>row.PB);const piProteins=piProteinPb.map(row=>row.Proteins);const piBarLayout=xTitle=>({{margin:{{l:12,r:30,t:16,b:52}},paper_bgcolor:'#ffffff',plot_bgcolor:'#ffffff',font:{{family:'Arial',color:'#172033'}},xaxis:{{title:xTitle,gridcolor:'#e2e8f0',rangemode:'tozero'}},yaxis:{{showticklabels:false,showgrid:false,zeroline:false,autorange:'reversed'}}}});Plotly.newPlot('pi-molecules',[{{type:'bar',orientation:'h',x:piMolecules,y:piMoleculeNames,marker:{{color:'#636EFA'}},hovertemplate:'PB: %{{y}}<br>Molecules: %{{x:,}}<extra></extra>'}}],piBarLayout('Molecules'),{{responsive:true}});Plotly.newPlot('pi-proteins',[{{type:'bar',orientation:'h',x:piProteins,y:piProteinNames,marker:{{color:'#EF553B'}},hovertemplate:'PB: %{{y}}<br>Proteins: %{{x:,}}<extra></extra>'}}],piBarLayout('Proteins'),{{responsive:true}});const piColors=['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A','#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52'];const piDistributionTraces=piDistribution.molecules.map((entry,index)=>({{type:'bar',orientation:'h',y:['Molecules','Proteins'],x:[entry.ratio,piDistribution.proteins[index].ratio],name:entry.label,marker:{{color:piColors[index]}},customdata:[[entry.count,entry.ratio],[piDistribution.proteins[index].count,piDistribution.proteins[index].ratio]],hovertemplate:'%{{y}} per PB<br>Count: '+entry.label+'<br>PBs: %{{customdata[0]:,}}<br>Ratio: %{{customdata[1]:.2f}}%<extra></extra>'}}));Plotly.newPlot('pi-distribution',piDistributionTraces,{{barmode:'stack',margin:{{l:80,r:24,t:20,b:92}},paper_bgcolor:'#ffffff',plot_bgcolor:'#ffffff',font:{{family:'Arial',color:'#172033'}},xaxis:{{title:'PB ratio (%)',range:[0,100],gridcolor:'#e2e8f0'}},yaxis:{{showgrid:false}},legend:{{orientation:'h',x:0,xanchor:'left',y:-.38,yanchor:'top',traceorder:'normal'}}}},{{responsive:true}});</script>"""


def _number(value):
    return int(float(value or 0))


def _barcode2_names(barcode2_segments):
    names = {}
    for segment in barcode2_segments:
        if getattr(segment, "whitelist", None):
            names.update(load_whitelist_names(segment.whitelist))
    return names


def _annotation_rate_and_details(config, sample, stats, barcode2_segments):
    whitelist_reads, missing_reads, details = _barcode2_annotation_gap(
        config, sample, stats, barcode2_segments)
    return ((whitelist_reads - missing_reads) / whitelist_reads if whitelist_reads else 0), details


def _annotation_warning_html(sample, rate, details, barcode2_names):
    rows = "".join(
        '<tr><td>{}</td><td>{:,}</td></tr>'.format(
            html.escape(barcode2_names.get(row["barcode2_sequence"], row["barcode2_sequence"])),
            int(row["reads_count"]),
        )
        for row in details
    )
    return (
        '<section class="panel"><h2>Warnings · {}</h2><p>Annotated is below 90% ({:.1%}). '
        'Top 10 FB_not_in_anno:</p><table><tr><th>FB_not_in_anno</th><th>Count</th></tr>{}</table></section>'
    ).format(html.escape(sample), rate, rows)


def _molecule_totals(config, sample):
    from .checkpoint import molecule_source
    preferred = molecule_source(config, sample)
    source = resolve_table(preferred, preferred.with_name(preferred.name.replace(".parquet", ".tsv.gz")))
    molecules = reads = 0
    for row in open_records(source):
        molecules += 1
        reads += int(row["reads_count"])
    return molecules, reads


def _analysis_rows(config, sample, dirs):
    if not getattr(config, "downsample_enabled", True):
        return [], []
    from .downsample import curve_rows
    return (curve_rows(list(open_tsv(dirs["downsample"] / f"{sample}.downsample_metrics.tsv"))),
            [selection_point(config, sample)])


def _pair_mapping_rows(config, samples):
    rows = [["Sample", "Total reads", "Valid reads", "Raw molecules", "Resolved molecules",
             "Supporting reads", "Retained PBs", "Retained CBs", "Retained molecules",
             "Retained supporting reads", "CB aliases", "Dominance basis", "Downsample status"]]
    for sample in samples:
        dirs = sample_dirs(config.output_dir, sample)
        payload = json.loads((dirs["matrix"] / f"{sample}.pb_cb_summary.json").read_text())
        stats = next((row for row in load_correction_stats(correction_qc_path(dirs["corrected"], sample))
                      if row["scope"] == "combined"), {})
        rows.append([sample, int(stats.get("total_reads", 0)), int(stats.get("valid_reads", 0)),
                     *[payload[key] for key in ("raw_molecules", "resolved_molecules", "supporting_reads",
                        "retained_pbs", "retained_cbs", "retained_molecules", "retained_supporting_reads", "cb_aliases")],
                     payload["dominance_basis"], "completed" if config.downsample_enabled else "not run"])
    return rows


def _pair_mapping_report(config, samples, sample_name=None):
    """PB-CB reports consume the mapping contract, independently of matrices."""
    from itertools import islice
    directory = (sample_dirs(config.output_dir, sample_name)["report"] if sample_name
                 else batch_report_dir(config.output_dir))
    directory.mkdir(parents=True, exist_ok=True)
    prefix = f"{sample_name}.report" if sample_name else "TagForge_batch_report"
    xlsx, report_html = directory / f"{prefix}.xlsx", directory / f"{prefix}.html"
    overview = _pair_mapping_rows(config, samples)
    if sample_name:
        sheets = [("Results", _result_rows(config, samples)), ("Overview", overview)]
        if config.downsample_enabled:
            sheets.append(("Downsample Definitions", [["Analysis", SAMPLING_NOTE]]))
            d = sample_dirs(config.output_dir, sample_name)["downsample"]
            for title, suffix in (("Downsample Metrics", "downsample_metrics.tsv"), ("Downsample Curve", "downsample_curve.tsv")):
                sheets.append((title, _rows(d / f"{sample_name}.{suffix}")))
            sheets.append(("Selected Downsample", _point_rows(config, sample_name)))
        for sample in samples:
            dirs = sample_dirs(config.output_dir, sample)
            mapping = list(islice(open_tsv(dirs["matrix"] / f"{sample}.pb_cb_map.tsv.gz"), 1000))
            sheets.append((sample, [list(mapping[0])] + [list(row.values()) for row in mapping]
                           if mapping else [[config.target_name("barcode1"), config.target_name("barcode2"), "molecule_count", "reads_count"]]))
        write_xlsx(xlsx, sheets)
    else:
        _write_batch_workbook(config, samples)
    cells = lambda row, tag: "".join(f"<{tag}>{html.escape(str(value))}</{tag}>" for value in row)
    table = '<tr>' + cells(overview[0], 'th') + '</tr>' + ''.join('<tr>' + cells(row, 'td') + '</tr>' for row in overview[1:])
    with atomic_text(report_html) as handle:
        handle.write(f'<!doctype html><html><head><meta charset="utf-8"><title>TagForge PB–CB report</title>'
                     '<style>body{font:15px Arial,sans-serif;margin:32px}table{border-collapse:collapse}th,td{padding:10px;border:1px solid #ddd}</style></head><body>'
                     f'<h1>PB–CB mapping</h1><p>Dominance uses read support; retained PB–CB–UMI molecules are merged after CB alias resolution. Threshold: {config.pb_cb_dominance_threshold} (strictly greater).</p>'
                     f'{_result_html(config, samples, directory)}'
                     f'<p>{html.escape(SAMPLING_NOTE) if config.downsample_enabled else "Downsample not run."}</p>'
                     f'<table>{table}</table><p>The workbook previews up to 1,000 retained mappings per sample. Full mappings remain in pb_cb_map.tsv.gz.</p></body></html>')
    return [xlsx, report_html]


def report_sample(config: TagForgeConfig, sample_name: str):
    if config.pb_cb_enabled:
        return _pair_mapping_report(config, [sample_name], sample_name)
    barcode1_name = config.target_name("barcode1")
    dirs = sample_dirs(config.output_dir, sample_name)
    extraction_stats = dirs["extracted"] / f"{sample_name}.extraction_stats.tsv"
    stats = correction_qc_path(dirs["corrected"], sample_name)
    metrics = dirs["downsample"] / f"{sample_name}.downsample_metrics.tsv"
    matrix = dirs["matrix"] / f"{sample_name}.raw_matrix_summary.json"
    xlsx = dirs["report"] / f"{sample_name}.report.xlsx"
    report_html = dirs["report"] / f"{sample_name}.report.html"
    stat_records = load_correction_stats(stats)
    metric_records, point_records = _analysis_rows(config, sample_name, dirs)
    b1_count, features = _matrix_summary(matrix)
    pi_seq_pb_records = _pi_seq_pb_summary(matrix) if config.pi_seq_enabled else None
    pi_qc = next(open_tsv(dirs["detail"] / f"{sample_name}.pi_seq_qc.tsv"), {}) if config.pi_seq_enabled else {}
    molecule_count, read_count = _molecule_totals(config, sample_name)
    combined = next((r for r in stat_records if r["scope"] == "combined"), {})
    barcode2_segments = [s for s in config.segments if s.target == "barcode2"]
    annotated_rate, missing_annotation_details = _annotation_rate_and_details(
        config, sample_name, stat_records, barcode2_segments)
    barcode2_names = _barcode2_names(barcode2_segments)
    summary = [
        ["Metric", "Value"], ["Sample", sample_name], ["Total Reads", combined.get("total_reads", 0)],
        ["Valid reads", combined.get("valid_reads", 0)], ["Combined valid rate", combined.get("valid_rate", 0)],
        ["Deduplicated molecules", molecule_count], ["Reads supporting molecules", read_count],
        [f"{barcode1_name} rows", b1_count], ["Features", len(features)], ["Annotated", annotated_rate],
        ["Downsample status", "completed" if config.downsample_enabled else "not run"],
    ]
    feature_rows = [["Feature", "Molecule count"]] + [[k, v] for k, v in features.most_common()]
    sheets = [("Overview", summary), ("Extraction Stats", _rows(extraction_stats)),
              ("Correction Stats", _rows(stats)), ("UMI Stats", [["Metric","Value"],["Molecules",molecule_count],["Supporting reads",read_count]]),
              ("Top Features", feature_rows[:1001])]
    if config.downsample_enabled:
        sheets.extend([("Downsample Metrics", _rows(metrics)), ("Selected Downsample", _point_rows(config, sample_name))])
        for title, suffix in (("Downsample Curve", "downsample_curve.tsv"),):
            path = dirs["downsample"] / f"{sample_name}.{suffix}"
            if path.is_file():
                sheets.append((title, _rows(path)))
        sheets.append(("Downsample Definitions", [["Metric", "Definition"], ["Analysis", SAMPLING_NOTE],
            ["sequencing_saturation / duplication_ratio (%)", "100 * (sampled reads - detected molecules) / sampled reads; zero when sampled reads = 0"],
            ["non_singleton_molecule_fraction (%)", "100 * molecules supported by >1 sampled read / detected molecules; zero when molecules = 0"],
            ["Expected counts", "Analytical expectations of reads, molecules and singletons; no expected fractions are reported"],
            ["Legacy optimal_* names", "Compatibility aliases for selected results under the recorded criterion"]]))
    if config.pi_seq_enabled:
        top_pbs = pi_seq_pb_records["top_molecules"]
        sheets.append(("PI-seq Matrix", [["PB", "Molecules", "Proteins"]] + [[row["PB"], row["Molecules"], row["Proteins"]] for row in top_pbs]))
    ratios = [float(r["downsample_ratio"]) for r in metric_records]
    saturation = [float(r["sequencing_saturation"]) for r in metric_records]
    non_singleton = [float(r["non_singleton_molecule_fraction"]) for r in metric_records]
    umi_types = [float(r["umi_types"]) for r in metric_records]
    singletons = [float(r["umi_detected_once"]) for r in metric_records]
    valid_rate = float(combined.get("valid_rate", 0) or 0)
    optimal = point_records[0] if point_records else {}
    optimal_saturation = float(optimal.get("sequencing_saturation", 0) or 0)
    max_non_singleton = float(optimal.get("max_non_singleton_molecule_fraction", 0) or 0)
    optimal_ratio = float(optimal.get("optimal_downsample_ratio", 0) or 0)
    feature_total = sum(features.values())
    flags = []
    if valid_rate < 0.5: flags.append("Combined barcode/UMI valid rate is below 50%.")
    if read_count and molecule_count / read_count < 0.1: flags.append("High duplication: fewer than 10% of supporting reads are unique molecules.")
    if annotated_rate < 0.9:
        flags.append(f"Annotated is below 90% ({annotated_rate:.1%}).")
    warning_rows = [["Warning", "FB_not_in_anno", "Count"]]
    warning_rows.extend([[flag, "", ""] for flag in flags if not flag.startswith("Annotated")])
    if annotated_rate < 0.9:
        warning_rows.append([f"Annotated is below 90% ({annotated_rate:.1%}).", "", ""])
        warning_rows.extend([["", barcode2_names.get(row["barcode2_sequence"], row["barcode2_sequence"]), row["reads_count"]]
                            for row in missing_annotation_details[:10]])
    sheets.insert(0, ("Results", _result_rows(config, [sample_name])))
    sheets.append(("Warnings", warning_rows))
    write_xlsx(xlsx, sheets)
    extraction_records = list(open_tsv(extraction_stats))
    barcode_segments = {segment.name: segment for segment in config.segments if segment.target in {"barcode1", "barcode2"}}
    correction_rates = {row["scope"]: float(row.get("valid_rate", 0) or 0) for row in stat_records}
    validation_rows = ''.join(
        f'<tr><td>{html.escape(barcode_segments[row["segment"]].target_name or barcode_segments[row["segment"]].target)}</td>'
        f'<td>{html.escape(row["segment"])}</td>'
        f'<td>{format(float(row["linker_success_rate"]), ".1%") if row["linker_attempted"] != "0" else "/"}</td>'
        f'<td>{correction_rates.get(row["segment"], 0):.1%}</td></tr>'
        for row in extraction_records if row["segment"] in barcode_segments
    )
    feature_table_rows = ''.join(
        f'<tr><td>{html.escape(name)}</td><td>{count:,}</td><td>{(count / feature_total if feature_total else 0):.2%}</td></tr>'
        for name, count in features.most_common()
    )
    top_features = features.most_common(20)
    warning_html = ""
    if flags:
        annotation_warning_rows = "".join(
            f'<tr><td>{html.escape(barcode2_names.get(row["barcode2_sequence"], row["barcode2_sequence"]))}</td><td>{int(row["reads_count"]):,}</td></tr>'
            for row in missing_annotation_details[:10]
        )
        annotation_warning = (
            f'<p>Annotated is below 90% ({annotated_rate:.1%}). Top 10 FB_not_in_anno:</p><table><tr><th>FB_not_in_anno</th><th>Count</th></tr>{annotation_warning_rows}</table>'
            if annotated_rate < 0.9 else "")
        warning_html = f'<section class="panel"><h2>Warnings</h2><ul>{"".join(f"<li>{html.escape(flag)}</li>" for flag in flags)}</ul>{annotation_warning}</section>'
    selection_cards = f'''<div class="card"><div class="sub">Sequencing Saturation / Duplication at selected ratio</div><div class="value">{optimal_saturation:.1f}%</div></div><div class="card"><div class="sub">Selected ratio non-singleton fraction</div><div class="value">{max_non_singleton:.1f}%</div></div><div class="card"><div class="sub">Selected downsample ratio</div><div class="value">{optimal_ratio:.4g}</div></div>''' if config.downsample_enabled else ''
    page = f'''<!doctype html><html><head><meta charset="utf-8"><title>TagForge · {html.escape(sample_name)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script><style>body{{margin:0;background:#f5f7fb;color:#172033;font:15px Arial,sans-serif}}main{{max-width:1180px;margin:auto;padding:36px}}h1{{font-size:36px;margin:0}}.section-title{{margin:36px 0 16px;padding-top:28px;border-top:1px solid #dce3ee}}.sub{{color:#64748b}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:14px;margin:28px 0}}.mini-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px;margin:14px 0}}.pi-column{{min-width:0}}.pi-metrics{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px;margin:0 0 14px}}.pi-top-panel h3{{margin:0 0 8px}}.pi-top-panel h4{{margin:8px 0 0}}.report-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(480px,1fr));gap:14px;align-items:stretch;margin:16px 0}}.report-grid .panel{{margin:0;display:flex;flex-direction:column;min-width:0;overflow:hidden}}#sat,#features{{width:100%;min-width:0}}.scroll-window{{height:520px;overflow:auto;border:1px solid #e6ebf2;border-radius:8px}}.scroll-window th{{position:sticky;top:0;z-index:1}}.icon-button{{display:inline-flex;align-items:center;justify-content:center;width:36px;height:36px;padding:0;border:1px solid #94a3b8;border-radius:8px;background:white;color:#172033;cursor:pointer;font-size:19px}}.card,.panel{{background:white;border:1px solid #dce3ee;border-radius:14px;padding:18px}}.value{{font-size:28px;font-weight:750;color:#155e75}}.panel{{margin:16px 0;overflow:auto}}table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;border-bottom:1px solid #e6ebf2;text-align:left}}th{{color:#475569;background:#f8fafc}}</style></head><body><main>
<style>.theme{{float:right}}body.dark{{background:#07111f;color:#e5edf8}}body.dark .card,body.dark .panel{{background:#0f1d30;border-color:#24344d}}body.dark .sub{{color:#94a3b8}}body.dark th{{background:#17263b;color:#93c5fd}}body.dark td,body.dark th,body.dark .scroll-window{{border-color:#24344d}}body.dark .value{{color:#5eead4}}body.dark .section-title{{border-color:#24344d}}body.dark .icon-button{{background:#0f1d30;color:#e5edf8;border-color:#64748b}}</style><button class="theme icon-button" aria-label="Switch to dark theme" title="Switch to dark theme">☾</button><p class="sub">TAGFORGE / SAMPLE REPORT</p><h1>{html.escape(sample_name)}</h1>{_result_html(config, [sample_name], dirs["report"])}<div class="cards">
<div class="card"><div class="sub">Total Reads</div><div class="value">{int(combined.get('total_reads', 0) or 0):,}</div></div><div class="card"><div class="sub">Features</div><div class="value">{len(features):,}</div></div><div class="card"><div class="sub">Valid Ratio</div><div class="value">{valid_rate:.1%}</div></div><div class="card"><div class="sub">Annotated</div><div class="value">{annotated_rate:.1%}</div></div><div class="card"><div class="sub">Molecules</div><div class="value">{molecule_count:,}</div></div><div class="card"><div class="sub">{html.escape(barcode1_name)}</div><div class="value">{b1_count:,}</div></div>{selection_cards}</div>
<div class="report-grid"><section class="panel"><h2>Barcode Validation</h2><table><tr><th>Barcode</th><th>Segment</th><th>Linker</th><th>Final</th></tr>{validation_rows}</table></section><section class="panel"><h2>Downsample</h2>{'<p>Not run. Counts use full-depth molecules.</p>' if not config.downsample_enabled else '<p>'+html.escape(SAMPLING_NOTE)+'</p>'}<div id="sat" style="height:430px" {'hidden' if not config.downsample_enabled else ''}></div></section></div>
<div class="report-grid"><section class="panel"><h2>Top 20 Feature</h2><div id="features" style="height:520px"></div></section><section class="panel"><h2>Feature Counts</h2><div class="scroll-window"><table><tr><th>Feature</th><th>Counts</th><th>Ratio</th></tr>{feature_table_rows}</table></div></section></div>
{warning_html}
{_pi_seq_matrix_section(pi_seq_pb_records, pi_qc, config.pi_seq_dominance_threshold) if config.pi_seq_enabled else ""}
<script>const featureData={{names:{json.dumps([name for name, _ in top_features])},counts:{json.dumps([count for _, count in top_features])}}};const featurePlot=document.getElementById('features');Plotly.newPlot(featurePlot,[{{type:'bar',orientation:'h',x:featureData.counts,y:featureData.names,marker:{{color:'#0f766e'}},hovertemplate:'%{{y}}<br>%{{x:,}} molecules<extra></extra>'}}],{{autosize:true,yaxis:{{autorange:'reversed',automargin:true}},margin:{{l:60,r:30,t:16,b:50}},paper_bgcolor:'#ffffff',plot_bgcolor:'#ffffff',font:{{family:'Arial',color:'#172033'}}}},{{responsive:true}});function toggleTheme(button){{const dark=document.body.classList.toggle('dark');button.textContent=dark?'☀':'☾';button.setAttribute('aria-label',dark?'Switch to light theme':'Switch to dark theme');button.title=button.getAttribute('aria-label');const grid=dark?'#24344d':'#e2e8f0';const surface=dark?'#0f1d30':'#ffffff';const font=dark?'#cbd5e1':'#172033';Plotly.relayout('sat',{{'paper_bgcolor':surface,'plot_bgcolor':surface,'font.color':font,'font.family':'Arial','xaxis.gridcolor':grid,'yaxis.gridcolor':grid,'yaxis2.gridcolor':grid,'yaxis3.gridcolor':grid}});Plotly.relayout('features',{{'paper_bgcolor':surface,'plot_bgcolor':surface,'font.color':font,'font.family':'Arial'}});if(document.getElementById('pi-molecules')){{const piTheme={{'paper_bgcolor':surface,'plot_bgcolor':surface,'font.color':font,'font.family':'Arial','xaxis.gridcolor':grid,'yaxis.gridcolor':grid}};Plotly.relayout('pi-molecules',piTheme);Plotly.relayout('pi-proteins',piTheme);Plotly.relayout('pi-distribution',{{...piTheme,'xaxis.gridcolor':grid,'yaxis.gridcolor':grid}});}}}}document.querySelector('.theme').onclick=function(){{toggleTheme(this)}};Plotly.newPlot('sat',[{{x:{json.dumps(ratios)},y:{json.dumps(saturation)},name:'Sequencing Saturation / Duplication',mode:'lines',yaxis:'y3',line:{{color:'#0f766e',width:3}}}},{{x:{json.dumps(ratios)},y:{json.dumps(non_singleton)},name:'Non-singleton Molecule Fraction',mode:'lines',yaxis:'y3',line:{{color:'#ea580c',width:3}}}},{{x:{json.dumps(ratios)},y:{json.dumps(umi_types)},name:'UMI Types',mode:'lines',yaxis:'y2',line:{{color:'#2563eb',width:3}}}},{{x:{json.dumps(ratios)},y:{json.dumps(singletons)},name:'UMI detected once',mode:'lines',yaxis:'y',line:{{color:'#9333ea',width:3}}}}],{{margin:{{l:60,r:60,t:12,b:70}},paper_bgcolor:'#ffffff',plot_bgcolor:'#ffffff',font:{{family:'Arial',color:'#172033'}},legend:{{orientation:'h',x:0.5,xanchor:'center',y:-0.15}},xaxis:{{title:'Downsample ratio',gridcolor:'#e2e8f0'}},yaxis:{{title:'UMI detected once',domain:[0.05,0.48],gridcolor:'#e2e8f0'}},yaxis2:{{title:'UMI Types',overlaying:'y',side:'right',domain:[0.05,0.48],gridcolor:'#e2e8f0'}},yaxis3:{{title:'',ticksuffix:'%',domain:[0.56,1],range:[0,100],gridcolor:'#e2e8f0'}}}},{{responsive:true}})</script></main></body></html>'''
    with atomic_text(report_html) as handle: handle.write(page)
    return [xlsx, report_html]


def _cluster_feature_order(names, bulk, samples):
    totals={s:sum(bulk[s].values()) or 1 for s in samples}
    return sorted(names,key=lambda f: tuple(round(bulk[s].get(f,0)/totals[s],6) for s in samples))


def batch_report(config, sample_names):
    if config.pb_cb_enabled:
        return _pair_mapping_report(config, sample_names)
    out=batch_report_dir(config.output_dir); out.mkdir(parents=True,exist_ok=True)
    b1=config.target_name('barcode1'); segs=[s for s in config.segments if s.target in {'barcode1','barcode2'}]
    barcode2_segments=[s for s in config.segments if s.target == 'barcode2']
    barcode2_names=_barcode2_names(barcode2_segments)
    bulk={}; rows=[]; downsample_curves={}; annotation_warnings=[]; pi_qcs={}; pi_pb_distributions={}
    for sample in sample_names:
        dirs=sample_dirs(config.output_dir,sample); stats=load_correction_stats(correction_qc_path(dirs['corrected'],sample))
        rates={r['scope']:float(r.get('valid_rate',0) or 0) for r in stats}; combined=next((r for r in stats if r['scope']=='combined'),{})
        annotated_rate, missing_annotation_details = _annotation_rate_and_details(config, sample, stats, barcode2_segments)
        if annotated_rate < 0.9: annotation_warnings.append((sample, annotated_rate, missing_annotation_details[:10]))
        matrix_path = dirs['matrix']/f'{sample}.raw_matrix_summary.json'
        b1_count, feats=_matrix_summary(matrix_path); bulk[sample]=feats
        if config.pi_seq_enabled:
            pb_records = _pi_seq_pb_summary(matrix_path)
            total_pbs = pb_records["total"] or 1
            pi_pb_distributions[sample] = {
                "Molecules": [count / total_pbs * 100 for count in pb_records["molecule_bins"]],
                "Proteins": [count / total_pbs * 100 for count in pb_records["protein_bins"]],
            }
        metrics=_analysis_rows(config, sample, dirs)[0]; downsample_curves[sample]={'x':[float(r['downsample_ratio']) for r in metrics],'saturation':[float(r['sequencing_saturation']) for r in metrics],'non_singleton':[float(r['non_singleton_molecule_fraction']) for r in metrics],'umi_types':[float(r['umi_types']) for r in metrics],'singletons':[float(r['umi_detected_once']) for r in metrics]}
        point_records = _analysis_rows(config, sample, dirs)[1]
        optimal = point_records[0] if point_records else {}
        pi_qc_path = dirs['detail'] / f'{sample}.pi_seq_qc.tsv'
        pi_qcs[sample] = next(open_tsv(pi_qc_path), {}) if config.pi_seq_enabled and pi_qc_path.is_file() else {}
        rows.append((sample,_number(combined.get('total_reads')),float(combined.get('valid_rate',0) or 0),sum(feats.values()),b1_count,len(feats),(float(optimal.get('max_non_singleton_molecule_fraction',0) or 0) if config.downsample_enabled else None),(float(optimal.get('sequencing_saturation',0) or 0) if config.downsample_enabled else None),(float(optimal.get('optimal_downsample_ratio',0) or 0) if config.downsample_enabled else None),[rates.get(s.name,0) for s in segs],annotated_rate))
    names=sorted(set().union(*[set(v) for v in bulk.values()]),key=lambda f:(-sum(bulk[s].get(f,0) for s in sample_names),f))
    raw=[[bulk[s].get(f,0) for s in sample_names] for f in names]; totals=[sum(bulk[s].values()) or 1 for s in sample_names]; ratio=[[v/totals[i]*100 for i,v in enumerate(row)] for row in raw]
    heat_names=names[:100]; cluster=_cluster_feature_order(heat_names,bulk,sample_names); cluster_raw=[[bulk[s].get(f,0) for s in sample_names] for f in cluster]; cluster_ratio=[[v/totals[i]*100 for i,v in enumerate(row)] for row in cluster_raw]
    over=''.join(
        '<tr><td>'+html.escape(r[0])+'</td><td>'+f'{r[1]:,}'+'</td><td>'+f'{r[2]:.1%}'+'</td>'
        + ''.join('<td>'+f'{rate:.1%}'+'</td>' for rate in r[9])
        + '<td>'+f'{r[3]:,}'+'</td><td>'+f'{r[4]:,}'+'</td><td>'+f'{r[5]:,}'+'</td>'
        + '<td>'+f'{r[10]:.1%}'+'</td><td>'+('Not run' if r[8] is None else f'{r[8]:.1%}')+'</td><td>'+('' if r[7] is None else f'{r[7]:.1f}%')+'</td><td>'+('' if r[6] is None else f'{r[6]:.1f}%')+'</td></tr>'
        for r in rows
    )
    counts=''.join('<tr><td>'+html.escape(f)+'</td>'+''.join('<td>'+f'{v:,}'+'</td>' for v in row)+'<td>'+f'{sum(row):,}'+'</td></tr>' for f,row in zip(names,raw))
    xlsx=out/'TagForge_batch_report.xlsx'; hp=out/'TagForge_batch_report.html'
    batch_warning_html = "".join(
        _annotation_warning_html(sample, rate, details, barcode2_names)
        for sample, rate, details in annotation_warnings
    )
    pi_matrix_rows = [['Sample', 'MP Ratio', 'Dominant MP Ratio', 'Removed Ratio', '≥2 Molecule PB Ratio', '≥2 Protein PB Ratio']]
    if config.pi_seq_enabled:
        for sample in sample_names:
            mp, dominant, removed = _pi_seq_ratios(pi_qcs[sample])
            distribution = pi_pb_distributions[sample]
            pi_matrix_rows.append([sample, mp, dominant, removed, 100 - distribution['Molecules'][0], 100 - distribution['Proteins'][0]])
    _write_batch_workbook(config, sample_names)
    pi_batch_html = ''
    if config.pi_seq_enabled:
        labels = [str(value) for value in range(1, 8)] + ['8+']
        pi_metrics_table = '<div class="overview"><table><tr><th>Sample</th><th>MP Ratio</th><th>Dominant MP Ratio</th><th>Removed Ratio</th><th>≥2 Molecule PB Ratio</th><th>≥2 Protein PB Ratio</th></tr>' + ''.join('<tr><td>{}</td><td>{:.1%}</td><td>{:.1%}</td><td>{:.1%}</td><td>{:.1f}%</td><td>{:.1f}%</td></tr>'.format(html.escape(row[0]), *row[1:]) for row in pi_matrix_rows[1:]) + '</table></div>'
        pi_batch_html = """<section class="panel"><h2>PI-seq Matrix</h2>""" + pi_metrics_table + """<div class="grid"><div id="pi-batch-molecules" style="height:360px"></div><div id="pi-batch-proteins" style="height:360px"></div></div></section><script>const piBatch=""" + json.dumps(pi_pb_distributions) + """;const piLabels=""" + json.dumps(labels) + """;const piColors=['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A','#19D3F3','#FF6692','#B6E880'];const piSamples=Object.keys(piBatch);const drawPiBatch=(id,field,title)=>Plotly.newPlot(id,piLabels.map((label,index)=>({type:'bar',orientation:'h',y:piSamples,x:piSamples.map(sample=>piBatch[sample][field][index]),name:label,marker:{color:piColors[index]},hovertemplate:'%{y}<br>Count: '+label+'<br>PB ratio: %{x:.2f}%<extra></extra>'})),{barmode:'stack',margin:{l:120,r:20,t:40,b:70},title,font:{family:'Arial',color:'#172033'},paper_bgcolor:'#ffffff',plot_bgcolor:'#ffffff',xaxis:{title:'PB ratio (%)',range:[0,100]},legend:{orientation:'h',x:0,xanchor:'left',y:-.25,traceorder:'normal'}},{responsive:true});drawPiBatch('pi-batch-molecules','Molecules','Molecules per PB');drawPiBatch('pi-batch-proteins','Proteins','Proteins per PB');</script>"""

    page=f'''<!doctype html><html><head><script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script><style>body{{margin:0;background:#f5f7fb;color:#172033;font:15px Arial,sans-serif}}main{{max-width:1200px;margin:auto;padding:36px}}h1{{margin:0 0 8px}}h2{{margin:0 0 14px}}.panel{{background:#fff;border:1px solid #dce3ee;border-radius:14px;padding:18px;margin:16px 0}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(480px,1fr));gap:18px;margin:20px 0}}.mini-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px}}.shared-legend{{display:flex;flex-wrap:wrap;justify-content:center;gap:14px;margin-top:14px}}.legend-item{{display:inline-flex;align-items:center;gap:6px}}.legend-swatch{{width:12px;height:12px;border-radius:50%}}.grid .panel{{margin:0}}table{{border-collapse:collapse;width:100%}}th,td{{padding:9px;border-bottom:1px solid #e6ebf2;text-align:left;white-space:nowrap}}th{{background:#f8fafc;color:#475569}}.controls{{display:flex;align-items:center;gap:8px;margin:0 0 14px}}.scroll{{height:520px;overflow:auto;border:1px solid #e6ebf2;border-radius:8px}}#ds{{margin-top:4px}}.overview{{max-height:300px;overflow:auto;border:1px solid #e6ebf2;border-radius:8px}}button{{font:inherit;padding:6px 10px;border:1px solid #94a3b8;border-radius:8px;background:#fff;cursor:pointer}}button.active{{background:#0f766e;color:#fff;border-color:#0f766e}}.theme{{float:right;width:36px;height:36px;padding:0;font-size:19px}}body.dark{{background:#07111f;color:#e5edf8}}body.dark .panel{{background:#0f1d30;border-color:#24344d}}body.dark .sub{{color:#94a3b8}}body.dark th{{background:#17263b;color:#93c5fd}}body.dark td,body.dark th,body.dark .scroll,body.dark .overview{{border-color:#24344d}}body.dark button{{background:#0f1d30;color:#e5edf8;border-color:#64748b}}body.dark button.active{{background:#0f766e}}</style></head><body><main><button class="theme">☾</button><p class="sub">TAGFORGE / BATCH RESULTS</p><h1>{len(sample_names)} samples</h1>{_result_html(config, sample_names, out)}<section class="panel"><h2>Overview</h2><div class="overview"><table><tr><th>Sample</th><th>Total Reads</th><th>Valid</th>{''.join('<th>'+html.escape(s.name)+'</th>' for s in segs)}<th>Molecules</th><th>{html.escape(b1)}</th><th>Features</th><th>Annotated</th><th>Selected DS</th><th>Read-level Sat. / Dup.</th><th>Selected non-singleton (%)</th></tr>{over}</table></div></section><section class="panel"><h2>Downsample</h2><p>{html.escape(SAMPLING_NOTE) if config.downsample_enabled else "Not run. Counts use full-depth molecules."}</p><div id="ds" style="height:400px"></div></section><div class="grid"><section class="panel"><h2>Feature Counts</h2><div class="controls"><button class="active mode" data-c="raw">Counts</button><button class="mode" data-c="ratio">Percent</button></div><div class="scroll"><table><thead><tr><th>Feature</th>{''.join('<th>'+html.escape(s)+'</th>' for s in sample_names)}<th id="aggregate-header">Total</th></tr></thead><tbody id="tbody">{counts}</tbody></table></div></section><section class="panel"><h2>Feature Heatmap</h2><div id="hm" style="height:600px"></div></section></div>{batch_warning_html}{pi_batch_html}<script>const curves={json.dumps(downsample_curves)};const plotlyColors=['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A','#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52'];const curveEntries=Object.entries(curves);const metrics=[['saturation','Read-level Saturation / Duplication (%)',[0,100],'x','y'],['non_singleton','Non-singleton Molecule Fraction (%)',[0,100],'x2','y2'],['umi_types','UMI Types',null,'x3','y3'],['singletons','UMI detected once',null,'x4','y4']];const dsTraces=metrics.flatMap(([key,title,range,xaxis,yaxis],metricIndex)=>curveEntries.map(([name,data],sampleIndex)=>({{x:data.x,y:data[key],name,mode:'lines',legendgroup:name,showlegend:metricIndex===0,line:{{width:3,color:plotlyColors[sampleIndex%plotlyColors.length]}},xaxis,yaxis}})));Plotly.newPlot('ds',dsTraces,{{margin:{{l:45,r:20,t:52,b:90}},legend:{{orientation:'h',y:-.18,groupclick:'togglegroup'}},font:{{family:'Arial'}},paper_bgcolor:'#fff',plot_bgcolor:'#fff',annotations:[{{text:'Read-level Saturation / Duplication (%)',x:.105,y:1.12,xref:'paper',yref:'paper',showarrow:false}},{{text:'Non-singleton Molecule Fraction (%)',x:.365,y:1.12,xref:'paper',yref:'paper',showarrow:false}},{{text:'UMI Types',x:.635,y:1.12,xref:'paper',yref:'paper',showarrow:false}},{{text:'UMI detected once',x:.905,y:1.12,xref:'paper',yref:'paper',showarrow:false}}],xaxis:{{domain:[0,.21],anchor:'y'}},yaxis:{{domain:[0,1],anchor:'x',range:[0,100]}},xaxis2:{{domain:[.26,.47],anchor:'y2'}},yaxis2:{{domain:[0,1],anchor:'x2',range:[0,100]}},xaxis3:{{domain:[.53,.74],anchor:'y3'}},yaxis3:{{domain:[0,1],anchor:'x3'}},xaxis4:{{domain:[.79,1],anchor:'y4'}},yaxis4:{{domain:[0,1],anchor:'x4'}}}},{{responsive:true}});const downsampleLayout=()=>{{const compact=document.getElementById('ds').clientWidth<1000;const domains=compact?[[0,.45],[.55,1],[0,.45],[.55,1]]:[[0,.21],[.26,.47],[.53,.74],[.79,1]];const yDomains=compact?[[.58,1],[.58,1],[0,.39],[0,.39]]:[[0,1],[0,1],[0,1],[0,1]];const titleY=compact?[1.04,1.04,.43,.43]:[1.04,1.04,1.04,1.04];const titles=['Read-level Saturation / Duplication (%)','Non-singleton Molecule Fraction (%)','UMI Types','UMI detected once'];const layout={{margin:{{l:55,r:25,t:compact?52:52,b:90}},annotations:titles.map((text,index)=>({{text,x:(domains[index][0]+domains[index][1])/2,y:titleY[index],xref:'paper',yref:'paper',xanchor:'center',yanchor:'bottom',showarrow:false,font:{{size:12}}}})),xaxis:{{domain:domains[0],anchor:'y',rangemode:'tozero'}},yaxis:{{domain:yDomains[0],anchor:'x',range:[0,100]}},xaxis2:{{domain:domains[1],anchor:'y2',rangemode:'tozero'}},yaxis2:{{domain:yDomains[1],anchor:'x2',range:[0,100]}},xaxis3:{{domain:domains[2],anchor:'y3',rangemode:'tozero'}},yaxis3:{{domain:yDomains[2],anchor:'x3',rangemode:'tozero'}},xaxis4:{{domain:domains[3],anchor:'y4',rangemode:'tozero'}},yaxis4:{{domain:yDomains[3],anchor:'x4',rangemode:'tozero'}}}};Plotly.relayout('ds',layout);}};downsampleLayout();window.addEventListener('resize',downsampleLayout);const n={json.dumps(names)},s={json.dumps(sample_names)},r={json.dumps(raw)},p={json.dumps(ratio)};let mode='raw';const hm=document.getElementById('hm');const draw=()=>{{const values=mode==='ratio'?p:r;Plotly.react(hm,[{{type:'heatmap',z:values.slice(0,100),x:s,y:n.slice(0,100),colorscale:[[0,'#ecfdf5'],[.5,'#14b8a6'],[1,'#0f766e']],colorbar:{{title:'Value',thickness:10,len:.55}},hovertemplate:mode==='ratio'?'%{{y}}<br>%{{x}}: %{{z:.2f}}%<extra></extra>':'%{{y}}<br>%{{x}}: %{{z:,}}<extra></extra>'}}],{{margin:{{l:140,r:50,t:76,b:50}},xaxis:{{side:'top'}},yaxis:{{autorange:'reversed'}},font:{{family:'Arial'}},paper_bgcolor:'#fff',plot_bgcolor:'#fff'}},{{responsive:true}})}};draw();document.querySelectorAll('[data-c]').forEach(button=>button.onclick=()=>{{mode=button.dataset.c;const ratioMode=mode==='ratio',values=ratioMode?p:r;document.getElementById('aggregate-header').textContent=ratioMode?'Average':'Total';document.getElementById('tbody').innerHTML=n.map((name,index)=>{{const row=values[index],aggregate=row.reduce((sum,value)=>sum+value,0)/(ratioMode?row.length:1);return '<tr><td>'+name+'</td>'+row.map(value=>'<td>'+(ratioMode?value.toFixed(2)+'%':value.toLocaleString())+'</td>').join('')+'<td>'+(ratioMode?aggregate.toFixed(2)+'%':aggregate.toLocaleString())+'</td></tr>'}}).join('');document.querySelectorAll('[data-c]').forEach(item=>item.classList.toggle('active',item===button));draw()}});const toggleBatchTheme=button=>{{const dark=document.body.classList.toggle('dark');button.textContent=dark?'☀':'☾';const surface=dark?'#0f1d30':'#ffffff',font=dark?'#cbd5e1':'#172033',grid=dark?'#24344d':'#e2e8f0';const theme={{'paper_bgcolor':surface,'plot_bgcolor':surface,'font.color':font,'font.family':'Arial','xaxis.gridcolor':grid,'yaxis.gridcolor':grid}};['ds','hm','pi-batch-molecules','pi-batch-proteins'].forEach(id=>{{if(document.getElementById(id))Plotly.relayout(id,theme)}})}};document.querySelector('.theme').onclick=function(){{toggleBatchTheme(this)}};</script></main></body></html>'''
    with atomic_text(hp) as h:h.write(page)
    return [xlsx,hp]


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
    if getattr(config, "pb_cb_enabled", False):
        stage, extra = "pair-map", [root / "04_matrix" / f"{sample}.pb_cb_summary.json"]
    else:
        stage = "downsample" if getattr(config, "downsample_enabled", True) else "matrix"
        extra = [root / "04_matrix" / f"{sample}.raw_matrix_summary.json"]
    if getattr(config, "downsample_enabled", True):
        selection = root / "05_downsample" / f"{sample}.selection.json"
        extra.append(selection if selection.exists() else selection.with_name(f"{sample}.optimal_saturation_point.tsv"))
    return stage, [
        correction_qc_path(root / "02_corrected", sample),
        *extra]


def _barcode2_annotation_gap(config, sample, stats, barcode2_segments):
    """Return whitelist reads, unannotated reads, and their barcode counts.

    Existing correction traces are scanned at most once per trace/annotation
    version.  The compact cache keeps subsequent summary refreshes cheap.
    """
    if getattr(config, "barcode2_sequence_only", False):
        by_scope = {row["scope"]: row for row in stats}
        valid_reads = int(by_scope.get(f"final_{config.target_name('barcode2')}", {}).get("valid_reads", 0) or 0)
        return valid_reads, 0, []
    root = config.output_dir / sample
    corrected = root / "02_corrected"
    qc_path = corrected / f"{sample}.barcode2_annotation_qc.tsv"
    detail_path = corrected / f"{sample}.barcode2_not_in_annotation.tsv.gz"
    current = corrected / f"{sample}.correction_qc.json"
    if current.exists():
        payload = json.loads(current.read_text())
        if payload.get("schema") != 1:
            raise ConfigError(f"Invalid correction QC: {current}")
        qc = payload["annotation"]
        return int(qc["whitelist_reads"]), int(qc["not_in_annotation_reads"]), list(open_tsv(detail_path))
    # New correction stages own these compact QC artifacts, even with trace
    # disabled. Their validity is covered by the correction manifest.
    if qc_path.is_file() and detail_path.is_file():
        qc = next(open_tsv(qc_path))
        return int(qc["whitelist_reads"]), int(qc["not_in_annotation_reads"]), list(open_tsv(detail_path))

    trace = corrected / f"{sample}.barcode_correction_trace.parquet"
    if not is_parquet_table(trace):
        trace = corrected / f"{sample}.barcode_correction_trace.tsv.gz"

    by_scope = {row["scope"]: row for row in stats}
    # A historical run without a trace can still report the rate for its usual
    # single barcode2 segment, but cannot recover the omitted sequences.
    if not trace.exists():
        if len(barcode2_segments) != 1:
            raise ConfigError("Legacy multi-segment annotation QC requires a correction trace; rerun correction")
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

    for row in open_records(trace):
        read_id = row["read_id"]
        if current_id is not None and read_id != current_id:
            consume(current); current = {}
        current_id = read_id
        if row["segment_name"] in segment_names and row["whitelist_hit"] == "true":
            current[row["segment_name"]] = row["corrected_sequence"]
    if current_id is not None:
        consume(current)
    details = [
        {"barcode2_sequence": sequence, "reads_count": count}
        for sequence, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    # Legacy compatibility reads do not mutate upstream completion artifacts.
    return whitelist_reads, sum(counts.values()), details


def _completed_summary_samples(config: TagForgeConfig):
    """Return samples in configuration order once their summary inputs exist."""
    completed = []
    for item in config.samples:
        stage, outputs = _summary_inputs(config, item.sample)
        memo = {}
        if stage_current(config, item.sample, stage, memo) and (
            not getattr(config, "downsample_enabled", True) or stage_current(config, item.sample, "downsample", memo)
        ) and all(path.exists() for path in outputs):
            completed.append(item.sample)
    return completed


def _downsample_table(config, samples):
    rows = [["Sample", "Status", "Selected ratio",
             "Supporting reads", "Molecules", "Singleton molecules",
             "Read-level saturation (%)", "Duplication (%)",
             "Non-singleton molecule fraction (%)", "Criterion"]]
    points = {}
    for sample in samples:
        enabled = getattr(config, "downsample_enabled", True)
        point = selection_point(config, sample) if enabled else {}
        points[sample] = point
        def number(key, cast=float):
            value = point.get(key)
            return cast(value) if value not in (None, "") else None
        rows.append([sample, "completed" if enabled else "not run",
                     number("optimal_downsample_ratio"),
                     number("reads_sampled", int), number("umi_types", int),
                     number("umi_detected_once", int), number("sequencing_saturation"),
                     number("duplication_ratio"), number("max_non_singleton_molecule_fraction"),
                     point.get("selection_criterion", "")])
    return rows, points


def _batch_workbook_sheets(config, samples):
    """One workbook contract for run, report and summary commands."""
    downsample, points = _downsample_table(config, samples)
    if getattr(config, "pb_cb_enabled", False):
        from itertools import islice
        sheets = [("Summary", _pair_mapping_rows(config, samples)), ("Downsample", downsample)]
        if getattr(config, "downsample_enabled", True):
            sheets.append(("Downsample Definitions", [["Analysis", SAMPLING_NOTE]]))
        for sample in samples:
            path = config.output_dir / sample / "04_matrix" / f"{sample}.pb_cb_map.tsv.gz"
            mapping = list(islice(open_tsv(path), 1000))
            sheets.append((sample, [list(mapping[0])] + [list(row.values()) for row in mapping]
                           if mapping else [[config.target_name("barcode1"), config.target_name("barcode2"),
                                             "molecule_count", "reads_count"]]))
        sheets.append(("Results", _result_rows(config, samples)))
        return sheets

    segments = [segment for segment in config.segments if segment.target in {"barcode1", "barcode2"}]
    barcode2_segments = [segment for segment in segments if segment.target == "barcode2"]
    names = _barcode2_names(barcode2_segments)
    pi_seq = getattr(config, "pi_seq_enabled", False)
    summary = [["Sample", "Reads", "Valid reads", "Valid",
                *[f"{segment.name} valid" for segment in segments],
                "Annotated", "Molecules",
                f"{config.target_name('barcode1')} count", "Features"] + (
                    ["MP", "Dominant MP", "Removed",
                     "MP groups", "Retained",
                     "PB >=2 mol", "PB >=2 feat"]
                    if pi_seq else [])]
    warnings = [["Sample", "Annotated rate (fraction)", "FB not in annotation", "Reads",
                 "Annotated below 90%"]]
    full_counts, selected_counts = {}, {}
    has_selected_counts = any(float(point.get("optimal_downsample_ratio", 1)) < 1
                              for point in points.values())
    for sample in samples:
        root = config.output_dir / sample
        stats = load_correction_stats(correction_qc_path(root / "02_corrected", sample))
        by_scope = {row["scope"]: row for row in stats}
        combined = by_scope.get("combined", {})
        annotated, details = _annotation_rate_and_details(config, sample, stats, barcode2_segments)
        matrix = root / "04_matrix" / f"{sample}.raw_matrix_summary.json"
        barcode_count, counts = _matrix_summary(matrix)
        full_counts[sample] = counts
        row = [sample, _number(combined.get("total_reads")), _number(combined.get("valid_reads")),
               float(combined.get("valid_rate", 0) or 0),
               *[float(by_scope.get(segment.name, {}).get("valid_rate", 0) or 0) for segment in segments],
               annotated, sum(counts.values()), barcode_count, len(counts)]
        if pi_seq:
            qc = next(open_tsv(root / "03_molecules" / f"{sample}.pi_seq_qc.tsv"), {})
            pb = _pi_seq_pb_summary(matrix)
            total = pb["total"]
            row += [*_pi_seq_ratios(qc), _number(qc.get("multi_pi_fb_umis")),
                    _number(qc.get("retained_molecules")),
                    (total - pb["molecule_bins"][0]) / total if total else 0.0,
                    (total - pb["protein_bins"][0]) / total if total else 0.0]
        summary.append(row)
        warnings.extend([sample, annotated,
                         names.get(detail["barcode2_sequence"], detail["barcode2_sequence"]),
                         int(detail["reads_count"]), "yes" if annotated < 0.9 else "no"]
                        for detail in details[:50])
        # Always validate enabled selections, including ratio 1 references.
        if getattr(config, "downsample_enabled", True):
            selected = _matrix_summary(selected_summary_path(config, sample))[1]
            if has_selected_counts:
                selected_counts[sample] = selected

    features = sorted(set().union(*(set(counts) for counts in full_counts.values())))
    def count_table(counts):
        return [["Sample"] + features] + [[sample] + [counts[sample].get(feature, 0) for feature in features]
                                          for sample in samples]
    sheets = [("Summary", summary), ("Counts_raw", count_table(full_counts))]
    if has_selected_counts:
        sheets.append(("Counts_final", count_table(selected_counts)))
    sheets += [("Downsample", downsample), ("Warnings", warnings)]
    sheets.append(("Results", _result_rows(config, samples)))
    return sheets


def _write_batch_workbook(config, samples):
    directory = batch_report_dir(config.output_dir)
    output = directory / "TagForge_batch_report.xlsx"
    write_xlsx(output, _batch_workbook_sheets(config, samples))
    # Retire duplicate generated workbooks only after successful publication.
    for name in ("00_summary.xlsx", "01_FB_not_in_anno.xlsx"):
        (directory / name).unlink(missing_ok=True)
    return output


def write_summary(config: TagForgeConfig):
    """Refresh the unified workbook from all current samples without recomputing."""
    with _summary_lock(config.output_dir / ".summary.lock"):
        samples = _completed_summary_samples(config)
        output = _write_batch_workbook(config, samples)
    return output, samples


def _point_rows(config, sample):
    point = selection_point(config, sample)
    return [list(point), list(point.values())]
