from __future__ import annotations

import html
import json
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from itertools import islice, repeat

from . import __version__
from .barcode_correct import (WhitelistCorrector, load_fb_annotation, load_whitelist,
                              evaluate_extracted_read, REJECTION_REASONS)
from .config import ConfigError, TagForgeConfig
from .external_tools import check_external_tools
from .extract import _extract_record, decode_method_payload, decode_segment_payload
from .fastq import paired_fastq
from .io_utils import atomic_text, sample_dirs, write_tsv
from .logging_utils import sample_logger
from .report_style import REPORT_CSS, DARK_CSS, THEME_JS
from .umi_correct import deduplicate_umis


def _quick_extract(record, segments, segment_columns=None):
    row, runtime = _extract_record(record, segments, segment_columns)
    sequence_stats = {
        "r1_length": len(record.r1_seq), "r2_length": len(record.r2_seq),
        "bases": len(record.r1_seq) + len(record.r2_seq),
        "n_bases": record.r1_seq.count("N") + record.r2_seq.count("N"),
    }
    return row, runtime, sequence_stats


def _rate(numerator, denominator):
    return numerator / denominator if denominator else 0.0


def _take_leading_records(records, size: int):
    selected = list(enumerate(islice(records, size), 1))
    return selected, len(selected)


def enforce_quick_test(summary):
    if summary["qc_failures"]:
        raise ConfigError("Configured quick-test threshold failed for " + summary["sample"]
                          + ": " + "; ".join(summary["qc_failures"]))


def quick_test_sample(config: TagForgeConfig, sample_name: str, max_reads: int):
    if max_reads < 1:
        raise ConfigError("quick-test --reads must be >= 1")
    versions = check_external_tools()
    sample = config.sample(sample_name)
    dirs = sample_dirs(config.output_dir, sample_name)
    stats_path = dirs["report"] / f"{sample_name}.quick_test.tsv"
    html_path = dirs["report"] / f"{sample_name}.quick_test.html"
    ids_path = dirs["report"] / f"{sample_name}.quick_test.sampled_read_ids.txt"
    read_qc_path = dirs["report"] / f"{sample_name}.quick_test.read_qc.tsv"
    log_path = dirs["logs"] / f"{sample_name}.quick_test.log"
    logger = sample_logger(f"{sample_name}.quick_test", log_path)
    logger.info(
        "quick_test_start\tversion=%s\tsampling=leading_reads\trequested_reads=%s\t"
        "workers=%s\tcutadapt=%s\tumi_tools=%s",
        __version__, max_reads, config.threads,
        versions.cutadapt, versions.umi_tools,
    )

    annotation = load_fb_annotation(config)
    correctors = {}
    for segment in config.segments:
        if segment.target == "umi":
            continue
        if segment.whitelist:
            correctors[segment.name] = WhitelistCorrector(
                load_whitelist(segment.whitelist), segment.correction, segment.length
            )
        elif segment.correction.enabled:
            raise ConfigError(f"Segment {segment.name}: correction enabled but whitelist is missing")
        else:
            correctors[segment.name] = WhitelistCorrector([], segment.correction, segment.length)

    segment_counts = {segment.name: Counter() for segment in config.segments}
    r1_lengths, r2_lengths = Counter(), Counter()
    totals = Counter()
    correction_counters = {s.name: Counter() for s in config.segments if s.target != "umi"}
    read_qc = []
    barcode1_top, feature_top = Counter(), Counter()
    umi_groups = defaultdict(Counter)
    segment_columns = {
        "barcode1": config.segment_column("barcode1"),
        "barcode2": config.segment_column("barcode2"),
        "umi": config.segment_column("umi"),
    }
    started = time.monotonic()
    sampling_started = time.monotonic()
    selected, reads_scanned = _take_leading_records(
        paired_fastq(sample.r1, sample.r2), max_reads
    )
    sampling_elapsed = time.monotonic() - sampling_started
    records = (record for _, record in selected)
    if config.quick_test_read_details:
        with atomic_text(ids_path) as handle:
            handle.write("read_index\tread_id\n")
            for index, record in selected:
                handle.write(f"{index}\t{record.read_id}\n")
    else:
        ids_path.unlink(missing_ok=True)

    if config.threads == 1:
        extracted = (
            _quick_extract(record, config.segments, segment_columns)
            for record in records
        )
        executor = None
        backend = "serial"
    else:
        executor = ThreadPoolExecutor(max_workers=config.threads)
        extracted = executor.map(
            _quick_extract, records, repeat(tuple(config.segments)),
            repeat(segment_columns),
        )
        backend = "thread"
    try:
        for row, runtime, seq_stats in extracted:
            totals["reads"] += 1
            totals["bases"] += seq_stats["bases"]
            totals["n_bases"] += seq_stats["n_bases"]
            r1_lengths[seq_stats["r1_length"]] += 1
            r2_lengths[seq_stats["r2_length"]] += 1
            segments_by_target = {
                target: [s for s in config.segments if s.target == target]
                for target in ("barcode1", "barcode2", "umi")
            }
            raw_by_target = {
                target: decode_segment_payload(row[segment_columns[target]], segments_by_target[target])
                for target in ("barcode1", "barcode2", "umi")
            }
            extraction_methods = decode_method_payload(
                row["methods"], config.segments
            )
            for segment in config.segments:
                count = segment_counts[segment.name]
                stat = runtime[segment.name]
                count["total"] += 1
                count["extraction_success"] += int(stat["success"])
                count[f"method_{stat['method']}"] += int(stat["success"])
                count["linker_attempted"] += int(segment.method in {"linker", "linker_fixed"})
                linker_success = stat["success"] and stat["method"] == "linker"
                count["linker_success"] += int(linker_success)
                if stat["linker_failure_reason"]:
                    count[f"failure_{stat['linker_failure_reason'].split(':', 1)[0]}"] += 1
                if stat["left_linker_matches"] > 1: count["multiple_left"] += 1
                if stat["right_linker_matches"] > 1: count["multiple_right"] += 1
                if stat["linker_candidate_pairs"] > 1: count["multiple_pairs"] += 1
                count["cutadapt_microseconds"] += round(stat["linker_elapsed_seconds"] * 1_000_000)
                if segment.target != "umi":
                    source_method = extraction_methods.get(segment.name, "unknown")
                    extra = (
                        segment.correction.max_shift
                        if source_method == "fixed" and segment.correction.enabled
                        and segment.correction.allow_shift else 0
                    )
                    anchor = min(extra, segment.start or 0) if source_method == "fixed" else 0
                    correction = correctors[segment.name].correct(
                        raw_by_target[segment.target].get(segment.name, ""), anchor
                    )
                    count["barcode_valid"] += int(correction.success)
                    count[f"correction_{correction.correction_type}"] += 1
            local_summary = Counter()
            valid_row, _ = evaluate_extracted_read(row, config, correctors, annotation,
                                                   correction_counters, local_summary)
            totals.update(local_summary)
            reason = next((r for r in REJECTION_REASONS if local_summary[f"failure_{r}"]), "")
            if config.quick_test_read_details:
                read_qc.append({"read_id": row["read_id"], "extraction_status": row["status"],
                            "eligible": str(valid_row is not None).lower(), "rejection_reason": reason})
            if valid_row is not None:
                barcode1 = valid_row[config.target_name("barcode1")]
                feature = valid_row[f"{config.target_name('barcode2')}_name"]
                umi = valid_row[config.target_name("umi")]
                barcode1_top[barcode1] += 1
                feature_top[feature] += 1
                umi_groups[(barcode1, feature)][umi] += 1
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)

    if not totals["reads"]:
        raise ValueError(f"No reads found for sample {sample_name}")
    molecules = 0
    for counts in umi_groups.values():
        assignments = deduplicate_umis(dict(counts), config.umi_method, config.umi_max_distance)
        molecules += len(set(assignments.values()))
    elapsed = time.monotonic() - started
    valid_reads = totals["combined_valid"]
    valid_rate = _rate(valid_reads, totals["reads"])
    n_base_rate = _rate(totals["n_bases"], totals["bases"])
    failures = []
    min_valid = config.quick_test_min_valid_rate
    max_n = config.quick_test_max_n_base_rate
    if min_valid is not None and valid_rate < min_valid:
        failures.append(f"valid_rate={valid_rate:.6f} < min_valid_rate={min_valid}")
    if max_n is not None and n_base_rate > max_n:
        failures.append(f"n_base_rate={n_base_rate:.6f} > max_n_base_rate={max_n}")
    qc_status = "failed" if failures else "passed" if min_valid is not None or max_n is not None else "diagnostic"
    metric_rows = [
        {"scope": "sample", "metric": "tagforge_version", "value": __version__},
        {"scope": "sample", "metric": "reads_examined", "value": totals["reads"]},
        {"scope": "sample", "metric": "reads_scanned", "value": reads_scanned},
        {"scope": "sample", "metric": "sampling_method", "value": "leading_reads"},
        {"scope": "sample", "metric": "qc_status", "value": qc_status},
        {"scope": "sample", "metric": "sampling_seconds", "value": f"{sampling_elapsed:.6f}"},
        {"scope": "sample", "metric": "scan_reads_per_second", "value": f"{reads_scanned / sampling_elapsed:.3f}" if sampling_elapsed else "0"},
        {"scope": "sample", "metric": "workers", "value": config.threads},
        {"scope": "sample", "metric": "parallel_backend", "value": backend},
        {"scope": "sample", "metric": "elapsed_seconds", "value": f"{elapsed:.6f}"},
        {"scope": "sample", "metric": "reads_per_second", "value": f"{totals['reads'] / elapsed:.3f}"},
        {"scope": "sample", "metric": "r1_length_min", "value": min(r1_lengths)},
        {"scope": "sample", "metric": "r1_length_max", "value": max(r1_lengths)},
        {"scope": "sample", "metric": "r1_length_mode", "value": r1_lengths.most_common(1)[0][0]},
        {"scope": "sample", "metric": "r2_length_min", "value": min(r2_lengths)},
        {"scope": "sample", "metric": "r2_length_max", "value": max(r2_lengths)},
        {"scope": "sample", "metric": "r2_length_mode", "value": r2_lengths.most_common(1)[0][0]},
        {"scope": "sample", "metric": "n_base_rate", "value": f"{_rate(totals['n_bases'], totals['bases']):.8f}"},
        {"scope": "sample", "metric": "combined_valid_reads", "value": valid_reads},
        {"scope": "sample", "metric": "combined_valid_rate", "value": f"{_rate(valid_reads, totals['reads']):.8f}"},
        {"scope": "sample", "metric": "deduplicated_molecules", "value": molecules},
        {"scope": "sample", "metric": "duplication_ratio", "value": f"{100 * _rate(valid_reads - molecules, valid_reads):.8f}"},
        {"scope": "sample", "metric": "duplication_ratio_unit", "value": "percent"},
        {"scope": "sample", "metric": "rate_unit", "value": "fraction"},
    ]
    metric_rows.extend({"scope": "rejection", "metric": reason, "value": totals[f"failure_{reason}"]}
                       for reason in REJECTION_REASONS)
    for segment in config.segments:
        count = segment_counts[segment.name]
        scope = f"segment:{segment.name}"
        metric_rows.extend([
            {"scope": scope, "metric": "extraction_success_rate", "value": f"{_rate(count['extraction_success'], count['total']):.8f}"},
            {"scope": scope, "metric": "linker_attempted", "value": count["linker_attempted"]},
            {"scope": scope, "metric": "linker_success", "value": count["linker_success"]},
            {"scope": scope, "metric": "linker_success_rate", "value": f"{_rate(count['linker_success'], count['linker_attempted']):.8f}" if count["linker_attempted"] else "NA"},
            {"scope": scope, "metric": "fixed_success", "value": count["method_fixed"]},
            {"scope": scope, "metric": "multiple_left_linker_reads", "value": count["multiple_left"]},
            {"scope": scope, "metric": "multiple_right_linker_reads", "value": count["multiple_right"]},
            {"scope": scope, "metric": "multiple_candidate_pair_reads", "value": count["multiple_pairs"]},
            {"scope": scope, "metric": "cutadapt_cpu_seconds", "value": f"{count['cutadapt_microseconds'] / 1_000_000:.6f}"},
        ])
        if segment.target != "umi":
            metric_rows.append({"scope": scope, "metric": "barcode_valid_rate", "value": f"{_rate(count['barcode_valid'], count['total']):.8f}"})
        for key in sorted(count):
            if key.startswith("failure_") or key.startswith("correction_"):
                metric_rows.append({"scope": scope, "metric": key, "value": count[key]})
    for rank, (barcode, count) in enumerate(barcode1_top.most_common(10), 1):
        metric_rows.append({"scope": f"top_{config.target_name('barcode1')}", "metric": f"{rank}:{barcode}", "value": count})
    for rank, (feature, count) in enumerate(feature_top.most_common(10), 1):
        metric_rows.append({"scope": "top_feature", "metric": f"{rank}:{feature}", "value": count})
    write_tsv(stats_path, ["scope", "metric", "value"], metric_rows)
    if config.quick_test_read_details:
        write_tsv(read_qc_path, ["read_id", "extraction_status", "eligible", "rejection_reason"], read_qc)
    else:
        read_qc_path.unlink(missing_ok=True)

    segment_rows = []
    for segment in config.segments:
        count = segment_counts[segment.name]
        linker_rate = (f"{_rate(count['linker_success'], count['linker_attempted']):.1%}"
                       if count['linker_attempted'] else "\\")
        final_rate = (f"{_rate(count['barcode_valid'], totals['reads']):.1%}"
                      if segment.target != "umi" else "\\")
        segment_rows.append(
            f"<tr><td>{html.escape(segment.target_name or segment.target)}</td>"
            f"<td>{html.escape(segment.name)}</td><td>{html.escape(segment.method)}</td>"
            f"<td>{_rate(count['extraction_success'], totals['reads']):.1%}</td>"
            f"<td>{linker_rate}</td><td>{final_rate}</td><td>{count['method_fixed']:,}</td></tr>"
        )
    feature_counts = feature_top.most_common()
    feature_total = sum(feature_top.values())
    feature_rows = "".join(
        f'<tr><td>{html.escape(name)}</td><td>{count:,}</td>'
        f'<td>{_rate(count, feature_total):.2%}</td></tr>'
        for name, count in feature_counts
    ) or '<tr><td colspan="3" class="sub">No eligible feature reads in this prefix.</td></tr>'
    top_features = feature_counts[:20]
    # Escape HTML delimiters as well as JSON quotes for annotation-derived names.
    feature_data = json.dumps({
        "names": [name for name, _ in top_features],
        "counts": [count for _, count in top_features],
    }).replace("<", "\\u003c")
    duplication = 100 * _rate(valid_reads - molecules, valid_reads)
    details = [
        ("QC status", qc_status.title()),
        ("Eligible reads", f"{valid_reads:,}"),
        ("N-base Ratio", f"{n_base_rate:.2%}"),
    ]
    if min_valid is not None:
        details.append(("Minimum valid ratio", f"{min_valid:.1%}"))
    if max_n is not None:
        details.append(("Maximum N-base ratio", f"{max_n:.2%}"))
    for read_name, lengths in (("R1", r1_lengths), ("R2", r2_lengths)):
        shortest, longest = min(lengths), max(lengths)
        length_summary = (
            f"{shortest} bp" if shortest == longest else
            f"{shortest}–{longest} bp (most common: {lengths.most_common(1)[0][0]} bp)"
        )
        details.append((f"{read_name} read length", length_summary))
    details.extend([
        ("Elapsed time", f"{elapsed:,.2f} s"),
        ("Speed", f"{totals['reads'] / elapsed:,.0f} reads/s"),
        ("TagForge version", __version__),
    ])
    detail_rows = "".join(f"<tr><td>{html.escape(label)}</td><td>{html.escape(value)}</td></tr>"
                          for label, value in details)
    warning_html = (
        '<section class="panel"><h2>Warnings</h2><ul>'
        + "".join(f"<li>{html.escape(failure)}</li>" for failure in failures)
        + '</ul></section>'
    ) if failures else ""
    feature_plot = (
        f"const featureData={feature_data};"
        "Plotly.newPlot('features',[{type:'bar',orientation:'h',x:featureData.counts,"
        "y:featureData.names,marker:{color:'#0f766e'},"
        "hovertemplate:'%{y}<br>%{x:,} eligible reads<extra></extra>'}],"
        "{autosize:true,yaxis:{autorange:'reversed',automargin:true},"
        "xaxis:{title:'Eligible reads',gridcolor:'#e2e8f0',rangemode:'tozero'},"
        "margin:{l:60,r:30,t:16,b:50},paper_bgcolor:'#ffffff',plot_bgcolor:'#ffffff',"
        "font:{family:'Arial',color:'#172033'}},{responsive:true});"
    ) if top_features else ""
    chart_content = ('<div id="features" style="height:520px"></div>' if top_features else
                     '<p class="sub">No eligible feature reads in this prefix.</p>')
    quick_css = (
        '.quick-report .cards{grid-template-columns:repeat(3,minmax(0,1fr))}'
        '.quick-report .report-grid{grid-template-columns:repeat(auto-fit,minmax(min(100%,480px),1fr))}'
        '.quick-report .table-wrap{overflow:auto}'
        '@media(max-width:600px){.quick-report{padding:20px}.quick-report h1{overflow-wrap:anywhere}'
        '.quick-report .cards .card{padding:12px}'
        '.quick-report .cards .value{font-size:clamp(18px,5vw,24px)}}'
    )
    page = f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>TagForge quick test · {html.escape(sample_name)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script><style>{REPORT_CSS}{DARK_CSS}{quick_css}</style></head><body><main class="quick-report">
<button class="theme icon-button" aria-label="Switch to dark theme" title="Switch to dark theme">☾</button><p class="sub">TAGFORGE / QUICK TEST</p><h1>{html.escape(sample_name)}</h1>
<p class="sub">Leading-read prefix diagnostic; not a random population sample. All metrics below describe the examined prefix.</p>
<div class="cards">
<div class="card"><div class="sub">Reads Examined</div><div class="value">{totals['reads']:,}</div></div>
<div class="card"><div class="sub">Features</div><div class="value">{len(feature_top):,}</div></div>
<div class="card"><div class="sub">Valid Ratio</div><div class="value">{valid_rate:.1%}</div></div>
<div class="card"><div class="sub">Molecules</div><div class="value">{molecules:,}</div></div>
<div class="card"><div class="sub">{html.escape(config.target_name('barcode1'))}</div><div class="value">{len(barcode1_top):,}</div></div>
<div class="card"><div class="sub">Duplication</div><div class="value">{duplication:.2f}%</div></div></div>
{warning_html}
<div class="report-grid"><section class="panel"><h2>Segment QC</h2><div class="table-wrap"><table><tr><th>Barcode</th><th>Segment</th><th>Mode</th><th>Extracted</th><th>Linker</th><th>Final</th><th>Fixed Reads</th></tr>{''.join(segment_rows)}</table></div></section>
<section class="panel"><h2>Quick Test Details</h2><table><tr><th>Metric</th><th>Value</th></tr>{detail_rows}</table></section></div>
<div class="report-grid"><section class="panel"><h2>Top 20 Feature</h2>{chart_content}</section>
<section class="panel"><h2>Feature Counts</h2><div class="scroll-window"><table><tr><th>Feature</th><th>Eligible Reads</th><th>Ratio</th></tr>{feature_rows}</table></div></section></div>
<script>{THEME_JS}document.querySelector('.theme').onclick=function(){{toggleTheme(this,{"['features']" if top_features else "[]"})}};{feature_plot}</script>
</main></body></html>'''
    with atomic_text(html_path) as handle:
        handle.write(page)
    for row in metric_rows:
        logger.info("quick_test_metric\tscope=%s\tmetric=%s\tvalue=%s", row["scope"], row["metric"], row["value"])
    logger.info(
        "quick_test_finish\treads_scanned=%s\treads_sampled=%s\tsampling_seconds=%.6f\t"
        "elapsed_seconds=%.6f\tstats=%s\thtml=%s\tread_ids=%s",
        reads_scanned, totals["reads"], sampling_elapsed, elapsed, stats_path, html_path, ids_path,
    )
    return {
        "sample": sample_name, "reads": totals["reads"], "valid_rate": valid_rate,
        "qc_status": qc_status, "qc_failures": failures, "read_qc": read_qc_path if config.quick_test_read_details else None,
        "molecules": molecules, "elapsed": elapsed, "reads_scanned": reads_scanned,
        "stats": stats_path, "html": html_path, "read_ids": ids_path if config.quick_test_read_details else None, "log": log_path,
    }
