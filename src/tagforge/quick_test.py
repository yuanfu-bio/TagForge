from __future__ import annotations

import html
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from itertools import islice, repeat

from . import __version__
from .barcode_correct import WhitelistCorrector, load_fb_annotation, load_whitelist
from .config import ConfigError, TagForgeConfig
from .external_tools import check_external_tools
from .extract import _extract_record, decode_method_payload, decode_segment_payload
from .fastq import paired_fastq
from .io_utils import atomic_text, sample_dirs, write_tsv
from .logging_utils import sample_logger
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


def quick_test_sample(config: TagForgeConfig, sample_name: str, max_reads: int):
    if max_reads < 1:
        raise ConfigError("quick-test --reads must be >= 1")
    versions = check_external_tools()
    sample = config.sample(sample_name)
    dirs = sample_dirs(config.output_dir, sample_name)
    stats_path = dirs["report"] / f"{sample_name}.quick_test.tsv"
    html_path = dirs["report"] / f"{sample_name}.quick_test.html"
    ids_path = dirs["report"] / f"{sample_name}.quick_test.sampled_read_ids.txt"
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
    with atomic_text(ids_path) as handle:
        handle.write("read_index\tread_id\n")
        for index, record in selected:
            handle.write(f"{index}\t{record.read_id}\n")
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
            corrected_by_target = {"barcode1": [], "barcode2": []}
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
                    corrected_by_target[segment.target].append(correction.corrected_sequence)
            barcode1 = "".join(corrected_by_target["barcode1"])
            barcode2 = "".join(corrected_by_target["barcode2"])
            umi = "".join(raw_by_target["umi"].get(s.name, "") for s in config.segments if s.target == "umi")
            feature = barcode2 if getattr(config, "barcode2_sequence_only", False) else annotation.get(barcode2, "")
            valid = bool(barcode1 and feature and umi and "N" not in umi)
            totals["combined_valid"] += int(valid)
            if valid:
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
    metric_rows = [
        {"scope": "sample", "metric": "tagforge_version", "value": __version__},
        {"scope": "sample", "metric": "reads_examined", "value": totals["reads"]},
        {"scope": "sample", "metric": "reads_scanned", "value": reads_scanned},
        {"scope": "sample", "metric": "sampling_method", "value": "leading_reads"},
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
        {"scope": "sample", "metric": "duplication_ratio", "value": f"{_rate(valid_reads - molecules, valid_reads):.8f}"},
    ]
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

    segment_html = "".join(
        f"<tr><td>{html.escape(segment.name)}</td><td>{html.escape(segment.method)}</td>"
        f"<td>{_rate(segment_counts[segment.name]['extraction_success'], totals['reads']):.1%}</td>"
        f"<td>{_rate(segment_counts[segment.name]['linker_success'], segment_counts[segment.name]['linker_attempted']):.1%}</td>"
        f"<td>{segment_counts[segment.name]['method_fixed']}</td></tr>"
        for segment in config.segments
    )
    page = f'''<!doctype html><html><head><meta charset="utf-8"><title>TagForge quick test · {html.escape(sample_name)}</title><style>
body{{margin:0;background:#07111f;color:#e5edf8;font:15px system-ui}}main{{max-width:1100px;margin:auto;padding:36px}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px}}.card,.panel{{background:#0f1d30;border:1px solid #24344d;border-radius:14px;padding:18px}}.value{{font-size:26px;color:#5eead4;font-weight:700}}.panel{{margin-top:18px}}table{{width:100%;border-collapse:collapse}}th,td{{padding:9px;text-align:left;border-bottom:1px solid #24344d}}</style></head><body><main>
<p>TAGFORGE {__version__} / QUICK TEST</p><h1>{html.escape(sample_name)}</h1><div class="cards">
<div class="card"><span>Reads</span><div class="value">{totals['reads']:,}</div></div><div class="card"><span>Valid</span><div class="value">{_rate(valid_reads, totals['reads']):.1%}</div></div><div class="card"><span>Molecules</span><div class="value">{molecules:,}</div></div><div class="card"><span>Speed</span><div class="value">{totals['reads']/elapsed:,.0f}/s</div></div></div>
<section class="panel"><h2>Segment QC</h2><table><tr><th>Segment</th><th>Mode</th><th>Extracted</th><th>Linker success</th><th>Fixed reads</th></tr>{segment_html}</table></section>
<section class="panel"><h2>Top features</h2><table>{''.join(f'<tr><td>{html.escape(k)}</td><td>{v}</td></tr>' for k,v in feature_top.most_common(10))}</table></section>
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
        "sample": sample_name, "reads": totals["reads"], "valid_rate": _rate(valid_reads, totals["reads"]),
        "molecules": molecules, "elapsed": elapsed, "reads_scanned": reads_scanned,
        "stats": stats_path, "html": html_path, "read_ids": ids_path, "log": log_path,
    }
