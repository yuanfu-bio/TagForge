from __future__ import annotations

import time
from pathlib import Path

from . import __version__
from .barcode_correct import correct_sample
from .config import TagForgeConfig
from .downsample import downsample_sample
from .extract import extract_sample
from .external_tools import check_external_tools
from .io_utils import sample_dirs, step_complete, touch_checkpoint
from .logging_utils import sample_logger
from .matrix import matrix_sample, pair_mapping_sample
from .reports import batch_report, report_sample, write_summary
from .umi_correct import dedup_sample


STEP_FUNCS = {"extract": extract_sample, "correct": correct_sample, "dedup": dedup_sample,
              "matrix": matrix_sample, "pair-map": pair_mapping_sample, "downsample": downsample_sample, "report": report_sample}


def expected_outputs(config: TagForgeConfig, sample: str, step: str):
    d = sample_dirs(config.output_dir, sample)
    return {
        "extract": [
            d["extracted"] / f"{sample}.extracted.tsv.gz",
            d["extracted"] / f"{sample}.extraction_stats.tsv",
        ],
        "correct": [d["detail"] / f"{sample}.valid_reads.tsv.gz", d["corrected"] / f"{sample}.barcode_correction_stats.tsv"] + ([d["corrected"] / f"{sample}.barcode_correction_trace.tsv.gz"] if config.trace_enabled else []),
        "dedup": [d["detail"] / f"{sample}.molecule_detail.tsv.gz"] + (
            [d["detail"] / f"{sample}.molecule_detail.rmMP.tsv.gz", d["detail"] / f"{sample}.pi_seq_qc.tsv"] if config.pi_seq_enabled else []),
        "matrix": [d["matrix"] / f"{sample}.raw_count_matrix.tsv.gz"],
        "pair-map": [d["matrix"] / f"{sample}.pb_cb_mapping.tsv.gz"] + (
            [d["matrix"] / f"{sample}.pb_cb_map.tsv.gz",
             d["matrix"] / f"{sample}.cb_pb_counts.tsv.gz",
             d["matrix"] / f"{sample}.cb_pb_count_distribution.tsv",
             d["matrix"] / f"{sample}.cb_observed_correction.tsv.gz"]
            if getattr(config, "pb_cb_enabled", False) else []),
        "downsample": [
            d["downsample"] / f"{sample}.downsample_metrics.tsv",
            d["downsample"] / f"{sample}.optimal_saturation_point.tsv",
            d["downsample"] / f"{sample}.downsample.html",
            d["detail"] / f"{sample}.optimal_saturation_molecule_detail.tsv.gz",
            d["matrix"] / f"{sample}.optimal_saturation_count_matrix.tsv.gz",
        ],
        "report": [d["report"] / f"{sample}.report.xlsx", d["report"] / f"{sample}.report.html"],
    }[step]


def run_step(config: TagForgeConfig, sample: str, step: str, overwrite: bool = False):
    dirs = sample_dirs(config.output_dir, sample)
    logger = sample_logger(sample, dirs["logs"] / f"{sample}.pipeline.log")
    logger.info("tagforge\tversion=%s\tstep=%s\telapsed_time=0", __version__, step)
    if step in {"extract", "correct", "dedup"}:
        versions = check_external_tools()
        backend = (
            f"cutadapt-python-api:{versions.cutadapt},workers={config.threads}"
            if step == "extract"
            else f"tagforge-barcode-correct,workers={config.barcode_workers or config.threads}"
            if step == "correct"
            else f"umi_tools-UMIClusterer:{versions.umi_tools},workers={config.umi_workers or config.threads}"
        )
        logger.info("%s\tbackend=%s\telapsed_time=0", step, backend)
    checkpoint = dirs["checkpoint"] / f"{step}.done"
    outputs = expected_outputs(config, sample, step)
    overwrite = overwrite or config.overwrite
    if step_complete(checkpoint, outputs, overwrite, __version__):
        logger.info("%s\tskipped_checkpoint\telapsed_time=0", step)
        return outputs, "skipped"
    start = time.monotonic(); logger.info("%s\tstart\telapsed_time=0", step)
    try:
        result = (
            extract_sample(config, sample, resume=not overwrite)
            if step == "extract" else
            correct_sample(config, sample, resume=not overwrite)
            if step == "correct" else STEP_FUNCS[step](config, sample)
        )
        if step == "extract" and isinstance(result, tuple) and isinstance(result[1], dict):
            summary = result[1]
            logger.info(
                "extract_summary\ttotal_reads=%s\textracted_reads=%s\tparallel_backend=%s\tworkers=%s\t"
                "wall_seconds=%.6f\tcutadapt_cpu_seconds=%.6f\treads_per_second=%.3f\telapsed_time=%.3f",
                summary.get("total_reads", 0), summary.get("extracted_reads", 0),
                summary.get("parallel_backend", "unknown"), config.threads,
                summary.get("wall_seconds", 0.0), summary.get("cutadapt_cpu_seconds", 0.0),
                summary.get("reads_per_second", 0.0), time.monotonic() - start,
            )
            for row in summary.get("segment_stats", []):
                linker_rate = row["linker_success_rate"]
                rescue_rate = row["fixed_rescue_rate"]
                logger.info(
                    "extract_segment\tsegment=%s\ttarget=%s\tread=%s\tmode=%s\t"
                    "linker_attempted=%s\tlinker_success=%s\tlinker_failed=%s\tlinker_success_rate=%s\t"
                    "left_linker_not_found=%s\tright_linker_not_found=%s\tunexpected_length=%s\t"
                    "multi_left_reads=%s\tmulti_right_reads=%s\tmulti_pair_reads=%s\t"
                    "candidate_pairs_total=%s\tselected_gap_min=%s\tselected_gap_max=%s\t"
                    "fixed_attempted=%s\tfixed_success=%s\tfixed_failed=%s\tfixed_rescue_rate=%s\t"
                    "final_success=%s\tfinal_success_rate=%.6f\tcutadapt_cpu_seconds=%s\telapsed_time=%.3f",
                    row["segment"], row["target"], row["read"], row["configured_mode"],
                    row["linker_attempted"], row["linker_success"], row["linker_failed"],
                    f"{float(linker_rate):.6f}" if linker_rate != "" else "NA",
                    row["left_linker_not_found"], row["right_linker_not_found"], row["unexpected_length"],
                    row["reads_with_multiple_left_matches"], row["reads_with_multiple_right_matches"],
                    row["reads_with_multiple_candidate_pairs"], row["linker_candidate_pairs_total"],
                    row["selected_gap_min"] if row["selected_gap_min"] != "" else "NA",
                    row["selected_gap_max"] if row["selected_gap_max"] != "" else "NA",
                    row["fixed_attempted"], row["fixed_success"], row["fixed_failed"],
                    f"{float(rescue_rate):.6f}" if rescue_rate != "" else "NA",
                    row["final_success"], float(row["final_success_rate"]), row["cutadapt_cpu_seconds"],
                    time.monotonic() - start,
                )
        if step == "correct" and isinstance(result, tuple) and isinstance(result[1], dict):
            summary = result[1]
            logger.info(
                "correction_summary\ttotal_reads=%s\textracted_reads=%s\t"
                "barcode1_valid=%s\tbarcode2_valid=%s\tcombined_valid=%s\t"
                "reads_per_second=%.3f\trequested_workers=%s\tworkers=%s\t"
                "chunk_size=%s\twall_seconds=%.3f\telapsed_time=%.3f",
                summary.get("total_reads", 0), summary.get("extracted_reads", 0),
                summary.get("barcode1_valid", 0), summary.get("barcode2_valid", 0),
                summary.get("combined_valid", 0), summary.get("reads_per_second", 0.0),
                summary.get("requested_workers", config.barcode_workers or config.threads),
                summary.get("workers", config.barcode_workers or config.threads),
                summary.get("chunk_size", config.chunk_size),
                summary.get("wall_seconds", 0.0), time.monotonic() - start,
            )
            for row in result[1].get("segment_method_qc", []):
                linker_rate = row["linker_barcode_valid_rate"]
                fixed_rate = row["fixed_barcode_valid_rate"]
                logger.info(
                    "barcode_method_qc\tsegment=%s\ttarget=%s\t"
                    "linker_extracted=%s\tlinker_barcode_valid=%s\tlinker_barcode_valid_rate=%s\t"
                    "fixed_extracted=%s\tfixed_barcode_valid=%s\tfixed_barcode_valid_rate=%s\t"
                    "elapsed_time=%.3f",
                    row["scope"], row["target_type"], row["linker_extracted_count"],
                    row["linker_barcode_valid_count"],
                    f"{float(linker_rate):.6f}" if linker_rate != "" else "NA",
                    row["fixed_extracted_count"], row["fixed_barcode_valid_count"],
                    f"{float(fixed_rate):.6f}" if fixed_rate != "" else "NA",
                    time.monotonic() - start,
                )
        if step == "dedup" and isinstance(result, tuple) and isinstance(result[1], dict):
            summary = result[1]
            logger.info(
                "dedup_summary\tvalid_reads=%s\tgroups=%s\ttotal_groups=%s\t"
                "raw_umis=%s\ttotal_raw_umis=%s\tmolecules=%s\t"
                "duplicates=%s\trequested_workers=%s\tworkers=%s\t"
                "requested_aggregation_workers=%s\taggregation_workers=%s\t"
                "umi_batch_size=%s\tbatches_submitted=%s\tbatches_completed=%s\t"
                "peak_batch_umis=%s\t"
                "aggregation_backend=%s\tsort_memory_mb=%s\tsqlite_cache_mb=%s\t"
                "aggregation_seconds=%.3f\tclustering_seconds=%.3f\t"
                "wall_seconds=%.3f\telapsed_time=%.3f",
                summary["valid_reads"], summary["groups"], summary["total_groups"],
                summary["raw_umis"], summary["total_raw_umis"], summary["molecules"],
                summary["duplicates"], summary["requested_workers"], summary["workers"],
                summary["requested_aggregation_workers"], summary["aggregation_workers"],
                summary["umi_batch_size"], summary["batches_submitted"],
                summary["batches_completed"], summary["peak_batch_umis"],
                summary.get("aggregation_backend", "sqlite"),
                summary.get("sort_memory_mb", "NA"),
                summary.get("sqlite_cache_mb", "NA"), summary["aggregation_seconds"],
                summary["clustering_seconds"], summary["wall_seconds"], time.monotonic() - start,
            )
        if not all(path.is_file() and path.stat().st_size > 0 for path in outputs):
            raise RuntimeError(f"Step {step} did not create all expected outputs")
        touch_checkpoint(checkpoint, __version__)
        logger.info("%s\tfinish\telapsed_time=%.3f", step, time.monotonic() - start)
        return result, "completed"
    except Exception:
        logger.exception("%s\tfailed\telapsed_time=%.3f", step, time.monotonic() - start)
        raise


def run_pipeline(config: TagForgeConfig, sample_names, overwrite: bool = False):
    versions = check_external_tools()
    for sample in sample_names:
        dirs = sample_dirs(config.output_dir, sample)
        logger = sample_logger(sample, dirs["logs"] / f"{sample}.pipeline.log")
        logger.info(
            "dependencies\tcutadapt=%s\tumi_tools=%s\telapsed_time=0",
            versions.cutadapt, versions.umi_tools,
        )
        steps = ("extract", "correct", "dedup", "pair-map") if getattr(config, "pb_cb_enabled", False) else (
            "extract", "correct", "dedup", "matrix", "downsample", "report")
        for step in steps:
            run_step(config, sample, step, overwrite)
            # Array tasks publish completed samples immediately; the writer
            # rescans all configured checkpoints under an interprocess lock.
            if step == "downsample":
                write_summary(config)
    if not getattr(config, "pb_cb_enabled", False):
        batch_report(config, sample_names)
