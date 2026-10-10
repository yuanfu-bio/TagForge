from __future__ import annotations

import time
from pathlib import Path

from . import __version__
from .barcode_correct import correct_sample
from .config import TagForgeConfig
from .checkpoint import (expected_outputs, invalidate, pending_publication, prepare_publication,
                         prepare_stage, publish_stage, require_dependencies, sample_lock,
                         stage_current, stage_request, stage_plan, manifest_path, migrate_manifest)
from .downsample import downsample_sample
from .extract import extract_sample
from .external_tools import check_external_tools
from .io_utils import sample_dirs, write_progress
from .logging_utils import sample_logger
from .matrix import matrix_sample, pair_mapping_sample
from .reports import _completed_summary_samples, _summary_lock, batch_report, report_sample
from .umi_correct import dedup_sample


STEP_FUNCS = {"extract": extract_sample, "correct": correct_sample, "dedup": dedup_sample,
              "matrix": matrix_sample, "pair-map": pair_mapping_sample, "downsample": downsample_sample, "report": report_sample}


def _canonical_result(value, source_root, destination_root):
    if isinstance(value, Path):
        try:
            return destination_root / value.relative_to(source_root)
        except ValueError:
            return value
    if isinstance(value, list):
        return [_canonical_result(item, source_root, destination_root) for item in value]
    if isinstance(value, tuple):
        return tuple(_canonical_result(item, source_root, destination_root) for item in value)
    return value


def run_step(config: TagForgeConfig, sample: str, step: str, overwrite: bool = False):
    with sample_lock(config, sample):
        return _run_step_locked(config, sample, step, overwrite)


def _run_step_locked(config: TagForgeConfig, sample: str, step: str, overwrite: bool = False):
    if step == "downsample" and not config.downsample_enabled:
        from .config import ConfigError
        raise ConfigError("downsample.enabled is false; enable it to run this stage")
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
    checkpoint = manifest_path(config, sample, step)
    outputs = expected_outputs(config, sample, step)
    overwrite = overwrite or config.overwrite or not config.resume_enabled
    if not overwrite and stage_current(config, sample, step):
        migrate_manifest(config, sample, step)
        logger.info("%s\tskipped_checkpoint\telapsed_time=0", step)
        if step == "extract":
            preview = dirs["extracted"] / f"{sample}.extracted.preview.tsv"
            if config.extraction_preview_reads > 0 and not preview.exists():
                from itertools import islice
                from .columnar import open_records
                from .extract import extract_fields
                from .io_utils import write_tsv
                write_tsv(preview, extract_fields(config), islice(open_records(outputs[0]), config.extraction_preview_reads))
            elif not config.extraction_preview_reads:
                preview.unlink(missing_ok=True)
        write_progress(dirs["logs"] / f"{sample}.progress.json", sample, step, {"status": "skipped"})
        return outputs, "skipped"
    require_dependencies(config, sample, step)
    request = stage_request(config, sample, step)
    if checkpoint.exists() or (dirs["checkpoint"] / f"{step}.done").exists():
        logger.info("%s\tcheckpoint_invalidated\treason=inputs_config_or_output_changed", step)
    invalidate(config, sample, step)
    staged, work = prepare_stage(config, sample, step, request, overwrite=overwrite)
    write_progress(dirs["logs"] / f"{sample}.progress.json", sample, step, {"status": "running"})
    start = time.monotonic(); logger.info("%s\tstart\telapsed_time=0", step)
    try:
        publication = pending_publication(work)
        if publication:
            result = outputs
            logger.info("%s\trecovering_publication\tgeneration=%s", step, publication["generation"])
        else:
            result = (
                extract_sample(staged, sample, resume=config.resume_enabled, keep_resume=True)
                if step == "extract" else
                correct_sample(staged, sample, resume=config.resume_enabled, keep_resume=True)
                if step == "correct" else STEP_FUNCS[step](staged, sample)
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
        publication = publication or prepare_publication(config, staged, sample, step, request, work)
        publish_stage(config, staged, sample, step, request, work, publication)
        logger.info("%s\tfinish\telapsed_time=%.3f", step, time.monotonic() - start)
        write_progress(dirs["logs"] / f"{sample}.progress.json", sample, step,
                       {"status": "completed", "elapsed_seconds": round(time.monotonic() - start, 3)})
        outputs = expected_outputs(config, sample, step)
        result = _canonical_result(result, staged.output_dir / sample, config.output_dir / sample)
        return result, "completed"
    except Exception:
        write_progress(dirs["logs"] / f"{sample}.progress.json", sample, step,
                       {"status": "failed", "elapsed_seconds": round(time.monotonic() - start, 3)})
        logger.exception("%s\tfailed\telapsed_time=%.3f", step, time.monotonic() - start)
        raise


def run_pipeline(config: TagForgeConfig, sample_names, overwrite: bool = False):
    versions = check_external_tools()
    for sample in sample_names:
        with sample_lock(config, sample):
            dirs = sample_dirs(config.output_dir, sample)
            if not config.quick_test_read_details:
                for suffix in ("quick_test.read_qc.tsv", "quick_test.sampled_read_ids.txt"):
                    (dirs["report"] / f"{sample}.{suffix}").unlink(missing_ok=True)
            logger = sample_logger(sample, dirs["logs"] / f"{sample}.pipeline.log")
            logger.info(
                "dependencies\tcutadapt=%s\tumi_tools=%s\telapsed_time=0",
                versions.cutadapt, versions.umi_tools,
            )
            steps = stage_plan(config)
            for step in steps:
                _run_step_locked(config, sample, step, overwrite)
    # Publish batch workbooks once per invocation, including one-sample array
    # tasks, rather than re-reading the entire batch after every sample.
    with _summary_lock(config.output_dir / ".summary.lock"):
        batch_report(config, _completed_summary_samples(config))
