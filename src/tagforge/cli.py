from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import ConfigError, load_config
from .external_tools import check_external_tools
from .pipeline import run_pipeline, run_step
from .sparse_matrix import export_matrix, pair_path
from .quick_test import quick_test_sample, enforce_quick_test
from .reports import _completed_summary_samples, _summary_lock, batch_report, write_summary
from .slurm import make_slurm


EXAMPLE_CONFIG = """# TagForge coordinates are 0-based; fixed segments use [start, start+length).
project:
  name: tagforge_example
  workdir: ..
  output_dir: 02_output
samples:
  auto:
    raw_dir: tests/fixtures/small_fastq/raw
    r1: "{sample}_raw_1.fq.gz"
    r2: "{sample}_raw_2.fq.gz"
barcode2_annotation:
  fb_info: tests/fixtures/small_fastq/config/FB_info.tsv
  id_column: FB_ID
  sequence_column: sequence
  name_column: antibody_name
correction_barcode:
  enabled: true
  max_shift: 1
  max_mismatch: 1
  allow_shift: true
  allow_mismatch: true
linker:
  max_mismatch: 0
correction_umi:
  method: directional
  max_distance: 1
# barcode1/barcode2/umi are output roles; name is the custom target name shown
# in extraction stats, correction trace, and QC logs.
barcode1:
  name: CELL
  segments:
    # `methods: [linker, fixed]` means linker-first with fixed fallback. The
    # fixed start is always relative to the original read.
    - segment: CELL
      read: R1
      methods: [linker, fixed]
      right_linker: AGGTC
      start: 0
      length: 8
      whitelist: tests/fixtures/small_fastq/config/CELL_WL.txt
barcode2:
  name: FB
  segments:
    - segment: FB
      read: R2
      method: fixed
      start: 6
      length: 6
      whitelist: tests/fixtures/small_fastq/config/FB_WL.txt
      correction:
        max_shift: 0
        allow_shift: false
umi:
  name: UMI
  segments:
    - segment: UMI
      read: R2
      method: fixed
      start: 0
      length: 6
      correction:
        enabled: false
downsample:
  enabled: true
  ratios: auto
  random_seed: 12345
quick_test:
  enabled: true
  reads: 10000
output:
  quick_test_read_details: false
performance:
  threads: 2
  # Barcode whitelist correction is parallelized by extracted-read chunks.
  # Omit barcode_workers to use threads.
  # barcode_workers: 2
  # Downsample ratio jobs use forked processes; omit to use threads.
  # downsample_workers: 2
  # UMI aggregation uses GNU sort by default; workers control sort parallelism.
  # umi_aggregation_backend: external_sort  # or sqlite for legacy compatibility
  # umi_aggregation_workers: 2
  # umi_sort_memory_mb: 512
  # UMI-tools groups are parallelized with processes. Omit umi_workers to use threads.
  # umi_workers: 2
  umi_batch_size: 5000
  # Connected components are not split to fit this batch target.
  # Optional budgets report the barcode/feature exceeding the limit.
  # umi_max_group_umis: 1000000
  # umi_max_component_umis: 100000
  umi_sqlite_cache_mb: 64
  chunk_size: 10000
  extraction_preview_reads: 0
  parquet_compression: snappy  # snappy (default) or zstd
  compression_level: 3  # gzip, and Parquet only when using zstd
resume:
  enabled: true
  overwrite: false
"""


def _common(parser):
    parser.add_argument("--config", required=True, help="YAML configuration file")
    parser.add_argument("--sample", action="append", help="Sample name (repeatable; default: all)")
    parser.add_argument(
        "--threads", type=int, default=None,
        help="Process workers for Cutadapt extraction and, unless overridden, UMI deduplication",
    )
    parser.add_argument("--overwrite", action="store_true", help="Ignore checkpoints and replace outputs")


def parser():
    root = argparse.ArgumentParser(prog="tagforge", description="Forge paired-end barcode libraries into count matrices")
    root.add_argument("--version", action="version", version=f"TagForge {__version__}")
    commands = root.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="Run quick-test followed by the complete pipeline"); _common(run)
    run.add_argument(
        "--skip-quick-test", "--no-quick-test", action="store_true",
        help="Skip the default small-subset QC before the formal pipeline",
    )
    for step in ("extract", "correct", "dedup", "matrix", "pair-map", "downsample", "report"):
        p = commands.add_parser(step, help=f"Run the {step} step"); _common(p)
    export = commands.add_parser("export-matrix", help="Export a sparse matrix on demand")
    export.add_argument("--config", required=True)
    export.add_argument("--sample", required=True)
    export.add_argument("--format", choices=("parquet", "mtx", "tsv"), default="parquet")
    export.add_argument("--out", required=True)
    migrate = commands.add_parser("migrate-parquet", help="Build sparse Parquet matrix artifacts from existing molecule detail")
    _common(migrate)
    summary = commands.add_parser("summary", help="Refresh the unified batch workbook from completed samples")
    summary.add_argument("--config", required=True, help="YAML configuration file")
    quick = commands.add_parser("quick-test", help="Inspect a small read subset without running the full pipeline")
    quick.add_argument("--config", required=True)
    quick.add_argument("--sample", action="append", help="Sample name (repeatable; default: all)")
    quick.add_argument("--reads", type=int, default=None, help="Paired reads to inspect per sample")
    quick.add_argument("--threads", type=int, default=None, help="Cutadapt quick-test workers")
    validate = commands.add_parser("validate-config", help="Validate configuration and inputs")
    validate.add_argument("--config", required=True)
    init = commands.add_parser("init-config", help="Write an example configuration")
    init.add_argument("--out", required=True); init.add_argument("--force", action="store_true")
    example = commands.add_parser("init-example", help="Generate a standalone deterministic synthetic project")
    example.add_argument("--out", required=True, help="Project directory")
    example.add_argument("--force", action="store_true")
    slurm = commands.add_parser("make-slurm", help="Generate Slurm scripts for one or many samples")
    slurm.add_argument("--config", required=True); slurm.add_argument("--out", required=True)
    slurm.add_argument("--threads", type=int, default=8); slurm.add_argument("--mem", default="16G"); slurm.add_argument("--time", default="24:00:00")
    slurm.add_argument("--partition"); slurm.add_argument("--account"); slurm.add_argument("--qos")
    slurm.add_argument("--constraint"); slurm.add_argument("--gres")
    slurm.add_argument("--nodes", type=int, default=1); slurm.add_argument("--ntasks", type=int, default=1)
    slurm.add_argument("--mail-user"); slurm.add_argument("--mail-type")
    slurm.add_argument("--conda-env", default="tagforge")
    slurm.add_argument("--skip-quick-test", action="store_true", help="Disable quick-test in generated jobs")
    slurm.add_argument("--extra-sbatch", action="append", default=[], help="Additional raw #SBATCH option; repeatable")
    slurm.add_argument("--mode", choices=("array", "per-sample"), default="array", help="Slurm output style")
    slurm.add_argument("--array-limit", type=int, default=None, help="Limit concurrent Slurm array tasks")
    return root


def _selected(config, requested):
    names = requested or [s.sample for s in config.samples]
    for name in names: config.sample(name)
    return names


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == "init-example":
            from .example_data import generate_example
            out = generate_example(Path(args.out), EXAMPLE_CONFIG, force=args.force)
            print(f"Wrote synthetic example: {out}")
            return 0
        if args.command == "init-config":
            out = Path(args.out)
            if out.exists() and not args.force: raise ConfigError(f"Refusing to overwrite existing file: {out} (use --force)")
            out.parent.mkdir(parents=True, exist_ok=True); out.write_text(EXAMPLE_CONFIG, encoding="utf-8")
            print(f"Wrote example configuration: {out}"); return 0
        config = load_config(args.config, check_files=True)
        if args.command == "validate-config":
            versions = check_external_tools()
            print(
                f"Configuration is valid: {config.path}\nSamples: {len(config.samples)}\n"
                f"Segments: {len(config.segments)}\ncutadapt: {versions.cutadapt}\n"
                f"umi_tools: {versions.umi_tools}"
            )
            return 0
        if args.command == "make-slurm":
            files = make_slurm(
                config, Path(args.out).resolve(), args.threads, args.mem, args.time,
                partition=args.partition, account=args.account, qos=args.qos,
                constraint=args.constraint, gres=args.gres, nodes=args.nodes,
                ntasks=args.ntasks, mail_user=args.mail_user, mail_type=args.mail_type,
                conda_env=args.conda_env, skip_quick_test=args.skip_quick_test,
                extra_sbatch=args.extra_sbatch, mode=args.mode,
                array_limit=args.array_limit,
            )
            if args.mode == "array":
                print(f"Generated Slurm array script: {files[-1]} with samples table {files[0]}")
            else:
                print(f"Generated {len(files)-1} Slurm jobs and {files[-1]}")
            return 0
        thread_override = getattr(args, "threads", None)
        if thread_override is not None:
            if thread_override < 1:
                raise ConfigError("--threads must be >= 1")
            config.threads = thread_override
        if args.command == "export-matrix":
            sample = config.sample(args.sample)
            del sample
            from .checkpoint import stage_current
            if not stage_current(config, args.sample, "matrix"):
                raise ConfigError(
                    f"Export requires a current matrix stage for {args.sample}; "
                    "rerun affected upstream stages and tagforge matrix first"
                )
            output = Path(args.out)
            export_matrix(pair_path(config.output_dir / args.sample / "04_matrix", args.sample), output, args.format, config.target_name("barcode1"), threads=config.threads, memory_mb=config.matrix_memory_mb, compression=config.parquet_compression, compression_level=config.compression_level)
            print(f"Exported {args.format} matrix: {output}")
            return 0
        if args.command == "migrate-parquet":
            for name in _selected(config, args.sample):
                run_step(config, name, "matrix", overwrite=True)
            print(f"Sparse matrix migration completed for: {', '.join(_selected(config, args.sample))}")
            return 0
        if args.command == "summary":
            output, samples = write_summary(config)
            print(f"TagForge summary completed for {len(samples)} samples: {output}")
            return 0
        names = _selected(config, args.sample)
        if args.command == "run":
            if config.quick_test_enabled and not args.skip_quick_test:
                for name in names:
                    summary = quick_test_sample(config, name, config.quick_test_reads)
                    print(
                        f"{name} quick-test: sampled={summary['reads']:,} reads, "
                        f"scanned={summary['reads_scanned']:,}, "
                        f"valid={summary['valid_rate']:.2%}, molecules={summary['molecules']:,}, "
                        f"elapsed={summary['elapsed']:.2f}s"
                    )
                    enforce_quick_test(summary)
            run_pipeline(config, names, args.overwrite)
        elif args.command == "quick-test":
            check_external_tools()
            read_limit = args.reads if args.reads is not None else config.quick_test_reads
            if read_limit < 1:
                raise ConfigError("--reads must be >= 1")
            for name in names:
                summary = quick_test_sample(config, name, read_limit)
                print(
                    f"{name}: sampled={summary['reads']:,} reads, scanned={summary['reads_scanned']:,}, "
                    f"valid={summary['valid_rate']:.2%}, "
                    f"molecules={summary['molecules']:,}, elapsed={summary['elapsed']:.2f}s\n"
                    f"  TSV: {summary['stats']}\n  HTML: {summary['html']}\n"
                    f"  Read IDs: {summary['read_ids'] or 'disabled'}\n  Log: {summary['log']}"
                )
                enforce_quick_test(summary)
        else:
            if args.command in {"extract", "dedup"}:
                check_external_tools()
            for name in names:
                if args.command == "report" and config.downsample_enabled:
                    run_step(config, name, "downsample", args.overwrite)
                run_step(config, name, args.command, args.overwrite)
                if args.command == "downsample":
                    write_summary(config)
            if args.command == "report":
                # Rescan every completed sample under the shared lock instead
                # of only the ones this invocation's --sample selected, so a
                # per-sample `tagforge report --sample X` (or a SLURM array
                # task) never clobbers the batch report with a single sample.
                with _summary_lock(config.output_dir / ".summary.lock"):
                    batch_report(config, _completed_summary_samples(config))
        print(f"TagForge {args.command} completed for: {', '.join(names)}")
        return 0
    except (ConfigError, FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"tagforge: error: {exc}", file=sys.stderr); return 2


if __name__ == "__main__":
    raise SystemExit(main())
