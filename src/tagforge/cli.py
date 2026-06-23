from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import ConfigError, load_config
from .external_tools import check_external_tools
from .pipeline import run_pipeline, run_step
from .quick_test import quick_test_sample
from .reports import batch_report
from .slurm import make_slurm


EXAMPLE_CONFIG = """# TagForge coordinates are 0-based; fixed segments use [start, start+length).
project:
  name: tagforge_example
  workdir: ..
  output_dir: 02_output
samples:
  auto:
    raw_dir: examples/small_fastq
    r1: "{sample}_raw_1.fq.gz"
    r2: "{sample}_raw_2.fq.gz"
barcode2_annotation:
  fb_info: examples/FB_info.tsv
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
      whitelist: examples/CELL_WL.txt
barcode2:
  name: FB
  segments:
    - segment: FB
      read: R2
      method: fixed
      start: 6
      length: 6
      whitelist: examples/FB_WL.txt
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
  repeats: 1
quick_test:
  enabled: true
  reads: 10000
performance:
  threads: 2
  # Barcode whitelist correction is parallelized by extracted-read chunks.
  # Omit barcode_workers to use threads.
  # barcode_workers: 2
  # UMI-tools groups are parallelized with processes. Omit umi_workers to use threads.
  # umi_workers: 2
  umi_batch_size: 5000
  umi_sqlite_cache_mb: 64
  chunk_size: 10000
  extraction_preview_reads: 1000
  compression_level: 3
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
    for step in ("extract", "correct", "dedup", "matrix", "downsample", "report"):
        p = commands.add_parser(step, help=f"Run the {step} step"); _common(p)
    quick = commands.add_parser("quick-test", help="Inspect a small read subset without running the full pipeline")
    quick.add_argument("--config", required=True)
    quick.add_argument("--sample", action="append", help="Sample name (repeatable; default: all)")
    quick.add_argument("--reads", type=int, default=None, help="Paired reads to inspect per sample")
    quick.add_argument("--threads", type=int, default=None, help="Cutadapt quick-test workers")
    validate = commands.add_parser("validate-config", help="Validate configuration and inputs")
    validate.add_argument("--config", required=True)
    init = commands.add_parser("init-config", help="Write an example configuration")
    init.add_argument("--out", required=True); init.add_argument("--force", action="store_true")
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
        names = _selected(config, args.sample)
        if args.threads is not None:
            if args.threads < 1:
                raise ConfigError("--threads must be >= 1")
            config.threads = args.threads
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
                    f"  Read IDs: {summary['read_ids']}\n  Log: {summary['log']}"
                )
        else:
            if args.command in {"extract", "dedup"}:
                check_external_tools()
            for name in names: run_step(config, name, args.command, args.overwrite)
            if args.command == "report": batch_report(config, names)
        print(f"TagForge {args.command} completed for: {', '.join(names)}")
        return 0
    except (ConfigError, FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"tagforge: error: {exc}", file=sys.stderr); return 2


if __name__ == "__main__":
    raise SystemExit(main())
