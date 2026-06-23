# Changelog

Every code-change release increments the patch component of the TagForge
version. The version is kept in package metadata, the CLI, and pipeline logs.

## 0.1.10 — 2026-06-23

- Add the new grouped configuration layout: `barcode1`, `barcode2`, and `umi`
  are top-level sequence modules with named `segments`, so segment entries no
  longer need a hard-coded `target`.
- Allow custom barcode target names such as `PB` and `FB` while preserving the
  output roles `barcode1` and `barcode2` for downstream matrix generation.
- Add top-level `correction_barcode` defaults with per-segment overrides.
- Add top-level `correction_umi` for UMI-tools method/distance settings and
  reserve `umi` for UMI segment definitions in the new layout.
- Add top-level `linker.max_mismatch` as the default linker mismatch setting,
  still overridable by individual segments.
- Add `samples.auto` to discover samples from one-level raw-data directories
  such as `01_raw/{sample}/{sample}_raw_1.fq.gz`.
- Make `make-slurm` default to one Slurm array script plus `samples.tsv`, with
  `--mode per-sample` retaining the previous script-per-sample output.

## 0.1.9 — 2026-06-23

- Parallelize barcode whitelist correction across extracted-read chunks using
  process workers while preserving ordered gzip-member commits for safe resume.
- Add `performance.barcode_workers` to tune barcode correction independently
  from extraction and UMI correction; when omitted it follows
  `performance.threads`.
- Report barcode correction backend, requested/effective workers, and chunk
  size in pipeline logs and correction summaries.

## 0.1.8 — 2026-06-23

- Add resumable barcode correction with atomic manifests and complete appended
  gzip members for both valid-read and correction-trace outputs.
- Restore correction counters and elapsed/input progress, truncate uncommitted
  output tails, and rapidly skip already committed extracted rows on resume.
- Add `correction_progress` pipeline log records and a live
  `correction_progress.tsv` containing throughput, input percentage, ETA,
  resume-skip percentage, valid-read counts, and temporary-output sizes.
- Recover safely if interruption occurs between final output renames and final
  correction-statistics/checkpoint creation.

## 0.1.7 — 2026-06-22

- Parallelize independent Barcode1–Barcode2 UMI correction groups across
  multiple processes while continuing to use UMI-tools `UMIClusterer`.
- Bound multiprocessing memory by batching complete groups by unique-UMI count
  and allowing at most one pending batch per worker.
- Add `performance.umi_workers`, `umi_batch_size`, and `umi_sqlite_cache_mb` so
  UMI CPU use and memory can be tuned independently from Cutadapt extraction.
- Use a bounded, disk-backed disposable SQLite aggregation cache and report
  worker count, phase timings, group counts, and peak batch size in the log.

## 0.1.6 — 2026-06-22

- Added resumable extraction using complete appended gzip members and an atomic
  manifest containing safe output bytes, completed reads, counters, elapsed
  time, version, and input/config fingerprint.
- Resume truncates any uncommitted gzip tail, rapidly skips completed FASTQ
  pairs without Cutadapt work, and continues from the next batch.
- `--overwrite` explicitly discards extraction resume state and starts over.
- Preserve committed input percentage across resume and report fast-forward
  progress separately, preventing percentage regression and invalid ETA values.

## 0.1.5 — 2026-06-22

- Changed parallel extraction to bounded batches so worker task submission and
  memory use no longer grow with total FASTQ size.
- Added live extraction progress with approximate input percentage, throughput,
  ETA, estimated finish time, and temporary-output size.
- Added an immediately readable extracted preview (1,000 rows by default) and
  periodic gzip flushing so users can assess results during long runs.
- Made example FASTQ gzip generation byte-reproducible.
- Fixed gzip progress to use physical compressed bytes rather than the
  uncompressed text position; non-final batches are capped below 100%.

## 0.1.4 — 2026-06-22

- Made fixed-position correction symmetric: retain `max_shift` bases on both
  sides and test signed left/right shifts around the configured start.
- Added separate left/right shift statistics while preserving signed distances
  in the correction trace.
- Replaced repeated raw-barcode/JSON fields with seven compact columns:
  comma-separated ordered segment values plus F/L/X method codes.
- Intentionally reject older extracted schemas and version checkpoints so an
  upgrade cannot silently reuse intermediates that lack the left shift margin.

## 0.1.3 — 2026-06-22

- Added `quick-test` for fast small-subset data inspection.
- Made `tagforge run` execute quick-test by default, with an explicit opt-out.
- Quick-test inspects the first 10,000 read pairs by default and stops without
  scanning the remainder, with examined read IDs and performance metrics.
- Added read-length, N-content, extraction, correction, UMI duplication, top
  Barcode1/feature, TSV, HTML, and dedicated-log summaries.

## 0.1.2 — 2026-06-22

- Enumerate all Cutadapt-verified left/right linker matches.
- Select the globally closest correctly oriented linker pair.
- Report multiple-linker hits, candidate-pair counts, and selected gap ranges.

## 0.1.1 — 2026-06-22

- Added strict linker-first/fixed-fallback extraction semantics.
- Added Cutadapt extraction timing, worker, throughput, success, failure,
  fallback-rescue, and barcode-validity statistics.
- Added configurable Slurm scheduler directives and Conda activation.
- Enabled parallel Cutadapt-backed extraction.
- Added automatic logarithmic downsampling ratios.
