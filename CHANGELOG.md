# Changelog

Every code-change release increments the patch component of the TagForge
version. The version is kept in package metadata, the CLI, and pipeline logs.

## 0.1.6 — 2026-06-22

- Added resumable extraction using complete appended gzip members and an atomic
  manifest containing safe output bytes, completed reads, counters, elapsed
  time, version, and input/config fingerprint.
- Resume truncates any uncommitted gzip tail, rapidly skips completed FASTQ
  pairs without Cutadapt work, and continues from the next batch.
- `--overwrite` explicitly discards extraction resume state and starts over.

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
