# Changelog

Every code-change release increments the patch component of the TagForge
version. The version is kept in package metadata, the CLI, and pipeline logs.

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
