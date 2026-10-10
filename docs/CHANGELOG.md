# Changelog

Major releases may change configuration and result contracts; minor and patch
releases contain compatible features and fixes. Package metadata, CLI and logs
use the version declared in `tagforge.__version__`.

## Unreleased

- Parallelise `extract`: the coordinator only reads raw FASTQ blocks; workers
  parse, extract, aggregate segment statistics and write their own Parquet
  fragments. Exact (`linker.max_mismatch: 0`) linker search uses `str.find`.
  Output tables and statistics are identical to the previous implementation;
  resume state written by earlier builds remains valid.
- Speed up `correct`: workers read and write whole fragments, distance-1
  whitelist correction enumerates neighbours instead of scanning the whitelist,
  and trace rows are built only when tracing is enabled.
- Add the `duckdb` UMI aggregation backend (now the default): DuckDB aggregates
  reads into hash buckets, workers cluster buckets independently and write
  fragments. PI-seq multi-PB filtering also runs in DuckDB. New
  `performance.umi_duckdb_memory_mb` (default 65536). `external_sort` and
  `sqlite` remain available. Molecule row order is no longer globally sorted.
- On ABL00074 PB-FB-A (373.7M valid reads), `dedup` including PI filtering fell
  from 8329 s to 765 s with identical molecule counts, PI QC and filtered tables.

## 2.0.0 — 2026-10-09 (local release preparation)

- Number processing directories from `01_extracted` to `06_report`; move valid
  reads into `02_corrected` and reserve `03_molecules` for molecule results.
  Nest checkpoints in `00_logs/checkpoints` with `01_extract.manifest.json`
  through `06_report.manifest.json`. Old output layouts require backup and rerun.
- Add main-result links to HTML reports and a Results sheet to Excel reports
  for standard matrices, PI-seq molecule details and PB–CB maps.

- Use the existing Python 3.11 TagForge environment for current validation.
  Remove historical test environments and archive experiment records under
  `docs/archive/testing/`.

- Remove the separate `08_tmp` output directory. Keep recoverable stage work
  under `00_logs/checkpoints/.work`, clean scratch after use and migrate interrupted
  stages from the former location without repeating committed batches.

- Reduce default PI-seq TSV outputs from 17 to 6 per sample: unify progress and
  selection JSON, combine correction QC and downsample curves, and make
  per-read debug exports opt-in. Preserve per-ratio sampling metrics and annotation detail.
- Merge batch report and summary workbooks into one `TagForge_batch_report.xlsx`,
  retaining full-depth/selected counts, QC, downsample metrics and annotation
  details; the summary command refreshes the same workbook.
- Consolidate all batch HTML/Excel reports under `{output_dir}/00_report/`,
  including summaries; upstream invalidation withdraws reports at that location.
- Default all Parquet outputs and matrix exports to Snappy for viewer
  compatibility; retain configurable ZSTD and invalidate storage checkpoints
  when the selected codec changes.

- Reject stale matrix exports after input/configuration changes; independent
  recomputation verifies the changed molecule counts before export.
- Stabilize stage fingerprints, dependency invalidation and recoverable artifact
  publication; explicit selected references prevent stale-count reuse.
- Unify complete-read QC and stable Parquet schemas; enforce configuration
  switches and independently report all three library types.
- Stream matrix aggregation/export, skip redundant exact-match searches, publish
  batch summaries once and diagnose oversized UMI scopes.
- Sample each downsampling ratio once, select by observed non-singleton fraction,
  retain analytical expected counts, and keep legacy linker behavior as the default after difference review.
- Validate installed wheels across Python 3.9/3.11/3.13, repair Excel/report
  compatibility and provide standalone public synthetic examples.
- Reject fractional/nonfinite configuration values and ignored UMI settings,
  enforce CB alias support ordering and validate migration thread overrides.
- Isolate integration fixtures and synchronize development package metadata.
  See MIGRATION.md and RELEASE_NOTES.md for breaking contracts, installation
  instructions and explicitly deferred research.


Historical candidate and 0.x releases are archived in
[docs/archive/CHANGELOG.md](archive/CHANGELOG.md).
