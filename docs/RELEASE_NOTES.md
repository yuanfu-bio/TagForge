# TagForge 2.0.0 release notes

TagForge 2.0 stabilizes result reuse, configuration validation and data interfaces
for barcode–feature, PI-seq and PB–CB libraries. Full-depth corrected molecule
counts remain the primary result. This file prepares the local release;
publishing, uploading and tagging are separate actions.

## Changes that affect existing projects

- The workflow no longer creates `08_tmp`. Recoverable stage work uses
  `00_logs/checkpoints/.work/`, removed when no unfinished stages remain; disposable
  sort and database files are cleaned after use.

- All Parquet tables and matrix exports default to Snappy compression for
  viewer compatibility. Set `performance.parquet_compression: zstd` to retain
  ZSTD. `compression_level` applies to gzip and ZSTD, not Snappy. Changing
  the Parquet codec invalidates affected checkpoints without changing counts.
- Processing directories are numbered `01_extracted` through `06_report`;
  valid reads move into `02_corrected` and molecule details into `03_molecules`.
  Stage completion records use numbered names such as `01_extract.manifest.json`
  in `00_logs/checkpoints`. Archive old-layout sample outputs and rerun.
  Valid unnumbered JSON manifests in the current directory migrate on reuse.
  Old version-only `.done` files require recomputation. Stage manifests now
  validate relevant settings, dependency generations and actual artifacts.
  Extraction/correction recover committed batches; other stages safely recompute.
- Current tables use stable Parquet schemas, including valid empty datasets and
  integer counts. Missing or incomplete current data cannot silently select an
  old TSV or another count source. Matrix exports support TSV, MTX and Parquet.
  Export also checks the current matrix and upstream completion manifests,
  rejecting existing counts after inputs or counting settings change.
- Quick-test and formal correction share complete-read eligibility. Exclusive
  rejection reasons reconcile input and eligible reads. Prefix QC is diagnostic
  by default; only explicitly configured thresholds block the workflow.
- Unknown, invalid and inapplicable configuration fields are rejected. Disabled
  downsampling uses full-depth reports; disabled resume recomputes all stages.
- `duplication_ratio` and `sequencing_saturation` use read support as their
  denominator and percent units. Quick-test duplication has changed from the
  previous 0–1 unit to 0–100; its other QC rates and thresholds stay in 0–1.
- Downsampling samples each ratio once and selects by the observed
  non-singleton fraction, preferring the lower ratio on ties. Observed metrics
  and separate analytical expected counts are retained. Historical `optimal_*` names
  remain aliases for selected diagnostics, not validated optimal sequencing depth.
- PB–CB alias resolution merges final PB–CB–UMI molecules and preserves read
  support. Dominance uses reads and a strict threshold. Alias parents require
  strictly higher support, preventing reciprocal equal-support aliases at ratio 1.

See [MIGRATION.md](MIGRATION.md) for configuration conversion, output contracts,
checkpoint invalidation and the limits of legacy result compatibility.

## Performance and reports

Anchored exact barcode hits skip further distance search. Matrix aggregation and
exports stream results with explicit DuckDB resources. Batch summaries publish
once per invocation. UMI components stream without splitting correction edges;
optional group/component budgets report which scope exceeded the limit.

One sequential shared-host benchmark measured four 50,000-pair PI-seq fixtures
at 336.94 → 260.59 seconds and a 1,000,000-pair prefix at 556.55 → 528.85 seconds.
The million-pair GNU time maximum RSS changed from 1,477,692 to 1,027,056 KiB.
Logical valid-read, molecule, filtered-molecule and matrix rows were identical.
These are single-run observations, with sampling and shared-memory limitations;
they do not establish a general speedup or biological accuracy.

HTML reports link each library’s main full-depth result; Excel reports add a
`Results` sheet listing the paths. Standard libraries use raw molecule pairs,
PI-seq uses PI-filtered molecule details and PB–CB uses the dominant mapping.

HTML and Excel show mean/SD curves and explicit metric definitions. Workbook
styles, freeze panes and sheet-name collisions are corrected. The installed
wheel can generate a standalone public synthetic project with `init-example`.
Real FASTQ and local experiment outputs are excluded from release distributions.
Batch HTML and the unified `TagForge_batch_report.xlsx` live under
`{output_dir}/00_report/`; sample reports remain in `06_report/`. The workbook
combines QC and annotation details with explicitly separate full-depth and
selected count sheets. Duplicate summary/annotation workbooks are retired. Progress, selection and
correction QC use single JSON records; observed downsample statistics and
analytical expectations share one curve TSV. Per-read preview/QC exports are
opt-in. Per-ratio sampling metrics and full unannotated feature detail are retained.
See MIGRATION.md for refreshing old report layouts.

## Supported and deferred behavior

Current development and release validation use the existing Linux/Python 3.11
TagForge environment. Historical cross-version checks remain archived evidence;
additional environments are no longer a release requirement. Exact historical
dependency versions are indexed in [RELEASE_READINESS.md](RELEASE_READINESS.md).
Follow README's tested pip/Conda
installation steps for the legacy UMI-tools source-build dependency constraint.

`legacy_shortest` remains the default linker strategy. Experimental
`length_position` is explicit opt-in: eight synthetic truth cases and all 865
changed reads in four fixtures were reviewed, including losses supported by
fixed-coordinate whitelist matches. That evidence does not justify changing
the default.

Full UMI reclustering and library filtering at each sampled depth, sequencing
depth recommendations, quality-aware correction, large-whitelist indexing,
complete solutions for huge UMI groups and new PB–CB dominance defaults remain
research or later optimization work, as scoped in the original 2.0 plan.
