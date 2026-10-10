# Migration to TagForge 2.0

This document describes TagForge 2.0.0. Candidate and final release validation
records identify the exact source snapshot, package and dependencies.

## Rebuild trustworthy completion state

Old version-only `.done` files are not trusted. Back up results you need to keep,
then run the full workflow with the new configuration. Do not edit `.done`
contents to reuse old results. Independent downstream commands require current,
complete upstream manifests and report the missing or stale stage.

Manifests identify stage settings, upstream generations and artifacts.
Whitelist/annotation fingerprints include file contents; FASTQ fingerprints use
path, size and nanosecond modification time by default. Set
`resume.input_fingerprint: sha256` when content hashing of FASTQ is required.
Parquet validation checks fragments, schema, footer and row metadata; it does
not checksum every data page. Recovery was tested for process interruption,
including publication boundaries, rather than sudden power loss.

Extraction/correction resume committed batches. Deduplication, matrices,
mapping and reports safely recompute their interrupted stage. Same-sample
writers are locked; separate samples can run concurrently.

## Configuration and read eligibility

Unknown keys, negative correction distances, invalid whitelist lengths and
conflicting output names are rejected. Correct spelling or remove obsolete
unused fields instead of relying on silent ignoring.
Coordinates, counts, seeds and resource integers reject fractional, boolean
and nonfinite values rather than truncating them. Exact integer-valued numbers
remain valid. CB alias parent/child ratios must be finite.
UMI segments reject barcode-correction activation and whitelists because they
do not participate in that correction step; use `correction_umi` for clustering.
An explicit segment `target` must agree with its containing role.

- `downsample.enabled: false` skips that stage. Reports use full-depth counts
  and leave selected-ratio metrics empty.
- `resume.enabled: false` disables both completion reuse and batch resume;
  every planned stage recomputes and publishes a fresh manifest.
- `output.correction_trace: false` disables the large trace, while compact
  full-feature annotation QC remains available.
- PB–CB requires `barcode2.sequence_only: true`, correction disabled for CB
  segments, and no CB whitelist. Dominance uses read support and retains support
  strictly greater than the configured threshold. Alias resolution merges
  duplicate final PB–CB–UMI molecules and sums their supporting reads.
  Alias parents must have strictly more reads than their children, including
  when `parent_child_ratio: 1`; equal-support CBs cannot alias each other.

Quick-test remains a prefix diagnostic. It uses the formal read eligibility
rule: all required segments extracted, barcodes accepted, feature annotation
accepted where required, and every UMI segment of configured length using ACGT.
Exclusive rejection reasons reconcile input reads with eligible reads.
Only explicitly configured `quick_test.min_valid_rate` and
`quick_test.max_n_base_rate` thresholds cause QC to stop a run.
Quick-test rates and threshold values use fractions (0–1).
Quick-test `duplication_ratio` now uses percent (0–100);
older quick-test files stored that metric as a fraction. The new TSV records
both unit conventions explicitly.

## Tables and matrix exports

Batch reports now live under `{output_dir}/00_report/` as
`TagForge_batch_report.html` and the unified `TagForge_batch_report.xlsx`.
`00_summary.xlsx` and `01_FB_not_in_anno.xlsx` are retired; their QC and
annotation details are included in `Summary`, `Downsample` and `Warnings`.
The unified workbook always writes both `Counts_raw` and `Counts_final` for
standard and PI-seq libraries. `Counts_final` uses the validated selected counts;
when the selected ratio is 1 or downsampling is disabled, its counts match
`Counts_raw`. Previous `count`/`counts` sheet consumers must
choose the intended source explicitly. Feature columns are ordered by name.
Rate columns in `Summary` use fractions; saturation and non-singleton metrics
in `Downsample` use percent. PB–CB retains mapping summaries/previews without
creating feature matrices.

Run `tagforge report --config ...` to refresh both formats; valid sample
checkpoints are reused. `tagforge summary` refreshes the same workbook while
preserving HTML. After successful workbook publication, duplicate workbooks
in `{output_dir}/00_report/` are removed. Back up old reports before upgrading.
Older layouts also placed summaries directly under the output root and the
batch report under `{workdir}/00_report/`; remove those old copies after
checking the new reports. For older report-only changes the sample data stayed in place; the current
pipeline directory change below requires backing up old sample outputs and rerunning.

Current tables use stable Parquet schemas, including integer count fields.
A dataset may be one file or a fragment directory; consumers must enumerate its
actual fragments. A zero-row dataset contains a valid schema-bearing fragment.
Missing/incomplete current Parquet is an error, even if an old TSV exists.
Legacy TSV is read explicitly when current Parquet is absent; it does not make
an old completion marker trustworthy. Legacy multi-segment annotation QC
without trace requires rerunning correction.

The primary full-depth matrix is `{sample}.raw_molecule_pairs.parquet`, with
`{sample}.raw_matrix_summary.json`. Use `export-matrix` for TSV, MTX or Parquet
exports. PI-seq uses full-depth PI-filtered molecule details as its main analysis
input. PB–CB reports consume mapping artifacts rather than count-matrix files.
`export-matrix` requires a current matrix manifest, including current upstream
stages. After changing inputs or counting parameters, rerun the affected stages
before exporting; an existing matrix file alone is not accepted as current.

## Selected ratio and diagnostic metrics

Downsampling reports two percentage metrics. `sequencing_saturation` is
`100*(reads-molecules)/reads`; the duplicate `duplication_ratio` column is
removed from per-ratio metrics, curves and selected-result views.
`non_singleton_molecule_fraction` is `100*(molecules-singletons)/molecules`.
Zero denominators yield zero. Units and denominators refer to the same eligible
read population and fixed corrected molecule groups.

Each ratio is sampled once. Remove `downsample.repeats` from older configs;
this setting is no longer accepted. The ratio maximizes the observed
non-singleton fraction; ties prefer the lower ratio. This
diagnostic does not establish recommended sequencing depth. Supporting reads
are thinned within full-depth corrected groups; UMI clustering and PI filtering
are not repeated for each depth.
PB–CB diagnostics use the deduplicated scopes before CB alias mapping; final
resolved and retained mapping counts are reported separately.

`downsample_metrics.tsv` retains one result per ratio. `downsample_curve.tsv`
stores observed metrics and analytical expected counts together, with separate
`observed_*` and `expected_*` fields. It never treats a ratio of expected counts
as an expected fraction. The former summary/expectations TSV files are retired.
`selection.json` alone records criterion, seed, metrics and artifacts;
selection schema 2 omits repeat numbers and repeat statistics. Rerun downsampling
and reports to replace older outputs; the two selection-point TSV aliases are retired. Reports derive selection
views from the JSON. Ratio 1 references full-depth artifacts; missing or changed
selected artifacts remain errors. Legacy `optimal_*` matrix/detail names retain
the selected-ratio meaning.

`correction_qc.json` combines barcode correction statistics and the annotation
summary; the full unannotated barcode table stays in TSV.gz. The corresponding
old stats/annotation QC TSVs are retired. Current stage manifests enforce the
new contract, so `tagforge run` reuses valid extraction and rebuilds correction
and its dependent stages. Do not edit completion manifests to migrate files.
Each sample now has one live `progress.json`; history remains in pipeline.log.
Extraction previews default to off (`performance.extraction_preview_reads: 0`),
and per-read quick-test QC/IDs require `output.quick_test_read_details: true`.
Disabling these outputs removes earlier generated debug copies when their
corresponding command runs; aggregate QC remains available.

Fixed seeds reproduce within tested runtime/algorithm versions. Python's
native binomial sampler and the older fallback can give different draws.

## Resource settings and linker policy

DuckDB matrix operations use `performance.threads`,
`performance.matrix_memory_mb` (default 16384, minimum 32) and local scratch.
UMI batch size is a target, not a hard cap on a connected component or full
barcode–feature group. Optional `performance.umi_max_group_umis` and
`umi_max_component_umis` budgets stop with a specific barcode/feature;
both default to null. Group/component/batch peaks are logged. Hamming-1
components are processed independently; connected correction edges are retained.

`legacy_shortest` remains the default linker policy. `linker.strategy:
length_position` and segment-level `linker_strategy` opt in to expected-length
and position ranking. The strategy affects extraction fingerprints. Synthetic
truth tests and four-fixture differences were audited, but the new policy can
change labels or reject previously accepted reads; a default change is deferred.

## Installation and public examples

Historical installed-wheel flows passed on Linux with Python 3.9, 3.11 and 3.13.
Their dependency versions and results remain archived. Current development and
release checks use `/data/home/yfu/miniconda3/envs/tagforge` (Python 3.11); the
old test environments have been removed. Package compatibility metadata is
unchanged, but other versions are outside the current validation scope.
UMI-tools 1.1.6's legacy source build needs compatible preinstalled build tools;
follow the tested pip or Conda commands in README.

`tagforge init-example --out DIRECTORY` creates a standalone public synthetic
project from the installed wheel. Real PI-seq FASTQ and experiment outputs are
excluded from distributions. The tiny synthetic source fixture is generated
deterministically and is the only FASTQ fixture allowed in the sdist.

## Parquet compression

Parquet output now defaults to Snappy, including fragment datasets and matrix
exports. To retain ZSTD, set `performance.parquet_compression: zstd`.
`performance.compression_level` still controls gzip and ZSTD; Snappy has no
level setting. Changing the codec invalidates affected stage checkpoints.
The table schema and molecule counts are unchanged. Fragment datasets remain
directories; read all fragments as one table.

## Pipeline directory layout and result entry points

Sample outputs follow processing order: `00_logs/`, `01_extracted/`,
`02_corrected/`, `03_molecules/`, `04_matrix/`, `05_downsample/` and `06_report/`.
Valid reads belong to `02_corrected/`; `03_molecules/` contains deduplicated
and PI-filtered molecules. Batch reports remain in `{output_dir}/00_report/`.
HTML reports link the primary result and Excel includes a `Results` sheet.
Standard libraries use `04_matrix/{sample}.raw_molecule_pairs.parquet`;
PI-seq uses `03_molecules/{sample}.molecule_detail.rmMP.parquet`;
PB–CB uses `04_matrix/{sample}.pb_cb_map.tsv.gz`. These are full-depth results.

Archive the entire old sample output and rerun to adopt this layout. The pipeline
rejects a sample directory containing the former `01_checkpoint/` directory.
Do not move manifests alone or edit fingerprints: input paths and selected
references must match the newly published outputs. Keep the backup until the
new data and checkpoints have been validated.

## Temporary work layout

The pipeline no longer creates `08_tmp`. Stage transactions use
`00_logs/checkpoints/.work/stage-{step}/`; completed stages remove their work,
while failed stages retain recovery data. Sort and SQLite scratch directories
are created inside staged `03_molecules/` and cleaned after use. Interrupted
`08_tmp/stage-{step}/` work is migrated when resuming a compatible sample layout.
Do not delete unfinished work before resuming.

## Stage manifest filenames

Completion records use `00_logs/checkpoints/{number}_{step}.manifest.json`:
`01_extract`, `02_correct`, `03_dedup`, `04_matrix` (or `04_pair-map`),
`05_downsample` and `06_report`. The number identifies the stage, not its
execution time. PB–CB mapping follows downsampling when enabled.
Earlier unnumbered JSON `.manifest.json` and JSON `.done` files in the current
checkpoint directory remain readable and migrate by renaming when a validated
stage is reused. A corrupt current manifest never falls back to an older file.
Version-only legacy markers require recomputation. This filename compatibility
does not migrate the former directory layout.
