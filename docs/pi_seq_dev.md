# PI-seq development fixture

The fixture contains four paired-end subsets from the raw reads in
`PI-seq/ABL00073_TagForge`. Each sample contains exactly 50,000 R1/R2 read
pairs:

- `100-1-N1` and `100-2-SLE1` were randomly sampled with seqkit 2.13.0
  `sample2 -n 50000 -2 -s 20261009`. The same seed was used for both mates to
  preserve pairing.
- `100-3-SLE2` and `100-5-N2` contain the first 50,000 records from each
  corresponding R1/R2 source FASTQ, extracted with `seqkit head -n 50000`.
  These are prefix subsets rather than random samples.

`00_config/config.yaml` keeps the source PI-seq barcode, linker, and UMI
coordinates. The feature annotation and all six required whitelists are copied
under `00_config`, so the configured fixture is self-contained. All four
samples are explicitly listed in the config. Downsampling uses five ratios to
keep a development run compact. Use the two random subsets for distributional
checks; the prefix subsets are useful for quick smoke tests but should not be
treated as representative for barcode or molecule-distribution analyses.

From the TagForge repository root, validate the fixture and run it with:

```bash
tagforge validate-config --config tests/fixtures/pi_seq_dev/00_config/config.yaml
tagforge run --config tests/fixtures/pi_seq_dev/00_config/config.yaml
```

Outputs are written to `tests/fixtures/pi_seq_dev/02_output/`.

The latest [half-width refresh](archive/testing/downsample-half-width-20261010-202545/summary.json)
changes the downsample chart and PDF width from 1000 to 500 pixels while keeping
the 540-pixel height, both panels, annotations, fonts and all data unchanged.

The [numerical annotation refresh](archive/testing/downsample-numeric-labels-20261010-202408/summary.json)
adds the sample name and full-depth (ratio 1) Saturation and non-singleton molecule
fraction percentages inside the upper panel, rounded to two decimal places.
Both panels and all curve data remain unchanged. Five relevant tests and the
Arial browser check pass; the [PDF preview](archive/testing/downsample-numeric-labels-20261010-202408/100-1-N1.downsample.pdf)
contains the same labels.

The [two-panel figure refresh](archive/testing/two-panel-downsample-figure-20261010-201817/summary.json)
keeps the original upper saturation/non-singleton panel and lower UMI Types/
UMI detected once panel, including their colours, axes and data. Only the
surrounding headings and small-print notes are removed; the sample name is
inside the chart and all text uses Arial. Five relevant tests and a browser
rendering check pass. The [PDF preview](archive/testing/two-panel-downsample-figure-20261010-201817/100-1-N1.downsample.pdf)
uses actual Arial. Numerical data and batch workbook cells are unchanged.

The prior [two-metric report refresh](archive/testing/two-downsample-metrics-20261010-195556/summary.json)
removes the duplicate downsample `duplication_ratio` field from metrics, curves
and selected-result views. Saturation and non-singleton fraction values are
unchanged; batch HTML uses `Sat.` and `Non-singleton (%)`. All 105 tests pass.

The current outputs use one seeded draw per downsampling ratio. Repeat-count
configuration, repeat identifiers, mean/SD columns and error bars have been
removed. Only the four downsample and four report stages were refreshed;
upstream artifacts and every observed per-ratio count are unchanged. The
[single-sampling audit](archive/testing/single-sampling-20261010-193032/summary.json) includes the
prior reports, source snapshot and 105 passing tests. Selection records use
schema 2, and all selected ratios remain 1.

The latest complete run on 2026-10-10 uses **TagForge 2.0.0** with the new
pipeline directory layout. All 24 stages were forcibly recomputed for 200,000
read pairs, including the four prefix quick-tests. The run took 256.50 seconds
on the shared host; this is a single observation, not a general speedup claim.
Parquet remains Snappy. Counts are 189,541 eligible reads, 185,505 deduplicated
molecules and 185,469 molecules after PI filtering.

Each sample now contains `00_logs/`, `01_extracted/`, `02_corrected/`,
`03_molecules/`, `04_matrix/`, `05_downsample/` and `06_report/`.
`00_logs/checkpoints/` holds `01_extract.manifest.json` through
`06_report.manifest.json`. Valid reads are in `02_corrected/`; the main PI-seq
analysis input is `03_molecules/{sample}.molecule_detail.rmMP.parquet`.
Sample HTML reports link this dataset, and Excel reports include its relative
path on the `Results` sheet. Batch HTML omits the main-result block and the
small-print downsample explanation. Read a Parquet fragment directory
as a single table for downstream analysis.

Batch reports are `02_output/00_report/TagForge_batch_report.html` and
`TagForge_batch_report.xlsx`. The default remains six TSV tables per sample,
with no per-read debug exports. The batch workbook always includes `Counts_raw`
and `Counts_final`. All selected ratios are 1, so the two sheets currently
contain identical counts.

The [complete-run baseline audit](archive/testing/pipeline-layout-20261010-183901/summary.json)
indexes the frozen runtime source, input checksums, complete outputs, old-output
backup, tests and run log. All 24 logical tables, QC, downsample curves and
selection values match the backup. Five Excel workbooks pass ZIP/XML structure, Arial styles and cell-value checks;
all prior statistical values are preserved and main-result paths resolve.
Runtime and throughput measurements change on rerun. Two samples have small
changes in quick-test segment linker diagnostics; sample-level counts and
correction QC are unchanged. This reproducibility issue is recorded for
investigation before release; its cause is not yet confirmed.
`02_output/latest_run.json` identifies the current result. The old layout was
backed up and the pipeline rerun rather than manually editing checkpoints.

Previous Snappy runs, report consolidation and compact-output migrations remain
historical evidence under `docs/archive/testing/`. See the [validation index](validation.md)
for dated records; frozen snapshots retain their original directory names.
