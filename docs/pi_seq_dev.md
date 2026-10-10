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
Sample and batch HTML reports link this dataset, and Excel reports include
its relative path on the `Results` sheet. Read a Parquet fragment directory
as a single table for downstream analysis.

Batch reports are `02_output/00_report/TagForge_batch_report.html` and
`TagForge_batch_report.xlsx`. The default remains six TSV tables per sample,
with no per-read debug exports. All selected ratios are 1, so the batch workbook
has no redundant `Counts_final` sheet.

The [current run audit](archive/testing/pipeline-layout-20261010-183901/summary.json)
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
