## 1. Project Goal

**Project Name: TagForge**

Build a production-ready Linux command-line pipeline for processing paired-end FASTQ sequencing data from antibody-conjugated oligonucleotide barcode libraries.

The pipeline should extract, correct, trace, deduplicate, downsample, and summarize Barcode1, Barcode2, and UMI information from raw paired-end FASTQ files, and generate clean count matrices, molecule-level detail files, saturation analysis results, logs, Excel reports, and HTML reports.

The project should be suitable for future GitHub release and routine use on Linux servers and Slurm HPC systems.

---

## 2. Biological Background and Terminology

The library contains multiple barcode and UMI components distributed across read1 and read2.

Common terms:

- **Barcode1**: Usually bead/cell/particle barcode, such as Cell Barcode, PI Barcode, PB Barcode, etc.
- **Barcode2**: Usually antibody feature barcode, also called FB.
- **UMI**: Unique Molecular Identifier.
- **Linker**: Fixed sequence adjacent to barcode segments, used for barcode extraction.
- **Whitelist**: Known valid barcode sequences used for barcode correction.
- **FB_info.tsv**: Feature barcode annotation table. It contains at least:
    - FB ID, e.g. `FB0288-V02`
    - FB sequence, e.g. `ATCACACGCA`
    - Antibody name, e.g. `CD4`
    - Optional additional metadata columns

For the final count matrix, Barcode2 should be represented by antibody name, such as `CD4`, not by raw FB sequence or FB ID.

---

## 3. Input Directory Structure

Example initial directory:

```
project/
├── 00_config/
│   ├── config.yaml
│   ├── submit.sh
│   ├── FB_info.tsv
│   ├── FB_WL.txt
│   ├── PB1_WL.txt
│   ├── PB2_WL.txt
│   └── ...
├── 01_raw/
│   ├── cDNA-M_M0-1/
│   │   ├── cDNA-M_M0-1_raw_1.fq.gz
│   │   └── cDNA-M_M0-1_raw_2.fq.gz
│   └── cDNA-M_M0-2/
│       ├── cDNA-M_M0-2_raw_1.fq.gz
│       └── cDNA-M_M0-2_raw_2.fq.gz
```

FASTQ input naming convention:

```
{sample}_raw_1.fq.gz
{sample}_raw_2.fq.gz
```

---

## 4. Library Structure

The library structure must be fully configurable.

Barcode1, Barcode2, and UMI may each be:

- a single segment
- multiple concatenated segments
- distributed across read1 and read2
- extracted by fixed position
- extracted using left/right linker sequences
- extracted using a combination of linker-based and fixed-position methods

Example:

```
Read1:
Barcode1_1 -- LINKER1 -- Barcode1_2 -- LINKER2 -- Barcode1_3 ...

Read2:
UMI1 -- Barcode2 -- UMI2 ...
```

Each barcode or UMI segment must be configurable with:

- segment name
- target type: `barcode1`, `barcode2`, or `umi`
- read source: `R1` or `R2`
- extraction method:
    - `fixed`
    - `linker`
- fixed start and length, if method is `fixed`
- left linker and/or right linker, if method is `linker`
- expected barcode length
- whitelist file, if applicable
- whether correction is enabled
- maximum forward shift
- maximum mismatch distance

---

## 5. Barcode Extraction Rules

Barcode extraction should happen in two stages:

### 5.1 Linker-based Extraction

For barcode segments with left and/or right linker specified, use **cutadapt** to locate linkers and extract the barcode between or adjacent to them.

Requirements:

- Use cutadapt where appropriate.
- Support linker mismatch parameters.
- Support paired-end input.
- Avoid writing excessive intermediate FASTQ files when possible.
- Preserve read-level information needed for downstream tracing.

### 5.2 Fixed-position Extraction

For barcode or UMI segments with fixed positions, extract directly from the raw read sequence according to configured start and length.

Position convention should be clearly defined:

- Use 0-based coordinates internally.
- Document whether YAML uses 0-based or 1-based coordinates.
- Prefer 0-based half-open intervals for clarity.

---

## 6. Barcode Correction Rules

Barcode correction should be applied to Barcode1 and Barcode2 segments according to their whitelists.

Correction should support both:

1. **Forward shift correction**
    - Example: extracted sequence may start several bases earlier than expected.
    - Parameter: `max_shift`
    - Default: `1`
    - Can be disabled.
2. **Mismatch correction**
    - Correct barcode sequence to nearest whitelist sequence if Hamming distance is within threshold.
    - Parameter: `max_mismatch`
    - Default: `1`
    - Can be disabled.

If a barcode consists of multiple segments, each segment may have its own whitelist and correction settings.

The pipeline must preserve both raw and corrected information.

For each segment, record:

- raw sequence
- shifted sequence, if shift correction was used
- corrected whitelist sequence
- correction status
- shift distance
- mismatch distance
- whether correction succeeded
- correction type:
    - `exact`
    - `mismatch_only`
    - `shift_only`
    - `shift_and_mismatch`
    - `failed`
    - `disabled`

After correction, concatenate corrected barcode segments to form final Barcode1 or Barcode2.

---

## 7. UMI Correction and Deduplication

UMI correction should use **UMI-tools** or an equivalent implementation of UMI-tools-compatible algorithms.

Requirements:

- Support common UMI-tools methods such as:
    - `unique`
    - `cluster`
    - `adjacency`
    - `directional`
- Default method should be `directional`.
- UMI correction/deduplication should be performed within each Barcode1-Barcode2 group.
- Preserve molecule-level detail before and after UMI correction.

For each final Barcode1-Barcode2 pair, retain:

```
Barcode1
Barcode2
raw_UMI
corrected_UMI
reads_count
```

This detail file is required for downstream molecule-level sampling.

---

## 8. Output Files

For each sample, create a structured output directory:

```
02_output/
└── {sample}/
    ├── 00_logs/
    ├── 01_checkpoint/
    ├── 02_extracted/
    ├── 03_corrected/
    ├── 04_matrix/
    ├── 05_detail/
    ├── 06_downsample/
    ├── 07_report/
    └── 08_tmp/
```

For all samples of this batch, create a structured output in 00_report, including excel and html:

```
Batch/
├── 00_report/
├── 00_config/
├── 01_raw/
├── 02_output/
```

### 8.1 Extracted Barcode Detail

File:

```
extracted/{sample}.extracted.tsv.gz
```

Compact columns (segment values follow configuration order):

```
read_id
barcode1_segments
barcode2_segments
umi_segments
methods
status
failure_reason
```

Segment fields are comma-separated strings. `methods` is a configuration-ordered
F/L/X code string. Combined raw values and per-read JSON keys are omitted to
avoid redundant storage. The file remains gzip-compressed because it can be
large.

### 8.2 Correction Trace

File:

```
corrected/{sample}.barcode_correction_trace.tsv.gz
```

Required columns:

```
read_id
segment_name
target_type
raw_sequence
shifted_sequence
corrected_sequence
whitelist_hit
correction_status
shift_distance
mismatch_distance
correction_type
```

This file is important for checking whether barcode correction is reliable.

### 8.3 Valid Reads Detail

File:

```
detail/{sample}.valid_reads.tsv.gz
```

Required columns:

```
read_id
barcode1
barcode2_sequence
barcode2_name
umi
correction_summary
```

Only reads with valid Barcode1, valid Barcode2, and valid UMI should be retained here.

### 8.4 Molecule Detail Before Downsampling

File:

```
detail/{sample}.molecule_detail.tsv.gz
```

Required columns:

```
barcode1
barcode2_name
corrected_umi
reads_count
raw_umi_count
```

Meaning:

- Each row represents one deduplicated molecule.
- `reads_count` is the number of reads supporting this molecule.
- `raw_umi_count` is the number of raw UMI sequences merged into this corrected UMI.

### 8.5 Raw Count Matrix

File:

```
matrix/{sample}.raw_count_matrix.tsv.gz
```

Rows:

```
Barcode1
```

Columns:

```
Antibody names from FB_info.tsv
```

Values:

```
Number of corrected UMIs for each Barcode1-antibody pair
```

### 8.6 Barcode Correction Statistics

File:

```
corrected/{sample}.barcode_correction_stats.tsv
```

Report for Barcode1 and Barcode2 separately, and for each barcode segment:

Required metrics:

```
total_reads
extracted_reads
valid_reads
invalid_reads
valid_rate
exact_count
mismatch_only_count
shift_only_count
shift_and_mismatch_count
failed_count
disabled_count
mismatch_0_count
mismatch_1_count
mismatch_2_count
...
shift_0_count
shift_1_count
shift_2_count
...
```

Also report:

- final Barcode1 valid rate
- final Barcode2 valid rate
- combined valid rate
- completely uncorrected barcode count
- number of barcodes corrected by mismatch
- number of barcodes corrected by shift only
- number of barcodes corrected by both shift and mismatch

### 8.7 Downsampling Saturation Metrics

File:

```
downsample/{sample}.downsample_metrics.tsv
```

Required columns:

```
sample
downsample_ratio
reads_sampled
umi_types
umi_detected_once
duplication_ratio
sequencing_saturation
```

Definitions should be implemented and documented clearly.

Suggested definitions:

```
UMI Types = number of unique corrected molecules detected
UMI detected once = number of corrected molecules supported by exactly one read
duplication_ratio = (dup) * 100 / (total)
seq_saturation = (1 - (single / n_duplicate_set)) * 100
```

If better definitions are used, document them clearly in README.

### 8.8 Optimal Saturation Matrix

When the library is over-sequenced, sequencing saturation may first increase and then decrease. The pipeline should identify the downsampling point where `Sequencing Saturation` is maximal.

At that ratio, generate:

```
matrix/{sample}.optimal_saturation_count_matrix.tsv.gz
detail/{sample}.optimal_saturation_molecule_detail.tsv.gz
```

The detail file should retain:

```
barcode1
barcode2_name
corrected_umi
reads_count_at_optimal_downsample
```

Also output:

```
downsample/{sample}.optimal_saturation_point.tsv
```

Required columns:

```
sample
optimal_downsample_ratio
max_sequencing_saturation
reads_sampled
umi_types
umi_detected_once
duplication_ratio
```

---

## 9. Reports

Generate both Excel and HTML reports.

### 9.1 Excel Report

File:

```
report/{sample}.report.xlsx
```

Suggested sheets:

- `Summary`
- `Barcode QC`
- `Correction Stats`
- `UMI Stats`
- `Downsample Metrics`
- `Optimal Saturation`
- `Top Barcode1`
- `Top Features`

### 9.2 HTML Report

File:

```
report/{sample}.report.html
```

The HTML report should include:

- sample summary
- read extraction statistics
- barcode correction statistics
- final valid rate
- UMI deduplication statistics
- count matrix summary
- downsampling saturation curves
- optimal saturation point
- warnings and QC flags

Use Plotly  plotting library.

### 9.3 Batch Summary

Batch
├── 00_report/

Summarize key results of all samples in the batch summary folder. In the excel, you should include at least two sheets, ”meta“ for important meta info, like sample name, barcode valid rartio, total reads, and so on. And “counts” for aggregated Barcode2 bulk Count Matrix, whose rows indicate sample, while columns indicate Barcode2 Info, if FB, then Antibody name, like CD4. More over, design a fancy html report.

---

## 10. Configuration File

Use YAML as the primary configuration format.

Example:

```yaml
project:
name: antibody_oligo_pipeline
workdir: /path/to/project
output_dir: 02_output
tmp_dir: tmp

samples:
-sample: cDNA-M_M0-1
r1: 01_raw/cDNA-M_M0-1/cDNA-M_M0-1_raw_1.fq.gz
r2: 01_raw/cDNA-M_M0-1/cDNA-M_M0-1_raw_2.fq.gz

-sample: cDNA-M_M0-2
r1: 01_raw/cDNA-M_M0-2/cDNA-M_M0-2_raw_1.fq.gz
r2: 01_raw/cDNA-M_M0-2/cDNA-M_M0-2_raw_2.fq.gz

barcode2_annotation:
fb_info: 00_config/FB_info.tsv
id_column: FB_ID
sequence_column: sequence
name_column: antibody_name

segments:
-name: PB1
target: barcode1
read: R1
method: linker
left_linker:null
right_linker: AGGTC
length:8
whitelist: 00_config/PB1_WL.txt
correction:
enabled:true
max_shift:1
max_mismatch:1
allow_shift:true
allow_mismatch:true

-name: PB2
target: barcode1
read: R1
method: linker
left_linker: AGGTC
right_linker: TTAAC
length:8
whitelist: 00_config/PB2_WL.txt
correction:
enabled:true
max_shift:1
max_mismatch:1
allow_shift:true
allow_mismatch:true

-name: FB
target: barcode2
read: R2
method: fixed
start:10
length:10
whitelist: 00_config/FB_WL.txt
correction:
enabled:true
max_shift:1
max_mismatch:1
allow_shift:true
allow_mismatch:true

-name: UMI1
target: umi
read: R2
method: fixed
start:0
length:10

-name: UMI2
target: umi
read: R2
method: fixed
start:20
length:8

umi:
correction_method: directional
max_distance:1

downsample:
enabled:true
ratios:[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
random_seed:12345
repeats:1
choose_optimal_by: max_sequencing_saturation

performance:
threads:8
chunk_size:100000
max_memory_gb:16
compression_level:3
keep_intermediate:false

resume:
enabled:true
overwrite:false

logging:
level: INFO
```

---

## 11. Command Line Interface

Implement a clear CLI.

Suggested package command name:

```bash
tagforge
```

### 11.1 Run Full Pipeline

```bash
tagforge run \
  --config 00_config/config.yaml \
  --sample cDNA-M_M0-1 \
  --threads 8
```

Run all samples in config:

```bash
tagforge run \
  --config 00_config/config.yaml \
  --threads 8
```

### 11.2 Run Specific Step

```bash
tagforge extract --config 00_config/config.yaml --sample cDNA-M_M0-1
tagforge correct --config 00_config/config.yaml --sample cDNA-M_M0-1
tagforge dedup --config 00_config/config.yaml --sample cDNA-M_M0-1
tagforge matrix --config 00_config/config.yaml --sample cDNA-M_M0-1
tagforge downsample --config 00_config/config.yaml --sample cDNA-M_M0-1
tagforge report --config 00_config/config.yaml --sample cDNA-M_M0-1
```

### 11.3 Validate Config

```bash
tagforge validate-config --config 00_config/config.yaml
```

### 11.4 Generate Example Config

```bash
tagforge init-config --out 00_config/config.example.yaml
```

### 11.5 Generate Slurm Scripts

```bash
tagforge make-slurm \
  --config 00_config/config.yaml \
  --out slurm_jobs \
  --threads 8 \
  --mem 16G \
  --time 24:00:00
```

This should generate one Slurm script per sample and optionally a master submission script.

---

## 12. Slurm Support

Generate scripts like:

```bash
#!/bin/bash
#SBATCH --job-name=tagforge_cDNA-M_M0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cDNA-M_M0-1.%j.out
#SBATCH --error=logs/cDNA-M_M0-1.%j.err

set -euo pipefail

tagforge run \
  --config 00_config/config.yaml \
  --sample cDNA-M_M0-1 \
  --threads 8
```

Also generate:

```bash
submit_all.sh
```

containing:

```bash
sbatch cDNA-M_M0-1.slurm
sbatch cDNA-M_M0-2.slurm
```

---

## 13. Checkpoint and Resume

The pipeline must support restart/resume.

Each step should write a checkpoint file only after successful completion.

Example:

```
checkpoint/extract.done
checkpoint/correct.done
checkpoint/dedup.done
checkpoint/matrix.done
checkpoint/downsample.done
checkpoint/report.done
```

Rules:

- If a checkpoint exists and output files exist, skip that step.
- If `-overwrite` is set, rerun the step and replace outputs.
- Write temporary outputs with `.tmp` suffix and atomically rename after success.
- Never mark a step as complete before all expected outputs are successfully written.

---

## 14. Performance Requirements

The data may be very large. Optimize for low memory usage.

Requirements:

- Process FASTQ files in streaming or chunk mode.
- Do not load all reads into memory.
- Use gzip streaming.
- Use chunked processing, configurable by `chunk_size`.
- Use efficient whitelist lookup.
- Use compact data structures for barcode correction.
- Use multiprocessing where safe and useful.
- Avoid storing huge uncompressed intermediate files.
- Write large detail files as `.tsv.gz/.parquet`.
- Use atomic writes for all important outputs.
- Allow users to disable large trace files if necessary, but default should preserve correction trace.

Suggested Python libraries:

- `click` or `typer` for CLI
- `pyyaml` for config
- `pydantic` for config validation
- `pandas` for summary-level tables only
- `polars` or streaming aggregation for large tables
- `gzip`, `xopen`, or `isal` for fast compressed IO
- `cutadapt` as an external dependency or Python module
- `umi_tools` for UMI correction
- `openpyxl` or `xlsxwriter` for Excel report
- `plotly` or `jinja2` for HTML report
- `pytest` for tests

Important:

- Do not use pandas to load huge read-level files into memory.
- For large detail files, use streaming aggregation or chunked processing.
- Count matrix generation should be memory-aware.
- If a high-cardinality aggregation may exceed memory, use SQLite, DuckDB, or disk-backed chunk aggregation.

---

## 15. Recommended Internal Architecture

Suggested repository structure:

```
tagforge/
├── pyproject.toml
├── README.md
├── LICENSE
├── .gitignore
├── configs/
│   └── config.example.yaml
├── examples/
│   ├── small_fastq/
│   ├── FB_info.tsv
│   └── tutorial.md
├── src/
│   └── tagforge/
│       ├── __init__.py
│       ├── cli.py
│       ├── config.py
│       ├── logging_utils.py
│       ├── checkpoint.py
│       ├── io.py
│       ├── fastq.py
│       ├── extract.py
│       ├── cutadapt_runner.py
│       ├── barcode_correct.py
│       ├── umi_correct.py
│       ├── aggregate.py
│       ├── matrix.py
│       ├── downsample.py
│       ├── report_excel.py
│       ├── report_html.py
│       ├── slurm.py
│       └── utils.py
└── tests/
    ├── test_config.py
    ├── test_extract.py
    ├── test_barcode_correct.py
    ├── test_umi_correct.py
    ├── test_downsample.py
    └── test_cli.py
```

---

## 16. Implementation Details

### 16.1 FASTQ Parser

Implement a fast paired-end FASTQ reader that streams records from gzipped files.

Each yielded record should contain:

```python
read_id
r1_seq
r1_qual
r2_seq
r2_qual
```

The parser should validate that read names match between R1 and R2.

### 16.2 Segment Extraction

Create a generic segment extractor.

Input:

```python
read_sequence
segment_config
```

Output:

```python
SegmentExtractionResult(
    segment_name,
    raw_sequence,
    success,
    failure_reason
)
```

For fixed extraction:

```python
raw_sequence = read_sequence[start:start + length]
```

For linker extraction:

- Use cutadapt if possible.
- Also provide a fallback Python implementation for simple exact-linker extraction.
- Linker mismatch handling should be configurable.

### 16.3 Barcode Correction

Implement correction as reusable logic.

For each raw barcode segment:

1. Try exact whitelist match.
2. Try forward shift correction up to `max_shift`.
3. Try mismatch correction up to `max_mismatch`.
4. Try shift plus mismatch correction.
5. If multiple whitelist candidates tie, mark as ambiguous and failed unless a deterministic tie-breaking rule is explicitly configured.

Correction result should include:

```python
raw_sequence
corrected_sequence
status
shift_distance
mismatch_distance
correction_type
ambiguous
```

For efficient whitelist correction:

- Precompute whitelist set for exact lookup.
- For max mismatch 1 or 2, precompute variant index if feasible.
- Use BK-tree or other approximate matching index for larger whitelist or mismatch thresholds.
- Cache correction results for repeated raw sequences.

### 16.4 UMI Correction

UMI correction should be done within each Barcode1-Barcode2 group.

Input records:

```
barcode1 barcode2 raw_umi read_count
```

Output:

```
barcode1 barcode2 corrected_umi reads_count raw_umi_count
```

Use UMI-tools compatible methods, preferably `directional`.

If using external UMI-tools command-line tools is difficult, implement directional UMI deduplication internally and document compatibility.

### 16.5 Count Matrix

Generate count matrix from molecule detail.

Rows:

```
barcode1
```

Columns:

```
barcode2_name
```

Values:

```
number of corrected UMIs
```

Use sparse or chunked aggregation internally.

### 16.6 Downsampling

Perform downsampling at read level or molecule-support level with reproducible random seed.

For each ratio:

1. Sample reads or read-support events.
2. Recompute UMI molecule counts.
3. Calculate:
    - Sequencing Saturation
    - Duplication Ratio
    - UMI detected once
    - UMI Types
4. Save metrics.
5. Identify the ratio with maximum Sequencing Saturation.
6. Generate the optimal saturation matrix and optimal molecule detail file.

Downsampling must be reproducible by sample name, ratio, repeat index, and random seed.

---

## 17. Logging

Each sample should have a log file:

```
logs/{sample}.pipeline.log
```

Log format should include:

```
timestamp
level
sample
step
message
elapsed_time
```

Log major events:

- config validation
- input file detection
- step start
- step finish
- number of reads processed
- valid read counts
- correction statistics
- output file paths
- skipped steps due to checkpoint
- warnings and errors

---

## 18. Error Handling

The pipeline should fail clearly and early when:

- FASTQ files are missing
- R1/R2 read IDs do not match
- config file is invalid
- whitelist files are missing
- FB_info.tsv is missing required columns
- duplicate FB sequences exist
- duplicate antibody names exist, unless allowed by config
- segment coordinates exceed read length
- cutadapt is required but unavailable
- output directory is not writable

Error messages should be human-readable and include the sample name and pipeline step.

---

## 19. README Requirements

Write a mature `README.md` including:

1. Project overview
2. Installation
3. Dependencies
4. Input file structure
5. Configuration file explanation
6. Library structure examples
7. Running one sample
8. Running all samples
9. Slurm usage
10. Output files
11. Report explanation
12. Resume/restart behavior
13. Performance tuning
14. Troubleshooting
15. FAQ
16. Minimal tutorial with example data

---

## 20. Installation

Support installation with `pip`.

Example:

```bash
git clone https://github.com/yourname/tagforge.git
cd tagforge
pip install -e .
```

Also support dependency installation:

```bash
pip install -e ".[dev]"
```

External tools:

```bash
cutadapt
umi_tools
```

The program should check whether these tools are available and provide clear installation hints if missing.

---

## 21. Tests

Add unit tests and small integration tests.

Minimum test coverage:

- YAML config validation
- FASTQ paired reader
- fixed-position extraction
- linker-based extraction
- barcode exact correction
- barcode mismatch correction
- barcode shift correction
- barcode shift plus mismatch correction
- ambiguous correction handling
- multi-segment barcode concatenation
- FB annotation mapping
- UMI deduplication
- count matrix generation
- downsampling metric calculation
- CLI smoke tests

Provide a tiny example dataset under `examples/` for testing and tutorial.

---

## 22. Acceptance Criteria

The implementation is acceptable when:

1. The package can be installed with `pip install -e .`.
2. `tagforge validate-config` validates the example config.
3. `tagforge run` can process the example FASTQ files end-to-end.
4. The pipeline supports paired-end gzipped FASTQ input.
5. Barcode1, Barcode2, and UMI can be extracted from configurable read segments.
6. Both linker-based and fixed-position extraction are supported.
7. Multi-segment Barcode1, Barcode2, and UMI are supported.
8. Barcode correction supports exact match, mismatch correction, shift correction, and shift plus mismatch correction.
9. Raw and corrected barcode information are retained.
10. UMI correction/deduplication works within each Barcode1-Barcode2 group.
11. Raw count matrix is generated using antibody names from `FB_info.tsv`.
12. Molecule detail file is generated.
13. Downsampling metrics are generated.
14. Optimal saturation matrix and detail files are generated.
15. Excel and HTML reports are generated.
16. Logs are written for each sample.
17. Checkpoints allow restart without rerunning completed steps.
18. Slurm scripts can be generated for all samples.
19. Large files are written compressed.
20. The implementation avoids loading the entire FASTQ dataset into memory.

---

## 23. Development Strategy for Codex

Please implement this project incrementally.

Recommended order:

1. Create repository structure and `pyproject.toml`.
2. Implement config schema and validation.
3. Implement logging and checkpoint utilities.
4. Implement paired FASTQ streaming reader.
5. Implement fixed-position extraction.
6. Implement simple linker extraction and cutadapt wrapper.
7. Implement whitelist loading and barcode correction.
8. Implement read-level extraction and correction pipeline.
9. Implement UMI aggregation and UMI correction.
10. Implement count matrix generation.
11. Implement downsampling and saturation metrics.
12. Implement Excel report.
13. Implement HTML report.
14. Implement Slurm script generation.
15. Add example data and tutorial.
16. Add tests.
17. Polish README and error messages.

Do not create a toy script only. Build a maintainable Python package with clear modules, tests, CLI, configuration validation, logs, checkpoints, and documentation.
