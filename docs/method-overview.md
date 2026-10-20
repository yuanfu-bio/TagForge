# TagForge: method overview and application scope

English | [简体中文](method-overview.zh-CN.md)

TagForge decodes known barcode designs into molecule-level evidence and counts.
Its shared model is **identity barcode × feature barcode × UMI**. Identity
barcodes distinguish the units being counted, such as cells or samples; feature
barcodes identify labels. The counting scope uses the corrected identity
barcode and the annotated feature name.

[Workflow figure](../figures/TagForge_core_principle.en.page-1.png) ·
[Core mechanisms](../figures/TagForge_core_principle.en.page-2.png) ·
[Two-page PDF](../figures/TagForge_core_principle.en.pdf) ·
[Self-contained HTML](../figures/TagForge_core_principle.en.html)

## From reads to molecules

1. **Describe the library.** Configure the role, read, length, position and/or
   linkers of each segment. Barcode and UMI segments can span R1 and R2.
2. **Extract segments.** Use fixed positions or linker matching. With
   `methods: [linker, fixed]`, a failed linker extraction falls back to fixed
   coordinates on the original read. A later whitelist failure does not trigger
   another extraction attempt.
3. **Correct identities.** Correct barcode segments against their own
   whitelists, then concatenate them. This uses segment-level design information
   without enumerating a full whitelist of every barcode combination.
4. **Define molecules.** Aggregate read support for each raw UMI within a
   corrected identity–feature pair, then use UMI-tools clustering. Each corrected
   UMI cluster contributes one molecule to the count matrix.
5. **Preserve evidence and report.** Molecule records retain the representative
   UMI, supporting read count and number of distinct raw UMIs. Sparse count
   tables, downsampling diagnostics and sample/batch reports use these results.

## Core decisions

### Barcode correction

An exact whitelist hit at the original window is accepted immediately. Otherwise,
TagForge searches the configured shifts and Hamming mismatch limits. Candidates
are ranked lexicographically by:

`(I[shift != 0] + I[mismatches > 0], abs(shift) + mismatches, mismatches, abs(shift))`

A unique best barcode identity is accepted. Different identities tied at the
best score are rejected as ambiguous; no legal candidate also causes rejection.
Multiple positions resolving to the same best identity do not create an identity
tie. Shifted-window search models positional offsets and substitutions, not
arbitrary indel alignment.

### UMI clustering

The default is UMI-tools `directional` with distance 1. A candidate directed
edge `u → v` requires distance at most 1 and `n(u) >= 2*n(v) - 1`, where `n` is
read support. UMI-tools determines the final clusters and representative UMIs.

In a group containing only two UMIs one base apart, support 8:1 merges into one
molecule, while support 8:7 remains two molecules. The figure compares independent
groups; all layouts, sequences and counts are illustrative.

The matrix counts clusters rather than supporting reads. For example, two
molecule records supported by 9 and 3 reads contribute **2 molecules**, not 12.

### Downsampling diagnostics

For a fixed full-depth molecule supported by `n` reads, sample support using
`Binomial(n, p)` at ratio `p`. Current diagnostics use one seeded draw per ratio
and provide analytical expectations for read, molecule and singleton counts.
They do not rerun extraction, UMI clustering or experiment-specific filtering at
each depth.

- Sequencing saturation: `1 - molecules / reads`.
- Non-singleton molecule fraction: `molecules with >=2 reads / molecules`.

The selected diagnostic ratio maximizes the observed non-singleton fraction,
with a lower ratio winning ties. Primary analysis uses full-depth results; the
selection is not a sequencing-depth recommendation.

## Engineering capabilities

- **Scale and optimize:** block input, Parquet storage, DuckDB aggregation with
  disk spill, sparse output, exact-match fast paths, neighbor lookup and cached
  correction.
- **Parallelize independent work:** each identity–feature pair remains in one
  hash bucket. Independent UMI components are processed without splitting their
  clustering boundaries.
- **Recover and validate:** reuse validated stages; extraction and correction
  can resume committed batches. Input/configuration fingerprints, upstream state
  and output checks govern reuse. FASTQ identity uses metadata by default, with
  full SHA256 available through configuration.
- **Operate in batches:** Slurm arrays, progress records, stage logs, segment QC,
  and sample/batch HTML and Excel reports.

These describe implemented designs, not a universal speed or memory advantage
over other tools. Resource needs still depend on data distribution and the
largest scopes processed by workers.

## Application scope

| Application | Counted relationship | Scope |
| --- | --- | --- |
| Antibody tags and combinatorial barcodes | Identity barcode × antibody/feature label | Primary application of the common pipeline. |
| Sample Tags | Cell × sample tag | Adaptation direction; downstream sample assignment remains a separate analysis. |
| CRISPR Guide Capture | Cell × sgRNA | Adaptation direction when guide identity is directly extractable; downstream perturbation assignment is separate. |
| Antigen Capture | Cell × antigen label | Adaptation direction; binding-specificity interpretation is separate. |
| Promoter/reporter barcodes | Identity barcode × design-linked tag | Adaptation direction for UMI-bearing designs; activity or expression interpretation is experiment-specific. |

Extensions require paired-end input supported by the pipeline, extractable
segments and UMIs, determinable feature identities, and whitelists wherever
whitelist correction is enabled. Validate each layout and dataset before use.
The shared count engine does not imply validated support for every named assay.

## Two dedicated libraries for PI-seq experiments

[Page 3](../figures/TagForge_core_principle.en.page-3.png) compares these extensions.
Both reuse the common core and apply experiment-specific assignment rules after
UMI clustering. PB is the identity barcode shared by the two libraries, FB is an
antibody feature barcode, and CB is a cell barcode.

| Library | Question | Assignment scope | Main result |
| --- | --- | --- | --- |
| PI-seq | Which PB carries this antibody-tagged molecule? | Compare PBs within `(FB_name, corrected_umi)` | Filtered PB × FB molecule details and counts. |
| PB–CB | Which cell barcode is supported by this PB? | Compare CB read support for each PB | PB → CB assignments and PB counts per cell. |

### PI-seq: resolve molecular PB identity

The development fixture has five 8-nt PB segments in R1, corrected individually
and concatenated, and `U1 (8 nt)–FB (12 nt)–U2 (8 nt)` in R2. UMI segments are
concatenated in configured order. Linker localization and fixed-position fallback
reuse the common core; the layout is configurable rather than universal.

After UMI clustering within `PB × FB_name`, the filter groups records by
`(FB_name, corrected_umi)` across PBs. Single-PB groups are retained directly.
For multi-PB groups, it keeps the highest-read PB only when its fraction of the
group's supporting reads is **strictly greater than**
`pi_seq.dominance_threshold` (default 0.8); otherwise the entire group is excluded.
Thus 9:1 support retains the 90% PB, while 8:2 support fails the default threshold.
Equal-read candidates are ordered by PB sequence before threshold evaluation.
Read support determines assignment; retained UMI molecules determine final counts.

The main output is `03_molecules/{sample}.molecule_detail.rmMP.parquet`, used for
PI-seq counts and downstream per-PB analysis. `pi_seq_qc.tsv` records multi-PB,
dominant-assignment and removal metrics.

### PB–CB: map PBs to cells

This library reads `PB × CB × UMI`. CB is retained as observed sequence:
`barcode2.sequence_only: true`, CB whitelist correction disabled, and no CB
whitelist. PB segments may still use configured whitelist correction.

Optional conservative observed-CB correction is enabled by default. Across the
sample, candidate parent CBs contain no N and have at least 5 supporting reads.
A child at Hamming distance 1, or with one matching N, is aliased only when it has
a unique qualifying parent with strictly higher read support and at least 10 times
its support. Ambiguous or ineligible sequences are not forcibly merged. These
settings are configurable under `pb_cb`.

After alias handling, records with the same `(PB, CB, corrected_umi)` are combined
again to prevent duplicate molecule counting after CB merging. For each PB,
TagForge aggregates molecule counts and supporting reads by CB. The leading CB's
**read support fraction** must be **strictly greater than**
`pb_cb.dominance_threshold` (default 0.8). Equal-read candidates are ordered by CB
sequence before threshold evaluation. At the default threshold, 9:1 support is
assigned and 8:2 support is unassigned.

Outputs in `04_matrix/` include raw pair support (`pb_cb_mapping.tsv.gz`), retained
assignments with molecule/read counts and dominance fractions (`pb_cb_map.tsv.gz`),
PB counts per CB and their distribution (`cb_pb_counts.tsv.gz`,
`cb_pb_count_distribution.tsv`), and CB aliases (`cb_observed_correction.tsv.gz`).

### Link the outputs through shared PBs

PI-seq supplies `PB × antibody feature` molecular evidence, while PB–CB supplies
`PB → cell` assignments. Outputs can be joined on shared PBs with consistent
identity definitions as input to cell-level feature integration. The join and
biological interpretation are separate downstream analyses. Multiple PBs may map
to the same CB; each retained PB has one CB assignment.

## Implementation references

- [Structure and configuration](../src/tagforge/config.py)
- [Extraction](../src/tagforge/extract.py) and [barcode correction](../src/tagforge/barcode_correct.py)
- [UMI clustering and aggregation](../src/tagforge/umi_correct.py)
- [PI-seq assignment filtering](../src/tagforge/pi_seq.py) and [PB–CB mapping](../src/tagforge/matrix.py)
- [Sparse counts](../src/tagforge/sparse_matrix.py) and [downsampling](../src/tagforge/downsample.py)
- [Checkpoints and dependency validation](../src/tagforge/checkpoint.py)
- [Slurm integration](../src/tagforge/slurm.py)

Linker matching uses Cutadapt and exact-match paths; UMI clustering uses
UMI-tools. The common workflow, data organization and execution logic are
implemented in TagForge.
