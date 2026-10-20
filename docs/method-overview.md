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

PI-seq multi-PB filtering and PB–CB mapping are dedicated experimental extensions.
PI-seq compares PB read support within `(feature, corrected UMI)` and, for multi-PB
groups, keeps a winner only above the configured dominance threshold. PB–CB
mapping selects a dominant CB for each PB using read support after the configured
observed-CB handling. These rules sit on top of the common processing core.

## Implementation references

- [Structure and configuration](../src/tagforge/config.py)
- [Extraction](../src/tagforge/extract.py) and [barcode correction](../src/tagforge/barcode_correct.py)
- [UMI clustering and aggregation](../src/tagforge/umi_correct.py)
- [Sparse counts](../src/tagforge/sparse_matrix.py) and [downsampling](../src/tagforge/downsample.py)
- [Checkpoints and dependency validation](../src/tagforge/checkpoint.py)
- [Slurm integration](../src/tagforge/slurm.py)

Linker matching uses Cutadapt and exact-match paths; UMI clustering uses
UMI-tools. The common workflow, data organization and execution logic are
implemented in TagForge.
