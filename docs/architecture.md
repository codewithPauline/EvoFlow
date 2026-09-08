# EvoFlow architecture

EvoFlow is organized as modular scientific software rather than a collection of analysis scripts. Configuration, genomic input access, metadata handling, validation, scientific engines, visualization, orchestration, and reporting are separated so that each layer can evolve without tightly coupling the rest of the workflow.

## Package layers

- `evoflow.config` — project configuration and YAML serialization.
- `evoflow.io` — VCF access, genotype extraction, metadata access, and input validation.
- `evoflow.modules` — scientific analysis engines and module-specific plotting code.
- `evoflow.core` — reserved for dependency resolution, workflow execution, manifests, provenance, checksums, restartability, and run lifecycle management.
- `evoflow.cli` — user-facing command-line entry points and current input-routing policy.

## Current analysis graph

```text
                        evoflow.yaml
                             |
                             v
                      EvoFlowConfig
                             |
                             v
                   strict input validation
                    /                \
                   v                  v
             VCF / VCF.gz       metadata CSV
                   |                  |
                   +--------+---------+
                            |
                            v
                    +---------------+
                    |   QC engine   |
                    +---------------+
                       |         |
                       |         +-------------------------------+
                       v                                         |
                  filtered.vcf                                  |
                    /       \                                    |
                   /         \                                   |
                  v           v                                  v
          +-------------+  +----------------+          QC tables / figures
          | LD pruning  |  | population     |
          +-------------+  | statistics     |
                  |        +-------+--------+
                  v                |       |
             ld_pruned.vcf         v       v
                  |           diversity    FST
                  v                |       |
          +---------------+        v       v
          |      PCA      |      tables + figures
          +---------------+
                  |
                  v
         scores/loadings/figures
```

The graph is intentionally not a single linear chain. Different statistical questions require different preprocessing choices:

- **PCA** prefers LD-pruned variants because correlated marker blocks can dominate ordination.
- **Diversity and FST** prefer the QC-filtered variant set and do not automatically inherit LD thinning.

This input routing is explicit in the CLI rather than hidden inside the scientific functions.

## Input layer

### `evoflow.io.vcf`

Responsibilities:

- stream uncompressed or gzip/bgzip-compatible text VCFs;
- read sample IDs from the `#CHROM` header;
- reject unsupported native BCF input with an explicit error;
- parse usable diploid biallelic genotype dosages; and
- extract per-record dosage vectors for downstream modules.

Shared genotype parsing prevents PCA, LD, diversity, and FST from independently implementing subtly different definitions of a usable genotype.

### `evoflow.io.metadata`

Responsibilities:

- read metadata CSV files; and
- map population labels to VCF sample indices while preserving VCF sample order.

### `evoflow.io.validation`

Responsibilities:

- verify configured input files exist;
- require a metadata `sample` column;
- require unique, non-empty sample IDs; and
- require exact VCF/metadata sample concordance.

## QC engine

`evoflow.modules.qc` performs one streaming VCF pass to calculate sample-, site-, and run-level metrics while applying configured variant thresholds. It writes:

- run JSON summary;
- sample QC table;
- site QC table;
- threshold-filtered VCF; and
- compact histogram state used for large per-site distributions.

Global depth and genotype-quality summaries use running sums/counts. MAF and missingness plots are generated from fixed-size histogram bins rather than complete in-memory vectors.

`evoflow.modules.qc_plots` is responsible only for visualization and consumes already-computed QC state.

## LD pruning engine

`evoflow.modules.ld` performs a greedy physical sliding-window pruning algorithm on informative diploid biallelic SNPs.

For each incoming SNP:

1. retained SNPs outside the configured physical window are removed from the active comparison set;
2. genotype-dosage `r²` is calculated against retained SNPs using pairwise-complete samples;
3. if any retained SNP reaches the configured threshold, the current SNP is removed; otherwise it is retained;
4. each removal records the blocking SNP and observed `r²`.

Memory use is controlled by the physical active window rather than total chromosome length.

The engine requires VCF positions to be sorted within chromosomes and writes pruning parameters into the output VCF header.

## PCA engine

`evoflow.modules.pca` operates on informative diploid biallelic SNPs.

For SNP allele frequency `p`, dosage `g` is standardized as:

```text
(g - 2p) / sqrt(2p(1-p))
```

Missing calls are mean-imputed to `2p`, which becomes zero after standardization.

Instead of materializing the full sample × SNP matrix `X`, the first VCF pass accumulates the sample Gram matrix:

```text
G = X X^T
```

This requires memory proportional primarily to the square of sample count rather than sample count × SNP count. Symmetric eigendecomposition of `G` produces sample principal-component scores. A second streaming pass projects SNP vectors onto the sample eigenvectors to calculate per-SNP loadings.

`evoflow.modules.pca_plots` generates the explained-variance and PC1–PC2 figures from machine-readable PCA tables.

## Population diversity engine

`evoflow.modules.diversity` groups samples using metadata-defined populations and streams diploid biallelic SNP records.

For each population/site combination with called genotypes it calculates:

- call rate;
- alternate-allele frequency;
- MAF;
- observed heterozygosity; and
- expected heterozygosity `2p(1-p)`.

Population summaries are arithmetic means across analyzed VCF variant records plus counts of observed and polymorphic sites.

This is deliberately not labeled genome-wide nucleotide diversity because the current VCF-based implementation does not model invariant callable bases or a genome-wide callable denominator.

`evoflow.modules.diversity_plots` renders population heterozygosity and population call-rate figures.

## Hudson FST engine

`evoflow.modules.fst` estimates pairwise Hudson FST for all metadata-defined population pairs.

At each usable biallelic SNP, allele counts are used to calculate:

- within-population pairwise diversity for each population;
- mean within-population diversity;
- between-population divergence;
- site numerator = between divergence - mean within diversity; and
- site denominator = between divergence.

The multi-site estimate is formed from summed components:

```text
FST = sum(numerators) / sum(denominators)
```

Per-site negative FST estimates are retained rather than silently clamped to zero. Pairwise outputs include the number of denominator-positive sites contributing to each aggregate estimate.

`evoflow.modules.fst_plots` consumes the symmetric FST matrix and renders a labeled heatmap.

## CLI input-routing policy

The current CLI keeps preprocessing decisions visible:

### `evoflow ld-prune`

```text
QC filtered.vcf if available -> otherwise configured VCF
```

### `evoflow pca`

```text
LD-pruned VCF if available -> QC-filtered VCF if available -> configured VCF
```

### `evoflow diversity` and `evoflow fst`

```text
QC-filtered VCF if available -> otherwise configured VCF
```

`--use-raw` explicitly bypasses generated intermediates for the relevant command.

## Deterministic output layout

Current scientific modules write:

```text
evoflow-results/
├── qc/
├── ld/
├── pca/
├── diversity/
└── fst/
```

Each directory contains machine-readable result files and, where appropriate, publication-oriented PNG/PDF figures.

## Scientific module contract

New modules should define:

1. scientific question and estimator;
2. required inputs;
3. model assumptions;
4. validation and failure rules;
5. configurable parameters and defaults;
6. deterministic output locations;
7. machine-readable tables/summaries;
8. scientific figures where appropriate;
9. provenance needed to reproduce the result; and
10. tests against known expected values and failure modes.

## Current architectural boundaries

The following are intentionally not hidden or treated as completed infrastructure:

- native BCF backend;
- genotype-level GQ filtering;
- Hardy-Weinberg filtering;
- generalized ploidy handling;
- additional diversity estimators using invariant/callable bases;
- FST uncertainty via blocks/jackknife;
- ancestry/structure inference;
- spatial, selection, and GEA engines;
- centralized workflow executor;
- run manifests/checksums/version capture; and
- automated scientific reporting.

## Next architectural step

The next major capability is the population-structure layer, followed by a formal orchestration layer that can resolve module dependencies and reproduce a full configured run from one command. The design goal is to preserve the current module independence while adding deterministic execution and provenance on top of it.
