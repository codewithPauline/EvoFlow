# EvoFlow

[![CI](https://github.com/codewithPauline/EvoFlow/actions/workflows/ci.yml/badge.svg)](https://github.com/codewithPauline/EvoFlow/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Development Status](https://img.shields.io/badge/status-active%20development-orange.svg)](#development-status)

> **From variants to evolutionary insight.**

**EvoFlow** is an open-source, species-agnostic command-line platform for reproducible population and landscape genomics. It is being built to connect genomic variants, sample metadata, population-genetic analyses, spatial information, scientific figures, and reproducible reporting through one transparent workflow.

**Author and lead developer:** [Pauline Owusu-Ansah](https://github.com/codewithPauline) (`@codewithPauline`)

---

## Why EvoFlow?

Population-genomic studies rarely consist of one analysis. A typical project may require variant and sample quality control, filtering, linkage-disequilibrium pruning, PCA, ancestry inference, diversity summaries, population differentiation, spatial analysis, selection scans, genotype-environment association, visualization, and reporting.

Those steps are often scattered across different command-line tools, R scripts, Python notebooks, file formats, and software environments. That fragmentation can make an analysis difficult to reproduce, audit, extend, transfer to another organism, or explain months later.

**EvoFlow is being developed as a transparent orchestration layer for that workflow.**

The project is guided by five principles:

1. **Reproducibility first** — inputs, parameters, software versions, intermediate files, and outputs should remain traceable.
2. **Species agnostic** — development examples may use salamander data, but the software is designed for population-genomic data from any organism that satisfies the current model assumptions.
3. **Modular analyses** — researchers should be able to run only the analyses they need.
4. **Transparent methods** — statistical assumptions and filtering decisions should remain visible rather than hidden behind an opaque interface.
5. **Scientific outputs** — the workflow should produce interpretable tables, figures, and machine-readable summaries, not merely successful command execution.

---

## Current workflow

```text
                 VCF / VCF.gz
                      +
                sample metadata
                      |
                      v
             +------------------+
             | Input validation |
             +------------------+
                      |
                      v
             +------------------+
             | Variant/sample QC|
             +------------------+
                      |
                      v
                 filtered.vcf
                 /          \
                /            \
               v              v
        +-------------+   +------------------+
        | LD pruning  |   | Diversity / FST  |
        +-------------+   +------------------+
               |
               v
          ld_pruned.vcf
               |
               v
        +-------------+
        |     PCA     |
        +-------------+
               |
               v
       tables + PNG/PDF figures
               |
               v
       future workflow modules
 structure -> spatial -> selection -> GEA -> report
```

A deliberate design decision is that **PCA prefers LD-pruned variants**, while **diversity and FST prefer QC-filtered variants without LD thinning**. LD pruning is useful for ordination/structure analyses because dense linked blocks can dominate those analyses; discarding linked loci is not automatically desirable for the currently implemented diversity and FST summaries.

---

## Development status

EvoFlow is in **active early development**, but it is no longer only a package scaffold. The repository now contains multiple functioning and tested population-genomic engines.

| Component | Purpose | Status |
|---|---|---|
| Python package | Installable `evoflow` package | ✅ Available |
| Command-line interface | User-facing EvoFlow commands | ✅ Available |
| YAML configuration | Reproducible project configuration | ✅ Available |
| VCF / metadata validation | Confirms input integrity and sample concordance | ✅ Available |
| Native VCF QC | Streaming sample/variant quality metrics | ✅ Available |
| QC filtering | MAF/missingness decisions + filtered VCF | ✅ Available |
| QC figures | PNG/PDF QC diagnostics | ✅ Available |
| DP / GQ summaries | Uses FORMAT `DP` and `GQ` when present | ✅ Available |
| LD pruning | Sliding-window genotype-dosage r² pruning | ✅ Available |
| PCA | Allele-frequency-standardized streaming PCA | ✅ Available |
| PCA loadings / figures | Scores, variance, loadings, scree, PC1–PC2 | ✅ Available |
| Population diversity | Population Ho, He, MAF, polymorphism, call-rate summaries | ✅ Available |
| Pairwise FST | Hudson FST with per-site and aggregate outputs | ✅ Available |
| FST heatmap | Population differentiation matrix in PNG/PDF | ✅ Available |
| Continuous integration | Install, lint, test on Python 3.10–3.12 | ✅ Passing |
| Population structure | Ancestry / clustering workflow | 🧭 Planned |
| Spatial genomics | IBD and spatial population-genomic analyses | 🧭 Planned |
| Selection scans | Outlier / selection analyses | 🧭 Planned |
| GEA | Genotype-environment association | 🧭 Planned |
| Workflow executor | Ordered module execution + provenance | 🧭 Planned |
| Automated report | Integrated methods, tables, figures, provenance | 🧭 Planned |

---

## Installation

EvoFlow is currently installed from source.

### Requirements

- Python 3.10 or newer
- Git

### macOS / Linux

```bash
git clone https://github.com/codewithPauline/EvoFlow.git
cd EvoFlow

python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

### Windows PowerShell

```powershell
git clone https://github.com/codewithPauline/EvoFlow.git
cd EvoFlow

python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

Confirm the installation:

```bash
evoflow version
```

---

## Quick start

### 1. Create a project configuration

```bash
evoflow init \
  --project salamander-demo \
  --vcf data/variants.vcf.gz \
  --metadata data/samples.csv
```

This creates `evoflow.yaml`.

### 2. Validate inputs

```bash
evoflow validate evoflow.yaml
```

### 3. Inspect the configured module plan

```bash
evoflow plan evoflow.yaml
```

### 4. Run quality control

```bash
evoflow qc evoflow.yaml --min-maf 0.05 --max-missing 0.20
```

### 5. LD-prune QC-filtered SNPs

```bash
evoflow ld-prune evoflow.yaml --r2 0.20 --window-bp 50000
```

### 6. Run PCA

```bash
evoflow pca evoflow.yaml --n-components 10
```

When available, PCA automatically prefers:

```text
LD-pruned VCF -> QC-filtered VCF -> configured/raw VCF
```

Use `--use-raw` to explicitly bypass generated intermediate VCFs.

### 7. Summarize population diversity

```bash
evoflow diversity evoflow.yaml
```

### 8. Estimate pairwise FST

```bash
evoflow fst evoflow.yaml
```

Both diversity and FST use `qc/filtered.vcf` when it exists, otherwise the configured VCF. Use `--use-raw` to force the configured VCF.

---

## Project configuration

Example:

```yaml
project: salamander-demo
vcf: data/variants.vcf.gz
metadata: data/samples.csv
output_dir: evoflow-results
modules:
  - qc
  - ld
  - pca
  - structure
  - diversity
  - fst
  - spatial
  - selection
  - gea
  - report
```

The module list describes the intended analysis plan. The currently executable scientific modules are `qc`, `ld`, `pca`, `diversity`, and `fst`.

---

# Input data

## Variant file

The native text parser currently supports:

- uncompressed VCF (`.vcf`)
- gzip/bgzip-compatible VCF text (`.vcf.gz`, `.vcf.bgz`, `.gz`, `.bgz`)

**Native BCF parsing is not yet implemented.** Convert BCF to VCF/VCF.gz before using the current native EvoFlow engines.

Current population-genomic engines are centered on diploid genotype (`GT`) data. Some QC metrics also use `DP` and `GQ` when those FORMAT fields are present.

## Sample metadata

The minimum metadata CSV contains a `sample` column:

```csv
sample,population,site,latitude,longitude
sample_01,OH,Site_A,39.5100,-84.7300
sample_02,IN,Site_B,39.1600,-86.5200
sample_03,KY,Site_C,38.0400,-84.5000
```

`population` is required by the current diversity and FST commands unless another metadata column is selected with `--population-column`.

---

# Strict input validation

Before scientific modules run, EvoFlow checks that:

- the configured variant file exists;
- the metadata CSV exists;
- metadata contains a `sample` column;
- metadata sample IDs are non-empty and unique;
- VCF sample IDs are unique; and
- VCF and metadata contain **exactly the same sample IDs**.

This strict behavior is intentional. Silent sample mismatches can invalidate downstream population-genomic analyses, so EvoFlow fails early rather than guessing.

---

# Quality-control engine

## Streaming approach

The QC engine streams VCF records rather than loading the entire variant file into memory. It maintains sample-level accumulators and fixed-size distribution bins for large per-site summaries.

## Metrics

The engine currently calculates:

- sample count;
- variant count;
- SNP count;
- biallelic SNP count;
- multiallelic variant count;
- transitions and transversions;
- Ti/Tv ratio when defined;
- called and missing genotype counts;
- overall genotype call rate;
- per-sample call rate and missingness;
- per-sample observed heterozygosity;
- per-site call rate and missingness;
- alternate-allele frequency;
- biallelic minor-allele frequency;
- per-sample/per-site mean depth when FORMAT `DP` is available; and
- per-sample/per-site mean genotype quality when FORMAT `GQ` is available.

## Filtering

```bash
evoflow qc evoflow.yaml --min-maf 0.05 --max-missing 0.20
```

Each record receives a `passes_qc` decision. Passing records are written to `filtered.vcf`, and thresholds are recorded in the VCF header:

```text
##evoflow_qc_min_maf=0.05
##evoflow_qc_max_missing=0.2
```

For multiallelic sites, EvoFlow currently leaves the biallelic-style MAF undefined rather than imposing an ambiguous definition.

## QC outputs

```text
evoflow-results/qc/
├── qc_summary.json
├── sample_qc.csv
├── variant_qc.csv
├── filtered.vcf
├── sample_call_rate.png
├── sample_call_rate.pdf
├── sample_heterozygosity.png
├── sample_heterozygosity.pdf
├── sample_depth.png
├── sample_depth.pdf
├── sample_gq.png
├── sample_gq.pdf
├── maf_distribution.png
├── maf_distribution.pdf
├── variant_missingness.png
└── variant_missingness.pdf
```

---

# Linkage-disequilibrium pruning

## Why it exists

Dense clusters of strongly linked SNPs can disproportionately influence PCA and ancestry-style analyses. EvoFlow therefore provides an explicit LD-pruning preprocessing stage instead of silently pruning inside PCA.

## Current method

```bash
evoflow ld-prune evoflow.yaml --r2 0.20 --window-bp 50000
```

The current implementation:

- uses informative diploid biallelic SNPs;
- scans a physical window along a position-sorted VCF;
- computes pairwise genotype-dosage `r²` using samples called at both SNPs;
- greedily retains the first acceptable SNP and removes later SNPs whose `r²` meets or exceeds the threshold with a retained SNP inside the window; and
- records the SNP responsible for each removal.

The main parameters are:

- `--r2` — pairwise r² threshold, default `0.20`;
- `--window-bp` — physical window, default `50,000 bp`;
- `--min-overlap` — minimum number of jointly called samples, default `3`.

The input VCF must be position-sorted within chromosomes.

## LD outputs

```text
evoflow-results/ld/
├── ld_pruned.vcf
├── ld_removed.csv
└── ld_summary.json
```

Pruning parameters are written into the output VCF header for provenance.

---

# PCA engine

## Statistical approach

PCA currently uses **informative diploid biallelic SNPs**. For a SNP with alternate-allele frequency `p`, genotype dosage `g` is standardized as:

```text
(g - 2p) / sqrt(2p(1-p))
```

Missing genotypes are mean-imputed to `2p`, which becomes zero after standardization. Monomorphic sites, indels, multiallelic sites, and records without usable diploid GT calls are excluded.

Instead of materializing the entire sample × SNP matrix, EvoFlow streams standardized variant vectors and accumulates the sample Gram matrix:

```text
X Xᵀ
```

Principal components are obtained through symmetric eigendecomposition. A second VCF pass calculates per-SNP loadings. Memory usage is therefore driven primarily by sample count rather than total SNP count.

## Input priority

```text
ld/ld_pruned.vcf exists  -> use LD-pruned variants
else qc/filtered.vcf     -> use QC-filtered variants
else                     -> use configured VCF
```

Use `--use-raw` to override this behavior.

## PCA outputs

```text
evoflow-results/pca/
├── pca_summary.json
├── pca_scores.csv
├── pca_variance.csv
├── pca_loadings.csv
├── pca_scree.png
├── pca_scree.pdf
├── pca_pc1_pc2.png
└── pca_pc1_pc2.pdf
```

`pca_scores.csv` preserves available metadata columns. If metadata contains `population`, the PC1–PC2 figure automatically groups points by population.

---

# Population diversity

## Scope

```bash
evoflow diversity evoflow.yaml
```

The current diversity engine reports **mean diversity summaries across analyzed diploid biallelic SNP records** for metadata-defined populations.

It calculates, by population and site:

- called sample count;
- call rate;
- alternate-allele frequency;
- MAF;
- observed heterozygosity (`Ho`); and
- expected heterozygosity (`He = 2p(1-p)`).

Population-level summaries include:

- number of samples;
- sites observed;
- polymorphic-site count;
- mean call rate;
- mean observed heterozygosity;
- mean expected heterozygosity; and
- mean MAF.

### Important interpretation note

These statistics are summaries across **variant records present in the analyzed VCF**. They should not be interpreted as genome-wide per-base nucleotide diversity (`π`) unless invariant callable sites and the appropriate denominator are represented. EvoFlow deliberately labels the current metric scope to avoid that common mistake.

## Diversity outputs

```text
evoflow-results/diversity/
├── diversity_summary.json
├── site_diversity.csv
├── population_diversity.csv
├── population_heterozygosity.png
├── population_heterozygosity.pdf
├── population_call_rate.png
└── population_call_rate.pdf
```

---

# Population differentiation: Hudson FST

## Method

```bash
evoflow fst evoflow.yaml
```

EvoFlow currently estimates **pairwise Hudson FST** between populations defined in metadata. For each usable biallelic SNP and population pair, the implementation calculates within-population diversity and between-population divergence, stores a per-site numerator/denominator, and computes the multi-site estimator as:

```text
FST = sum(site numerators) / sum(site denominators)
```

Per-site FST estimates can be negative because of sampling variation; EvoFlow does not silently clamp them to zero. The pairwise aggregate is reported from the summed numerator and denominator components.

At least two populations are required.

## FST outputs

```text
evoflow-results/fst/
├── fst_summary.json
├── site_fst.csv
├── pairwise_fst.csv
├── fst_matrix.csv
├── fst_heatmap.png
└── fst_heatmap.pdf
```

`fst_matrix.csv` provides a symmetric population-by-population matrix suitable for later reporting and visualization.

---

## Analysis modules

| Module | Role | Current state |
|---|---|---|
| `qc` | Variant and sample quality control | ✅ Implemented |
| `ld` | Linkage-disequilibrium pruning | ✅ Implemented |
| `pca` | Principal component analysis | ✅ Implemented |
| `diversity` | Population diversity summaries | ✅ Implemented foundation |
| `fst` | Pairwise Hudson FST | ✅ Implemented foundation |
| `structure` | Ancestry / clustering | 🧭 Planned |
| `spatial` | Spatial population genomics | 🧭 Planned |
| `selection` | Selection / outlier scans | 🧭 Planned |
| `gea` | Genotype-environment association | 🧭 Planned |
| `report` | Reproducible scientific reporting | 🧭 Planned |

---

# Output organization

Current modules write deterministic, analysis-specific directories:

```text
evoflow-results/
├── qc/
├── ld/
├── pca/
├── diversity/
└── fst/
```

As the orchestration layer matures, these directories will become part of a formal run manifest and provenance system.

---

# Repository architecture

```text
EvoFlow/
├── .github/
│   └── workflows/
│       └── ci.yml
├── docs/
│   └── architecture.md
├── examples/
│   └── evoflow.example.yaml
├── src/
│   └── evoflow/
│       ├── __init__.py
│       ├── cli.py
│       ├── config.py
│       ├── core/
│       ├── io/
│       │   ├── metadata.py
│       │   ├── validation.py
│       │   └── vcf.py
│       └── modules/
│           ├── diversity.py
│           ├── diversity_plots.py
│           ├── fst.py
│           ├── fst_plots.py
│           ├── ld.py
│           ├── pca.py
│           ├── pca_plots.py
│           ├── qc.py
│           ├── qc_plots.py
│           └── registry.py
├── tests/
│   ├── test_config.py
│   ├── test_ld.py
│   ├── test_pca.py
│   ├── test_population_stats.py
│   ├── test_qc.py
│   ├── test_registry.py
│   └── test_validation.py
├── CONTRIBUTING.md
├── LICENSE
├── README.md
└── pyproject.toml
```

See [`docs/architecture.md`](docs/architecture.md) for implementation-level design notes.

---

# Testing and continuous integration

Every push to `main` is checked with GitHub Actions on:

- Python 3.10
- Python 3.11
- Python 3.12

Each environment performs:

```text
Install -> Ruff linting -> Pytest
```

The test suite includes hand-checkable scientific expectations, not only import/smoke tests. Current coverage includes:

- VCF and compressed-VCF parsing;
- strict metadata/VCF sample concordance;
- expected QC values;
- DP and GQ summaries;
- filtered VCF contents and provenance headers;
- QC PNG/PDF generation;
- known LD patterns with `r² = 1` and `r² = 0`;
- physical-window LD behavior;
- PCA informative-SNP selection and explained variance;
- PCA metadata propagation, loadings, and figures;
- expected population Ho/He/MAF summaries; and
- Hudson FST components, aggregate pairwise FST, matrix symmetry, and heatmap generation.

The CI badge at the top of this README reflects the current `main` branch.

---

# Current limitations

EvoFlow is **not yet a complete end-to-end population-genomics pipeline**. Important limitations are documented explicitly:

- native BCF parsing is not yet available;
- current scientific genotype engines assume diploid biallelic dosage data;
- QC does not yet filter genotypes by `GQ`, although GQ is summarized;
- Hardy-Weinberg filtering is not yet implemented;
- multiallelic MAF is not currently assigned;
- filtered/pruned VCF outputs are currently uncompressed;
- LD pruning is a greedy physical-window implementation rather than a PLINK-compatible sliding variant-count algorithm;
- PCA mean-imputes missing genotypes using variant allele frequency;
- diversity metrics currently summarize analyzed variant records and are not genome-wide per-base `π`;
- current FST is pairwise Hudson FST only; confidence intervals/block-jackknife support is not yet implemented;
- population-structure/ancestry inference is not yet implemented;
- spatial, selection, GEA, and reporting modules remain under development; and
- full run manifests, checksums, external-tool version capture, restartability, and dependency resolution are not yet implemented.

---

# Roadmap

## Phase 0 — Software foundation

- [x] Repository architecture
- [x] Installable Python package
- [x] CLI
- [x] YAML project configuration
- [x] Analysis registry
- [x] Automated tests
- [x] Python 3.10–3.12 CI
- [x] Documentation foundation

## Phase 1 — Genomic quality control

- [x] Stream VCF / VCF.gz
- [x] Variant and sample counts
- [x] Exact VCF/metadata sample validation
- [x] Sample/site missingness
- [x] Genotype call rate
- [x] Observed heterozygosity
- [x] Biallelic AF / MAF
- [x] DP summaries
- [x] GQ summaries
- [x] Ti/Tv
- [x] MAF/missingness filtering
- [x] Filtered VCF
- [x] QC tables
- [x] PNG/PDF QC figures
- [ ] Genotype-level GQ filtering
- [ ] Compressed filtered VCF
- [ ] Native BCF backend

## Phase 2 — Population structure and differentiation

- [x] LD pruning
- [x] Allele-frequency-standardized PCA
- [x] Streaming Gram-matrix PCA
- [x] PCA scores / variance / SNP loadings
- [x] PCA figures
- [x] Population Ho / He / MAF summaries
- [x] Pairwise Hudson FST
- [x] FST matrix / heatmap
- [ ] Population-structure / ancestry inference
- [ ] Additional diversity estimators
- [ ] FST uncertainty / block jackknife
- [ ] Hardy-Weinberg filtering

## Phase 3 — Spatial and adaptive genomics

- [ ] Isolation by distance
- [ ] Spatial population-genomic summaries
- [ ] Selection / outlier scans
- [ ] Genotype-environment association
- [ ] Geographic visualization interfaces

## Phase 4 — Workflow orchestration

- [ ] Module dependency resolution
- [ ] Ordered execution engine
- [ ] Run manifests and checksums
- [ ] Parameter provenance
- [ ] External-tool version capture
- [ ] Restartable / resumable execution

## Phase 5 — Reporting and distribution

- [ ] Automated HTML scientific report
- [ ] Integrated methods summary
- [ ] Example biological dataset
- [ ] Extended user documentation
- [ ] Containerized release
- [ ] Versioned software release
- [ ] DOI / archival release workflow

---

# Scientific scope

EvoFlow is being developed for evolutionary biologists, population geneticists, conservation genomicists, landscape genomicists, molecular ecologists, and bioinformaticians working with questions such as:

- How are individuals genetically structured across populations?
- Which SNPs dominate major axes of genomic variation?
- How much heterozygosity and polymorphism is represented within populations?
- Which populations are most differentiated?
- Is genomic differentiation associated with geography?
- Are there barriers or corridors to gene flow?
- Which loci show unusual differentiation?
- Are genotypes associated with environmental variation?
- Can the full computational path be reproduced from one documented project configuration?

EvoFlow does not replace biological interpretation. Its purpose is to make the computational path to that interpretation more coherent, transparent, testable, and reproducible.

---

# Contributing

EvoFlow is under active development. Issues, test cases, documentation improvements, scientifically justified module proposals, and code contributions are welcome.

See [`CONTRIBUTING.md`](CONTRIBUTING.md).

Substantial analysis modules should clearly define:

1. required inputs;
2. statistical assumptions;
3. configurable parameters;
4. output files;
5. dependency requirements;
6. validation and failure rules;
7. provenance information; and
8. tests against known expected results where feasible.

---

# Citation

EvoFlow does not yet have a versioned archival release or DOI.

Until one is available, researchers referring to the project should cite the repository and record the EvoFlow version or Git commit used. Formal citation instructions will be added with the first archival release.

---

# Author

**Pauline Owusu-Ansah**  
Creator, author, and lead developer of EvoFlow  
GitHub: [@codewithPauline](https://github.com/codewithPauline)

---

# License

**EvoFlow is open-source software released under the MIT License.**

**Copyright © 2026 Pauline Owusu-Ansah.**

The MIT License permits use, copying, modification, merging, publication, distribution, sublicensing, and sale of copies of the software, provided that the copyright and permission notice are retained in copies or substantial portions of the software.

The software is provided **“as is”**, without warranty of any kind, as described in the full license text.

See the complete license: **[`LICENSE`](LICENSE)**.

---

<p align="center">
  <strong>EvoFlow</strong><br>
  <em>From variants to evolutionary insight.</em>
</p>
