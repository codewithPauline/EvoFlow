# EvoFlow

[![CI](https://github.com/codewithPauline/EvoFlow/actions/workflows/ci.yml/badge.svg)](https://github.com/codewithPauline/EvoFlow/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Development Status](https://img.shields.io/badge/status-active%20development-orange.svg)](#development-status)

> **From variants to evolutionary insight.**

**EvoFlow** is an open-source, species-agnostic workflow platform for reproducible population and landscape genomics. It is being developed to connect genomic variants, sample metadata, population-genetic analyses, spatial information, and scientific reporting through one transparent command-line workflow.

**Author and lead developer:** [Pauline Owusu-Ansah](https://github.com/codewithPauline) (`@codewithPauline`)

---

## Why EvoFlow?

Population-genomic studies rarely involve a single analysis. A typical project may require variant and sample quality control, filtering, PCA, ancestry inference, diversity statistics, FST, spatial analyses, selection scans, genotype-environment association, visualization, and reproducible reporting.

In practice, those steps are often distributed across different command-line tools, R scripts, Python notebooks, file formats, environments, and plotting workflows. That fragmentation can make analyses difficult to reproduce, audit, extend, and transfer between projects or organisms.

**EvoFlow is being built as a reproducible orchestration layer for that workflow.**

The project is guided by five principles:

1. **Reproducibility first** — inputs, parameters, software versions, and outputs should remain traceable.
2. **Species agnostic** — development examples may use salamander data, but the software is intended for any organism.
3. **Modular analyses** — researchers should be able to run only the analyses they need.
4. **Transparent methods** — statistical assumptions and parameters should remain visible rather than hidden behind an opaque interface.
5. **Scientific outputs** — the end product should be interpretable tables, figures, provenance, and reports, not merely successful commands.

---

## Conceptual workflow

```text
        VCF / VCF.gz + sample metadata
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
             filtered VCF
                    |
       +------------+-------------+
       |            |             |
       v            v             v
      PCA      Population      Diversity
               structure
       |            |             |
       +------------+-------------+
                    |
                    v
            Differentiation / FST
                    |
                    v
            Spatial population genomics
                    |
          +---------+----------+
          |                    |
          v                    v
   Selection scans            GEA
          |                    |
          +---------+----------+
                    |
                    v
        Figures + tables + provenance
                    |
                    v
          Reproducible scientific report
```

The complete workflow is the long-term target. EvoFlow deliberately distinguishes functionality that works today from modules that remain under development.

---

## Development status

EvoFlow is in **active early development**. The software foundation is installable and tested, and two scientific components are now functional: native VCF quality control and population-genomic PCA.

| Component | Purpose | Status |
|---|---|---|
| Python package | Installable `evoflow` package | ✅ Available |
| Command-line interface | User-facing EvoFlow commands | ✅ Available |
| YAML configuration | Reproducible project configuration | ✅ Available |
| VCF / metadata validation | Confirms files exist and sample IDs match exactly | ✅ Available |
| Native VCF QC engine | Streams VCF/VCF.gz and computes genomic QC metrics | ✅ Available |
| QC filtering | MAF/missingness decisions and filtered VCF output | ✅ Available |
| QC tables | Run-, sample-, and variant-level summaries | ✅ Available |
| QC figures | PNG/PDF diagnostic plots | ✅ Available |
| GQ / depth summaries | Uses FORMAT `GQ` and `DP` when present | ✅ Available |
| PCA engine | Allele-frequency-standardized population-genomic PCA | ✅ Available |
| PCA tables | Scores, loadings, variance, summary | ✅ Available |
| PCA figures | Scree and PC1–PC2 plots in PNG/PDF | ✅ Available |
| Continuous integration | Install, lint, and test on Python 3.10–3.12 | ✅ Passing |
| Population structure | Ancestry / clustering workflow | 🧭 Planned |
| Diversity statistics | Population diversity summaries | 🧭 Planned |
| FST | Population differentiation | 🧭 Planned |
| Spatial genomics | IBD and spatial population-genomic analyses | 🧭 Planned |
| Selection scans | Outlier / selection analyses | 🧭 Planned |
| GEA | Genotype-environment association | 🧭 Planned |
| Workflow execution engine | Ordered module execution and provenance | 🧭 Planned |
| Automated report | Integrated figures, tables, methods, and provenance | 🧭 Planned |

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

### 2. Validate genomic and metadata inputs

```bash
evoflow validate evoflow.yaml
```

Validation confirms that the configured files exist and that VCF and metadata sample identities match exactly.

### 3. Inspect the configured workflow

```bash
evoflow plan evoflow.yaml
```

### 4. Run quality control

```bash
evoflow qc evoflow.yaml --min-maf 0.05 --max-missing 0.20
```

### 5. Run PCA

```bash
evoflow pca evoflow.yaml --n-components 10
```

If `evoflow-results/qc/filtered.vcf` exists, PCA uses it automatically. To force PCA to use the original configured VCF:

```bash
evoflow pca evoflow.yaml --use-raw
```

---

## Project configuration

A project is described with YAML:

```yaml
project: salamander-demo
vcf: data/variants.vcf.gz
metadata: data/samples.csv
output_dir: evoflow-results
modules:
  - qc
  - pca
  - structure
  - diversity
  - fst
  - spatial
  - selection
  - gea
  - report
```

At present, `qc` and `pca` are executable. The remaining module names define EvoFlow's planned analysis vocabulary.

---

## Input validation

Before analysis, EvoFlow checks that:

- the configured VCF exists;
- the metadata CSV exists;
- metadata contains a `sample` column;
- metadata sample IDs are non-empty and unique;
- VCF sample IDs are unique; and
- VCF and metadata contain **exactly the same sample IDs**.

This strict behavior is intentional. Silent sample mismatches can invalidate downstream population-genetic analyses, so EvoFlow fails early instead of guessing.

---

# Quality-control engine

## Native streaming VCF QC

The QC engine reads `.vcf` and `.vcf.gz` files without loading the complete variant dataset into memory.

It currently calculates:

- sample count;
- variant count;
- SNP count;
- biallelic SNP count;
- multiallelic variant count;
- transition and transversion counts;
- Ti/Tv ratio when defined;
- called and missing genotype counts;
- overall genotype call rate;
- per-sample call rate and missingness;
- per-sample observed heterozygosity;
- per-site call rate and missingness;
- alternate-allele frequency;
- biallelic minor-allele frequency;
- per-sample and per-site mean depth when FORMAT `DP` is available; and
- per-sample and per-site mean genotype quality when FORMAT `GQ` is available.

## QC filtering

Variants can be evaluated with explicit MAF and missingness thresholds:

```bash
evoflow qc evoflow.yaml --min-maf 0.05 --max-missing 0.20
```

Each variant receives a `passes_qc` decision. Passing records are also written to `filtered.vcf`.

The output VCF records the selected thresholds in its header:

```text
##evoflow_qc_min_maf=0.05
##evoflow_qc_max_missing=0.2
```

For multiallelic sites, the current implementation does not assign a biallelic-style MAF. Those records are evaluated by missingness but not by the MAF threshold.

## QC outputs

```text
evoflow-results/
└── qc/
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

The variant-distribution plots use fixed-size streaming histogram bins so large VCFs do not require retaining millions of QC values in memory.

---

# PCA engine

## Statistical approach

EvoFlow's PCA module currently operates on **informative diploid biallelic SNPs**.

For each SNP with alternate-allele frequency `p`, genotype dosage `g` is standardized as:

```text
(g - 2p) / sqrt(2p(1-p))
```

Missing genotypes are mean-imputed to `2p`, which becomes zero after standardization. Monomorphic variants, indels, multiallelic variants, sites without usable diploid GT calls, and other non-informative records are excluded from PCA.

Instead of retaining the complete sample × SNP matrix, EvoFlow streams standardized SNP vectors and accumulates the sample Gram matrix:

```text
X Xᵀ
```

The principal components are then recovered through symmetric eigendecomposition. A second streaming pass calculates per-SNP loadings. This keeps PCA memory usage driven primarily by the number of samples rather than the number of SNPs.

## PCA input selection

By default:

```text
QC filtered.vcf exists  -> use filtered.vcf
otherwise                -> use configured VCF
```

Use `--use-raw` to explicitly bypass an existing QC-filtered VCF.

## PCA outputs

```text
evoflow-results/
└── pca/
    ├── pca_summary.json
    ├── pca_scores.csv
    ├── pca_variance.csv
    ├── pca_loadings.csv
    ├── pca_scree.png
    ├── pca_scree.pdf
    ├── pca_pc1_pc2.png
    └── pca_pc1_pc2.pdf
```

### `pca_scores.csv`

Contains one row per sample, preserves available metadata columns, and appends PC coordinates.

If the metadata contains a `population` column, the PC1–PC2 plot automatically groups samples by population.

### `pca_variance.csv`

Contains:

- component name;
- eigenvalue;
- explained variance;
- explained-variance ratio; and
- cumulative explained variance.

### `pca_loadings.csv`

Contains chromosome, position, reference/alternate alleles, and loadings for each reported component.

### `pca_summary.json`

Records the input VCF, sample count, informative SNP count, number of components, explained-variance ratios, and PCA method.

---

## Input data

### Variant file

The native parser currently supports:

- uncompressed VCF (`.vcf`)
- gzip/bgzip-compressed VCF (`.vcf.gz`, `.vcf.bgz`, `.gz`, `.bgz`)

**BCF is not yet parsed natively.** Convert BCF to VCF/VCF.gz before using the current native engines.

### Sample metadata

The minimum metadata file is CSV with a `sample` column:

```csv
sample,population,site,latitude,longitude
sample_01,OH,Site_A,39.5100,-84.7300
sample_02,IN,Site_B,39.1600,-86.5200
sample_03,KY,Site_C,38.0400,-84.5000
```

Only `sample` is mandatory today. Other columns are preserved in PCA score output and will support later population, spatial, and landscape-genomic modules.

---

## Analysis modules

| Module | Intended role | Current state |
|---|---|---|
| `qc` | Variant and sample quality control | ✅ Implemented |
| `pca` | Principal component analysis | ✅ Implemented foundation |
| `structure` | Population structure / ancestry inference | 🧭 Planned |
| `diversity` | Population diversity statistics | 🧭 Planned |
| `fst` | Population differentiation | 🧭 Planned |
| `spatial` | Spatial population-genomic analyses | 🧭 Planned |
| `selection` | Selection and outlier scans | 🧭 Planned |
| `gea` | Genotype-environment association | 🧭 Planned |
| `report` | Reproducible analysis reporting | 🧭 Planned |

---

## Architecture

EvoFlow uses a `src`-based Python package layout and separates configuration, input handling, scientific modules, plotting, orchestration, and the CLI.

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
│       │   ├── validation.py
│       │   └── vcf.py
│       └── modules/
│           ├── pca.py
│           ├── pca_plots.py
│           ├── qc.py
│           ├── qc_plots.py
│           └── registry.py
├── tests/
│   ├── test_config.py
│   ├── test_pca.py
│   ├── test_qc.py
│   ├── test_registry.py
│   └── test_validation.py
├── CONTRIBUTING.md
├── LICENSE
├── README.md
└── pyproject.toml
```

See [`docs/architecture.md`](docs/architecture.md) for the evolving internal design.

---

## Testing and continuous integration

Every push to `main` is checked with GitHub Actions on:

- Python 3.10
- Python 3.11
- Python 3.12

Each environment performs:

```text
Install  ->  Ruff linting  ->  Pytest
```

The test suite currently covers:

- VCF parsing;
- compressed VCF support;
- strict metadata/VCF sample matching;
- expected QC statistics;
- DP and GQ summaries;
- threshold-filtered VCF output and provenance headers;
- QC diagnostic PNG/PDF files;
- PCA informative-SNP selection;
- PCA explained variance and component structure;
- metadata propagation into PCA scores;
- PCA loadings; and
- PCA figure generation.

The CI badge at the top of this README reflects the current state of `main`.

---

## Current limitations

EvoFlow is **not yet a complete end-to-end population-genomics pipeline**. Important current limitations include:

- no native BCF parser;
- no LD-pruning module yet;
- no Hardy-Weinberg filtering yet;
- no genotype-quality (`GQ`) filtering threshold yet, although GQ summaries are calculated;
- multiallelic MAF is not currently assigned;
- filtered VCF output is currently uncompressed;
- PCA currently assumes diploid biallelic genotype dosages;
- PCA mean-imputes missing genotypes at each SNP;
- PCA does not yet perform LD pruning automatically, so users should interpret analyses of dense linked SNP datasets accordingly;
- ancestry/structure, diversity, FST, spatial, selection, GEA, and reporting modules remain under development; and
- full provenance manifests and restartable workflow execution are not yet implemented.

These limitations are documented intentionally so current functionality is not confused with the roadmap.

---

## Roadmap

### Phase 0 — Software foundation

- [x] Repository architecture
- [x] Python package structure
- [x] Command-line interface
- [x] YAML project configuration
- [x] Analysis-module registry
- [x] Automated tests
- [x] Continuous integration
- [x] Project documentation foundation

### Phase 1 — Genomic quality control

- [x] Stream and summarize VCF content
- [x] Support VCF and VCF.gz
- [x] Variant and sample counts
- [x] Exact metadata / VCF sample validation
- [x] Per-sample and per-site missingness
- [x] Genotype call rate
- [x] Observed heterozygosity
- [x] Biallelic allele-frequency / MAF summaries
- [x] Depth summaries from FORMAT/DP
- [x] Genotype-quality summaries from FORMAT/GQ
- [x] Transition / transversion summaries
- [x] Configurable MAF and missingness QC decisions
- [x] Run-, sample-, and variant-level QC tables
- [x] Threshold-filtered VCF output
- [x] Publication-oriented QC figures in PNG/PDF
- [ ] GQ-based filtering thresholds
- [ ] Compressed filtered VCF output
- [ ] Native BCF support

### Phase 2 — Population structure and diversity

- [x] Allele-frequency-standardized PCA engine
- [x] Streaming sample Gram-matrix PCA
- [x] PCA sample scores and explained variance
- [x] SNP loadings
- [x] PCA scree and PC1–PC2 figures
- [x] Metadata-aware population grouping in PCA plots
- [ ] LD pruning
- [ ] Population-structure / ancestry workflow
- [ ] Diversity statistics
- [ ] Pairwise and global FST
- [ ] Standardized population-genetic figures

### Phase 3 — Spatial and adaptive genomics

- [ ] Isolation-by-distance analysis
- [ ] Spatial population-genomic summaries
- [ ] Selection / outlier scans
- [ ] Genotype-environment association
- [ ] Geographic visualization interfaces

### Phase 4 — Reproducible workflow execution

- [ ] Module dependency resolution
- [ ] Ordered execution engine
- [ ] Provenance manifests
- [ ] Restartable / resumable runs
- [ ] Standardized result directories
- [ ] External-tool version tracking

### Phase 5 — Reporting and distribution

- [ ] Automated HTML scientific report
- [ ] Integrated methods summary
- [ ] Example biological dataset
- [ ] Extended user documentation
- [ ] Containerized release
- [ ] Versioned software release
- [ ] DOI / archival release workflow

---

## Scientific scope

EvoFlow is being developed for evolutionary biologists, population geneticists, conservation genomicists, landscape genomicists, molecular ecologists, and bioinformaticians working with questions such as:

- How are individuals genetically structured across populations or geography?
- Which populations are most differentiated?
- How much genomic diversity exists within populations?
- Is genetic differentiation associated with geographic distance?
- Are there potential barriers or corridors to gene flow?
- Which genomic regions show unusual differentiation?
- Are genotypes associated with environmental variation?
- Can these analyses be reproduced from one documented project configuration?

EvoFlow will not replace biological interpretation. Its purpose is to make the computational path to that interpretation more coherent, transparent, testable, and reproducible.

---

## Contributing

EvoFlow is under active development. Issues, test cases, documentation improvements, scientifically justified module proposals, and code contributions are welcome.

See [`CONTRIBUTING.md`](CONTRIBUTING.md).

Substantial analysis modules should clearly define:

1. required inputs;
2. output files;
3. configurable parameters;
4. software dependencies;
5. validation rules;
6. provenance information; and
7. tests.

---

## Citation

EvoFlow does not yet have a versioned archival release or DOI.

Until one is available, researchers referring to the project should cite the repository and record the EvoFlow version or Git commit used. Formal citation instructions will be added with the first archival release.

---

## Author

**Pauline Owusu-Ansah**  
Creator, author, and lead developer of EvoFlow  
GitHub: [@codewithPauline](https://github.com/codewithPauline)

---

## License

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
