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

Population-genomic studies usually require many independent steps: input validation, variant filtering, sample quality control, PCA, ancestry inference, diversity statistics, FST, spatial analyses, selection scans, genotype-environment association, plotting, and reporting.

Those steps are often spread across different command-line tools, R scripts, Python notebooks, file formats, and software environments. That fragmentation makes analyses harder to reproduce, audit, extend, and transfer between organisms or projects.

**EvoFlow is being built as a reproducible orchestration layer for that workflow.**

The project is designed around five principles:

1. **Reproducibility first** — inputs, configuration, parameters, software versions, and outputs should remain traceable.
2. **Species agnostic** — development examples may use salamander data, but the software is intended for any organism.
3. **Modular analyses** — researchers should be able to run only the analyses they need.
4. **Transparent methods** — EvoFlow should expose methods and parameters rather than hide them behind an opaque interface.
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

The full workflow is the long-term target. The repository clearly separates **implemented** functionality from **planned** modules.

---

## Development status

EvoFlow is in **active early development**. The software foundation is installable and tested, and the first biological analysis engine — native VCF quality control — is now implemented.

| Component | Purpose | Status |
|---|---|---|
| Python package | Installable `evoflow` package | ✅ Available |
| Command-line interface | User-facing EvoFlow commands | ✅ Available |
| YAML configuration | Reproducible project configuration | ✅ Available |
| VCF / metadata validation | Confirms files exist and sample IDs match exactly | ✅ Available |
| Analysis registry | Defines supported EvoFlow modules | ✅ Available |
| Native VCF QC engine | Streams VCF/VCF.gz and computes genomic QC metrics | ✅ Available |
| QC tables | Run-, sample-, and variant-level CSV/JSON summaries | ✅ Available |
| MAF / missingness thresholds | Configurable per-variant QC decisions | ✅ Available |
| Filtered VCF | Writes variants that pass configured QC thresholds | ✅ Available |
| Continuous integration | Install, lint, and test on Python 3.10–3.12 | ✅ Passing |
| QC figures | Publication-quality diagnostic plots | 🚧 In development |
| PCA | Population-genomic dimensionality reduction | 🧭 Planned |
| Population structure | Ancestry / clustering workflow | 🧭 Planned |
| Diversity statistics | Population diversity summaries | 🧭 Planned |
| FST | Population differentiation | 🧭 Planned |
| Spatial genomics | IBD and spatial population-genomic analyses | 🧭 Planned |
| Selection scans | Outlier / selection analyses | 🧭 Planned |
| GEA | Genotype-environment association | 🧭 Planned |
| Workflow execution engine | Ordered module execution and provenance | 🧭 Planned |
| Automated report | Integrated figures, tables, methods, and provenance | 🧭 Planned |

---

## What works today

### 1. Project configuration

EvoFlow projects are described with YAML:

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

At this stage, `qc` is executable. The remaining module names define the planned workflow vocabulary and are being implemented incrementally.

### 2. Strict input validation

Before QC runs, EvoFlow checks that:

- the configured VCF exists;
- the metadata CSV exists;
- the metadata contains a `sample` column;
- metadata sample IDs are non-empty and unique;
- VCF sample IDs are unique; and
- the VCF and metadata contain **exactly the same sample IDs**.

That last check is intentionally strict. Silent sample mismatches can invalidate downstream population-genomic analyses, so EvoFlow fails early instead of guessing.

### 3. Native streaming VCF QC

The first QC engine reads `.vcf` and `.vcf.gz` files without loading the entire variant dataset into memory.

It currently calculates:

- number of samples;
- number of variants;
- SNP count;
- biallelic SNP count;
- multiallelic variant count;
- transition count;
- transversion count;
- Ti/Tv ratio when defined;
- total genotype calls;
- called and missing genotypes;
- overall genotype call rate;
- per-sample call rate;
- per-sample missingness;
- per-sample observed heterozygosity;
- per-sample mean depth when `DP` is available;
- per-site call rate and missingness;
- alternate-allele frequency;
- biallelic minor-allele frequency; and
- per-site mean depth when `DP` is available.

### 4. Configurable QC thresholds and filtered VCF output

The QC command can filter variants using minor-allele-frequency and missingness thresholds:

```bash
evoflow qc evoflow.yaml --min-maf 0.05 --max-missing 0.20
```

Each variant receives a `passes_qc` decision in `variant_qc.csv`. Variants that pass are also written to `filtered.vcf`.

The filtered VCF preserves the original VCF metadata/header and records the EvoFlow thresholds in additional header lines:

```text
##evoflow_qc_min_maf=0.05
##evoflow_qc_max_missing=0.2
```

For multiallelic sites, the current engine does not assign a biallelic-style MAF; those sites are therefore evaluated by missingness but not by the MAF threshold.

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

### 2. Validate the inputs

```bash
evoflow validate evoflow.yaml
```

Successful validation confirms that the genomic and metadata sample identities match.

### 3. Inspect the analysis plan

```bash
evoflow plan evoflow.yaml
```

### 4. Run VCF quality control

```bash
evoflow qc evoflow.yaml
```

Or apply explicit QC thresholds:

```bash
evoflow qc evoflow.yaml --min-maf 0.05 --max-missing 0.20
```

---

## QC output

The QC command writes:

```text
evoflow-results/
└── qc/
    ├── qc_summary.json
    ├── sample_qc.csv
    ├── variant_qc.csv
    └── filtered.vcf
```

### `qc_summary.json`

Run-level metrics, including sample and variant counts, retained-variant count, SNP/multiallelic counts, transitions/transversions, overall genotype call rate, mean depth when available, and the QC thresholds used.

### `sample_qc.csv`

One row per sample:

```text
sample
called_genotypes
missing_genotypes
call_rate
missing_rate
heterozygosity
mean_depth
```

### `variant_qc.csv`

One row per variant:

```text
chrom
pos
ref
alt
is_snp
is_biallelic
call_rate
missing_rate
alt_allele_frequency
maf
mean_depth
passes_qc
```

For multiallelic sites, `maf` is intentionally left blank in the current engine rather than applying an ambiguous biallelic definition.

### `filtered.vcf`

A VCF containing only records that pass the selected MAF and missingness criteria. The EvoFlow filter thresholds are recorded in the VCF header.

---

## Input data

### Variant file

The native QC engine currently supports:

- uncompressed VCF (`.vcf`)
- gzip/bgzip-compressed VCF (`.vcf.gz`, `.vcf.bgz`, `.gz`, `.bgz`)

**BCF is not yet parsed natively.** Convert BCF to VCF/VCF.gz before running the current QC engine. Native BCF support can be added later through an appropriate binary variant backend.

The parser uses genotype (`GT`) information when present and depth (`DP`) when available in FORMAT fields.

### Sample metadata

The minimum metadata file is CSV with a `sample` column:

```csv
sample,population,site,latitude,longitude
sample_01,OH,Site_A,39.5100,-84.7300
sample_02,IN,Site_B,39.1600,-86.5200
sample_03,KY,Site_C,38.0400,-84.5000
```

Only `sample` is required today. Columns such as `population`, `site`, `latitude`, `longitude`, and environmental covariates will support later population, spatial, and landscape-genomic modules.

---

## Analysis modules

| Module | Intended role | Current state |
|---|---|---|
| `qc` | Variant and sample quality control | ✅ Implemented foundation |
| `pca` | Principal component analysis | 🧭 Planned |
| `structure` | Population structure / ancestry inference | 🧭 Planned |
| `diversity` | Population diversity statistics | 🧭 Planned |
| `fst` | Population differentiation | 🧭 Planned |
| `spatial` | Spatial population-genomic analyses | 🧭 Planned |
| `selection` | Selection and outlier scans | 🧭 Planned |
| `gea` | Genotype-environment association | 🧭 Planned |
| `report` | Reproducible analysis reporting | 🧭 Planned |

---

## Architecture

EvoFlow uses a `src`-based Python package layout and separates configuration, input handling, analysis modules, orchestration, and the command-line interface.

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
│           ├── qc.py
│           └── registry.py
├── tests/
│   ├── test_config.py
│   ├── test_qc.py
│   ├── test_registry.py
│   └── test_validation.py
├── CONTRIBUTING.md
├── LICENSE
├── README.md
└── pyproject.toml
```

### Architectural layers

- **`evoflow.config`** — project configuration and serialization.
- **`evoflow.io`** — VCF access and biological metadata validation.
- **`evoflow.modules`** — analysis implementations and registry.
- **`evoflow.core`** — reserved for orchestration, execution, provenance, and restartable run management.
- **`evoflow.cli`** — user-facing commands.

The native QC implementation is intentionally streaming so that basic summaries and filtering do not require retaining an entire VCF in Python memory.

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

The tests include VCF parsing, gzip support, expected QC statistics, filtered VCF output, module validation, configuration handling, and strict VCF/metadata sample matching.

The CI badge at the top of this README reflects the current state of `main`.

---

## Current limitations

EvoFlow is usable for its implemented QC foundation, but it is **not yet a complete end-to-end population-genomics pipeline**.

Current limitations include:

- no native BCF parser yet;
- no genotype-quality (`GQ`) filtering yet;
- multiallelic MAF is not currently assigned;
- QC plots are not yet generated;
- filtered VCF output is currently written uncompressed;
- PCA, structure, diversity, FST, spatial, selection, GEA, and report modules remain under development; and
- provenance manifests and restartable workflow execution are not yet implemented.

These limitations are documented deliberately so users can distinguish current functionality from the roadmap.

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
- [x] Variant counts
- [x] Sample counts
- [x] Exact metadata / VCF sample validation
- [x] Per-sample missingness
- [x] Per-site missingness
- [x] Genotype call rate
- [x] Observed heterozygosity
- [x] Biallelic allele-frequency / MAF summaries
- [x] Depth summaries when FORMAT/DP is present
- [x] Transition / transversion summaries
- [x] Configurable MAF and missingness QC decisions
- [x] Run-, sample-, and variant-level QC tables
- [x] Threshold-filtered VCF output
- [ ] Genotype-quality (`GQ`) summaries and thresholds
- [ ] Publication-quality QC figures
- [ ] Compressed filtered VCF output
- [ ] Native BCF support

### Phase 2 — Population structure and diversity

- [ ] PCA engine
- [ ] Population-structure workflow
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
- How much genomic diversity is present within populations?
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
