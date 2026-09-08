# EvoFlow

[![CI](https://github.com/codewithPauline/EvoFlow/actions/workflows/ci.yml/badge.svg)](https://github.com/codewithPauline/EvoFlow/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Development Status](https://img.shields.io/badge/status-active%20development-orange.svg)](#development-status)

> **From variants to evolutionary insight.**

**EvoFlow** is an open-source, species-agnostic workflow platform for reproducible population and landscape genomics. It is being designed to connect variant data, biological metadata, spatial information, population-genetic analyses, and publication-ready reporting through one transparent and consistent interface.

**Author and lead developer:** [Pauline Owusu-Ansah](https://github.com/codewithPauline) (`@codewithPauline`)

---

## Overview

Population-genomic projects rarely consist of a single analysis. A typical study may require variant filtering, sample quality control, PCA, ancestry inference, diversity statistics, FST, spatial analyses, selection scans, genotype-environment association, plotting, and finally the assembly of all parameters and results into a reproducible report.

In practice, those steps are often distributed across different command-line programs, R scripts, Python notebooks, file formats, environments, and plotting workflows. That fragmentation can make analyses difficult to reproduce, audit, extend, or transfer to another organism.

**EvoFlow is being built to solve that workflow problem.**

Rather than replacing established population-genomic methods, EvoFlow is intended to act as a transparent orchestration layer that:

- validates genomic and metadata inputs before analysis;
- organizes analyses into explicit modules;
- records the configuration used for each project;
- provides one command-line interface for the workflow;
- keeps underlying methods and parameters visible to the researcher;
- standardizes tables, figures, and output organization; and
- ultimately produces a reproducible analysis report linking variants to evolutionary interpretation.

---

## Conceptual workflow

```text
            Variant data
             VCF / BCF
                 |
                 v
       +--------------------+
       |  Input validation  | <---- sample metadata
       +--------------------+ <---- coordinates / environment
                 |
                 v
       +--------------------+
       | Variant & sample QC|
       +--------------------+
                 |
       +---------+----------+----------------+
       |                    |                |
       v                    v                v
      PCA            Population          Diversity
                      structure
       |                    |                |
       +----------+---------+----------------+
                  |
                  v
          Population differentiation
                  FST
                  |
                  v
        Spatial population genomics
                  |
         +--------+--------+
         |                 |
         v                 v
   Selection scans       GEA
                         genotype-
                         environment
                         association
         |                 |
         +--------+--------+
                  |
                  v
       Figures + tables + provenance
                  |
                  v
        Reproducible scientific report
```

---

## Development status

EvoFlow is in **active early development**. The software foundation is working, tested, and installable, while the biological analysis engines are being implemented incrementally.

The distinction below is intentional: the README should describe what EvoFlow can do **today** without presenting planned modules as completed software.

| Component | Purpose | Status |
|---|---|---|
| Python package | Installable `evoflow` package | ✅ Available |
| Command-line interface | User-facing EvoFlow commands | ✅ Available |
| YAML configuration | Reproducible project configuration | ✅ Available |
| Variant-file path validation | Confirms configured VCF/BCF input exists | ✅ Available |
| Metadata validation | Confirms metadata exists, contains a `sample` column, and contains samples | ✅ Available |
| Analysis registry | Defines supported EvoFlow analysis modules | ✅ Available |
| Continuous integration | Install, lint, and test across Python 3.10, 3.11, and 3.12 | ✅ Passing |
| Variant/sample QC engine | Genomic QC metrics, filtering, summaries, and plots | 🚧 Next milestone |
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

## Analysis modules

EvoFlow currently defines the following modular analysis vocabulary:

| Module | Intended role |
|---|---|
| `qc` | Variant and sample quality control |
| `pca` | Principal component analysis |
| `structure` | Population structure and ancestry inference |
| `diversity` | Population diversity statistics |
| `fst` | Population differentiation using FST |
| `spatial` | Spatial population-genomic analyses |
| `selection` | Selection and outlier scans |
| `gea` | Genotype-environment association |
| `report` | Reproducible analysis reporting |

Modules are being implemented independently so that users will eventually be able to run only the analyses required for a particular study while retaining one shared project configuration and output structure.

---

## Design principles

EvoFlow is being developed around several principles that are especially important in evolutionary genomics.

### 1. Reproducibility first

A computational result is only useful if another researcher can understand how it was produced. EvoFlow is being designed so that project inputs, configuration, module settings, software versions, and analysis provenance can be associated with each run.

### 2. Species agnostic

EvoFlow is not an amphibian-specific workflow. Salamander datasets may be used as development examples because they provide realistic evolutionary-genomic use cases, but the software architecture is intended for genomic datasets from any organism.

### 3. Modular rather than monolithic

Researchers should not have to execute an entire pipeline simply to obtain one analysis. QC, PCA, structure, diversity, FST, spatial analyses, selection, GEA, and reporting are treated as separate modules with explicit roles.

### 4. Transparent methods

EvoFlow is intended to orchestrate established methods rather than hide them behind an opaque interface. Parameters, dependencies, inputs, and outputs should remain inspectable.

### 5. Scientific outputs, not just software logs

The long-term goal is not merely successful command execution. EvoFlow should generate interpretable tables, publication-quality figures, reproducible methods information, and organized outputs that support biological inference.

---

## Architecture

The project uses a `src`-based Python package layout and separates configuration, input handling, analysis modules, orchestration, and the command-line interface.

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
│       │   └── validation.py
│       └── modules/
│           └── registry.py
├── tests/
├── CONTRIBUTING.md
├── LICENSE
├── README.md
└── pyproject.toml
```

### Architectural layers

- **`evoflow.config`** — project configuration and serialization.
- **`evoflow.io`** — validation of genomic files and associated biological metadata.
- **`evoflow.modules`** — analysis-module definitions and future implementations.
- **`evoflow.core`** — reserved for orchestration, execution, provenance, and run management as those capabilities are implemented.
- **`evoflow.cli`** — the user-facing command-line interface.

This separation is intended to keep the codebase maintainable as the number of genomic analyses grows.

---

## Installation

EvoFlow is currently intended for **development installation from source**.

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

### 1. Create an EvoFlow project configuration

```bash
evoflow init \
  --project salamander-demo \
  --vcf data/variants.vcf \
  --metadata data/samples.csv
```

This creates an `evoflow.yaml` configuration file.

### 2. Validate the project inputs

```bash
evoflow validate evoflow.yaml
```

The current validator checks that:

- the configured variant file exists;
- the configured metadata file exists;
- the metadata contains a `sample` column; and
- the metadata contains at least one sample.

### 3. Inspect the configured analysis plan

```bash
evoflow plan evoflow.yaml
```

The `plan` command reads the requested module list and prints the ordered analyses without executing them.

---

## Project configuration

An EvoFlow project is described with YAML.

Example:

```yaml
project: salamander-demo
vcf: data/variants.vcf
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

At this stage, the module list defines the intended workflow plan. Analysis execution will be connected to the registry as individual engines are implemented.

---

## Input data

### Variant data

EvoFlow is designed around standard population-genomic variant files such as VCF/BCF. Current validation confirms that the configured file path exists; deeper VCF inspection is part of the QC milestone.

### Sample metadata

The current minimum metadata requirement is a CSV file containing a `sample` column.

```csv
sample,population,site,latitude,longitude
sample_01,OH,Site_A,39.5100,-84.7300
sample_02,IN,Site_B,39.1600,-86.5200
sample_03,KY,Site_C,38.0400,-84.5000
```

Only `sample` is currently enforced by the validator. Columns such as `population`, `site`, `latitude`, `longitude`, and environmental covariates are intended to support later population, spatial, and landscape-genomic modules.

A central design requirement is that sample identifiers in metadata remain traceable to the individuals represented in the genomic dataset.

---

## Planned results structure

As the execution engine is implemented, EvoFlow is intended to organize results consistently rather than scattering output across tool-specific folders.

A target structure is:

```text
evoflow-results/
├── 00_provenance/
├── 01_qc/
├── 02_pca/
├── 03_structure/
├── 04_diversity/
├── 05_fst/
├── 06_spatial/
├── 07_selection/
├── 08_gea/
└── report/
```

This is a design target, not yet a promise of current generated output.

---

## Reproducibility and provenance

A major objective of EvoFlow is to make a genomic workflow auditable from input to interpretation.

The execution layer is being designed to eventually record information such as:

- input file paths and checksums;
- project configuration;
- selected modules;
- analysis parameters;
- EvoFlow version;
- external-tool versions;
- command history;
- generated files; and
- run timestamps.

These provenance features will be introduced alongside the workflow execution engine.

---

## Testing and continuous integration

Every change pushed to the repository is checked with GitHub Actions.

The current CI matrix tests EvoFlow on:

- Python 3.10
- Python 3.11
- Python 3.12

Each environment performs:

```text
Install  ->  Ruff linting  ->  Pytest
```

The CI badge at the top of this README reflects the current state of the `main` branch.

---

## Roadmap

### Phase 0 — Software foundation

- [x] Repository architecture
- [x] Python package structure
- [x] Command-line interface
- [x] YAML project configuration
- [x] Input-path and metadata validation
- [x] Modular analysis registry
- [x] Automated tests
- [x] Continuous integration
- [x] Project documentation foundation

### Phase 1 — Genomic quality control

- [ ] Read and summarize VCF content
- [ ] Variant counts
- [ ] Sample counts
- [ ] Per-sample missingness
- [ ] Per-site missingness
- [ ] Minor-allele-frequency summaries
- [ ] Depth / genotype-quality summaries where available
- [ ] Configurable QC filtering
- [ ] QC tables
- [ ] Publication-quality QC figures

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
- [ ] User documentation
- [ ] Containerized release
- [ ] Versioned software release
- [ ] DOI / archival release workflow

---

## Intended users

EvoFlow is being developed for:

- evolutionary biologists;
- population geneticists;
- conservation genomicists;
- landscape genomicists;
- molecular ecologists;
- bioinformaticians; and
- researchers who need a reproducible path from variant data to population-level evolutionary inference.

---

## Scientific scope

EvoFlow is intended for research workflows involving questions such as:

- How are individuals genetically structured across populations or geography?
- Which populations are most differentiated?
- How much genomic diversity is present within populations?
- Is genetic differentiation associated with geographic distance?
- Are there potential barriers or corridors to gene flow?
- Which genomic regions show unusual differentiation?
- Are genotypes associated with environmental variation?
- Can all of these analyses be reproduced from one documented project configuration?

The software will not replace biological interpretation. Its purpose is to make the computational path to that interpretation more coherent, transparent, and reproducible.

---

## Contributing

EvoFlow is under active development. Issues, suggestions, documentation improvements, test cases, and scientifically justified module proposals are welcome as the project matures.

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the current contribution guidance.

For substantial new analysis modules, the intended standard is that the implementation should clearly define:

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

Until one is available, researchers referring to the project should cite the repository and include the software version or commit used in their analysis. Formal citation instructions will be added with the first archival release.

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
