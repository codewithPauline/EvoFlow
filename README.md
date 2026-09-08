# EvoFlow

[![CI](https://github.com/codewithPauline/EvoFlow/actions/workflows/ci.yml/badge.svg)](https://github.com/codewithPauline/EvoFlow/actions/workflows/ci.yml)

> **From variants to evolutionary insight.**

EvoFlow is an open-source workflow platform for reproducible population and landscape genomics. It is designed to take researchers from variant data and sample metadata to quality control, population structure, differentiation, diversity, spatial analyses, selection scans, and publication-ready outputs through one consistent interface.

## Why EvoFlow?

Population-genomic studies often require researchers to stitch together multiple command-line programs, R scripts, file formats, and plotting workflows. EvoFlow aims to provide a transparent orchestration layer that keeps those steps reproducible while preserving access to the underlying tools and parameters.

## Planned workflow

```text
VCF / BCF + metadata + coordinates
            |
            v
     Input validation
            |
            v
     Variant & sample QC
            |
   +--------+---------+------------------+
   |                  |                  |
   v                  v                  v
  PCA          Population structure   Diversity
   |                  |                  |
   +---------+--------+------------------+
             |
             v
      Differentiation / FST
             |
             v
      Spatial population genomics
             |
             v
      Selection / GEA modules
             |
             v
   Figures + tables + reproducible report
```

## Current status

EvoFlow is in active early development. The first milestone establishes:

- a Python package and command-line interface;
- validated project configuration;
- VCF and metadata input checks;
- a modular analysis registry;
- reproducible run directories and manifests;
- automated tests and continuous integration across Python 3.10, 3.11, and 3.12.

## Installation (development)

```bash
git clone https://github.com/codewithPauline/EvoFlow.git
cd EvoFlow
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install -e ".[dev]"
```

## Quick start

Create a project configuration:

```bash
evoflow init --project salamander-demo --vcf data/variants.vcf --metadata data/samples.csv
```

Validate inputs:

```bash
evoflow validate evoflow.yaml
```

Inspect the analysis plan:

```bash
evoflow plan evoflow.yaml
```

## Input metadata

A minimal metadata table should include:

```csv
sample,population
sample_01,OH
sample_02,IN
sample_03,KY
```

Optional columns such as `latitude`, `longitude`, `site`, and other biological covariates can support spatial and landscape-genomic modules.

## Design principles

1. **Reproducibility first** — every run records inputs, configuration, software version, and module settings.
2. **Modular analyses** — users can run only the analyses they need.
3. **Tool transparency** — EvoFlow orchestrates established population-genomic tools rather than hiding them.
4. **Species agnostic** — example datasets may come from salamanders, but EvoFlow is designed for any organism.
5. **Publication-ready outputs** — tables, figures, and reports should be ready for scientific interpretation and refinement.

## Roadmap

- [x] Project architecture
- [x] CLI scaffold
- [x] Configuration and validation layer
- [ ] VCF QC module
- [ ] PCA module
- [ ] Population-structure module
- [ ] FST and diversity statistics
- [ ] Isolation-by-distance and spatial summaries
- [ ] Selection scans
- [ ] Genotype-environment association
- [ ] Workflow execution engine
- [ ] HTML report generation
- [ ] Example biological dataset
- [ ] Containerized release

## Intended audience

EvoFlow is being developed for evolutionary biologists, population geneticists, conservation genomicists, and bioinformatics researchers who want a reproducible path from genomic variants to evolutionary interpretation.

## License

MIT License. See `LICENSE`.
