# EvoFlow architecture

EvoFlow is organized as a modular scientific-software system rather than a collection of analysis scripts. Configuration, input access, validation, analysis engines, visualization, orchestration, and reporting are separated so that each layer can evolve without tightly coupling the rest of the workflow.

## Package layers

- `evoflow.config` — project configuration and YAML serialization.
- `evoflow.io` — genomic file access and biological metadata validation.
- `evoflow.modules` — analysis engines, module-specific outputs, and visualization code.
- `evoflow.core` — reserved for workflow orchestration, run manifests, provenance, dependency resolution, and restartable execution.
- `evoflow.cli` — the user-facing command-line interface.

## Current QC data flow

```text
          evoflow.yaml
               |
               v
       EvoFlowConfig
               |
               v
    strict input validation
     /                   \
    v                     v
VCF / VCF.gz        metadata CSV
    |                     |
    +----------+----------+
               |
               v
       streaming VCF parser
               |
      +--------+---------+
      |        |         |
      v        v         v
 sample QC  variant QC  run QC
      |        |         |
      +--------+---------+
               |
       +-------+--------+
       |                |
       v                v
 threshold filter    diagnostics
       |                |
       v                v
 filtered.vcf       PNG + PDF
       |
       +----------------------+
               |
               v
      structured QC outputs
```

The QC engine streams VCF records instead of materializing the full file in memory. Global depth and genotype-quality statistics use running sums/counts, while MAF and missingness distributions use fixed-size histogram bins. This keeps the basic QC memory footprint driven primarily by sample-level state rather than total variant count.

## QC responsibilities

`evoflow.io.vcf` is responsible for text VCF access and sample-header extraction. It supports uncompressed VCF and gzip/bgzip-compatible text streams. Native BCF support is intentionally not claimed yet.

`evoflow.io.validation` verifies that required files exist and that VCF and metadata sample identifiers are non-empty, unique, and exactly concordant.

`evoflow.modules.qc` performs streaming genotype/variant summaries, threshold decisions, filtered VCF writing, and machine-readable QC outputs.

`evoflow.modules.qc_plots` converts the streaming summaries into diagnostic figures without requiring the complete per-variant metric arrays to be retained in memory.

## Module contract

As new modules are implemented, each should define:

1. required inputs and validation rules;
2. configurable parameters and defaults;
3. deterministic output locations;
4. machine-readable result tables;
5. scientific figures where appropriate;
6. dependency and software-version information;
7. provenance needed to reproduce the analysis; and
8. tests covering both expected results and failure modes.

## Next architectural step

The next major layer is a PCA engine that consumes QC-filtered biallelic genotypes, creates a sample-by-variant dosage matrix with documented missing-data handling, computes principal components, joins sample metadata, and writes both numerical scores and reproducible visualizations. Later, `evoflow.core` will coordinate these modules as an ordered workflow rather than requiring separate CLI invocations.
