# EvoFlow architecture

EvoFlow separates configuration, input validation, analysis modules, orchestration, and reporting so each layer can evolve independently.

## Layers

- `evoflow.config`: project configuration and serialization.
- `evoflow.io`: VCF, metadata, coordinate, and environmental-data validation.
- `evoflow.modules`: analysis-module registry and module implementations.
- `evoflow.core`: orchestration, run manifests, provenance, and execution.
- `evoflow.cli`: user-facing command-line interface.

Analysis modules should expose explicit inputs, outputs, parameters, dependency checks, and provenance metadata.
