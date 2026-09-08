from __future__ import annotations

MODULES = {
    "qc": "Variant and sample quality control",
    "pca": "Principal component analysis",
    "structure": "Population structure / ancestry inference",
    "diversity": "Population diversity statistics",
    "fst": "Population differentiation (FST)",
    "spatial": "Spatial population-genomic analyses",
    "selection": "Selection and outlier scans",
    "gea": "Genotype-environment association",
    "report": "Reproducible analysis report",
}


def validate_modules(modules: list[str]) -> None:
    unknown = sorted(set(modules) - MODULES.keys())
    if unknown:
        raise ValueError(f"Unknown EvoFlow modules: {', '.join(unknown)}")
