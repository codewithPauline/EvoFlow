from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from evoflow.io.vcf import open_vcf_text, read_vcf_samples
from evoflow.modules.pca_plots import generate_pca_figures


@dataclass(slots=True)
class PCAResult:
    samples: int
    variants_used: int
    components: int
    explained_variance_ratio: list[float]
    input_vcf: str
    method: str = "allele-frequency-standardized streaming Gram-matrix PCA"


def _parse_diploid_biallelic_dosage(gt: str) -> float | None:
    if not gt or gt in {".", "./.", ".|."}:
        return None

    if "/" in gt:
        parts = gt.split("/")
    elif "|" in gt:
        parts = gt.split("|")
    else:
        return None

    if len(parts) != 2 or any(part == "." for part in parts):
        return None

    try:
        alleles = [int(part) for part in parts]
    except ValueError:
        return None

    if any(allele not in {0, 1} for allele in alleles):
        return None
    return float(sum(alleles))


def _standardized_variant(fields: list[str], expected_samples: int) -> np.ndarray | None:
    if len(fields) < 10:
        return None

    ref = fields[3]
    alt_field = fields[4]
    alts = [] if alt_field == "." else alt_field.split(",")
    if len(ref) != 1 or len(alts) != 1 or len(alts[0]) != 1:
        return None

    format_keys = fields[8].split(":")
    if "GT" not in format_keys:
        return None
    gt_index = format_keys.index("GT")

    sample_fields = fields[9:]
    if len(sample_fields) != expected_samples:
        raise ValueError(
            f"VCF record {fields[0]}:{fields[1]} has {len(sample_fields)} sample columns; "
            f"expected {expected_samples}."
        )

    dosages = np.full(expected_samples, np.nan, dtype=float)
    for index, sample_value in enumerate(sample_fields):
        parts = sample_value.split(":")
        gt = parts[gt_index] if gt_index < len(parts) else "."
        dosage = _parse_diploid_biallelic_dosage(gt)
        if dosage is not None:
            dosages[index] = dosage

    called = dosages[~np.isnan(dosages)]
    if called.size == 0:
        return None

    p = float(called.sum() / (2.0 * called.size))
    if p <= 0.0 or p >= 1.0:
        return None

    expected_dosage = 2.0 * p
    scale = math.sqrt(2.0 * p * (1.0 - p))
    dosages[np.isnan(dosages)] = expected_dosage
    return (dosages - expected_dosage) / scale


def _iter_standardized_variants(vcf: Path, expected_samples: int):
    with open_vcf_text(vcf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError("Encountered a malformed VCF record with fewer than 8 columns.")
            standardized = _standardized_variant(fields, expected_samples)
            if standardized is not None:
                yield fields[0], fields[1], fields[3], fields[4], standardized


def _read_metadata(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = {row["sample"]: row for row in reader}
    return fieldnames, rows


def run_pca(
    vcf: str | Path,
    metadata: str | Path,
    output_dir: str | Path,
    *,
    n_components: int = 10,
) -> PCAResult:
    """Run allele-frequency-standardized PCA using a streaming sample Gram matrix."""
    if n_components < 1:
        raise ValueError("n_components must be at least 1.")

    vcf_path = Path(vcf)
    metadata_path = Path(metadata)
    sample_names = read_vcf_samples(vcf_path)
    if len(sample_names) < 2:
        raise ValueError("PCA requires at least two samples.")

    gram = np.zeros((len(sample_names), len(sample_names)), dtype=float)
    variants_used = 0
    for _chrom, _pos, _ref, _alt, standardized in _iter_standardized_variants(
        vcf_path, len(sample_names)
    ):
        gram += np.outer(standardized, standardized)
        variants_used += 1

    if variants_used < 2:
        raise ValueError("PCA requires at least two informative biallelic SNPs.")

    all_eigenvalues, all_eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(all_eigenvalues)[::-1]
    all_eigenvalues = np.clip(all_eigenvalues[order], 0.0, None)
    all_eigenvectors = all_eigenvectors[:, order]

    positive = all_eigenvalues > np.finfo(float).eps
    positive_count = int(positive.sum())
    if positive_count == 0:
        raise ValueError("PCA could not find any non-zero genomic variance.")

    components = min(n_components, positive_count, len(sample_names) - 1, variants_used)
    eigenvalues = all_eigenvalues[:components]
    eigenvectors = all_eigenvectors[:, :components]
    singular_values = np.sqrt(eigenvalues)
    scores = eigenvectors * singular_values

    total_variance = float(all_eigenvalues.sum())
    explained_ratio = eigenvalues / total_variance if total_variance else np.zeros_like(eigenvalues)
    explained_variance = eigenvalues / (len(sample_names) - 1)
    cumulative = np.cumsum(explained_ratio)

    pca_dir = Path(output_dir) / "pca"
    pca_dir.mkdir(parents=True, exist_ok=True)

    metadata_fields, metadata_rows = _read_metadata(metadata_path)
    extra_fields = [field for field in metadata_fields if field != "sample"]
    score_fields = ["sample", *extra_fields, *[f"PC{i}" for i in range(1, components + 1)]]
    with (pca_dir / "pca_scores.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=score_fields)
        writer.writeheader()
        for sample_index, sample in enumerate(sample_names):
            row: dict[str, str | float] = {"sample": sample}
            metadata_row = metadata_rows.get(sample, {})
            for field in extra_fields:
                row[field] = metadata_row.get(field, "")
            for component_index in range(components):
                row[f"PC{component_index + 1}"] = round(
                    float(scores[sample_index, component_index]), 8
                )
            writer.writerow(row)

    with (pca_dir / "pca_variance.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "component",
            "eigenvalue",
            "explained_variance",
            "explained_variance_ratio",
            "cumulative_explained_variance",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(components):
            writer.writerow(
                {
                    "component": f"PC{index + 1}",
                    "eigenvalue": round(float(eigenvalues[index]), 8),
                    "explained_variance": round(float(explained_variance[index]), 8),
                    "explained_variance_ratio": round(float(explained_ratio[index]), 8),
                    "cumulative_explained_variance": round(float(cumulative[index]), 8),
                }
            )

    loading_fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        *[f"PC{i}" for i in range(1, components + 1)],
    ]
    with (pca_dir / "pca_loadings.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=loading_fields)
        writer.writeheader()
        for chrom, pos, ref, alt, standardized in _iter_standardized_variants(
            vcf_path, len(sample_names)
        ):
            row: dict[str, str | float] = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
            }
            for component_index in range(components):
                singular = singular_values[component_index]
                loading = (
                    float(np.dot(eigenvectors[:, component_index], standardized) / singular)
                    if singular > 0
                    else 0.0
                )
                row[f"PC{component_index + 1}"] = round(loading, 8)
            writer.writerow(row)

    result = PCAResult(
        samples=len(sample_names),
        variants_used=variants_used,
        components=components,
        explained_variance_ratio=[round(float(value), 8) for value in explained_ratio],
        input_vcf=str(vcf_path),
    )
    (pca_dir / "pca_summary.json").write_text(
        json.dumps(asdict(result), indent=2) + "\n", encoding="utf-8"
    )
    generate_pca_figures(pca_dir)
    return result
