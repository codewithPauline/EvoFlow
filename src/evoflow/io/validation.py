from __future__ import annotations

import csv

from evoflow.config import EvoFlowConfig
from evoflow.io.vcf import read_vcf_samples


def _preview(values: list[str], limit: int = 5) -> str:
    shown = ", ".join(values[:limit])
    return shown if len(values) <= limit else f"{shown}, ..."


def validate_config(config: EvoFlowConfig) -> list[str]:
    """Return validation messages or raise for fatal input problems."""
    messages: list[str] = []

    if not config.vcf.exists():
        raise FileNotFoundError(f"VCF not found: {config.vcf}")
    messages.append(f"Variant file found: {config.vcf}")

    if not config.metadata.exists():
        raise FileNotFoundError(f"Metadata file not found: {config.metadata}")
    messages.append(f"Metadata file found: {config.metadata}")

    with config.metadata.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        if "sample" not in fields:
            raise ValueError("Metadata must contain a 'sample' column.")
        rows = list(reader)
        if not rows:
            raise ValueError("Metadata contains no samples.")

    metadata_samples = [str(row.get("sample", "")).strip() for row in rows]
    if any(not sample for sample in metadata_samples):
        raise ValueError("Metadata contains one or more blank sample IDs.")
    if len(set(metadata_samples)) != len(metadata_samples):
        raise ValueError("Metadata sample IDs must be unique.")
    messages.append(f"Metadata contains {len(metadata_samples)} samples")

    vcf_samples = read_vcf_samples(config.vcf)
    if len(set(vcf_samples)) != len(vcf_samples):
        raise ValueError("VCF sample IDs must be unique.")
    messages.append(f"VCF contains {len(vcf_samples)} samples")

    metadata_set = set(metadata_samples)
    vcf_set = set(vcf_samples)
    missing_metadata = sorted(vcf_set - metadata_set)
    extra_metadata = sorted(metadata_set - vcf_set)
    if missing_metadata or extra_metadata:
        details: list[str] = []
        if missing_metadata:
            details.append(f"missing from metadata: {_preview(missing_metadata)}")
        if extra_metadata:
            details.append(f"not present in VCF: {_preview(extra_metadata)}")
        raise ValueError("VCF/metadata sample mismatch (" + "; ".join(details) + ").")

    messages.append("VCF and metadata sample IDs match exactly")
    return messages
