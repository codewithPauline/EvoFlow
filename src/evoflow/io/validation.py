from __future__ import annotations

import csv

from evoflow.config import EvoFlowConfig


def validate_config(config: EvoFlowConfig) -> list[str]:
    """Return validation messages or raise for fatal input problems."""
    messages: list[str] = []

    if not config.vcf.exists():
        raise FileNotFoundError(f"VCF/BCF not found: {config.vcf}")
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

    messages.append(f"Metadata contains {len(rows)} samples")
    return messages
