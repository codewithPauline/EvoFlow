from __future__ import annotations

from pathlib import Path

import pytest

from evoflow.config import EvoFlowConfig
from evoflow.io.validation import validate_config


VCF_TEXT = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0
"""


def _config(tmp_path: Path, metadata_text: str) -> EvoFlowConfig:
    vcf = tmp_path / "tiny.vcf"
    metadata = tmp_path / "samples.csv"
    vcf.write_text(VCF_TEXT, encoding="utf-8")
    metadata.write_text(metadata_text, encoding="utf-8")
    return EvoFlowConfig(project="test", vcf=vcf, metadata=metadata)


def test_validation_requires_exact_sample_match(tmp_path: Path) -> None:
    config = _config(tmp_path, "sample,population\nS1,A\nS2,B\n")
    messages = validate_config(config)
    assert "VCF and metadata sample IDs match exactly" in messages


def test_validation_rejects_sample_mismatch(tmp_path: Path) -> None:
    config = _config(tmp_path, "sample,population\nS1,A\nS3,B\n")
    with pytest.raises(ValueError, match="VCF/metadata sample mismatch"):
        validate_config(config)


def test_validation_rejects_duplicate_metadata_ids(tmp_path: Path) -> None:
    config = _config(tmp_path, "sample,population\nS1,A\nS1,B\n")
    with pytest.raises(ValueError, match="must be unique"):
        validate_config(config)
