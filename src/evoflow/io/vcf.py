from __future__ import annotations

import gzip
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO


@contextmanager
def open_vcf_text(path: str | Path) -> Iterator[TextIO]:
    """Open a text VCF or bgzip/gzip-compressed VCF for streaming reads."""
    vcf_path = Path(path)
    lower_name = vcf_path.name.lower()

    if lower_name.endswith(".bcf"):
        raise ValueError(
            "BCF parsing is not available in the native EvoFlow engines yet. "
            "Convert BCF to VCF/VCF.gz before running EvoFlow."
        )

    if lower_name.endswith((".vcf.gz", ".vcf.bgz", ".gz", ".bgz")):
        with gzip.open(vcf_path, mode="rt", encoding="utf-8") as handle:
            yield handle
        return

    with vcf_path.open(mode="rt", encoding="utf-8") as handle:
        yield handle


def read_vcf_samples(path: str | Path) -> list[str]:
    """Return sample IDs from the #CHROM header of a VCF."""
    with open_vcf_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                return fields[9:] if len(fields) > 9 else []

    raise ValueError("VCF is missing the required #CHROM header line.")


def parse_diploid_biallelic_dosage(gt: str) -> float | None:
    """Return alternate-allele dosage 0/1/2 for a usable diploid biallelic GT."""
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


def extract_biallelic_snp_dosages(
    fields: list[str], expected_samples: int
) -> list[float | None] | None:
    """Extract diploid dosages from one biallelic SNP VCF record."""
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

    dosages: list[float | None] = []
    for sample_value in sample_fields:
        parts = sample_value.split(":")
        gt = parts[gt_index] if gt_index < len(parts) else "."
        dosages.append(parse_diploid_biallelic_dosage(gt))
    return dosages
