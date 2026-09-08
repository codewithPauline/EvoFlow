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
            "BCF parsing is not available in the native QC engine yet. "
            "Convert BCF to VCF/VCF.gz before running EvoFlow QC."
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
