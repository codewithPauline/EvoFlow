from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from evoflow.io.vcf import extract_biallelic_snp_dosages, open_vcf_text, read_vcf_samples


@dataclass(slots=True)
class LDPruneResult:
    samples: int
    records_seen: int
    informative_snps: int
    retained_snps: int
    removed_for_ld: int
    skipped_noninformative: int
    r2_threshold: float
    window_bp: int
    min_overlap: int
    input_vcf: str


@dataclass(slots=True)
class _RetainedVariant:
    chrom: str
    pos: int
    dosages: list[float | None]


def _is_informative(dosages: list[float | None]) -> bool:
    called = np.array([value for value in dosages if value is not None], dtype=float)
    return called.size >= 2 and float(np.var(called)) > 0.0


def pairwise_r2(
    left: list[float | None],
    right: list[float | None],
    *,
    min_overlap: int = 3,
) -> float | None:
    """Calculate genotype-dosage r² using samples called at both variants."""
    pairs = [
        (left_value, right_value)
        for left_value, right_value in zip(left, right, strict=True)
        if left_value is not None and right_value is not None
    ]
    if len(pairs) < min_overlap:
        return None

    x = np.array([pair[0] for pair in pairs], dtype=float)
    y = np.array([pair[1] for pair in pairs], dtype=float)
    if float(np.var(x)) == 0.0 or float(np.var(y)) == 0.0:
        return None

    correlation = float(np.corrcoef(x, y)[0, 1])
    if not np.isfinite(correlation):
        return None
    return correlation * correlation


def run_ld_prune(
    vcf: str | Path,
    output_dir: str | Path,
    *,
    r2_threshold: float = 0.2,
    window_bp: int = 50_000,
    min_overlap: int = 3,
) -> LDPruneResult:
    """Greedily prune correlated biallelic SNPs within a physical sliding window."""
    if not 0.0 <= r2_threshold <= 1.0:
        raise ValueError("r2_threshold must be between 0 and 1.")
    if window_bp < 1:
        raise ValueError("window_bp must be at least 1.")
    if min_overlap < 2:
        raise ValueError("min_overlap must be at least 2.")

    vcf_path = Path(vcf)
    sample_names = read_vcf_samples(vcf_path)
    if len(sample_names) < min_overlap:
        raise ValueError(
            f"LD pruning requires at least {min_overlap} samples for the configured overlap."
        )

    ld_dir = Path(output_dir) / "ld"
    ld_dir.mkdir(parents=True, exist_ok=True)
    output_vcf = ld_dir / "ld_pruned.vcf"
    removed_csv = ld_dir / "ld_removed.csv"

    records_seen = 0
    informative_snps = 0
    retained_snps = 0
    removed_for_ld = 0
    skipped_noninformative = 0
    retained_window: list[_RetainedVariant] = []
    current_chrom: str | None = None
    previous_pos: int | None = None

    with (
        open_vcf_text(vcf_path) as source,
        output_vcf.open("w", encoding="utf-8") as destination,
        removed_csv.open("w", newline="", encoding="utf-8") as removed_handle,
    ):
        removed_writer = csv.DictWriter(
            removed_handle,
            fieldnames=["chrom", "pos", "blocked_by_pos", "r2"],
        )
        removed_writer.writeheader()

        for line in source:
            if line.startswith("##"):
                destination.write(line)
                continue
            if line.startswith("#CHROM"):
                destination.write(f"##evoflow_ld_r2_threshold={r2_threshold}\n")
                destination.write(f"##evoflow_ld_window_bp={window_bp}\n")
                destination.write(f"##evoflow_ld_min_overlap={min_overlap}\n")
                destination.write(line)
                continue
            if line.startswith("#") or not line.strip():
                destination.write(line)
                continue

            records_seen += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError("Encountered a malformed VCF record with fewer than 8 columns.")

            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"VCF record has a non-integer position: {fields[1]}") from exc

            if chrom != current_chrom:
                current_chrom = chrom
                previous_pos = None
                retained_window = []
            elif previous_pos is not None and pos < previous_pos:
                raise ValueError(
                    f"VCF must be position-sorted within each chromosome for LD pruning: "
                    f"{chrom}:{pos} follows {previous_pos}."
                )
            previous_pos = pos

            dosages = extract_biallelic_snp_dosages(fields, len(sample_names))
            if dosages is None or not _is_informative(dosages):
                skipped_noninformative += 1
                continue

            informative_snps += 1
            retained_window = [
                variant for variant in retained_window if pos - variant.pos <= window_bp
            ]

            blocker: _RetainedVariant | None = None
            blocker_r2: float | None = None
            for retained in retained_window:
                r2 = pairwise_r2(retained.dosages, dosages, min_overlap=min_overlap)
                if r2 is not None and r2 >= r2_threshold:
                    blocker = retained
                    blocker_r2 = r2
                    break

            if blocker is not None:
                removed_for_ld += 1
                removed_writer.writerow(
                    {
                        "chrom": chrom,
                        "pos": pos,
                        "blocked_by_pos": blocker.pos,
                        "r2": round(float(blocker_r2), 8),
                    }
                )
                continue

            destination.write(line)
            retained_snps += 1
            retained_window.append(_RetainedVariant(chrom=chrom, pos=pos, dosages=dosages))

    result = LDPruneResult(
        samples=len(sample_names),
        records_seen=records_seen,
        informative_snps=informative_snps,
        retained_snps=retained_snps,
        removed_for_ld=removed_for_ld,
        skipped_noninformative=skipped_noninformative,
        r2_threshold=r2_threshold,
        window_bp=window_bp,
        min_overlap=min_overlap,
        input_vcf=str(vcf_path),
    )
    (ld_dir / "ld_summary.json").write_text(
        json.dumps(asdict(result), indent=2) + "\n", encoding="utf-8"
    )
    return result
