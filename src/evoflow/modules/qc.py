from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from evoflow.io.vcf import open_vcf_text, read_vcf_samples


@dataclass(slots=True)
class SampleQC:
    sample: str
    called_genotypes: int = 0
    missing_genotypes: int = 0
    heterozygous_genotypes: int = 0
    depth_sum: float = 0.0
    depth_observations: int = 0

    def to_row(self) -> dict[str, str | int | float]:
        total = self.called_genotypes + self.missing_genotypes
        call_rate = self.called_genotypes / total if total else 0.0
        heterozygosity = (
            self.heterozygous_genotypes / self.called_genotypes
            if self.called_genotypes
            else 0.0
        )
        mean_depth = self.depth_sum / self.depth_observations if self.depth_observations else 0.0
        return {
            "sample": self.sample,
            "called_genotypes": self.called_genotypes,
            "missing_genotypes": self.missing_genotypes,
            "call_rate": round(call_rate, 6),
            "missing_rate": round(1.0 - call_rate, 6),
            "heterozygosity": round(heterozygosity, 6),
            "mean_depth": round(mean_depth, 6),
        }


@dataclass(slots=True)
class QCSummary:
    samples: int
    variants: int
    retained_variants: int
    snps: int
    biallelic_snps: int
    multiallelic_variants: int
    transitions: int
    transversions: int
    ti_tv_ratio: float | None
    genotype_calls: int
    called_genotypes: int
    missing_genotypes: int
    overall_call_rate: float
    mean_depth: float | None
    min_maf: float
    max_missing: float


_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def _parse_genotype(gt: str) -> list[int] | None:
    if not gt or gt in {".", "./.", ".|."}:
        return None
    separator = "/" if "/" in gt else "|" if "|" in gt else None
    parts = gt.split(separator) if separator else [gt]
    if any(part == "." for part in parts):
        return None
    try:
        return [int(part) for part in parts]
    except ValueError:
        return None


def _safe_float(value: str) -> float | None:
    if not value or value == ".":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def run_qc(
    vcf: str | Path,
    output_dir: str | Path,
    *,
    min_maf: float = 0.0,
    max_missing: float = 1.0,
) -> QCSummary:
    """Stream a VCF and write variant-, sample-, and run-level QC summaries."""
    if not 0.0 <= min_maf <= 0.5:
        raise ValueError("min_maf must be between 0.0 and 0.5.")
    if not 0.0 <= max_missing <= 1.0:
        raise ValueError("max_missing must be between 0.0 and 1.0.")

    vcf_path = Path(vcf)
    qc_dir = Path(output_dir) / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    sample_names = read_vcf_samples(vcf_path)
    sample_stats = [SampleQC(sample=name) for name in sample_names]

    variants = 0
    retained_variants = 0
    snps = 0
    biallelic_snps = 0
    multiallelic_variants = 0
    transitions = 0
    transversions = 0
    called_genotypes = 0
    missing_genotypes = 0
    depth_sum = 0.0
    depth_observations = 0

    variant_path = qc_dir / "variant_qc.csv"
    with variant_path.open("w", newline="", encoding="utf-8") as variant_handle:
        writer = csv.DictWriter(
            variant_handle,
            fieldnames=[
                "chrom",
                "pos",
                "ref",
                "alt",
                "is_snp",
                "is_biallelic",
                "call_rate",
                "missing_rate",
                "alt_allele_frequency",
                "maf",
                "mean_depth",
                "passes_qc",
            ],
        )
        writer.writeheader()

        with open_vcf_text(vcf_path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ValueError("Encountered a malformed VCF record with fewer than 8 columns.")

                chrom, pos, _variant_id, ref, alt_field = fields[:5]
                alts = [] if alt_field == "." else alt_field.split(",")
                is_snp = bool(alts) and len(ref) == 1 and all(len(alt) == 1 for alt in alts)
                is_biallelic = len(alts) == 1

                variants += 1
                if is_snp:
                    snps += 1
                if is_snp and is_biallelic:
                    biallelic_snps += 1
                    pair = (ref.upper(), alts[0].upper())
                    if pair in _TRANSITIONS:
                        transitions += 1
                    else:
                        transversions += 1
                if len(alts) > 1:
                    multiallelic_variants += 1

                format_keys = fields[8].split(":") if len(fields) > 8 else []
                sample_fields = fields[9:] if len(fields) > 9 else []
                gt_index = format_keys.index("GT") if "GT" in format_keys else None
                dp_index = format_keys.index("DP") if "DP" in format_keys else None

                site_called = 0
                site_missing = 0
                allele_counts = [0] * (len(alts) + 1)
                site_depth_sum = 0.0
                site_depth_observations = 0

                for index, sample_value in enumerate(sample_fields):
                    parts = sample_value.split(":")
                    gt = parts[gt_index] if gt_index is not None and gt_index < len(parts) else "."
                    genotype = _parse_genotype(gt)
                    sample_stat = sample_stats[index] if index < len(sample_stats) else None

                    if genotype is None:
                        site_missing += 1
                        missing_genotypes += 1
                        if sample_stat is not None:
                            sample_stat.missing_genotypes += 1
                    else:
                        site_called += 1
                        called_genotypes += 1
                        if sample_stat is not None:
                            sample_stat.called_genotypes += 1
                            if len(set(genotype)) > 1:
                                sample_stat.heterozygous_genotypes += 1
                        for allele in genotype:
                            if 0 <= allele < len(allele_counts):
                                allele_counts[allele] += 1

                    if dp_index is not None and dp_index < len(parts):
                        depth = _safe_float(parts[dp_index])
                        if depth is not None:
                            site_depth_sum += depth
                            site_depth_observations += 1
                            depth_sum += depth
                            depth_observations += 1
                            if sample_stat is not None:
                                sample_stat.depth_sum += depth
                                sample_stat.depth_observations += 1

                expected_samples = len(sample_names)
                if expected_samples and len(sample_fields) != expected_samples:
                    raise ValueError(
                        f"VCF record {chrom}:{pos} has {len(sample_fields)} sample columns; "
                        f"expected {expected_samples}."
                    )

                site_total = site_called + site_missing
                call_rate = site_called / site_total if site_total else 0.0
                missing_rate = 1.0 - call_rate
                total_alleles = sum(allele_counts)
                alt_count = sum(allele_counts[1:])
                alt_frequency = alt_count / total_alleles if total_alleles else None

                maf: float | None = None
                if is_biallelic and total_alleles:
                    ref_frequency = allele_counts[0] / total_alleles
                    maf = min(ref_frequency, 1.0 - ref_frequency)

                passes_missing = missing_rate <= max_missing
                passes_maf = maf is None or maf >= min_maf
                passes_qc = passes_missing and passes_maf
                if passes_qc:
                    retained_variants += 1

                site_mean_depth = (
                    site_depth_sum / site_depth_observations if site_depth_observations else None
                )
                writer.writerow(
                    {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt_field,
                        "is_snp": is_snp,
                        "is_biallelic": is_biallelic,
                        "call_rate": round(call_rate, 6),
                        "missing_rate": round(missing_rate, 6),
                        "alt_allele_frequency": (
                            "" if alt_frequency is None else round(alt_frequency, 6)
                        ),
                        "maf": "" if maf is None else round(maf, 6),
                        "mean_depth": (
                            "" if site_mean_depth is None else round(site_mean_depth, 6)
                        ),
                        "passes_qc": passes_qc,
                    }
                )

    sample_path = qc_dir / "sample_qc.csv"
    with sample_path.open("w", newline="", encoding="utf-8") as sample_handle:
        fieldnames = [
            "sample",
            "called_genotypes",
            "missing_genotypes",
            "call_rate",
            "missing_rate",
            "heterozygosity",
            "mean_depth",
        ]
        writer = csv.DictWriter(sample_handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(stat.to_row() for stat in sample_stats)

    genotype_calls = called_genotypes + missing_genotypes
    overall_call_rate = called_genotypes / genotype_calls if genotype_calls else 0.0
    ti_tv_ratio = transitions / transversions if transversions else None
    mean_depth = depth_sum / depth_observations if depth_observations else None

    summary = QCSummary(
        samples=len(sample_names),
        variants=variants,
        retained_variants=retained_variants,
        snps=snps,
        biallelic_snps=biallelic_snps,
        multiallelic_variants=multiallelic_variants,
        transitions=transitions,
        transversions=transversions,
        ti_tv_ratio=None if ti_tv_ratio is None else round(ti_tv_ratio, 6),
        genotype_calls=genotype_calls,
        called_genotypes=called_genotypes,
        missing_genotypes=missing_genotypes,
        overall_call_rate=round(overall_call_rate, 6),
        mean_depth=None if mean_depth is None else round(mean_depth, 6),
        min_maf=min_maf,
        max_missing=max_missing,
    )

    summary_path = qc_dir / "qc_summary.json"
    summary_path.write_text(json.dumps(asdict(summary), indent=2) + "\n", encoding="utf-8")
    return summary
