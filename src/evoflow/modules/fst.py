from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path

from evoflow.io.metadata import population_indices
from evoflow.io.vcf import extract_biallelic_snp_dosages, open_vcf_text, read_vcf_samples
from evoflow.modules.fst_plots import generate_fst_heatmap


@dataclass(slots=True)
class FSTResult:
    populations: int
    population_pairs: int
    samples: int
    records_seen: int
    biallelic_snps: int
    pair_site_records: int
    input_vcf: str
    estimator: str = "Hudson FST (Hudson 1992; Bhatia et al. 2013 formulation)"


@dataclass(slots=True)
class _PairAccumulator:
    numerator_sum: float = 0.0
    denominator_sum: float = 0.0
    sites_used: int = 0


def _allele_counts(
    dosages: list[float | None], indices: list[int]
) -> tuple[int, int] | None:
    called = [dosages[index] for index in indices if dosages[index] is not None]
    if not called:
        return None
    alt = round(sum(called))
    total = 2 * len(called)
    ref = total - alt
    return ref, alt


def hudson_fst_components(
    ac1: tuple[int, int], ac2: tuple[int, int]
) -> tuple[float, float] | None:
    """Return per-site Hudson FST numerator and denominator from allele counts."""
    ref1, alt1 = ac1
    ref2, alt2 = ac2
    an1 = ref1 + alt1
    an2 = ref2 + alt2
    if an1 <= 1 or an2 <= 1:
        return None

    within1 = 2.0 * ref1 * alt1 / (an1 * (an1 - 1))
    within2 = 2.0 * ref2 * alt2 / (an2 * (an2 - 1))
    within = (within1 + within2) / 2.0
    between = (ref1 * alt2 + alt1 * ref2) / (an1 * an2)
    return between - within, between


def run_fst(
    vcf: str | Path,
    metadata: str | Path,
    output_dir: str | Path,
    *,
    population_column: str = "population",
) -> FSTResult:
    """Estimate pairwise Hudson FST across diploid biallelic SNP records."""
    vcf_path = Path(vcf)
    sample_names = read_vcf_samples(vcf_path)
    populations = population_indices(
        metadata,
        sample_names,
        population_column=population_column,
    )
    if len(populations) < 2:
        raise ValueError("FST requires at least two populations.")

    pairs = list(combinations(sorted(populations), 2))
    accumulators = {pair: _PairAccumulator() for pair in pairs}

    fst_dir = Path(output_dir) / "fst"
    fst_dir.mkdir(parents=True, exist_ok=True)
    site_path = fst_dir / "site_fst.csv"

    records_seen = 0
    biallelic_snps = 0
    pair_site_records = 0

    with site_path.open("w", newline="", encoding="utf-8") as site_handle:
        fieldnames = [
            "chrom",
            "pos",
            "population_1",
            "population_2",
            "numerator",
            "denominator",
            "fst",
        ]
        writer = csv.DictWriter(site_handle, fieldnames=fieldnames)
        writer.writeheader()

        with open_vcf_text(vcf_path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                records_seen += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ValueError(
                        "Encountered a malformed VCF record with fewer than 8 columns."
                    )

                dosages = extract_biallelic_snp_dosages(fields, len(sample_names))
                if dosages is None:
                    continue
                biallelic_snps += 1

                allele_counts = {
                    population: _allele_counts(dosages, indices)
                    for population, indices in populations.items()
                }
                for pair in pairs:
                    ac1 = allele_counts[pair[0]]
                    ac2 = allele_counts[pair[1]]
                    if ac1 is None or ac2 is None:
                        continue
                    components = hudson_fst_components(ac1, ac2)
                    if components is None:
                        continue
                    numerator, denominator = components
                    fst_value = numerator / denominator if denominator > 0.0 else None

                    writer.writerow(
                        {
                            "chrom": fields[0],
                            "pos": fields[1],
                            "population_1": pair[0],
                            "population_2": pair[1],
                            "numerator": round(numerator, 10),
                            "denominator": round(denominator, 10),
                            "fst": "" if fst_value is None else round(fst_value, 10),
                        }
                    )
                    pair_site_records += 1

                    if denominator > 0.0:
                        accumulator = accumulators[pair]
                        accumulator.numerator_sum += numerator
                        accumulator.denominator_sum += denominator
                        accumulator.sites_used += 1

    pairwise_path = fst_dir / "pairwise_fst.csv"
    pairwise_values: dict[tuple[str, str], float | None] = {}
    with pairwise_path.open("w", newline="", encoding="utf-8") as pairwise_handle:
        fieldnames = [
            "population_1",
            "population_2",
            "sites_used",
            "numerator_sum",
            "denominator_sum",
            "fst",
        ]
        writer = csv.DictWriter(pairwise_handle, fieldnames=fieldnames)
        writer.writeheader()
        for pair in pairs:
            accumulator = accumulators[pair]
            fst_value = (
                accumulator.numerator_sum / accumulator.denominator_sum
                if accumulator.denominator_sum > 0.0
                else None
            )
            pairwise_values[pair] = fst_value
            writer.writerow(
                {
                    "population_1": pair[0],
                    "population_2": pair[1],
                    "sites_used": accumulator.sites_used,
                    "numerator_sum": round(accumulator.numerator_sum, 10),
                    "denominator_sum": round(accumulator.denominator_sum, 10),
                    "fst": "" if fst_value is None else round(fst_value, 10),
                }
            )

    population_names = sorted(populations)
    matrix_path = fst_dir / "fst_matrix.csv"
    with matrix_path.open("w", newline="", encoding="utf-8") as matrix_handle:
        writer = csv.writer(matrix_handle)
        writer.writerow(["population", *population_names])
        for population_1 in population_names:
            row: list[str | float] = [population_1]
            for population_2 in population_names:
                if population_1 == population_2:
                    row.append(0.0)
                    continue
                pair = tuple(sorted((population_1, population_2)))
                value = pairwise_values[pair]
                row.append("NA" if value is None else round(value, 10))
            writer.writerow(row)

    result = FSTResult(
        populations=len(populations),
        population_pairs=len(pairs),
        samples=len(sample_names),
        records_seen=records_seen,
        biallelic_snps=biallelic_snps,
        pair_site_records=pair_site_records,
        input_vcf=str(vcf_path),
    )
    (fst_dir / "fst_summary.json").write_text(
        json.dumps(asdict(result), indent=2) + "\n", encoding="utf-8"
    )
    generate_fst_heatmap(fst_dir)
    return result
