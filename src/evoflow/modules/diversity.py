from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from evoflow.io.metadata import population_indices
from evoflow.io.vcf import extract_biallelic_snp_dosages, open_vcf_text, read_vcf_samples
from evoflow.modules.diversity_plots import generate_diversity_figures


@dataclass(slots=True)
class DiversityResult:
    populations: int
    samples: int
    records_seen: int
    biallelic_snps: int
    site_population_records: int
    input_vcf: str
    population_column: str
    metric_scope: str = "mean diversity across analyzed biallelic SNP records"


@dataclass(slots=True)
class _PopulationAccumulator:
    samples: int
    sites: int = 0
    polymorphic_sites: int = 0
    call_rate_sum: float = 0.0
    observed_het_sum: float = 0.0
    expected_het_sum: float = 0.0
    maf_sum: float = 0.0


def run_diversity(
    vcf: str | Path,
    metadata: str | Path,
    output_dir: str | Path,
    *,
    population_column: str = "population",
) -> DiversityResult:
    """Summarize population diversity across analyzed diploid biallelic SNPs."""
    vcf_path = Path(vcf)
    sample_names = read_vcf_samples(vcf_path)
    populations = population_indices(
        metadata,
        sample_names,
        population_column=population_column,
    )

    diversity_dir = Path(output_dir) / "diversity"
    diversity_dir.mkdir(parents=True, exist_ok=True)
    site_path = diversity_dir / "site_diversity.csv"
    population_path = diversity_dir / "population_diversity.csv"

    accumulators = {
        population: _PopulationAccumulator(samples=len(indices))
        for population, indices in populations.items()
    }
    records_seen = 0
    biallelic_snps = 0
    site_population_records = 0

    with site_path.open("w", newline="", encoding="utf-8") as site_handle:
        fieldnames = [
            "chrom",
            "pos",
            "population",
            "called_samples",
            "population_samples",
            "call_rate",
            "alt_allele_frequency",
            "maf",
            "observed_heterozygosity",
            "expected_heterozygosity",
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

                for population, indices in populations.items():
                    called = [dosages[index] for index in indices if dosages[index] is not None]
                    if not called:
                        continue

                    called_count = len(called)
                    population_count = len(indices)
                    call_rate = called_count / population_count
                    alt_frequency = float(sum(called) / (2.0 * called_count))
                    maf = min(alt_frequency, 1.0 - alt_frequency)
                    observed_heterozygosity = (
                        sum(1 for dosage in called if dosage == 1.0) / called_count
                    )
                    expected_heterozygosity = 2.0 * alt_frequency * (1.0 - alt_frequency)

                    writer.writerow(
                        {
                            "chrom": fields[0],
                            "pos": fields[1],
                            "population": population,
                            "called_samples": called_count,
                            "population_samples": population_count,
                            "call_rate": round(call_rate, 8),
                            "alt_allele_frequency": round(alt_frequency, 8),
                            "maf": round(maf, 8),
                            "observed_heterozygosity": round(observed_heterozygosity, 8),
                            "expected_heterozygosity": round(expected_heterozygosity, 8),
                        }
                    )
                    site_population_records += 1

                    accumulator = accumulators[population]
                    accumulator.sites += 1
                    accumulator.polymorphic_sites += int(0.0 < alt_frequency < 1.0)
                    accumulator.call_rate_sum += call_rate
                    accumulator.observed_het_sum += observed_heterozygosity
                    accumulator.expected_het_sum += expected_heterozygosity
                    accumulator.maf_sum += maf

    with population_path.open("w", newline="", encoding="utf-8") as population_handle:
        fieldnames = [
            "population",
            "samples",
            "sites_observed",
            "polymorphic_sites",
            "mean_call_rate",
            "mean_observed_heterozygosity",
            "mean_expected_heterozygosity",
            "mean_maf",
        ]
        writer = csv.DictWriter(population_handle, fieldnames=fieldnames)
        writer.writeheader()
        for population in sorted(accumulators):
            accumulator = accumulators[population]
            sites = accumulator.sites
            writer.writerow(
                {
                    "population": population,
                    "samples": accumulator.samples,
                    "sites_observed": sites,
                    "polymorphic_sites": accumulator.polymorphic_sites,
                    "mean_call_rate": round(accumulator.call_rate_sum / sites, 8)
                    if sites
                    else 0.0,
                    "mean_observed_heterozygosity": round(
                        accumulator.observed_het_sum / sites, 8
                    )
                    if sites
                    else 0.0,
                    "mean_expected_heterozygosity": round(
                        accumulator.expected_het_sum / sites, 8
                    )
                    if sites
                    else 0.0,
                    "mean_maf": round(accumulator.maf_sum / sites, 8) if sites else 0.0,
                }
            )

    result = DiversityResult(
        populations=len(populations),
        samples=len(sample_names),
        records_seen=records_seen,
        biallelic_snps=biallelic_snps,
        site_population_records=site_population_records,
        input_vcf=str(vcf_path),
        population_column=population_column,
    )
    (diversity_dir / "diversity_summary.json").write_text(
        json.dumps(asdict(result), indent=2) + "\n", encoding="utf-8"
    )
    generate_diversity_figures(diversity_dir)
    return result
