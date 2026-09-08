import csv
import json
from pathlib import Path

import pytest

from evoflow.modules.diversity import run_diversity
from evoflow.modules.fst import hudson_fst_components, run_fst

VCF_TEXT = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1
chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1
"""

METADATA = """sample,population
S1,A
S2,A
S3,B
S4,B
"""


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    vcf = tmp_path / "population.vcf"
    metadata = tmp_path / "samples.csv"
    vcf.write_text(VCF_TEXT, encoding="utf-8")
    metadata.write_text(METADATA, encoding="utf-8")
    return vcf, metadata


def test_diversity_summaries_match_expected_values(tmp_path: Path) -> None:
    vcf, metadata = _write_inputs(tmp_path)
    result = run_diversity(vcf, metadata, tmp_path / "results")

    assert result.populations == 2
    assert result.biallelic_snps == 2
    assert result.site_population_records == 4

    diversity_dir = tmp_path / "results" / "diversity"
    with (diversity_dir / "population_diversity.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = {row["population"]: row for row in csv.DictReader(handle)}

    for population in ["A", "B"]:
        assert int(rows[population]["sites_observed"]) == 2
        assert int(rows[population]["polymorphic_sites"]) == 1
        assert float(rows[population]["mean_call_rate"]) == pytest.approx(1.0)
        assert float(rows[population]["mean_observed_heterozygosity"]) == pytest.approx(0.5)
        assert float(rows[population]["mean_expected_heterozygosity"]) == pytest.approx(0.25)
        assert float(rows[population]["mean_maf"]) == pytest.approx(0.25)

    summary = json.loads(
        (diversity_dir / "diversity_summary.json").read_text(encoding="utf-8")
    )
    assert "analyzed biallelic SNP records" in summary["metric_scope"]

    for stem in ["population_heterozygosity", "population_call_rate"]:
        for suffix in [".png", ".pdf"]:
            figure = diversity_dir / f"{stem}{suffix}"
            assert figure.exists()
            assert figure.stat().st_size > 0


def test_hudson_fst_components_and_aggregate(tmp_path: Path) -> None:
    fixed = hudson_fst_components((4, 0), (0, 4))
    assert fixed is not None
    assert fixed[0] == pytest.approx(1.0)
    assert fixed[1] == pytest.approx(1.0)

    shared = hudson_fst_components((2, 2), (2, 2))
    assert shared is not None
    assert shared[0] == pytest.approx(-1.0 / 6.0)
    assert shared[1] == pytest.approx(0.5)

    vcf, metadata = _write_inputs(tmp_path)
    result = run_fst(vcf, metadata, tmp_path / "results")
    assert result.populations == 2
    assert result.population_pairs == 1
    assert result.biallelic_snps == 2

    fst_dir = tmp_path / "results" / "fst"
    with (fst_dir / "pairwise_fst.csv").open(newline="", encoding="utf-8") as handle:
        pairwise = list(csv.DictReader(handle))
    assert len(pairwise) == 1
    assert pairwise[0]["population_1"] == "A"
    assert pairwise[0]["population_2"] == "B"
    assert int(pairwise[0]["sites_used"]) == 2
    assert float(pairwise[0]["fst"]) == pytest.approx(5.0 / 9.0)

    with (fst_dir / "site_fst.csv").open(newline="", encoding="utf-8") as handle:
        sites = list(csv.DictReader(handle))
    assert float(sites[0]["fst"]) == pytest.approx(1.0)
    assert float(sites[1]["fst"]) == pytest.approx(-1.0 / 3.0)

    with (fst_dir / "fst_matrix.csv").open(newline="", encoding="utf-8") as handle:
        matrix = list(csv.reader(handle))
    assert matrix[0] == ["population", "A", "B"]
    assert float(matrix[1][2]) == pytest.approx(5.0 / 9.0)
    assert float(matrix[2][1]) == pytest.approx(5.0 / 9.0)

    for suffix in [".png", ".pdf"]:
        figure = fst_dir / f"fst_heatmap{suffix}"
        assert figure.exists()
        assert figure.stat().st_size > 0
