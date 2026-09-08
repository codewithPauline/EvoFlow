import csv
import json
from pathlib import Path

import pytest

from evoflow.modules.ld import pairwise_r2, run_ld_prune

VCF_TEXT = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6
chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1
chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1
chr1\t30\t.\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\t0/1\t1/1
chr1\t40\t.\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0
chr1\t1000\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1
"""


def test_pairwise_r2_handles_identical_and_independent_patterns() -> None:
    x = [0.0, 0.0, 0.0, 2.0, 2.0, 2.0]
    identical = [0.0, 0.0, 0.0, 2.0, 2.0, 2.0]
    independent = [0.0, 1.0, 2.0, 0.0, 1.0, 2.0]

    assert pairwise_r2(x, identical, min_overlap=4) == pytest.approx(1.0)
    assert pairwise_r2(x, independent, min_overlap=4) == pytest.approx(0.0)


def test_ld_prune_removes_linked_snp_within_window(tmp_path: Path) -> None:
    vcf = tmp_path / "ld.vcf"
    vcf.write_text(VCF_TEXT, encoding="utf-8")

    result = run_ld_prune(
        vcf,
        tmp_path / "results",
        r2_threshold=0.8,
        window_bp=100,
        min_overlap=4,
    )

    assert result.samples == 6
    assert result.records_seen == 5
    assert result.informative_snps == 4
    assert result.retained_snps == 3
    assert result.removed_for_ld == 1
    assert result.skipped_noninformative == 1

    ld_dir = tmp_path / "results" / "ld"
    summary = json.loads((ld_dir / "ld_summary.json").read_text(encoding="utf-8"))
    assert summary["r2_threshold"] == pytest.approx(0.8)
    assert summary["window_bp"] == 100

    pruned_text = (ld_dir / "ld_pruned.vcf").read_text(encoding="utf-8")
    assert "##evoflow_ld_r2_threshold=0.8" in pruned_text
    assert "##evoflow_ld_window_bp=100" in pruned_text
    records = [line for line in pruned_text.splitlines() if not line.startswith("#")]
    assert [record.split("\t")[1] for record in records] == ["10", "30", "1000"]

    with (ld_dir / "ld_removed.csv").open(newline="", encoding="utf-8") as handle:
        removed = list(csv.DictReader(handle))
    assert len(removed) == 1
    assert removed[0]["pos"] == "20"
    assert removed[0]["blocked_by_pos"] == "10"
    assert float(removed[0]["r2"]) == pytest.approx(1.0)
