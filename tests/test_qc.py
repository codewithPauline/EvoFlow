import csv
import gzip
import json
from pathlib import Path

from evoflow.io.vcf import read_vcf_samples
from evoflow.modules.qc import run_qc

VCF_TEXT = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:10\t0/0:12
chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT:DP\t./.:.\t0/1:8
chr1\t30\t.\tG\tA,C\t.\tPASS\t.\tGT:DP\t1/2:20\t0/2:18
chr1\t40\t.\tA\tAT\t.\tPASS\t.\tGT:DP\t0/1:7\t0/0:6
"""


def test_run_qc_writes_expected_metrics(tmp_path: Path) -> None:
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(VCF_TEXT, encoding="utf-8")

    summary = run_qc(vcf, tmp_path / "results", min_maf=0.3, max_missing=0.4)

    assert summary.samples == 2
    assert summary.variants == 4
    assert summary.retained_variants == 1
    assert summary.snps == 3
    assert summary.biallelic_snps == 2
    assert summary.multiallelic_variants == 1
    assert summary.transitions == 2
    assert summary.transversions == 0
    assert summary.ti_tv_ratio is None
    assert summary.called_genotypes == 7
    assert summary.missing_genotypes == 1
    assert summary.overall_call_rate == 0.875
    assert summary.mean_depth == 11.571429

    qc_dir = tmp_path / "results" / "qc"
    stored_summary = json.loads((qc_dir / "qc_summary.json").read_text(encoding="utf-8"))
    assert stored_summary["retained_variants"] == 1

    with (qc_dir / "sample_qc.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["sample"] == "S1"
    assert float(rows[0]["call_rate"]) == 0.75
    assert float(rows[0]["heterozygosity"]) == 1.0
    assert rows[1]["sample"] == "S2"
    assert float(rows[1]["call_rate"]) == 1.0
    assert float(rows[1]["heterozygosity"]) == 0.5

    with (qc_dir / "variant_qc.csv").open(newline="", encoding="utf-8") as handle:
        variant_rows = list(csv.DictReader(handle))
    assert variant_rows[0]["maf"] == "0.25"
    assert variant_rows[0]["passes_qc"] == "False"
    assert variant_rows[2]["maf"] == ""
    assert variant_rows[2]["passes_qc"] == "True"


def test_gzipped_vcf_header_is_supported(tmp_path: Path) -> None:
    vcf = tmp_path / "tiny.vcf.gz"
    with gzip.open(vcf, "wt", encoding="utf-8") as handle:
        handle.write(VCF_TEXT)

    assert read_vcf_samples(vcf) == ["S1", "S2"]
    summary = run_qc(vcf, tmp_path / "gz-results")
    assert summary.samples == 2
    assert summary.variants == 4
