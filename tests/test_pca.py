import csv
import json
from pathlib import Path

import pytest

from evoflow.modules.pca import run_pca

VCF_TEXT = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1
chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1
chr1\t30\t.\tG\tA\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\t1/1
chr1\t40\t.\tT\tC\t.\tPASS\t.\tGT\t1/1\t0/0\t1/1\t0/0
chr1\t50\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0
chr1\t60\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1
"""

METADATA = """sample,population,site
S1,A,west
S2,A,west
S3,B,east
S4,B,east
"""


def test_streaming_pca_outputs_expected_structure(tmp_path: Path) -> None:
    vcf = tmp_path / "pca.vcf"
    metadata = tmp_path / "samples.csv"
    vcf.write_text(VCF_TEXT, encoding="utf-8")
    metadata.write_text(METADATA, encoding="utf-8")

    result = run_pca(vcf, metadata, tmp_path / "results", n_components=3)

    assert result.samples == 4
    assert result.variants_used == 4
    assert result.components == 2
    assert sum(result.explained_variance_ratio) == pytest.approx(1.0)
    assert result.explained_variance_ratio[0] == pytest.approx(0.5)
    assert result.explained_variance_ratio[1] == pytest.approx(0.5)

    pca_dir = tmp_path / "results" / "pca"
    summary = json.loads((pca_dir / "pca_summary.json").read_text(encoding="utf-8"))
    assert summary["variants_used"] == 4
    assert "streaming Gram-matrix PCA" in summary["method"]

    with (pca_dir / "pca_scores.csv").open(newline="", encoding="utf-8") as handle:
        scores = list(csv.DictReader(handle))
    assert [row["sample"] for row in scores] == ["S1", "S2", "S3", "S4"]
    assert [row["population"] for row in scores] == ["A", "A", "B", "B"]
    assert [row["site"] for row in scores] == ["west", "west", "east", "east"]

    with (pca_dir / "pca_loadings.csv").open(newline="", encoding="utf-8") as handle:
        loadings = list(csv.DictReader(handle))
    assert len(loadings) == 4
    assert [row["pos"] for row in loadings] == ["10", "20", "30", "40"]

    with (pca_dir / "pca_variance.csv").open(newline="", encoding="utf-8") as handle:
        variance = list(csv.DictReader(handle))
    assert len(variance) == 2
    assert float(variance[-1]["cumulative_explained_variance"]) == pytest.approx(1.0)

    for stem in ["pca_scree", "pca_pc1_pc2"]:
        for suffix in [".png", ".pdf"]:
            figure = pca_dir / f"{stem}{suffix}"
            assert figure.exists()
            assert figure.stat().st_size > 0
