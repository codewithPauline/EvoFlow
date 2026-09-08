from pathlib import Path

import yaml

from evoflow.config import EvoFlowConfig


def test_config_round_trip(tmp_path: Path) -> None:
    config_file = tmp_path / "evoflow.yaml"
    config_file.write_text(
        yaml.safe_dump(
            {
                "project": "demo",
                "vcf": "data/demo.vcf",
                "metadata": "data/samples.csv",
                "modules": ["qc", "pca"],
            }
        )
    )
    config = EvoFlowConfig.from_yaml(config_file)
    assert config.project == "demo"
    assert config.modules == ["qc", "pca"]
