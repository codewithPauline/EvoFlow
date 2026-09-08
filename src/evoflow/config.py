from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(slots=True)
class EvoFlowConfig:
    project: str
    vcf: Path
    metadata: Path
    output_dir: Path = Path("evoflow-results")
    modules: list[str] = field(default_factory=lambda: ["qc"])

    @classmethod
    def from_yaml(cls, path: str | Path) -> EvoFlowConfig:
        config_path = Path(path)
        data: dict[str, Any] = yaml.safe_load(config_path.read_text()) or {}
        required = {"project", "vcf", "metadata"}
        missing = sorted(required - data.keys())
        if missing:
            raise ValueError(f"Missing required configuration fields: {', '.join(missing)}")
        return cls(
            project=str(data["project"]),
            vcf=Path(data["vcf"]),
            metadata=Path(data["metadata"]),
            output_dir=Path(data.get("output_dir", "evoflow-results")),
            modules=list(data.get("modules", ["qc"])),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "project": self.project,
            "vcf": str(self.vcf),
            "metadata": str(self.metadata),
            "output_dir": str(self.output_dir),
            "modules": self.modules,
        }
