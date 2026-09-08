from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
import yaml

from evoflow import __version__
from evoflow.config import EvoFlowConfig
from evoflow.io.validation import validate_config
from evoflow.modules.registry import MODULES, validate_modules

app = typer.Typer(help="EvoFlow — from variants to evolutionary insight.")


@app.command()
def version() -> None:
    """Show the installed EvoFlow version."""
    typer.echo(f"EvoFlow {__version__}")


@app.command()
def init(
    project: Annotated[str, typer.Option(help="Project name.")],
    vcf: Annotated[Path, typer.Option(help="Input VCF/BCF path.")],
    metadata: Annotated[Path, typer.Option(help="Sample metadata CSV path.")],
    output: Annotated[Path, typer.Option(help="Configuration file to create.")] = Path(
        "evoflow.yaml"
    ),
) -> None:
    """Create an EvoFlow project configuration."""
    config = EvoFlowConfig(project=project, vcf=vcf, metadata=metadata)
    output.write_text(yaml.safe_dump(config.to_dict(), sort_keys=False), encoding="utf-8")
    typer.echo(f"Created {output}")


@app.command()
def validate(config: Annotated[Path, typer.Argument(exists=True)]) -> None:
    """Validate EvoFlow configuration and input files."""
    cfg = EvoFlowConfig.from_yaml(config)
    validate_modules(cfg.modules)
    for message in validate_config(cfg):
        typer.echo(f"✓ {message}")
    typer.echo("✓ EvoFlow inputs are valid")


@app.command()
def plan(config: Annotated[Path, typer.Argument(exists=True)]) -> None:
    """Print the ordered analysis plan without running it."""
    cfg = EvoFlowConfig.from_yaml(config)
    validate_modules(cfg.modules)
    typer.echo(f"Project: {cfg.project}")
    for index, module in enumerate(cfg.modules, start=1):
        typer.echo(f"{index}. {module}: {MODULES[module]}")


if __name__ == "__main__":
    app()
