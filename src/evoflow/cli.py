from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
import yaml

from evoflow import __version__
from evoflow.config import EvoFlowConfig
from evoflow.io.validation import validate_config
from evoflow.modules.ld import run_ld_prune
from evoflow.modules.pca import run_pca
from evoflow.modules.qc import run_qc
from evoflow.modules.registry import MODULES, validate_modules

app = typer.Typer(help="EvoFlow — from variants to evolutionary insight.")


@app.command()
def version() -> None:
    """Show the installed EvoFlow version."""
    typer.echo(f"EvoFlow {__version__}")


@app.command()
def init(
    project: Annotated[str, typer.Option(help="Project name.")],
    vcf: Annotated[Path, typer.Option(help="Input VCF or VCF.gz path.")],
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


@app.command()
def qc(
    config: Annotated[Path, typer.Argument(exists=True)],
    min_maf: Annotated[
        float,
        typer.Option(help="Minimum minor allele frequency for the QC pass flag."),
    ] = 0.0,
    max_missing: Annotated[
        float,
        typer.Option(help="Maximum allowed missing-genotype proportion per variant."),
    ] = 1.0,
) -> None:
    """Run native variant and sample QC on a VCF/VCF.gz file."""
    cfg = EvoFlowConfig.from_yaml(config)
    for message in validate_config(cfg):
        typer.echo(f"✓ {message}")

    summary = run_qc(
        cfg.vcf,
        cfg.output_dir,
        min_maf=min_maf,
        max_missing=max_missing,
    )
    typer.echo(f"✓ QC complete: {summary.samples} samples, {summary.variants} variants")
    typer.echo(
        f"✓ {summary.retained_variants} variants pass MAF/missingness thresholds "
        f"(call rate {summary.overall_call_rate:.3f})"
    )
    typer.echo(f"✓ Results: {cfg.output_dir / 'qc'}")


@app.command(name="ld-prune")
def ld_prune(
    config: Annotated[Path, typer.Argument(exists=True)],
    r2_threshold: Annotated[
        float,
        typer.Option("--r2", help="Maximum allowed pairwise genotype-dosage r²."),
    ] = 0.2,
    window_bp: Annotated[
        int,
        typer.Option(help="Physical sliding-window size in base pairs."),
    ] = 50_000,
    min_overlap: Annotated[
        int,
        typer.Option(help="Minimum samples called at both SNPs for an r² calculation."),
    ] = 3,
    use_raw: Annotated[
        bool,
        typer.Option("--use-raw", help="Use the configured VCF instead of QC filtered.vcf."),
    ] = False,
) -> None:
    """Greedily prune linked biallelic SNPs in a physical sliding window."""
    cfg = EvoFlowConfig.from_yaml(config)
    for message in validate_config(cfg):
        typer.echo(f"✓ {message}")

    filtered_vcf = cfg.output_dir / "qc" / "filtered.vcf"
    input_vcf = cfg.vcf if use_raw or not filtered_vcf.exists() else filtered_vcf
    if input_vcf == filtered_vcf:
        typer.echo(f"✓ LD input: QC-filtered variants ({filtered_vcf})")
    else:
        typer.echo(f"✓ LD input: configured variant file ({cfg.vcf})")

    result = run_ld_prune(
        input_vcf,
        cfg.output_dir,
        r2_threshold=r2_threshold,
        window_bp=window_bp,
        min_overlap=min_overlap,
    )
    typer.echo(
        f"✓ LD pruning complete: {result.retained_snps}/{result.informative_snps} "
        f"informative SNPs retained"
    )
    typer.echo(f"✓ Results: {cfg.output_dir / 'ld'}")


@app.command()
def pca(
    config: Annotated[Path, typer.Argument(exists=True)],
    n_components: Annotated[
        int,
        typer.Option(help="Maximum number of principal components to report."),
    ] = 10,
    use_raw: Annotated[
        bool,
        typer.Option(
            "--use-raw",
            help="Use the configured VCF instead of LD-pruned or QC-filtered variants.",
        ),
    ] = False,
) -> None:
    """Run allele-frequency-standardized population-genomic PCA."""
    cfg = EvoFlowConfig.from_yaml(config)
    for message in validate_config(cfg):
        typer.echo(f"✓ {message}")

    ld_vcf = cfg.output_dir / "ld" / "ld_pruned.vcf"
    filtered_vcf = cfg.output_dir / "qc" / "filtered.vcf"
    if use_raw:
        input_vcf = cfg.vcf
    elif ld_vcf.exists():
        input_vcf = ld_vcf
    elif filtered_vcf.exists():
        input_vcf = filtered_vcf
    else:
        input_vcf = cfg.vcf

    if input_vcf == ld_vcf:
        typer.echo(f"✓ PCA input: LD-pruned variants ({ld_vcf})")
    elif input_vcf == filtered_vcf:
        typer.echo(f"✓ PCA input: QC-filtered variants ({filtered_vcf})")
    else:
        typer.echo(f"✓ PCA input: configured variant file ({cfg.vcf})")

    result = run_pca(
        input_vcf,
        cfg.metadata,
        cfg.output_dir,
        n_components=n_components,
    )
    typer.echo(
        f"✓ PCA complete: {result.samples} samples, {result.variants_used} informative SNPs, "
        f"{result.components} components"
    )
    typer.echo(f"✓ Results: {cfg.output_dir / 'pca'}")


if __name__ == "__main__":
    app()
