from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def _sample_metric_plot(
    sample_rows: list[dict[str, str | int | float]],
    metric: str,
    *,
    ylabel: str,
    title: str,
    stem: str,
    output_dir: Path,
    unit_interval: bool = False,
) -> None:
    if not sample_rows:
        return

    values = [float(row[metric]) for row in sample_rows]
    labels = [str(row["sample"]) for row in sample_rows]
    positions = list(range(1, len(values) + 1))

    width = max(6.5, min(14.0, 0.22 * len(values) + 5.0))
    fig, ax = plt.subplots(figsize=(width, 4.8))
    ax.scatter(positions, values, s=26)
    ax.set_title(title)
    ax.set_xlabel("Sample")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.25)

    if unit_interval:
        ax.set_ylim(0.0, 1.05)

    if len(labels) <= 40:
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
    else:
        ax.set_xlabel("Sample index")

    fig.tight_layout()
    _save_figure(fig, output_dir, stem)


def _histogram_plot(
    counts: list[int],
    *,
    maximum: float,
    xlabel: str,
    title: str,
    stem: str,
    output_dir: Path,
) -> None:
    if not counts or sum(counts) == 0:
        return

    bin_width = maximum / len(counts)
    centers = [(index + 0.5) * bin_width for index in range(len(counts))]

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.bar(centers, counts, width=bin_width * 0.9)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Variant count")
    ax.set_xlim(0.0, maximum)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    _save_figure(fig, output_dir, stem)


def write_qc_plots(
    output_dir: str | Path,
    sample_rows: list[dict[str, str | int | float]],
    *,
    maf_histogram: list[int],
    missingness_histogram: list[int],
    has_depth: bool,
    has_gq: bool,
) -> None:
    """Write PNG and PDF QC diagnostics from streaming summaries."""
    qc_dir = Path(output_dir)

    _sample_metric_plot(
        sample_rows,
        "call_rate",
        ylabel="Genotype call rate",
        title="Per-sample genotype call rate",
        stem="sample_call_rate",
        output_dir=qc_dir,
        unit_interval=True,
    )
    _sample_metric_plot(
        sample_rows,
        "heterozygosity",
        ylabel="Observed heterozygosity",
        title="Per-sample observed heterozygosity",
        stem="sample_heterozygosity",
        output_dir=qc_dir,
        unit_interval=True,
    )

    if has_depth:
        _sample_metric_plot(
            sample_rows,
            "mean_depth",
            ylabel="Mean depth (DP)",
            title="Per-sample mean read depth",
            stem="sample_mean_depth",
            output_dir=qc_dir,
        )

    if has_gq:
        _sample_metric_plot(
            sample_rows,
            "mean_gq",
            ylabel="Mean genotype quality (GQ)",
            title="Per-sample mean genotype quality",
            stem="sample_mean_gq",
            output_dir=qc_dir,
        )

    _histogram_plot(
        maf_histogram,
        maximum=0.5,
        xlabel="Minor allele frequency",
        title="Biallelic minor allele frequency distribution",
        stem="variant_maf_distribution",
        output_dir=qc_dir,
    )
    _histogram_plot(
        missingness_histogram,
        maximum=1.0,
        xlabel="Missing genotype proportion",
        title="Variant missingness distribution",
        stem="variant_missingness_distribution",
        output_dir=qc_dir,
    )
