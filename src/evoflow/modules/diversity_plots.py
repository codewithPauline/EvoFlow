from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _save(fig: plt.Figure, base: Path) -> None:
    fig.tight_layout()
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def generate_diversity_figures(diversity_dir: str | Path) -> list[Path]:
    """Generate population-level diversity diagnostic figures."""
    output_dir = Path(diversity_dir)
    summary_path = output_dir / "population_diversity.csv"
    with summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return []

    populations = [row["population"] for row in rows]
    x = np.arange(len(populations), dtype=float)
    outputs: list[Path] = []

    fig, ax = plt.subplots(figsize=(max(7.0, 0.7 * len(populations)), 4.8))
    width = 0.38
    ax.bar(
        x - width / 2,
        [float(row["mean_observed_heterozygosity"]) for row in rows],
        width,
        label="Observed",
    )
    ax.bar(
        x + width / 2,
        [float(row["mean_expected_heterozygosity"]) for row in rows],
        width,
        label="Expected",
    )
    ax.set_xticks(x, populations, rotation=45, ha="right")
    ax.set_ylabel("Mean heterozygosity across analyzed SNPs")
    ax.set_title("Population heterozygosity")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    heterozygosity_base = output_dir / "population_heterozygosity"
    _save(fig, heterozygosity_base)
    outputs.extend(
        [heterozygosity_base.with_suffix(".png"), heterozygosity_base.with_suffix(".pdf")]
    )

    fig, ax = plt.subplots(figsize=(max(7.0, 0.7 * len(populations)), 4.8))
    ax.bar(x, [float(row["mean_call_rate"]) for row in rows])
    ax.set_xticks(x, populations, rotation=45, ha="right")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Mean genotype call rate")
    ax.set_title("Population call rate across analyzed SNPs")
    ax.spines[["top", "right"]].set_visible(False)
    call_rate_base = output_dir / "population_call_rate"
    _save(fig, call_rate_base)
    outputs.extend([call_rate_base.with_suffix(".png"), call_rate_base.with_suffix(".pdf")])
    return outputs
