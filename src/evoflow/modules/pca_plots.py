from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


def _save_figure(fig: plt.Figure, output_base: Path) -> None:
    fig.tight_layout()
    fig.savefig(output_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def generate_pca_figures(pca_dir: str | Path) -> list[Path]:
    """Generate PCA scatter and scree plots from EvoFlow PCA tables."""
    output_dir = Path(pca_dir)

    with (output_dir / "pca_variance.csv").open(newline="", encoding="utf-8") as handle:
        variance_rows = list(csv.DictReader(handle))
    with (output_dir / "pca_scores.csv").open(newline="", encoding="utf-8") as handle:
        score_rows = list(csv.DictReader(handle))

    if not variance_rows or not score_rows:
        return []

    ratios = [float(row["explained_variance_ratio"]) for row in variance_rows]
    labels = [row["component"] for row in variance_rows]
    outputs: list[Path] = []

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.bar(labels, [ratio * 100.0 for ratio in ratios])
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Explained variance (%)")
    ax.set_title("PCA explained variance")
    ax.spines[["top", "right"]].set_visible(False)
    scree_base = output_dir / "pca_scree"
    _save_figure(fig, scree_base)
    outputs.extend([scree_base.with_suffix(".png"), scree_base.with_suffix(".pdf")])

    if len(variance_rows) < 2 or "PC2" not in score_rows[0]:
        return outputs

    pc1_pct = ratios[0] * 100.0
    pc2_pct = ratios[1] * 100.0
    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    population_available = "population" in score_rows[0] and any(
        row.get("population", "") for row in score_rows
    )
    if population_available:
        populations = sorted({row.get("population", "") or "Unassigned" for row in score_rows})
        for population in populations:
            group = [
                row
                for row in score_rows
                if (row.get("population", "") or "Unassigned") == population
            ]
            ax.scatter(
                [float(row["PC1"]) for row in group],
                [float(row["PC2"]) for row in group],
                label=population,
                alpha=0.85,
            )
        ax.legend(title="Population", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    else:
        ax.scatter(
            [float(row["PC1"]) for row in score_rows],
            [float(row["PC2"]) for row in score_rows],
            alpha=0.85,
        )

    ax.axhline(0, linewidth=0.7, alpha=0.35)
    ax.axvline(0, linewidth=0.7, alpha=0.35)
    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}%)")
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}%)")
    ax.set_title("Population-genomic PCA")
    ax.spines[["top", "right"]].set_visible(False)
    scatter_base = output_dir / "pca_pc1_pc2"
    _save_figure(fig, scatter_base)
    outputs.extend([scatter_base.with_suffix(".png"), scatter_base.with_suffix(".pdf")])
    return outputs
