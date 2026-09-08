from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def generate_fst_heatmap(fst_dir: str | Path) -> list[Path]:
    """Render the EvoFlow pairwise FST matrix as PNG and PDF."""
    output_dir = Path(fst_dir)
    matrix_path = output_dir / "fst_matrix.csv"
    with matrix_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 2 or len(rows[0]) < 2:
        return []

    populations = rows[0][1:]
    matrix = np.array(
        [
            [float(value) if value not in {"", "NA"} else np.nan for value in row[1:]]
            for row in rows[1:]
        ],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(max(6.0, 0.75 * len(populations)), 5.5))
    image = ax.imshow(matrix, aspect="equal")
    ax.set_xticks(np.arange(len(populations)), populations, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(populations)), populations)
    ax.set_title("Pairwise Hudson FST")
    fig.colorbar(image, ax=ax, label="FST")

    for row_index in range(matrix.shape[0]):
        for col_index in range(matrix.shape[1]):
            value = matrix[row_index, col_index]
            if np.isfinite(value):
                ax.text(col_index, row_index, f"{value:.3f}", ha="center", va="center")

    fig.tight_layout()
    base = output_dir / "fst_heatmap"
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return [base.with_suffix(".png"), base.with_suffix(".pdf")]
