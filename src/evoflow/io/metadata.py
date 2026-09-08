from __future__ import annotations

import csv
from pathlib import Path


def read_metadata(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a metadata CSV and return field names and rows."""
    metadata_path = Path(path)
    with metadata_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    return fieldnames, rows


def population_indices(
    path: str | Path,
    sample_names: list[str],
    *,
    population_column: str = "population",
) -> dict[str, list[int]]:
    """Map population labels to VCF sample indices using metadata sample IDs."""
    fieldnames, rows = read_metadata(path)
    if "sample" not in fieldnames:
        raise ValueError("Metadata must contain a 'sample' column.")
    if population_column not in fieldnames:
        raise ValueError(
            f"Metadata must contain a '{population_column}' column for population analyses."
        )

    row_by_sample = {row["sample"]: row for row in rows}
    populations: dict[str, list[int]] = {}
    for index, sample in enumerate(sample_names):
        if sample not in row_by_sample:
            raise ValueError(f"Sample '{sample}' is missing from metadata.")
        population = row_by_sample[sample].get(population_column, "").strip()
        if not population:
            raise ValueError(
                f"Sample '{sample}' has an empty '{population_column}' value in metadata."
            )
        populations.setdefault(population, []).append(index)

    if len(populations) < 1:
        raise ValueError("No populations were found in metadata.")
    return populations
