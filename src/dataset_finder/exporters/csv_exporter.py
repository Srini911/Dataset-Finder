"""CSV export support for normalized dataset records."""

from __future__ import annotations

import csv
from pathlib import Path

from dataset_finder.models import DatasetRecord


def export_csv(
    records: list[DatasetRecord],
    output_path: str | Path,
) -> Path:
    """Write normalized dataset records to a CSV file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(DatasetRecord.__dataclass_fields__)

    try:
        with path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as output_file:
            writer = csv.DictWriter(
                output_file,
                fieldnames=fieldnames,
                extrasaction="ignore",
            )
            writer.writeheader()

            for record in records:
                writer.writerow(record.to_dict())
    except OSError:
        raise

    return path
