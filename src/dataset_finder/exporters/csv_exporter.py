"""CSV export support for Dataset Finder records."""

from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from dataset_finder.models import DatasetRecord

FIELD_NAMES = (
    "uid",
    "accession",
    "title",
    "organism",
    "study_type",
    "sample_count",
    "publication_date",
    "url",
)


def export_csv(
    records: list[DatasetRecord],
    output_path: str | Path,
) -> Path:
    """Write normalized dataset records to a UTF-8 CSV file."""
    path = Path(output_path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELD_NAMES)
        writer.writeheader()

        for record in records:
            writer.writerow(asdict(record))

    return path.resolve()
