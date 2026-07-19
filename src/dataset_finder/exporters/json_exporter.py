"""JSON export support for Dataset Finder records."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from dataset_finder.models import DatasetRecord


def export_json(
    records: list[DatasetRecord],
    output_path: str | Path,
) -> Path:
    """Write normalized dataset records to a formatted JSON file."""
    path = Path(output_path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)

    payload = [asdict(record) for record in records]

    with path.open("w", encoding="utf-8") as handle:
        json.dump(
            payload,
            handle,
            ensure_ascii=False,
            indent=2,
        )
        handle.write("\n")

    return path.resolve()
