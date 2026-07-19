"""Tests for CSV and JSON dataset exporters."""

from __future__ import annotations

import csv
import json

from dataset_finder.exporters import export_csv, export_json
from dataset_finder.models import DatasetRecord


def make_record() -> DatasetRecord:
    """Create a representative normalized dataset record."""
    return DatasetRecord(
        uid="12345",
        accession="GSE12345",
        title="Example RNA-seq dataset",
        organism="Homo sapiens",
        study_type="Expression profiling by high throughput sequencing",
        sample_count=8,
        publication_date="2025-01-15",
        url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12345",
    )


def test_export_csv_writes_normalized_records(tmp_path) -> None:
    """CSV export should include a header and normalized record values."""
    output_path = tmp_path / "results.csv"

    resolved_path = export_csv([make_record()], output_path)

    assert resolved_path == output_path.resolve()
    assert output_path.exists()

    with output_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    assert rows[0]["accession"] == "GSE12345"
    assert rows[0]["organism"] == "Homo sapiens"
    assert rows[0]["sample_count"] == "8"


def test_export_json_writes_formatted_records(tmp_path) -> None:
    """JSON export should contain a list of normalized record objects."""
    output_path = tmp_path / "results.json"

    resolved_path = export_json([make_record()], output_path)

    assert resolved_path == output_path.resolve()
    assert output_path.exists()

    payload = json.loads(output_path.read_text(encoding="utf-8"))

    assert len(payload) == 1
    assert payload[0]["accession"] == "GSE12345"
    assert payload[0]["organism"] == "Homo sapiens"
    assert payload[0]["sample_count"] == 8


def test_exporters_create_parent_directories(tmp_path) -> None:
    """Exporters should create missing parent directories."""
    csv_path = tmp_path / "nested" / "csv" / "results.csv"
    json_path = tmp_path / "nested" / "json" / "results.json"

    export_csv([], csv_path)
    export_json([], json_path)

    assert csv_path.exists()
    assert json_path.exists()
