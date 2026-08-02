"""Tests for Excel workbook exports."""

from openpyxl import load_workbook

from dataset_finder.batch import BatchSearchIssue, BatchSearchResult
from dataset_finder.exporters.excel_exporter import export_excel
from dataset_finder.models import DatasetRecord


def test_export_excel_creates_expected_sheets(tmp_path) -> None:
    record = DatasetRecord(
        uid="123",
        accession="GSE123",
        title="Drosophila brain RNA-seq",
        organism="Drosophila melanogaster",
        study_type="RNA-seq",
        sample_count=6,
        publication_date="2026-01-01",
        url="https://example.org/GSE123",
        database="GEO",
        gene="bru1",
        gene_set="RBP",
        technique="RNA_seq",
    )

    result = BatchSearchResult(
        records=(record,),
        issues=(
            BatchSearchIssue(
                gene="bad",
                database="geo",
                message="Test error",
            ),
        ),
        genes=("bru1", "bad"),
        gene_set="RBP",
        database="geo",
    )

    output_path = export_excel(
        result,
        tmp_path / "screening.xlsx",
    )

    workbook = load_workbook(
        output_path,
        read_only=True,
    )

    assert "README" in workbook.sheetnames
    assert "Gene_Summary" in workbook.sheetnames
    assert "All_Datasets" in workbook.sheetnames
    assert "RNA_seq" in workbook.sheetnames
    assert "Errors" in workbook.sheetnames

    all_datasets = workbook["All_Datasets"]
    headers = [
        cell.value
        for cell in next(all_datasets.iter_rows())
    ]

    assert "Gene" in headers
    assert "Technique" in headers
    assert "Dataset URL" in headers
