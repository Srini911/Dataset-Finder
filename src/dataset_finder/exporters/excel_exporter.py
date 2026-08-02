"""Formatted Excel workbook export for batch dataset searches."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from dataset_finder.assay_classifier import TECHNIQUE_ORDER
from dataset_finder.batch import BatchSearchResult
from dataset_finder.models import DatasetRecord

RECORD_COLUMNS = [
    "gene",
    "gene_set",
    "official_symbol",
    "flybase_id",
    "synonyms",
    "technique",
    "technique_subtype",
    "database",
    "accession",
    "project_accession",
    "experiment_accessions",
    "sample_accessions",
    "biosample_accessions",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "platform",
    "title",
    "description",
    "organism",
    "study_type",
    "tissue",
    "cell_type",
    "developmental_stage",
    "sex",
    "genotype",
    "strain",
    "treatment",
    "control_status",
    "disease",
    "time_point",
    "perturbation",
    "control",
    "sample_count",
    "publication",
    "publication_date",
    "url",
    "flybase_url",
    "flyatlas_url",
    "flyatlas_brain_male_fpkm",
    "flyatlas_brain_female_fpkm",
    "flyatlas_brain_larval_fpkm",
    "flyatlas_head_male_fpkm",
    "flyatlas_head_female_fpkm",
    "flyatlas_top_male_tissue",
    "flyatlas_top_male_fpkm",
    "flyatlas_top_female_tissue",
    "flyatlas_top_female_fpkm",
    "flyatlas_top_larval_tissue",
    "flyatlas_top_larval_fpkm",
    "evidence_text",
    "match_type",
    "confidence",
    "search_date",
    "uid",
]


DISPLAY_NAMES = {
    "gene": "Gene",
    "gene_set": "Gene Set",
    "official_symbol": "Official Symbol",
    "flybase_id": "FlyBase ID",
    "synonyms": "Synonyms",
    "technique": "Technique",
    "technique_subtype": "Technique Subtype",
    "database": "Database",
    "accession": "Accession",
    "project_accession": "Project Accession",
    "experiment_accessions": "Experiment Accessions",
    "sample_accessions": "Run Accessions",
    "biosample_accessions": "BioSample Accessions",
    "library_strategy": "Library Strategy",
    "library_source": "Library Source",
    "library_selection": "Library Selection",
    "library_layout": "Library Layout",
    "platform": "Platform",
    "title": "Title",
    "description": "Description",
    "organism": "Organism",
    "study_type": "Study Type",
    "tissue": "Tissue",
    "cell_type": "Cell Type",
    "developmental_stage": "Developmental Stage",
    "sex": "Sex",
    "genotype": "Genotype",
    "strain": "Strain",
    "treatment": "Treatment",
    "control_status": "Control Status",
    "disease": "Disease",
    "time_point": "Time Point",
    "perturbation": "Perturbation",
    "control": "Control",
    "sample_count": "Sample Count",
    "publication": "Publication",
    "publication_date": "Publication Date",
    "url": "Dataset URL",
    "flybase_url": "FlyBase URL",
    "flyatlas_url": "FlyAtlas URL",
    "flyatlas_brain_male_fpkm": "FlyAtlas Brain Male FPKM",
    "flyatlas_brain_female_fpkm": "FlyAtlas Brain Female FPKM",
    "flyatlas_brain_larval_fpkm": "FlyAtlas Brain Larval FPKM",
    "flyatlas_head_male_fpkm": "FlyAtlas Head Male FPKM",
    "flyatlas_head_female_fpkm": "FlyAtlas Head Female FPKM",
    "flyatlas_top_male_tissue": "FlyAtlas Top Male Tissue",
    "flyatlas_top_male_fpkm": "FlyAtlas Top Male FPKM",
    "flyatlas_top_female_tissue": "FlyAtlas Top Female Tissue",
    "flyatlas_top_female_fpkm": "FlyAtlas Top Female FPKM",
    "flyatlas_top_larval_tissue": "FlyAtlas Top Larval Tissue",
    "flyatlas_top_larval_fpkm": "FlyAtlas Top Larval FPKM",
    "evidence_text": "Evidence Text",
    "match_type": "Match Type",
    "confidence": "Confidence",
    "search_date": "Search Date",
    "uid": "UID",
}


def _record_rows(records: tuple[DatasetRecord, ...]) -> list[dict[str, object]]:
    """Convert records to workbook-friendly dictionaries."""
    rows: list[dict[str, object]] = []

    for record in records:
        row = record.to_dict()
        rows.append(
            {
                column: row.get(column, "")
                for column in RECORD_COLUMNS
            }
        )

    return rows


def _records_dataframe(
    records: tuple[DatasetRecord, ...],
) -> pd.DataFrame:
    """Create a consistently ordered records dataframe."""
    dataframe = pd.DataFrame(
        _record_rows(records),
        columns=RECORD_COLUMNS,
    )
    return dataframe.rename(columns=DISPLAY_NAMES)


def _safe_sheet_name(value: str) -> str:
    """Return a valid Excel sheet name."""
    replacements = {
        "[": "(",
        "]": ")",
        ":": "-",
        "*": "-",
        "?": "",
        "/": "-",
        "\\": "-",
    }

    for original, replacement in replacements.items():
        value = value.replace(original, replacement)

    return value[:31] or "Sheet"


def _format_dataframe_sheet(
    *,
    worksheet,
    dataframe: pd.DataFrame,
    workbook,
) -> None:
    """Apply consistent formatting to a dataframe worksheet."""
    header_format = workbook.add_format(
        {
            "bold": True,
            "font_color": "#FFFFFF",
            "bg_color": "#1F4E78",
            "border": 1,
            "align": "center",
            "valign": "vcenter",
        }
    )
    text_format = workbook.add_format(
        {
            "valign": "top",
            "text_wrap": True,
        }
    )
    center_format = workbook.add_format(
        {
            "align": "center",
            "valign": "top",
        }
    )
    link_format = workbook.add_format(
        {
            "font_color": "#0563C1",
            "underline": True,
            "valign": "top",
        }
    )

    worksheet.freeze_panes(1, 0)
    worksheet.set_row(0, 30)

    if len(dataframe.columns):
        worksheet.autofilter(
            0,
            0,
            len(dataframe),
            len(dataframe.columns) - 1,
        )

    centered_columns = {
        "Gene",
        "Gene Set",
        "Technique",
        "Database",
        "Accession",
        "Sample Count",
        "Sex",
        "Confidence",
        "Search Date",
    }

    url_columns = {
        "Dataset URL",
        "FlyBase URL",
        "FlyAtlas URL",
    }

    for column_index, column_name in enumerate(dataframe.columns):
        worksheet.write(
            0,
            column_index,
            column_name,
            header_format,
        )

        values = dataframe[column_name].fillna("").astype(str)
        content_width = max(
            [len(column_name)]
            + [len(value) for value in values.head(200)]
        )

        if column_name in {"Title", "Description", "Evidence Text"}:
            width = min(max(content_width + 2, 28), 55)
        elif column_name in url_columns:
            width = 18
        elif column_name in {"Synonyms", "Sample Accessions"}:
            width = min(max(content_width + 2, 20), 35)
        else:
            width = min(max(content_width + 2, 12), 28)

        cell_format = (
            center_format
            if column_name in centered_columns
            else text_format
        )

        worksheet.set_column(
            column_index,
            column_index,
            width,
            cell_format,
        )

        if column_name in url_columns:
            for row_index, url in enumerate(values, start=1):
                if url.startswith(("http://", "https://")):
                    worksheet.write_url(
                        row_index,
                        column_index,
                        url,
                        link_format,
                        string="Open",
                    )


def _write_readme_sheet(
    *,
    writer: pd.ExcelWriter,
    result: BatchSearchResult,
) -> None:
    """Write workbook metadata and usage notes."""
    summary_rows = [
        ("Workbook", "Dataset Finder multi-gene screening report"),
        ("Genes searched", len(result.genes)),
        ("Gene set", result.gene_set or "Custom"),
        ("Databases requested", result.database),
        ("Datasets found", len(result.records)),
        ("Errors recorded", len(result.issues)),
        (
            "Workbook structure",
            "All_Datasets contains every result. Technique sheets contain classified subsets.",
        ),
        (
            "Interpretation warning",
            "Search matches should be reviewed before biological conclusions are made.",
        ),
    ]

    dataframe = pd.DataFrame(
        summary_rows,
        columns=["Field", "Value"],
    )
    dataframe.to_excel(
        writer,
        sheet_name="README",
        index=False,
    )

    workbook = writer.book
    worksheet = writer.sheets["README"]

    header_format = workbook.add_format(
        {
            "bold": True,
            "font_color": "#FFFFFF",
            "bg_color": "#1F4E78",
            "border": 1,
        }
    )
    wrap_format = workbook.add_format(
        {
            "text_wrap": True,
            "valign": "top",
        }
    )

    worksheet.set_column("A:A", 25)
    worksheet.set_column("B:B", 85, wrap_format)
    worksheet.set_row(0, 28)
    worksheet.freeze_panes(1, 0)

    for column_index, column_name in enumerate(dataframe.columns):
        worksheet.write(
            0,
            column_index,
            column_name,
            header_format,
        )


def _write_gene_annotations(
    *,
    writer: pd.ExcelWriter,
    result: BatchSearchResult,
) -> None:
    """Write FlyBase and FlyAtlas annotation rows for every submitted gene."""
    records_by_gene = {
        record.gene.casefold(): record
        for record in result.records
    }

    rows: list[dict[str, object]] = []

    for gene_annotation in result.gene_annotations:
        annotation = gene_annotation.gene
        expression = gene_annotation.flyatlas

        record = records_by_gene.get(
            annotation.submitted_symbol.casefold()
        )

        rows.append(
            {
                "Gene": annotation.submitted_symbol,
                "Gene Set": result.gene_set,
                "Official Symbol": annotation.official_symbol,
                "FlyBase ID": annotation.flybase_id,
                "Current Full Name": annotation.current_fullname,
                "Synonyms": "; ".join(annotation.synonyms),
                "Resolution Type": annotation.match_type,
                "Ambiguous Alias": (
                    "Yes"
                    if annotation.ambiguous
                    else "No"
                ),
                "Accepted Datasets": sum(
                    1
                    for dataset in result.records
                    if dataset.gene.casefold()
                    == annotation.submitted_symbol.casefold()
                ),
                "Dataset Match Type": (
                    record.match_type
                    if record
                    else "No accepted datasets"
                ),
                "Dataset Confidence": (
                    record.confidence
                    if record
                    else ""
                ),
                "FlyBase URL": annotation.flybase_url,
                "FlyAtlas URL": annotation.flyatlas_url,
                "Brain Male FPKM": (
                    expression.brain_male_fpkm
                ),
                "Brain Female FPKM": (
                    expression.brain_female_fpkm
                ),
                "Brain Larval FPKM": (
                    expression.brain_larval_fpkm
                ),
                "Head Male FPKM": (
                    expression.head_male_fpkm
                ),
                "Head Female FPKM": (
                    expression.head_female_fpkm
                ),
                "Top Male Tissue": (
                    expression.top_male_tissue
                ),
                "Top Male FPKM": (
                    expression.top_male_fpkm
                ),
                "Top Female Tissue": (
                    expression.top_female_tissue
                ),
                "Top Female FPKM": (
                    expression.top_female_fpkm
                ),
                "Top Larval Tissue": (
                    expression.top_larval_tissue
                ),
                "Top Larval FPKM": (
                    expression.top_larval_fpkm
                ),
            }
        )

    dataframe = pd.DataFrame(rows)

    dataframe.to_excel(
        writer,
        sheet_name="Gene_Annotations",
        index=False,
    )

    _format_dataframe_sheet(
        worksheet=writer.sheets["Gene_Annotations"],
        dataframe=dataframe,
        workbook=writer.book,
    )



def _write_gene_summary(
    *,
    writer: pd.ExcelWriter,
    result: BatchSearchResult,
) -> None:
    """Write per-gene and per-technique counts."""
    records_df = _records_dataframe(result.records)

    rows: list[dict[str, object]] = []

    for gene in result.genes:
        gene_records = records_df[
            records_df["Gene"].astype(str).str.casefold()
            == gene.casefold()
        ]

        row: dict[str, object] = {
            "Gene": gene,
            "Total Datasets": len(gene_records),
        }

        for database_name in (
            "GEO",
            "ENCODE",
            "SRA",
            "BioProject",
            "BioStudies",
        ):
            row[database_name] = int(
                (
                    gene_records["Database"]
                    == database_name
                ).sum()
            )

        for technique in TECHNIQUE_ORDER:
            row[technique] = int(
                (
                    gene_records["Technique"]
                    == technique
                ).sum()
            )

        rows.append(row)

    dataframe = pd.DataFrame(rows)
    dataframe.to_excel(
        writer,
        sheet_name="Gene_Summary",
        index=False,
    )
    _format_dataframe_sheet(
        worksheet=writer.sheets["Gene_Summary"],
        dataframe=dataframe,
        workbook=writer.book,
    )


def export_excel(
    result: BatchSearchResult,
    output_path: str | Path,
) -> Path:
    """Export a complete multi-sheet Excel screening workbook."""
    path = Path(output_path)

    if path.suffix.casefold() != ".xlsx":
        path = path.with_suffix(".xlsx")

    path.parent.mkdir(parents=True, exist_ok=True)

    all_records_df = _records_dataframe(result.records)

    with pd.ExcelWriter(
        path,
        engine="xlsxwriter",
        engine_kwargs={
            "options": {
                "strings_to_urls": False,
                "constant_memory": False,
            }
        },
    ) as writer:
        _write_readme_sheet(
            writer=writer,
            result=result,
        )

        _write_gene_summary(
            writer=writer,
            result=result,
        )

        _write_gene_annotations(
            writer=writer,
            result=result,
        )

        all_records_df.to_excel(
            writer,
            sheet_name="All_Datasets",
            index=False,
        )
        _format_dataframe_sheet(
            worksheet=writer.sheets["All_Datasets"],
            dataframe=all_records_df,
            workbook=writer.book,
        )

        for technique in TECHNIQUE_ORDER:
            technique_df = all_records_df[
                all_records_df["Technique"] == technique
            ].copy()

            sheet_name = _safe_sheet_name(technique)

            technique_df.to_excel(
                writer,
                sheet_name=sheet_name,
                index=False,
            )
            _format_dataframe_sheet(
                worksheet=writer.sheets[sheet_name],
                dataframe=technique_df,
                workbook=writer.book,
            )

        status_df = pd.DataFrame(
            [
                {
                    "Gene": status.gene,
                    "Database": status.database,
                    "Status": (
                        "Success"
                        if status.success
                        else "Failed"
                    ),
                    "Results": status.result_count,
                    "Error": status.error,
                }
                for status in result.database_statuses
            ],
            columns=[
                "Gene",
                "Database",
                "Status",
                "Results",
                "Error",
            ],
        )
        status_df.to_excel(
            writer,
            sheet_name="Database_Status",
            index=False,
        )
        _format_dataframe_sheet(
            worksheet=writer.sheets["Database_Status"],
            dataframe=status_df,
            workbook=writer.book,
        )

        issue_df = pd.DataFrame(
            [
                {
                    "Gene": issue.gene,
                    "Database": issue.database,
                    "Error": issue.message,
                }
                for issue in result.issues
            ],
            columns=["Gene", "Database", "Error"],
        )
        issue_df.to_excel(
            writer,
            sheet_name="Errors",
            index=False,
        )
        _format_dataframe_sheet(
            worksheet=writer.sheets["Errors"],
            dataframe=issue_df,
            workbook=writer.book,
        )

    return path
