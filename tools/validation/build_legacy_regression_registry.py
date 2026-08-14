"""Build a regression registry from validated legacy dataset audits."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

VERIFIED = "VERIFIED_LEGACY_POSITIVE"


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""

    return str(value).strip()


def first_present(
    row: pd.Series,
    *columns: str,
) -> str:
    for column in columns:
        if column not in row.index:
            continue

        value = clean(row[column])

        if value:
            return value

    return ""


def normalize_technique(value: object) -> str:
    text = clean(value).casefold()

    mapping = {
        "rna_seq": "RNA_seq",
        "rna-seq": "RNA_seq",
        "rnaseq": "RNA_seq",
        "chip_seq": "ChIP_seq",
        "chip-seq": "ChIP_seq",
        "chipseq": "ChIP_seq",
        "cut_run": "CUT_RUN",
        "cut&run": "CUT_RUN",
        "cutnrun": "CUT_RUN",
        "cut_tag": "CUT_TAG",
        "cut&tag": "CUT_TAG",
        "cutntag": "CUT_TAG",
        "eclip": "eCLIP",
        "iclip": "iCLIP",
        "par_clip": "PAR_CLIP",
        "par-clip": "PAR_CLIP",
        "hits_clip": "HITS_CLIP",
        "hits-clip": "HITS_CLIP",
        "clip": "CLIP",
        "clip_seq": "CLIP",
        "clip-seq": "CLIP",
    }

    return mapping.get(
        text,
        clean(value),
    )


def build_rows(
    *,
    audit_path: Path,
    gene_set: str,
) -> list[dict[str, object]]:
    dataframe = pd.read_excel(
        audit_path,
        sheet_name="All_Legacy_Calls",
    )

    if "Legacy Verdict" not in dataframe.columns:
        raise ValueError(
            f"{audit_path}: missing 'Legacy Verdict' column"
        )

    dataframe = dataframe[
        dataframe["Legacy Verdict"]
        .fillna("")
        .astype(str)
        .str.strip()
        == VERIFIED
    ].copy()

    rows: list[dict[str, object]] = []

    for _, row in dataframe.iterrows():
        gene = first_present(
            row,
            "Gene",
        )

        technique = normalize_technique(
            first_present(
                row,
                "Resolved Technique",
                "Legacy Technique",
            )
        )

        study = first_present(
            row,
            "Resolved Study",
            "study_accession",
            "Study Accession",
        )

        experiment = first_present(
            row,
            "Resolved Experiment",
            "experiment_accession",
            "Experiment Accession",
        )

        old_accession = first_present(
            row,
            "Old Accession",
            "Old Stored Accession",
        )

        title = first_present(
            row,
            "Study Title",
            "study_title",
        )

        strategy = first_present(
            row,
            "Library Strategy",
            "library_strategy",
        )

        if not gene or not technique:
            continue

        rows.append(
            {
                "Gene Set": gene_set,
                "Gene": gene,
                "Technique": technique,
                "Expected Study Accession": study,
                "Expected Experiment Accession": experiment,
                "Legacy Stored Accession": old_accession,
                "Library Strategy": strategy,
                "Study Title": title,
                "Legacy Verdict": VERIFIED,
                "Source Audit": audit_path.name,
            }
        )

    return rows


def deduplicate(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    keys = [
        "Gene Set",
        "Gene",
        "Technique",
        "Expected Study Accession",
        "Expected Experiment Accession",
        "Legacy Stored Accession",
    ]

    return (
        dataframe
        .drop_duplicates(
            subset=keys,
            keep="first",
        )
        .sort_values(
            [
                "Gene Set",
                "Gene",
                "Technique",
                "Expected Study Accession",
            ],
            kind="stable",
        )
        .reset_index(drop=True)
    )


def write_registry(
    *,
    dataframe: pd.DataFrame,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary = (
        dataframe.groupby(
            [
                "Gene Set",
                "Technique",
            ],
            dropna=False,
        )
        .size()
        .reset_index(
            name="Validated Targets"
        )
    )

    gene_summary = (
        dataframe.groupby(
            "Gene Set",
            dropna=False,
        )["Gene"]
        .nunique()
        .reset_index(
            name="Genes With Validated Targets"
        )
    )

    with pd.ExcelWriter(
        output_path,
        engine="xlsxwriter",
    ) as writer:
        dataframe.to_excel(
            writer,
            sheet_name="Validated_Targets",
            index=False,
        )

        summary.to_excel(
            writer,
            sheet_name="Technique_Summary",
            index=False,
        )

        gene_summary.to_excel(
            writer,
            sheet_name="Gene_Summary",
            index=False,
        )

        workbook = writer.book

        header = workbook.add_format(
            {
                "bold": True,
                "font_color": "#FFFFFF",
                "bg_color": "#1F4E78",
                "border": 1,
                "align": "center",
                "valign": "vcenter",
            }
        )

        for sheet_name, frame in (
            ("Validated_Targets", dataframe),
            ("Technique_Summary", summary),
            ("Gene_Summary", gene_summary),
        ):
            worksheet = writer.sheets[
                sheet_name
            ]

            worksheet.freeze_panes(
                1,
                0,
            )

            worksheet.autofilter(
                0,
                0,
                len(frame),
                max(len(frame.columns) - 1, 0),
            )

            for index, column in enumerate(
                frame.columns
            ):
                worksheet.write(
                    0,
                    index,
                    column,
                    header,
                )

                if column == "Study Title":
                    width = 55
                elif "Accession" in column:
                    width = 24
                elif column in {
                    "Gene",
                    "Technique",
                    "Gene Set",
                }:
                    width = 18
                else:
                    width = 25

                worksheet.set_column(
                    index,
                    index,
                    width,
                )

    csv_path = output_path.with_suffix(
        ".csv"
    )

    dataframe.to_csv(
        csv_path,
        index=False,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build validated legacy dataset "
            "regression targets."
        )
    )

    parser.add_argument(
        "--rbp-audit",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--tf-audit",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    return parser


def main() -> int:
    args = build_parser().parse_args()

    rows = [
        *build_rows(
            audit_path=args.rbp_audit,
            gene_set="RBP",
        ),
        *build_rows(
            audit_path=args.tf_audit,
            gene_set="TF",
        ),
    ]

    registry = deduplicate(
        pd.DataFrame(rows)
    )

    write_registry(
        dataframe=registry,
        output_path=args.output,
    )

    print(
        f"Validated regression targets: "
        f"{len(registry)}"
    )

    print(
        "Unique genes: "
        f"{registry['Gene'].nunique()}"
    )

    print()
    print(
        registry.groupby(
            [
                "Gene Set",
                "Technique",
            ]
        )
        .size()
        .to_string()
    )

    print()
    print(
        f"Workbook: {args.output}"
    )

    print(
        "CSV: "
        f"{args.output.with_suffix('.csv')}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
