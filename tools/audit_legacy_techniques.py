"""Audit legacy technique calls against resolved repository metadata."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

LEGACY_TECHNIQUES = {
    "RNA-Seq": "RNA_seq",
    "ChIP-Seq": "ChIP_seq",
    "CUT&RUN": "CUT_RUN",
    "CUT&Tag": "CUT_TAG",
    "eCLIP": "eCLIP",
}

LIBRARY_TECHNIQUE_MAP = {
    "rna-seq": "RNA_seq",
    "chip-seq": "ChIP_seq",
    "cut&run": "CUT_RUN",
    "cut and run": "CUT_RUN",
    "cut&tag": "CUT_TAG",
    "cut and tag": "CUT_TAG",
    "eclip": "eCLIP",
}


def normalize(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def positive(value: object) -> bool:
    text = normalize(value).casefold()
    return text == "yes" or text.startswith("yes")


def normalize_technique(value: object) -> str:
    text = normalize(value).casefold()

    for key, technique in LIBRARY_TECHNIQUE_MAP.items():
        if key in text:
            return technique

    return ""


def audit(
    *,
    legacy_path: Path,
    resolved_path: Path,
    output_path: Path,
) -> None:
    legacy = pd.read_excel(legacy_path)
    resolved = pd.read_excel(
        resolved_path,
        sheet_name="Resolved_Old_Results",
    )

    resolved["gene_key"] = (
        resolved["Gene"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
    )

    rows: list[dict[str, object]] = []

    for _, legacy_row in legacy.iterrows():
        gene = normalize(
            legacy_row.get("Gene Symbol")
        )

        if not gene:
            continue

        gene_key = gene.casefold()

        matches = resolved[
            resolved["gene_key"] == gene_key
        ]

        resolved_row = (
            matches.iloc[0]
            if not matches.empty
            else None
        )

        for legacy_column, expected_technique in (
            LEGACY_TECHNIQUES.items()
        ):
            if not positive(
                legacy_row.get(legacy_column)
            ):
                continue

            actual_technique = ""
            resolution_status = ""
            study_accession = ""
            experiment_accession = ""
            old_accession = normalize(
                legacy_row.get("Accession")
            )
            study_title = ""
            library_strategy = ""

            if resolved_row is not None:
                resolution_status = normalize(
                    resolved_row.get(
                        "resolution_status"
                    )
                )
                study_accession = normalize(
                    resolved_row.get(
                        "study_accession"
                    )
                )
                experiment_accession = normalize(
                    resolved_row.get(
                        "experiment_accession"
                    )
                )
                study_title = normalize(
                    resolved_row.get(
                        "study_title"
                    )
                )
                library_strategy = normalize(
                    resolved_row.get(
                        "library_strategy"
                    )
                )

                actual_technique = (
                    normalize_technique(
                        library_strategy
                    )
                    or normalize_technique(
                        study_title
                    )
                    or normalize_technique(
                        resolved_row.get(
                            "raw_title"
                        )
                    )
                )

            if not resolution_status:
                verdict = "UNRESOLVED"

            elif (
                resolution_status.casefold()
                != "resolved"
            ):
                verdict = "UNRESOLVED"

            elif not actual_technique:
                verdict = "NEEDS_MANUAL_REVIEW"

            elif actual_technique == expected_technique:
                verdict = "VERIFIED_LEGACY_POSITIVE"

            else:
                verdict = "LIKELY_LEGACY_FALSE_POSITIVE"

            rows.append(
                {
                    "Gene": gene,
                    "Legacy Technique": (
                        expected_technique
                    ),
                    "Old Accession": old_accession,
                    "Resolution Status": (
                        resolution_status
                    ),
                    "Resolved Study": (
                        study_accession
                    ),
                    "Resolved Experiment": (
                        experiment_accession
                    ),
                    "Library Strategy": (
                        library_strategy
                    ),
                    "Resolved Technique": (
                        actual_technique
                    ),
                    "Study Title": study_title,
                    "Legacy Verdict": verdict,
                }
            )

    audit_df = pd.DataFrame(rows)

    summary = (
        audit_df["Legacy Verdict"]
        .value_counts()
        .rename_axis("Legacy Verdict")
        .reset_index(name="Count")
    )

    with pd.ExcelWriter(
        output_path,
        engine="xlsxwriter",
    ) as writer:
        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

        audit_df.to_excel(
            writer,
            sheet_name="All_Legacy_Calls",
            index=False,
        )

        for verdict in (
            "VERIFIED_LEGACY_POSITIVE",
            "LIKELY_LEGACY_FALSE_POSITIVE",
            "UNRESOLVED",
            "NEEDS_MANUAL_REVIEW",
        ):
            subset = audit_df[
                audit_df["Legacy Verdict"]
                == verdict
            ]

            subset.to_excel(
                writer,
                sheet_name=verdict[:31],
                index=False,
            )

        workbook = writer.book

        header_format = workbook.add_format(
            {
                "bold": True,
                "font_color": "white",
                "bg_color": "#1F4E78",
                "border": 1,
            }
        )

        for sheet_name in writer.sheets:
            worksheet = writer.sheets[
                sheet_name
            ]
            dataframe = (
                summary
                if sheet_name == "Summary"
                else None
            )

            worksheet.freeze_panes(1, 0)

            if dataframe is not None:
                columns = dataframe.columns
            else:
                columns = audit_df.columns

            for index, column in enumerate(
                columns
            ):
                worksheet.write(
                    0,
                    index,
                    column,
                    header_format,
                )
                worksheet.set_column(
                    index,
                    index,
                    24,
                )

    print("Legacy technique audit completed.")
    print()
    print(
        audit_df["Legacy Verdict"]
        .value_counts()
        .to_string()
    )
    print()
    print(f"Report: {output_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate legacy technique calls "
            "against resolved repository metadata."
        )
    )

    parser.add_argument(
        "--legacy",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--resolved",
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

    audit(
        legacy_path=args.legacy,
        resolved_path=args.resolved,
        output_path=args.output,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
