"""Compare verified legacy technique calls with Dataset Finder results."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

COMPATIBLE_TECHNIQUES = {
    "RNA_seq": {
        "RNA_seq",
        "scRNA_seq",
        "snRNA_seq",
    },
    "ChIP_seq": {
        "ChIP_seq",
    },
    "CUT_RUN": {
        "CUT_RUN",
    },
    "CUT_TAG": {
        "CUT_TAG",
    },
    "eCLIP": {
        "eCLIP",
        "CLIP",
        "iCLIP",
        "PAR_CLIP",
        "HITS_CLIP",
    },
}


def normalize(value: object) -> str:
    if value is None or pd.isna(value):
        return ""

    return str(value).strip()


def verified_legacy_pairs(
    audit_df: pd.DataFrame,
) -> set[tuple[str, str]]:
    """Return only repository-verified legacy positives."""
    verified = audit_df[
        audit_df["Legacy Verdict"]
        == "VERIFIED_LEGACY_POSITIVE"
    ]

    return {
        (
            normalize(row["Gene"]).casefold(),
            normalize(row["Legacy Technique"]),
        )
        for _, row in verified.iterrows()
        if normalize(row["Gene"])
        and normalize(row["Legacy Technique"])
    }


def new_pairs(
    new_df: pd.DataFrame,
) -> set[tuple[str, str]]:
    """Return verified gene-technique pairs from the new pipeline."""
    pairs: set[tuple[str, str]] = set()

    for _, row in new_df.iterrows():
        gene = normalize(row.get("Gene"))

        if not gene:
            continue

        technique = normalize(
            row.get("Technique")
        )

        technique_match = normalize(
            row.get("Technique Match")
        ).casefold()

        if technique_match == "mismatch":
            continue

        for expected, compatible in (
            COMPATIBLE_TECHNIQUES.items()
        ):
            if technique in compatible:
                pairs.add(
                    (
                        gene.casefold(),
                        expected,
                    )
                )

    return pairs


def accession_lookup(
    new_df: pd.DataFrame,
) -> dict[tuple[str, str], set[str]]:
    """Map gene-technique pairs to stable accessions."""
    lookup: dict[
        tuple[str, str],
        set[str],
    ] = {}

    for _, row in new_df.iterrows():
        gene = normalize(row.get("Gene"))
        technique = normalize(
            row.get("Technique")
        )
        accession = normalize(
            row.get("Accession")
        )

        if not gene:
            continue

        if normalize(
            row.get("Technique Match")
        ).casefold() == "mismatch":
            continue

        for expected, compatible in (
            COMPATIBLE_TECHNIQUES.items()
        ):
            if technique not in compatible:
                continue

            key = (
                gene.casefold(),
                expected,
            )

            lookup.setdefault(
                key,
                set(),
            )

            if accession:
                lookup[key].add(accession)

    return lookup


def compare(
    *,
    audit_path: Path,
    new_path: Path,
    output_path: Path,
) -> None:
    """Compare verified legacy positives with new results."""
    audit_df = pd.read_excel(
        audit_path,
        sheet_name="All_Legacy_Calls",
    )

    new_df = pd.read_excel(
        new_path,
        sheet_name="All_Datasets",
    )

    legacy_pairs = verified_legacy_pairs(
        audit_df
    )

    discovered_pairs = new_pairs(
        new_df
    )

    accessions = accession_lookup(
        new_df
    )

    rows: list[dict[str, object]] = []

    for gene, technique in sorted(
        legacy_pairs
    ):
        recovered = (
            gene,
            technique,
        ) in discovered_pairs

        rows.append(
            {
                "Gene": gene,
                "Verified Legacy Technique": technique,
                "Recovered By New Pipeline": (
                    "Yes"
                    if recovered
                    else "No"
                ),
                "New Accessions": "; ".join(
                    sorted(
                        accessions.get(
                            (gene, technique),
                            set(),
                        )
                    )
                ),
            }
        )

    comparison_df = pd.DataFrame(
        rows
    )

    recovered_count = int(
        (
            comparison_df[
                "Recovered By New Pipeline"
            ]
            == "Yes"
        ).sum()
    )

    total = len(
        comparison_df
    )

    missing_count = (
        total
        - recovered_count
    )

    recovery_rate = (
        recovered_count
        / total
        * 100
        if total
        else 0.0
    )

    missing_df = comparison_df[
        comparison_df[
            "Recovered By New Pipeline"
        ]
        == "No"
    ].copy()

    summary_df = pd.DataFrame(
        [
            {
                "Verified Legacy Positives": total,
                "Recovered": recovered_count,
                "Missing": missing_count,
                "Recovery Rate (%)": round(
                    recovery_rate,
                    2,
                ),
            }
        ]
    )

    with pd.ExcelWriter(
        output_path,
        engine="xlsxwriter",
    ) as writer:
        summary_df.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

        comparison_df.to_excel(
            writer,
            sheet_name="Verified_Legacy_Pairs",
            index=False,
        )

        missing_df.to_excel(
            writer,
            sheet_name="True_Missing",
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

        green_format = workbook.add_format(
            {
                "bg_color": "#E2F0D9",
            }
        )

        red_format = workbook.add_format(
            {
                "bg_color": "#FCE4D6",
            }
        )

        for sheet_name, dataframe in (
            ("Summary", summary_df),
            (
                "Verified_Legacy_Pairs",
                comparison_df,
            ),
            (
                "True_Missing",
                missing_df,
            ),
        ):
            worksheet = writer.sheets[
                sheet_name
            ]

            worksheet.freeze_panes(
                1,
                0,
            )

            for column_index, column in enumerate(
                dataframe.columns
            ):
                worksheet.write(
                    0,
                    column_index,
                    column,
                    header_format,
                )

                worksheet.set_column(
                    column_index,
                    column_index,
                    28,
                )

        if not comparison_df.empty:
            status_column = (
                comparison_df.columns.get_loc(
                    "Recovered By New Pipeline"
                )
            )

            worksheet = writer.sheets[
                "Verified_Legacy_Pairs"
            ]

            worksheet.conditional_format(
                1,
                status_column,
                len(comparison_df),
                status_column,
                {
                    "type": "text",
                    "criteria": "containing",
                    "value": "Yes",
                    "format": green_format,
                },
            )

            worksheet.conditional_format(
                1,
                status_column,
                len(comparison_df),
                status_column,
                {
                    "type": "text",
                    "criteria": "containing",
                    "value": "No",
                    "format": red_format,
                },
            )

    print(
        f"Verified legacy positives: {total}"
    )
    print(
        f"Recovered: {recovered_count}"
    )
    print(
        f"Missing: {missing_count}"
    )
    print(
        "Verified recovery rate: "
        f"{recovery_rate:.2f}%"
    )
    print(
        f"Report: {output_path}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare verified legacy positives "
            "with Dataset Finder results."
        )
    )

    parser.add_argument(
        "--audit",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--new",
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

    compare(
        audit_path=args.audit,
        new_path=args.new,
        output_path=args.output,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
