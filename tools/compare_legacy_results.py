"""Compare legacy technique screens with Dataset Finder results."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

LEGACY_TECHNIQUES = {
    "RNA-Seq": "RNA_seq",
    "ChIP-Seq": "ChIP_seq",
    "CUT&RUN": "CUT_RUN",
    "CUT&Tag": "CUT_TAG",
    "eCLIP": "eCLIP",
}

NEW_COMPATIBLE_TECHNIQUES = {
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


def _normalize(value: object) -> str:
    if value is None:
        return ""

    return str(value).strip()


def _positive(value: object) -> bool:
    text = _normalize(value).casefold()

    return (
        text == "yes"
        or text.startswith("yes ")
        or text.startswith("yes(")
    )


def _legacy_gene_column(dataframe: pd.DataFrame) -> str:
    candidates = (
        "Gene Symbol",
        "Gene",
        "Symbol",
    )

    for candidate in candidates:
        if candidate in dataframe.columns:
            return candidate

    raise ValueError(
        "Could not identify the legacy gene-symbol column."
    )


def legacy_positive_pairs(
    dataframe: pd.DataFrame,
) -> set[tuple[str, str]]:
    """Return positive legacy gene-technique pairs."""
    gene_column = _legacy_gene_column(dataframe)
    pairs: set[tuple[str, str]] = set()

    for _, row in dataframe.iterrows():
        gene = _normalize(row.get(gene_column))

        if not gene:
            continue

        for legacy_column, normalized_technique in (
            LEGACY_TECHNIQUES.items()
        ):
            if legacy_column not in dataframe.columns:
                continue

            if _positive(row.get(legacy_column)):
                pairs.add(
                    (
                        gene.casefold(),
                        normalized_technique,
                    )
                )

    return pairs


def _split_techniques(value: object) -> set[str]:
    text = _normalize(value)

    if not text:
        return set()

    return {
        part.strip()
        for part in re.split(r"[;,|]", text)
        if part.strip()
    }


def new_positive_pairs(
    dataframe: pd.DataFrame,
) -> set[tuple[str, str]]:
    """Return verified gene-technique pairs from Dataset Finder."""
    required = {
        "Gene",
        "Technique",
    }

    missing = required.difference(dataframe.columns)

    if missing:
        raise ValueError(
            "New workbook is missing columns: "
            + ", ".join(sorted(missing))
        )

    pairs: set[tuple[str, str]] = set()

    for _, row in dataframe.iterrows():
        gene = _normalize(row.get("Gene"))

        if not gene:
            continue

        verified = _split_techniques(
            row.get("Technique")
        )

        technique_match = _normalize(
            row.get("Technique Match")
        )

        if technique_match.casefold() == "mismatch":
            continue

        for expected, compatible in (
            NEW_COMPATIBLE_TECHNIQUES.items()
        ):
            if verified.intersection(compatible):
                pairs.add(
                    (
                        gene.casefold(),
                        expected,
                    )
                )

    return pairs


def accession_lookup(
    dataframe: pd.DataFrame,
) -> dict[tuple[str, str], set[str]]:
    """Map gene-technique pairs to discovered accessions."""
    lookup: dict[
        tuple[str, str],
        set[str],
    ] = {}

    if "Gene" not in dataframe.columns:
        return lookup

    for _, row in dataframe.iterrows():
        gene = _normalize(row.get("Gene"))

        if not gene:
            continue

        technique = _normalize(
            row.get("Technique")
        )
        accession = _normalize(
            row.get("Accession")
        )
        technique_match = _normalize(
            row.get("Technique Match")
        )

        if technique_match.casefold() == "mismatch":
            continue

        for expected, compatible in (
            NEW_COMPATIBLE_TECHNIQUES.items()
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
    legacy_path: Path,
    new_path: Path,
    output_path: Path,
) -> None:
    """Create an Excel regression report."""
    legacy_df = pd.read_excel(legacy_path)
    new_df = pd.read_excel(
        new_path,
        sheet_name="All_Datasets",
    )

    legacy_pairs = legacy_positive_pairs(
        legacy_df
    )
    discovered_pairs = new_positive_pairs(
        new_df
    )
    accessions = accession_lookup(new_df)

    rows: list[dict[str, object]] = []

    for gene, technique in sorted(legacy_pairs):
        recovered = (
            gene,
            technique,
        ) in discovered_pairs

        rows.append(
            {
                "Gene": gene,
                "Legacy Technique": technique,
                "Legacy Positive": "Yes",
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

    comparison_df = pd.DataFrame(rows)

    recovered_count = int(
        (
            comparison_df[
                "Recovered By New Pipeline"
            ]
            == "Yes"
        ).sum()
    )

    total = len(comparison_df)
    missed = total - recovered_count

    recovery_rate = (
        recovered_count / total * 100
        if total
        else 0.0
    )

    summary_df = pd.DataFrame(
        [
            {
                "Legacy Positive Pairs": total,
                "Recovered": recovered_count,
                "Missed": missed,
                "Recovery Rate (%)": round(
                    recovery_rate,
                    2,
                ),
            }
        ]
    )

    missing_df = comparison_df[
        comparison_df[
            "Recovered By New Pipeline"
        ]
        == "No"
    ].copy()

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
            sheet_name="All_Legacy_Pairs",
            index=False,
        )
        missing_df.to_excel(
            writer,
            sheet_name="Missing_From_New",
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

        missing_format = workbook.add_format(
            {
                "bg_color": "#FCE4D6",
            }
        )

        recovered_format = workbook.add_format(
            {
                "bg_color": "#E2F0D9",
            }
        )

        for sheet_name, dataframe in (
            ("Summary", summary_df),
            ("All_Legacy_Pairs", comparison_df),
            ("Missing_From_New", missing_df),
        ):
            worksheet = writer.sheets[
                sheet_name
            ]

            worksheet.freeze_panes(1, 0)

            for column_index, column in enumerate(
                dataframe.columns
            ):
                worksheet.write(
                    0,
                    column_index,
                    column,
                    header_format,
                )

                width = min(
                    max(
                        len(column) + 2,
                        18,
                    ),
                    55,
                )

                worksheet.set_column(
                    column_index,
                    column_index,
                    width,
                )

        if not comparison_df.empty:
            recovered_column = (
                comparison_df.columns.get_loc(
                    "Recovered By New Pipeline"
                )
            )

            worksheet = writer.sheets[
                "All_Legacy_Pairs"
            ]

            worksheet.conditional_format(
                1,
                recovered_column,
                len(comparison_df),
                recovered_column,
                {
                    "type": "text",
                    "criteria": "containing",
                    "value": "Yes",
                    "format": recovered_format,
                },
            )

            worksheet.conditional_format(
                1,
                recovered_column,
                len(comparison_df),
                recovered_column,
                {
                    "type": "text",
                    "criteria": "containing",
                    "value": "No",
                    "format": missing_format,
                },
            )

    print(f"Legacy positive pairs: {total}")
    print(f"Recovered: {recovered_count}")
    print(f"Missing: {missed}")
    print(
        "Recovery rate: "
        f"{recovery_rate:.2f}%"
    )
    print(f"Report: {output_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare legacy technique-screen positives "
            "with Dataset Finder results."
        )
    )

    parser.add_argument(
        "--legacy",
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
        legacy_path=args.legacy,
        new_path=args.new,
        output_path=args.output,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
