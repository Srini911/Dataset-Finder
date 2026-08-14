"""Compare new Dataset Finder workbooks against the validated legacy registry."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_new_workbooks(paths: list[Path]) -> pd.DataFrame:
    frames = []

    for path in paths:
        frame = pd.read_excel(
            path,
            sheet_name="All_Datasets",
        )

        frame["Source Workbook"] = path.name
        frames.append(frame)

    return pd.concat(
        frames,
        ignore_index=True,
    )


def normalize_text(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
    )


def compare(
    *,
    registry_path: Path,
    new_paths: list[Path],
    output_path: Path,
) -> None:
    registry = pd.read_excel(
        registry_path,
        sheet_name="Validated_Targets",
    )

    new = load_new_workbooks(
        new_paths
    )

    registry["gene_key"] = (
        normalize_text(registry["Gene"])
        .str.casefold()
    )

    registry["technique_key"] = normalize_text(
        registry["Technique"]
    )

    registry["study_key"] = (
        normalize_text(
            registry["Expected Study Accession"]
        )
        .str.upper()
    )

    new["gene_key"] = (
        normalize_text(new["Gene"])
        .str.casefold()
    )

    new["technique_key"] = normalize_text(
        new["Technique"]
    )

    new["accession_key"] = (
        normalize_text(new["Accession"])
        .str.upper()
    )

    related_accession_columns = (
        "Project Accession",
        "Related Accessions",
        "Related GEO Accessions",
        "Related SRA / ENA Studies",
        "Related BioProjects",
        "Related Study Accessions",
    )

    available_related_columns = [
        column
        for column in related_accession_columns
        if column in new.columns
    ]

    if available_related_columns:
        new["related_accession_key"] = (
            new[
                available_related_columns
            ]
            .fillna("")
            .astype(str)
            .agg(
                " ; ".join,
                axis=1,
            )
            .str.upper()
        )
    else:
        new["related_accession_key"] = ""

    rows = []

    for _, target in registry.iterrows():
        direct = new[
            (new["gene_key"] == target["gene_key"])
            & (
                new["technique_key"]
                == target["technique_key"]
            )
            & (
                new["accession_key"]
                == target["study_key"]
            )
        ]

        related = new[
            (new["gene_key"] == target["gene_key"])
            & (
                new["technique_key"]
                == target["technique_key"]
            )
            & (
                new["related_accession_key"]
                .str.contains(
                    target["study_key"],
                    regex=False,
                )
            )
        ]

        combined = pd.concat(
            [direct, related],
            ignore_index=True,
        ).drop_duplicates()

        recovered = not combined.empty

        rows.append(
            {
                "Gene Set": target["Gene Set"],
                "Gene": target["Gene"],
                "Technique": target["Technique"],
                "Expected Study Accession": (
                    target[
                        "Expected Study Accession"
                    ]
                ),
                "Recovered": (
                    "Yes"
                    if recovered
                    else "No"
                ),
                "Recovered Accessions": (
                    "; ".join(
                        sorted(
                            set(
                                combined[
                                    "Accession"
                                ]
                                .fillna("")
                                .astype(str)
                            )
                        )
                    )
                    if recovered
                    else ""
                ),
                "Source Workbook": (
                    "; ".join(
                        sorted(
                            set(
                                combined[
                                    "Source Workbook"
                                ]
                                .fillna("")
                                .astype(str)
                            )
                        )
                    )
                    if recovered
                    else ""
                ),
                "Study Title": target["Study Title"],
            }
        )

    comparison = pd.DataFrame(
        rows
    )

    recovered_count = (
        comparison["Recovered"]
        == "Yes"
    ).sum()

    total = len(comparison)

    missing = comparison[
        comparison["Recovered"]
        == "No"
    ].copy()

    summary = (
        comparison.groupby(
            [
                "Gene Set",
                "Technique",
                "Recovered",
            ]
        )
        .size()
        .reset_index(
            name="Count"
        )
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with pd.ExcelWriter(
        output_path,
        engine="xlsxwriter",
    ) as writer:
        comparison.to_excel(
            writer,
            sheet_name="All_Targets",
            index=False,
        )

        missing.to_excel(
            writer,
            sheet_name="True_Missing",
            index=False,
        )

        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

    print(
        f"Validated targets: {total}"
    )

    print(
        f"Recovered: {recovered_count}"
    )

    print(
        f"Missing: {total - recovered_count}"
    )

    print(
        "Recovery rate: "
        f"{recovered_count / total * 100:.2f}%"
    )

    print()
    print(
        comparison.groupby(
            [
                "Gene Set",
                "Recovered",
            ]
        )
        .size()
        .to_string()
    )

    print()
    print(
        f"Report: {output_path}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--registry",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--new",
        required=True,
        nargs="+",
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
        registry_path=args.registry,
        new_paths=args.new,
        output_path=args.output,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
