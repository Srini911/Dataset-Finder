"""Check Dataset Finder workbooks against known validated dataset targets."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

KNOWN_DATASETS = {
    ("CG7804", "RNA_seq"): {"SRP110269"},
    ("CG7804", "ChIP_seq"): {"SRP110266"},
    ("Hrb98DE", "RNA_seq"): {
        "SRP186005",
        "SRP001537",
    },
    ("bru1", "RNA_seq"): {
        "SRP377648",
        "SRP050336",
    },
    ("snf", "RNA_seq"): {"SRP055034"},
    ("snf", "ChIP_seq"): {"SRP131779"},
}


def load_workbooks(paths: list[Path]) -> pd.DataFrame:
    frames = []

    for path in paths:
        frame = pd.read_excel(
            path,
            sheet_name="All_Datasets",
        )

        frame["Source Workbook"] = path.name
        frames.append(frame)

    if not frames:
        return pd.DataFrame()

    return pd.concat(
        frames,
        ignore_index=True,
    )


def check(workbook_paths: list[Path]) -> int:
    dataframe = load_workbooks(
        workbook_paths
    )

    dataframe["gene_key"] = (
        dataframe["Gene"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
    )

    dataframe["accession_key"] = (
        dataframe["Accession"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

    dataframe["technique_key"] = (
        dataframe["Technique"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    rows = []
    total = 0
    recovered = 0

    for (
        gene,
        technique,
    ), accessions in KNOWN_DATASETS.items():
        for accession in sorted(accessions):
            total += 1

            matches = dataframe[
                (
                    dataframe["gene_key"]
                    == gene.casefold()
                )
                & (
                    dataframe["technique_key"]
                    == technique
                )
                & (
                    dataframe["accession_key"]
                    == accession.upper()
                )
            ]

            found = not matches.empty

            if found:
                recovered += 1

            rows.append(
                {
                    "Gene": gene,
                    "Technique": technique,
                    "Expected Accession": accession,
                    "Recovered": (
                        "Yes"
                        if found
                        else "No"
                    ),
                    "Source Workbook": (
                        "; ".join(
                            sorted(
                                set(
                                    matches[
                                        "Source Workbook"
                                    ].astype(str)
                                )
                            )
                        )
                        if found
                        else ""
                    ),
                }
            )

    report = pd.DataFrame(rows)

    print(report.to_string(index=False))
    print()
    print(
        f"Recovered known datasets: "
        f"{recovered}/{total}"
    )
    print(
        "Recovery rate: "
        f"{recovered / total * 100:.2f}%"
    )

    return 0 if recovered == total else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--workbook",
        required=True,
        type=Path,
        nargs="+",
    )

    return parser


def main() -> int:
    args = build_parser().parse_args()

    return check(
        args.workbook
    )


if __name__ == "__main__":
    raise SystemExit(main())
