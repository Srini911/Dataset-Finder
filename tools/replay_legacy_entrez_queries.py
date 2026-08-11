"""Replay legacy-style Entrez queries for remaining regression misses."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd
import requests

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SPECIES = "Drosophila melanogaster"

TECHNIQUES = {
    "RNA_seq": [
        "RNA-seq",
        "RNA seq",
        "RNAseq",
        "transcriptome sequencing",
        "transcriptomic sequencing",
    ],
    "ChIP_seq": [
        "ChIP-seq",
        "ChIP seq",
        "ChIPseq",
        "chromatin immunoprecipitation sequencing",
    ],
}

session = requests.Session()


def esearch(term: str, retmax: int = 100) -> list[str]:
    response = session.get(
        BASE + "esearch.fcgi",
        params={
            "db": "sra",
            "term": term,
            "retmode": "json",
            "retmax": retmax,
        },
        timeout=60,
    )
    response.raise_for_status()
    time.sleep(0.34)

    return [
        str(value)
        for value in response.json()
        .get("esearchresult", {})
        .get("idlist", [])
    ]


def esummary(ids: list[str]) -> list[dict]:
    if not ids:
        return []

    response = session.get(
        BASE + "esummary.fcgi",
        params={
            "db": "sra",
            "id": ",".join(ids),
            "retmode": "json",
        },
        timeout=60,
    )
    response.raise_for_status()
    time.sleep(0.34)

    result = response.json().get(
        "result",
        {},
    )

    return [
        result[uid]
        for uid in result.get("uids", [])
        if uid in result
    ]


def build_query(
    gene: str,
    technique: str,
) -> str:
    technique_query = " OR ".join(
        f'"{term}"[All Fields]'
        for term in TECHNIQUES[technique]
    )

    return (
        f'("{gene}"[All Fields]) AND '
        f"({technique_query}) AND "
        f'"{SPECIES}"[Organism]'
    )


def summary_contains_study(
    summary: dict,
    study: str,
) -> bool:
    text = " ".join(
        str(value)
        for value in summary.values()
        if value is not None
    ).upper()

    return study.upper() in text


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--missing-report",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    missing = pd.read_excel(
        args.missing_report,
        sheet_name="True_Missing",
    )

    rows = []

    cache: dict[
        tuple[str, str],
        tuple[list[str], list[dict]],
    ] = {}

    total = len(missing)

    for number, (_, row) in enumerate(
        missing.iterrows(),
        start=1,
    ):
        gene = str(row["Gene"]).strip()
        technique = str(
            row["Technique"]
        ).strip()

        study = str(
            row["Expected Study Accession"]
        ).strip()

        key = (
            gene.casefold(),
            technique,
        )

        print(
            f"[{number}/{total}] "
            f"{gene} | {technique} | {study}"
        )

        if technique not in TECHNIQUES:
            rows.append(
                {
                    "Gene Set": row["Gene Set"],
                    "Gene": gene,
                    "Technique": technique,
                    "Expected Study": study,
                    "Legacy Query Reproduced": "Unsupported",
                    "Candidate Count": 0,
                    "Search Query": "",
                    "Study Title": row["Study Title"],
                }
            )
            continue

        if key not in cache:
            query = build_query(
                gene,
                technique,
            )

            ids = esearch(
                query,
                retmax=100,
            )

            summaries = esummary(ids)

            cache[key] = (
                ids,
                summaries,
            )
        else:
            query = build_query(
                gene,
                technique,
            )

            ids, summaries = cache[key]

        found = any(
            summary_contains_study(
                summary,
                study,
            )
            for summary in summaries
        )

        rows.append(
            {
                "Gene Set": row["Gene Set"],
                "Gene": gene,
                "Technique": technique,
                "Expected Study": study,
                "Legacy Query Reproduced": (
                    "Yes"
                    if found
                    else "No"
                ),
                "Candidate Count": len(ids),
                "Search Query": query,
                "Study Title": row["Study Title"],
            }
        )

    result = pd.DataFrame(rows)

    with pd.ExcelWriter(
        args.output,
        engine="xlsxwriter",
    ) as writer:
        result.to_excel(
            writer,
            sheet_name="Replay",
            index=False,
        )

        summary = (
            result.groupby(
                [
                    "Gene Set",
                    "Legacy Query Reproduced",
                ]
            )
            .size()
            .reset_index(name="Count")
        )

        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

    print("\n===== LEGACY QUERY REPLAY =====\n")

    print(
        result[
            "Legacy Query Reproduced"
        ]
        .value_counts()
        .to_string()
    )

    print("\nBY GENE SET\n")

    print(
        pd.crosstab(
            result["Gene Set"],
            result["Legacy Query Reproduced"],
        ).to_string()
    )

    reproduced = result[
        result["Legacy Query Reproduced"]
        == "Yes"
    ]

    print(
        "\nOLD ENTrez QUERY REALLY FINDS "
        "THE EXPECTED STUDY\n"
    )

    if reproduced.empty:
        print("None")
    else:
        print(
            reproduced[
                [
                    "Gene Set",
                    "Gene",
                    "Technique",
                    "Expected Study",
                    "Candidate Count",
                    "Study Title",
                ]
            ].to_string(index=False)
        )

    print("\nReport:", args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
