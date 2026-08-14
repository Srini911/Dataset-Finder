"""Direct Drosophila RBP/TF public dataset screening."""

from __future__ import annotations

import argparse
import os
import re
import time
import urllib.parse as up
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from dataset_finder.builtin_gene_sets import load_builtin_gene_set
from dataset_finder.flybase_resolver import FlyBaseResolver

EMAIL = "srinivas@umb.edu"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SPECIES = "Drosophila melanogaster"
MAX_RESULTS = 100
REQUEST_DELAY = 0.34


TECHNIQUES = {
    "CUT_RUN": [
        "CUT&RUN",
        "CUT and RUN",
        "CUT-AND-RUN",
        "CUT N RUN",
        "CUT-RUN",
        "cleavage under targets and release using nuclease",
    ],
    "CUT_TAG": [
        "CUT&Tag",
        "CUT and Tag",
        "CUT-TAG",
        "CUT N TAG",
        "CUT-AND-TAG",
        "cleavage under targets and tagmentation",
    ],
    "ChIP_seq": [
        "ChIP-seq",
        "ChIP seq",
        "ChIPseq",
        "chromatin immunoprecipitation sequencing",
    ],
    "CLIP": [
        "CLIP",
        "CLIP-seq",
        "CLIP seq",
        "HITS-CLIP",
        "iCLIP",
        "PAR-CLIP",
        "eCLIP",
        "crosslinking immunoprecipitation",
    ],
    "RNA_seq": [
        "RNA-seq",
        "RNA seq",
        "RNAseq",
        "transcriptome sequencing",
        "transcriptomic sequencing",
    ],
}


session = requests.Session()
session.headers.update(
    {
        "User-Agent": "Dataset-Finder-direct-screen/1.0",
        "Accept": "application/json",
    }
)


def request_json(
    endpoint: str,
    params: dict[str, Any],
    attempts: int = 4,
) -> dict:
    query = {
        "tool": "dataset_finder_direct_screen",
        "email": EMAIL,
        **params,
    }

    if NCBI_API_KEY:
        query["api_key"] = NCBI_API_KEY

    last_error: Exception | None = None

    for attempt in range(1, attempts + 1):
        try:
            response = session.get(
                NCBI_BASE + endpoint,
                params=query,
                timeout=60,
            )
            response.raise_for_status()
            time.sleep(REQUEST_DELAY)
            return response.json()

        except Exception as exc:
            last_error = exc

            if attempt < attempts:
                time.sleep(
                    min(
                        2**attempt,
                        10,
                    )
                )

    raise RuntimeError(
        f"NCBI request failed after {attempts} attempts: "
        f"{last_error}"
    )


def esearch(
    database: str,
    term: str,
) -> list[str]:
    payload = request_json(
        "esearch.fcgi",
        {
            "db": database,
            "term": term,
            "retmode": "json",
            "retmax": MAX_RESULTS,
        },
    )

    return [
        str(value)
        for value in (
            payload
            .get("esearchresult", {})
            .get("idlist", [])
        )
    ]


def esummary(
    database: str,
    identifiers: list[str],
) -> list[dict[str, Any]]:
    if not identifiers:
        return []

    entries: list[dict[str, Any]] = []

    for start in range(
        0,
        len(identifiers),
        100,
    ):
        batch = identifiers[
            start : start + 100
        ]

        payload = request_json(
            "esummary.fcgi",
            {
                "db": database,
                "id": ",".join(batch),
                "retmode": "json",
            },
        )

        result = payload.get(
            "result",
            {},
        )

        for uid in result.get(
            "uids",
            [],
        ):
            entry = result.get(uid)

            if isinstance(
                entry,
                dict,
            ):
                entry = dict(entry)
                entry["_uid"] = str(uid)
                entries.append(entry)

    return entries


def clean(value: object) -> str:
    if value is None:
        return ""

    return str(value).strip()


def unique_terms(
    values: list[str],
) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()

    for value in values:
        value = clean(value)

        if not value:
            continue

        key = value.casefold()

        if key in seen:
            continue

        seen.add(key)
        output.append(value)

    return output


def gene_search_terms(
    gene: str,
) -> tuple[list[str], list[str], dict[str, str]]:
    resolver = FlyBaseResolver()
    resolved = resolver.resolve(gene)

    canonical = unique_terms(
        [
            resolved.flybase_id,
            resolved.annotation_id,
            resolved.current_fullname,
        ]
    )

    safe_synonyms = []

    for synonym in resolved.synonyms:
        compact = re.sub(
            r"[^A-Za-z0-9]+",
            "",
            synonym,
        )

        if (
            len(compact) >= 4
            and not synonym.startswith("NEST:")
            and not synonym.startswith("BcDNA:")
            and not synonym.startswith("ms(")
        ):
            safe_synonyms.append(
                synonym
            )

    canonical.extend(
        safe_synonyms[:10]
    )

    submitted_compact = re.sub(
        r"[^A-Za-z0-9]+",
        "",
        gene,
    )

    official_compact = re.sub(
        r"[^A-Za-z0-9]+",
        "",
        resolved.official_symbol,
    )

    if (
        resolved.official_symbol
        and len(official_compact) >= 4
    ):
        canonical.append(
            resolved.official_symbol
        )

    if len(submitted_compact) >= 4:
        canonical.append(gene)

    canonical = unique_terms(
        canonical
    )

    fallback = unique_terms(
        [
            gene,
            resolved.official_symbol,
        ]
    )

    risky_symbol = (
        len(submitted_compact) <= 3
        or resolved.ambiguous
    )

    metadata = {
        "Official Symbol":
            resolved.official_symbol,
        "FlyBase ID":
            resolved.flybase_id,
        "Annotation ID":
            resolved.annotation_id,
        "Full Name":
            resolved.current_fullname,
        "Resolver Ambiguous":
            (
                "Yes"
                if resolved.ambiguous
                else "No"
            ),
        "Risky Symbol":
            (
                "Yes"
                if risky_symbol
                else "No"
            ),
    }

    return (
        canonical,
        fallback,
        metadata,
    )


def build_search_term(
    gene_terms: list[str],
    technique_terms: list[str],
) -> str:
    gene_query = " OR ".join(
        f'"{term}"[All Fields]'
        for term in gene_terms
    )

    technique_query = " OR ".join(
        f'"{term}"[All Fields]'
        for term in technique_terms
    )

    return (
        f"(({gene_query}) AND "
        f"({technique_query})) AND "
        f'"{SPECIES}"[Organism]'
    )


def summary_text(
    entry: dict[str, Any],
) -> str:
    return " ".join(
        clean(value)
        for value in entry.values()
        if value
    )


def extract_accession(
    database: str,
    entry: dict[str, Any],
) -> str:
    text = summary_text(
        entry
    )

    if database == "sra":
        patterns = (
            r"\bSRP\d+\b",
            r"\bERP\d+\b",
            r"\bDRP\d+\b",
        )
    else:
        patterns = (
            r"\bGSE\d+\b",
            r"\bGDS\d+\b",
            r"\bGSM\d+\b",
        )

    for pattern in patterns:
        match = re.search(
            pattern,
            text,
            flags=re.IGNORECASE,
        )

        if match:
            return match.group(
                0
            ).upper()

    for key in (
        "accession",
        "acc",
        "gse",
    ):
        value = clean(
            entry.get(key)
        )

        if value:
            return value

    return clean(
        entry.get("_uid")
    )


def dataset_link(
    database: str,
    accession: str,
) -> str:
    if database == "SRA":
        return (
            "https://www.ncbi.nlm.nih.gov/sra/"
            f"?term={up.quote(accession)}"
        )

    return (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={up.quote(accession)}"
    )


def collect_route(
    *,
    gene: str,
    gene_set: str,
    technique: str,
    technique_terms: list[str],
    database: str,
    database_label: str,
    gene_terms: list[str],
    route: str,
    metadata: dict[str, str],
) -> list[dict[str, object]]:
    if not gene_terms:
        return []

    query = build_search_term(
        gene_terms,
        technique_terms,
    )

    identifiers = esearch(
        database,
        query,
    )

    entries = esummary(
        database,
        identifiers,
    )

    rows: list[dict[str, object]] = []

    for entry in entries:
        accession = extract_accession(
            database,
            entry,
        )

        rows.append(
            {
                "Gene": gene,
                "Gene Set":
                    gene_set.upper(),
                **metadata,
                "Technique":
                    technique,
                "Database":
                    database_label,
                "Accession":
                    accession,
                "Title":
                    clean(
                        entry.get(
                            "title"
                        )
                    ),
                "Search Route":
                    route,
                "Gene Query Used":
                    "; ".join(
                        gene_terms
                    ),
                "Search Query Used":
                    query,
                "Link":
                    dataset_link(
                        database_label,
                        accession,
                    ),
                "Status":
                    "OK",
                "Error":
                    "",
            }
        )

    return rows


def collect_gene(
    *,
    gene: str,
    gene_set: str,
) -> list[dict[str, object]]:
    (
        canonical_terms,
        fallback_terms,
        metadata,
    ) = gene_search_terms(
        gene
    )

    rows: list[
        dict[str, object]
    ] = []

    for technique, synonyms in (
        TECHNIQUES.items()
    ):
        for database, label in (
            ("sra", "SRA"),
            ("gds", "GEO"),
        ):
            canonical_rows = []

            try:
                canonical_rows = collect_route(
                    gene=gene,
                    gene_set=gene_set,
                    technique=technique,
                    technique_terms=synonyms,
                    database=database,
                    database_label=label,
                    gene_terms=canonical_terms,
                    route="Canonical",
                    metadata=metadata,
                )

                rows.extend(
                    canonical_rows
                )

                fallback_rows = collect_route(
                    gene=gene,
                    gene_set=gene_set,
                    technique=technique,
                    technique_terms=synonyms,
                    database=database,
                    database_label=label,
                    gene_terms=fallback_terms,
                    route="Legacy symbol fallback",
                    metadata=metadata,
                )

                rows.extend(
                    fallback_rows
                )

            except Exception as exc:
                rows.append(
                    {
                        "Gene":
                            gene,
                        "Gene Set":
                            gene_set.upper(),
                        **metadata,
                        "Technique":
                            technique,
                        "Database":
                            label,
                        "Accession":
                            "",
                        "Title":
                            "",
                        "Search Route":
                            "",
                        "Gene Query Used":
                            "",
                        "Search Query Used":
                            "",
                        "Link":
                            "",
                        "Status":
                            "ERROR",
                        "Error":
                            str(exc),
                    }
                )

    return rows


def load_validated_legacy_pairs() -> set[tuple[str, str, str]]:
    registry_path = (
        Path.home()
        / "Downloads"
        / "Validated_Legacy_Dataset_Registry.xlsx"
    )

    if not registry_path.exists():
        return set()

    registry = pd.read_excel(
        registry_path,
        sheet_name="Validated_Targets",
    )

    pairs: set[tuple[str, str, str]] = set()

    for _, row in registry.iterrows():
        gene = clean(row.get("Gene")).casefold()
        technique = clean(row.get("Technique"))
        accession = clean(
            row.get("Expected Study Accession")
        ).upper()

        if gene and technique and accession:
            pairs.add(
                (
                    gene,
                    technique,
                    accession,
                )
            )

    return pairs


def write_workbook(
    dataframe: pd.DataFrame,
    output: Path,
) -> None:
    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    ok = dataframe[
        dataframe["Status"] == "OK"
    ].copy()

    ok = ok.drop_duplicates(
        subset=[
            "Gene",
            "Technique",
            "Database",
            "Accession",
            "Search Route",
        ]
    )

    validated_legacy_pairs = (
        load_validated_legacy_pairs()
    )

    def is_validated_legacy(
        row: pd.Series,
    ) -> bool:
        return (
            clean(row["Gene"]).casefold(),
            clean(row["Technique"]),
            clean(row["Accession"]).upper(),
        ) in validated_legacy_pairs

    ok["Validated Legacy"] = ok.apply(
        lambda row: (
            "Yes"
            if is_validated_legacy(row)
            else "No"
        ),
        axis=1,
    )

    risky_fallback = ok[
        (ok["Search Route"] == "Legacy symbol fallback")
        & (ok["Risky Symbol"] == "Yes")
        & (ok["Validated Legacy"] != "Yes")
    ].copy()

    accepted = ok[
        ~(
            (ok["Search Route"] == "Legacy symbol fallback")
            & (ok["Risky Symbol"] == "Yes")
            & (ok["Validated Legacy"] != "Yes")
        )
    ].copy()

    accepted = accepted.drop_duplicates(
        subset=[
            "Gene",
            "Technique",
            "Database",
            "Accession",
        ]
    )

    summary = (
        accepted.groupby(
            [
                "Gene",
                "Technique",
                "Database",
            ]
        )
        .size()
        .reset_index(
            name="Dataset Count"
        )
    )

    with pd.ExcelWriter(
        output,
        engine="xlsxwriter",
        engine_kwargs={
            "options": {
                "strings_to_urls":
                    False,
            }
        },
    ) as writer:
        accepted.to_excel(
            writer,
            sheet_name="All_Datasets",
            index=False,
        )

        risky_fallback.to_excel(
            writer,
            sheet_name="Legacy_Candidates",
            index=False,
        )

        ok.to_excel(
            writer,
            sheet_name="All_Discovery",
            index=False,
        )

        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

        for technique in TECHNIQUES:
            subset = accepted[
                accepted["Technique"]
                == technique
            ].copy()

            subset.to_excel(
                writer,
                sheet_name=technique,
                index=False,
            )

        errors = dataframe[
            dataframe["Status"]
            == "ERROR"
        ].copy()

        errors.to_excel(
            writer,
            sheet_name="Errors",
            index=False,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--gene-set",
        choices=[
            "rbp",
            "tf",
        ],
        default="tf",
    )

    parser.add_argument(
        "--genes",
        nargs="*",
        default=[],
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    return parser


def main() -> int:
    args = build_parser().parse_args()

    if args.genes:
        genes = args.genes
    else:
        genes = load_builtin_gene_set(
            args.gene_set
        )

    print(
        f"Gene set: {args.gene_set.upper()}"
    )
    print(
        f"Genes: {len(genes)}"
    )
    print(
        "Technique families:",
        len(TECHNIQUES),
    )
    print(
        "Databases: SRA + GEO"
    )
    print()

    rows: list[
        dict[str, object]
    ] = []

    for index, gene in enumerate(
        genes,
        start=1,
    ):
        print(
            f"[{index}/{len(genes)}] "
            f"{gene}",
            flush=True,
        )

        rows.extend(
            collect_gene(
                gene=gene,
                gene_set=args.gene_set,
            )
        )

    dataframe = pd.DataFrame(
        rows
    )

    write_workbook(
        dataframe,
        args.output,
    )

    ok = dataframe[
        dataframe["Status"] == "OK"
    ]

    print()
    print("=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(
        "Raw rows:",
        len(ok),
    )
    print(
        "Genes with hits:",
        ok["Gene"].nunique(),
    )
    print(
        "Errors:",
        (
            dataframe["Status"]
            == "ERROR"
        ).sum(),
    )
    print(
        "Output:",
        args.output,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
