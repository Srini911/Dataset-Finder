import os
import re
import time
import urllib.parse as up
from pathlib import Path

import pandas as pd
import requests

from dataset_finder.builtin_gene_sets import load_builtin_gene_set

SPECIES = "Drosophila melanogaster"
MAX_RESULTS = 100
REQUEST_DELAY = 0.34

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")

OUTPUT = Path(
    "outputs/final_screening/"
    "Drosophila_TF_Technique_Datasets.xlsx"
)

CHECKPOINT = Path(
    "outputs/final_screening/"
    "Drosophila_TF_Technique_Datasets.checkpoint.csv"
)

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
        "User-Agent": "Dataset-Finder-fast-screen/1.0",
        "Accept": "application/json",
    }
)


def request_json(endpoint, params, attempts=4):
    query = {
        "tool": "dataset_finder_fast_screen",
        "email": "srinivas@umb.edu",
        **params,
    }

    if NCBI_API_KEY:
        query["api_key"] = NCBI_API_KEY

    last_error = None

    for attempt in range(1, attempts + 1):
        try:
            response = session.get(
                NCBI_BASE + endpoint,
                params=query,
                timeout=45,
            )
            response.raise_for_status()

            time.sleep(REQUEST_DELAY)

            return response.json()

        except Exception as exc:
            last_error = exc

            if attempt < attempts:
                time.sleep(min(2 ** attempt, 8))

    raise RuntimeError(str(last_error))


def esearch(database, term):
    payload = request_json(
        "esearch.fcgi",
        {
            "db": database,
            "term": term,
            "retmode": "json",
            "retmax": MAX_RESULTS,
        },
    )

    return payload.get(
        "esearchresult",
        {},
    ).get(
        "idlist",
        [],
    )


def esummary(database, ids):
    if not ids:
        return []

    payload = request_json(
        "esummary.fcgi",
        {
            "db": database,
            "id": ",".join(ids),
            "retmode": "json",
        },
    )

    result = payload.get("result", {})
    uids = result.get("uids", [])

    return [
        result[uid]
        for uid in uids
        if isinstance(
            result.get(uid),
            dict,
        )
    ]


def build_query(gene, terms):
    technique = " OR ".join(
        f'"{term}"[All Fields]'
        for term in terms
    )

    return (
        f'("{gene}"[All Fields]) AND '
        f"({technique}) AND "
        f'"{SPECIES}"[Organism]'
    )


def extract_accession(database, uid, entry):
    text = " ".join(
        str(value)
        for value in entry.values()
        if value
    )

    if database == "sra":
        patterns = [
            r"\bSRP\d+\b",
            r"\bERP\d+\b",
            r"\bDRP\d+\b",
        ]
    else:
        patterns = [
            r"\bGSE\d+\b",
            r"\bGDS\d+\b",
            r"\bGSM\d+\b",
        ]

    for pattern in patterns:
        match = re.search(
            pattern,
            text,
            flags=re.IGNORECASE,
        )

        if match:
            return match.group(0).upper()

    return str(
        entry.get("accession")
        or entry.get("acc")
        or uid
    )


def dataset_link(database, accession):
    if database == "SRA":
        return (
            "https://www.ncbi.nlm.nih.gov/sra/"
            f"?term={up.quote(accession)}"
        )

    return (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={up.quote(accession)}"
    )


genes = load_builtin_gene_set("tf")

print("TF genes:", len(genes))
print("Technique families:", len(TECHNIQUES))
print("Databases: SRA + GEO")
print()

rows = []

for index, gene in enumerate(genes, start=1):
    print(
        f"[{index}/{len(genes)}] {gene}",
        flush=True,
    )

    for technique, synonyms in TECHNIQUES.items():
        query = build_query(
            gene,
            synonyms,
        )

        for database, label in (
            ("sra", "SRA"),
            ("gds", "GEO"),
        ):
            try:
                ids = esearch(
                    database,
                    query,
                )

                entries = esummary(
                    database,
                    ids,
                )

            except Exception as exc:
                rows.append(
                    {
                        "Gene": gene,
                        "Technique": technique,
                        "Database": label,
                        "Accession": "",
                        "Title": "",
                        "Query": query,
                        "Status": "ERROR",
                        "Error": str(exc),
                        "Link": "",
                    }
                )

                continue

            for uid, entry in zip(
                ids,
                entries,
                strict=False,
            ):
                accession = extract_accession(
                    database,
                    uid,
                    entry,
                )

                rows.append(
                    {
                        "Gene": gene,
                        "Technique": technique,
                        "Database": label,
                        "Accession": accession,
                        "Title": str(
                            entry.get(
                                "title",
                                "",
                            )
                        ),
                        "Query": query,
                        "Status": "OK",
                        "Error": "",
                        "Link": dataset_link(
                            label,
                            accession,
                        ),
                    }
                )

    pd.DataFrame(rows).to_csv(
        CHECKPOINT,
        index=False,
    )


df = pd.DataFrame(rows)

if not df.empty:
    df = df.drop_duplicates(
        subset=[
            "Gene",
            "Technique",
            "Database",
            "Accession",
        ]
    )

with pd.ExcelWriter(
    OUTPUT,
    engine="xlsxwriter",
    engine_kwargs={
        "options": {
            "strings_to_urls": False,
        }
    },
) as writer:

    df.to_excel(
        writer,
        sheet_name="All_Datasets",
        index=False,
    )

    for technique in TECHNIQUES:
        subset = df[
            df["Technique"] == technique
        ].copy()

        subset.to_excel(
            writer,
            sheet_name=technique,
            index=False,
        )

    summary = (
        df[
            df["Status"] == "OK"
        ]
        .groupby(
            [
                "Gene",
                "Technique",
                "Database",
            ]
        )
        .size()
        .reset_index(name="Dataset Count")
    )

    summary.to_excel(
        writer,
        sheet_name="Summary",
        index=False,
    )

print()
print("COMPLETE")
print("Rows:", len(df))
print("Genes:", df["Gene"].nunique())
print("Output:", OUTPUT)
