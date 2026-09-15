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

NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")
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
        **params,
    }

    if NCBI_EMAIL:
        query["email"] = NCBI_EMAIL

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



def xml_attribute(
    text: str,
    tag: str,
    attribute: str,
) -> str:
    if not text:
        return ""

    match = re.search(
        rf"<{tag}\b[^>]*\b{attribute}=[\"']([^\"']+)[\"']",
        text,
        flags=re.IGNORECASE,
    )

    return clean(match.group(1)) if match else ""


def xml_text(
    text: str,
    tag: str,
) -> str:
    if not text:
        return ""

    match = re.search(
        rf"<{tag}\b[^>]*>(.*?)</{tag}>",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )

    if not match:
        return ""

    value = re.sub(
        r"<[^>]+>",
        " ",
        match.group(1),
    )

    return re.sub(r"\s+", " ", value).strip()


def extract_all_accessions(
    text: str,
    prefixes: tuple[str, ...],
) -> list[str]:
    if not text:
        return []

    prefix_pattern = "|".join(
        re.escape(prefix)
        for prefix in prefixes
    )

    values = re.findall(
        rf"\b(?:{prefix_pattern})\d+\b",
        text,
        flags=re.IGNORECASE,
    )

    return unique_terms(
        [value.upper() for value in values]
    )


def sra_metadata(
    entry: dict[str, Any],
) -> dict[str, str]:
    expxml = clean(entry.get("expxml"))
    runs = clean(entry.get("runs"))

    experiment_accession = xml_attribute(
        expxml,
        "Experiment",
        "acc",
    )

    sample_accession = xml_attribute(
        expxml,
        "Sample",
        "acc",
    )

    study_accession = xml_attribute(
        expxml,
        "Study",
        "acc",
    )

    study_title = xml_attribute(
        expxml,
        "Study",
        "name",
    )

    experiment_title = xml_text(
        expxml,
        "Title",
    )

    bioproject = xml_text(
        expxml,
        "Bioproject",
    )

    biosample = xml_text(
        expxml,
        "Biosample",
    )

    library_strategy = xml_text(
        expxml,
        "LIBRARY_STRATEGY",
    )

    library_source = xml_text(
        expxml,
        "LIBRARY_SOURCE",
    )

    library_selection = xml_text(
        expxml,
        "LIBRARY_SELECTION",
    )

    if re.search(
        r"<PAIRED\b",
        expxml,
        flags=re.IGNORECASE,
    ):
        library_layout = "PAIRED"
    elif re.search(
        r"<SINGLE\b",
        expxml,
        flags=re.IGNORECASE,
    ):
        library_layout = "SINGLE"
    else:
        library_layout = ""

    platform_match = re.search(
        r'instrument_model=[\"\']([^\"\']+)[\"\']',
        expxml,
        flags=re.IGNORECASE,
    )

    platform = (
        clean(platform_match.group(1))
        if platform_match
        else ""
    )

    run_accessions = extract_all_accessions(
        runs,
        ("SRR", "ERR", "DRR"),
    )

    biosample_accessions = unique_terms(
        [
            sample_accession,
            biosample,
            *extract_all_accessions(
                expxml,
                ("SAMN", "SAMEA", "SAMD"),
            ),
        ]
    )

    evidence = re.sub(
        r"<[^>]+>",
        " ",
        expxml,
    )
    evidence = re.sub(
        r"\s+",
        " ",
        evidence,
    ).strip()

    return {
        "Project Accession":
            bioproject or study_accession,
        "Study Accession":
            study_accession,
        "Experiment Accession":
            experiment_accession,
        "Run Accessions":
            "; ".join(run_accessions),
        "BioSample Accessions":
            "; ".join(biosample_accessions),
        "Sample Accession":
            sample_accession,
        "Experiment Title":
            experiment_title,
        "Study Title":
            study_title,
        "Library Strategy":
            library_strategy,
        "Library Source":
            library_source,
        "Library Selection":
            library_selection,
        "Library Layout":
            library_layout,
        "Platform":
            platform,
        "Metadata Evidence":
            evidence,
    }


def contains_gene(
    text: str,
    gene: str,
) -> bool:
    text = clean(text)

    if not text or not gene:
        return False

    return bool(
        re.search(
            rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])",
            text,
            flags=re.IGNORECASE,
        )
    )


def technique_matches_sample(
    technique: str,
    sample_metadata: dict[str, str],
) -> bool:
    text = " ".join(
        [
            clean(sample_metadata.get("Experiment Title")),
            clean(sample_metadata.get("Library Strategy")),
            clean(sample_metadata.get("Metadata Evidence")),
        ]
    ).casefold()

    strategy = clean(
        sample_metadata.get("Library Strategy")
    ).casefold()

    if technique == "RNA_seq":
        return (
            strategy in {"rna-seq", "rna seq", "transcriptomic"}
            or "rna-seq" in text
            or "rnaseq" in text
        )

    if technique == "ChIP_seq":
        return (
            "chip-seq" in text
            or "chip seq" in text
            or "chipseq" in text
        )

    if technique == "CUT_RUN":
        return (
            "cut&run" in text
            or "cut and run" in text
            or "cut-run" in text
            or "cut n run" in text
        )

    if technique == "CUT_TAG":
        return (
            "cut&tag" in text
            or "cut and tag" in text
            or "cut-tag" in text
            or "cut n tag" in text
        )

    if technique == "CLIP":
        return any(
            value in text
            for value in (
                "clip-seq",
                "clip seq",
                "eclip",
                "iclip",
                "par-clip",
                "hits-clip",
                "rip-seq",
            )
        )

    return True


def classify_group(
    gene: str,
    experiment_title: str,
    metadata_evidence: str,
) -> tuple[str, str, str]:
    title = clean(experiment_title)
    evidence = clean(metadata_evidence)

    title_lower = title.casefold()
    evidence_lower = evidence.casefold()

    title_normalized = re.sub(
        r"[_-]+",
        " ",
        title_lower,
    )
    evidence_normalized = re.sub(
        r"[_-]+",
        " ",
        evidence_lower,
    )

    target_in_title = contains_gene(
        title,
        gene,
    )
    target_in_evidence = contains_gene(
        evidence,
        gene,
    )

    gene_pattern = re.escape(gene.casefold())

    target_kd = bool(
        re.search(
            rf"\bkd\s+{gene_pattern}\b",
            evidence_normalized,
            flags=re.IGNORECASE,
        )
        or re.search(
            rf"\b{gene_pattern}\s+kd\b",
            evidence_normalized,
            flags=re.IGNORECASE,
        )
        or re.search(
            rf"\bknock\s*down\s+{gene_pattern}\b",
            evidence_normalized,
            flags=re.IGNORECASE,
        )
        or re.search(
            rf"\b{gene_pattern}\s+knock\s*down\b",
            evidence_normalized,
            flags=re.IGNORECASE,
        )
    )

    target_rnai = bool(
        target_in_evidence
        and re.search(
            r"\brnai\b",
            evidence_normalized,
            flags=re.IGNORECASE,
        )
    )

    negative_control_patterns = {
        "GFP RNAi control":
            r"\bgfp\s*rnai\b",
        "IgG control":
            r"\bigg\b",
        "Input control":
            r"\binput\b",
        "wild type":
            r"\bwild\s*type\b|\bwildtype\b|\bwt\b",
        "mock":
            r"\bmock\b",
        "vehicle":
            r"\bvehicle\b",
        "untreated":
            r"\buntreated\b",
        "control":
            r"\bcontrol\b",
    }

    perturbation_patterns = {
        "RNAi":
            r"\brnai\b",
        "knockdown":
            r"\bknock\s*down\b|\bkd\b",
        "knockout":
            r"\bknock\s*out\b|\bko\b",
        "mutant":
            r"\bmutant\b",
        "depletion":
            r"\bdeplet(?:ed|ion)\b",
        "overexpression":
            r"\boverexpress(?:ion|ed)?\b",
        "CRISPR":
            r"\bcrispr\b",
        "treated":
            r"\btreated\b",
    }

    control_hits = [
        label
        for label, pattern in negative_control_patterns.items()
        if re.search(
            pattern,
            title_normalized,
            flags=re.IGNORECASE,
        )
    ]

    perturbation_hits = [
        label
        for label, pattern in perturbation_patterns.items()
        if re.search(
            pattern,
            title_normalized,
            flags=re.IGNORECASE,
        )
    ]

    if target_in_title and perturbation_hits:
        return (
            "Experiment",
            "High",
            (
                f"Target gene {gene} present; "
                + "; ".join(perturbation_hits)
            ),
        )

    if target_in_title and (
        "rip" in title_lower
        or "clip" in title_lower
        or "chip" in title_lower
    ):
        return (
            "Experiment",
            "High",
            f"Target gene {gene} present in target-enrichment sample",
        )

    if control_hits and not target_in_title:
        return (
            "Control",
            "High",
            "; ".join(control_hits),
        )

    if perturbation_hits and not target_in_title:
        return (
            "Other",
            "High",
            (
                "Non-target perturbation: "
                + "; ".join(perturbation_hits)
            ),
        )

    if target_kd:
        return (
            "Experiment",
            "High",
            f"Target-specific knockdown evidence for {gene}",
        )

    if target_rnai:
        return (
            "Experiment",
            "High",
            f"Target-specific RNAi evidence for {gene}",
        )

    if re.search(
        r"\bgfp(?:\s+rnai)?\b",
        evidence_normalized,
        flags=re.IGNORECASE,
    ):
        return (
            "Control",
            "High",
            "GFP control evidence in sample metadata",
        )

    generic_perturbation = [
        label
        for label, pattern in perturbation_patterns.items()
        if re.search(
            pattern,
            evidence_normalized,
            flags=re.IGNORECASE,
        )
    ]

    if generic_perturbation and not target_in_evidence:
        return (
            "Other",
            "High",
            (
                "Non-target perturbation in metadata: "
                + "; ".join(generic_perturbation)
            ),
        )

    if target_in_title:
        return (
            "Experiment",
            "Medium",
            f"Target gene {gene} present in sample title",
        )

    return (
        "Unclear",
        "Low",
        "No target-specific control/experiment evidence",
    )


def classify_sex(
    experiment_title: str,
    metadata_evidence: str,
) -> tuple[str, str, str]:
    title = clean(experiment_title)
    evidence = clean(metadata_evidence)

    title_lower = title.casefold()
    evidence_lower = evidence.casefold()

    male_pattern = r"\bmales?\b"
    female_pattern = r"\bfemales?\b"

    title_male = bool(
        re.search(male_pattern, title_lower)
    )
    title_female = bool(
        re.search(female_pattern, title_lower)
    )

    if title_male and title_female:
        return (
            "Mixed",
            "High",
            "Male and female explicitly present in experiment title",
        )

    if title_male:
        return (
            "Male",
            "High",
            "Male explicitly present in experiment title",
        )

    if title_female:
        return (
            "Female",
            "High",
            "Female explicitly present in experiment title",
        )

    evidence_male = bool(
        re.search(male_pattern, evidence_lower)
    )
    evidence_female = bool(
        re.search(female_pattern, evidence_lower)
    )

    if evidence_male and evidence_female:
        return (
            "Unclear",
            "Low",
            "Conflicting male and female terms in metadata evidence",
        )

    if evidence_male:
        return (
            "Male",
            "Medium",
            "Male explicitly present in metadata evidence",
        )

    if evidence_female:
        return (
            "Female",
            "Medium",
            "Female explicitly present in metadata evidence",
        )

    return (
        "Unclear",
        "Low",
        "No explicit sex information found",
    )


def expand_sra_study(
    study_accession: str,
) -> list[dict[str, Any]]:
    if not study_accession:
        return []

    term = (
        f'"{study_accession}"[All Fields] AND '
        f'"{SPECIES}"[Organism]'
    )

    identifiers = esearch(
        "sra",
        term,
    )

    return esummary(
        "sra",
        identifiers,
    )

def empty_sample_metadata() -> dict[str, str]:
    return {
        "Project Accession": "",
        "Study Accession": "",
        "Experiment Accession": "",
        "Run Accessions": "",
        "BioSample Accessions": "",
        "Sample Accession": "",
        "Experiment Title": "",
        "Study Title": "",
        "Library Strategy": "",
        "Library Source": "",
        "Library Selection": "",
        "Library Layout": "",
        "Platform": "",
        "Metadata Evidence": "",
        "Group": "Unclear",
        "Group Confidence": "Low",
        "Group Evidence":
            "Sample-level metadata unavailable from this summary record",
        "Sex": "Unclear",
        "Sex Confidence": "Low",
        "Sex Evidence":
            "Sample-level sex metadata unavailable from this summary record",
    }

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

    if database == "sra":
        studies: list[str] = []

        for entry in entries:
            sm = sra_metadata(entry)

            study = clean(
                sm.get("Study Accession")
            )

            if study and study not in studies:
                studies.append(study)

        expanded_entries: list[dict[str, Any]] = []
        seen_experiments: set[str] = set()

        for study in studies:
            for entry in expand_sra_study(study):
                sm = sra_metadata(entry)

                experiment = clean(
                    sm.get("Experiment Accession")
                )

                identity = (
                    experiment
                    or clean(entry.get("_uid"))
                )

                if identity in seen_experiments:
                    continue

                seen_experiments.add(identity)
                expanded_entries.append(entry)

        if expanded_entries:
            entries = expanded_entries

    rows: list[dict[str, object]] = []

    for entry in entries:
        accession = extract_accession(
            database,
            entry,
        )

        if database == "sra":
            sample_metadata = sra_metadata(
                entry
            )

            if not technique_matches_sample(
                technique,
                sample_metadata,
            ):
                continue

            experiment_title = sample_metadata[
                "Experiment Title"
            ]

            (
                group,
                group_confidence,
                group_evidence,
            ) = classify_group(
                gene,
                experiment_title,
                sample_metadata[
                    "Metadata Evidence"
                ],
            )

            sample_metadata["Group"] = group
            sample_metadata[
                "Group Confidence"
            ] = group_confidence
            sample_metadata[
                "Group Evidence"
            ] = group_evidence

            (
                sex,
                sex_confidence,
                sex_evidence,
            ) = classify_sex(
                experiment_title,
                sample_metadata[
                    "Metadata Evidence"
                ],
            )

            sample_metadata["Sex"] = sex
            sample_metadata[
                "Sex Confidence"
            ] = sex_confidence
            sample_metadata[
                "Sex Evidence"
            ] = sex_evidence

            title = (
                sample_metadata["Study Title"]
                or experiment_title
                or clean(entry.get("title"))
            )

        else:
            sample_metadata = empty_sample_metadata()
            title = clean(
                entry.get("title")
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
                **sample_metadata,
                "Title":
                    title,
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


def load_validated_legacy_pairs(
    registry_path: Path | None,
) -> set[tuple[str, str, str]]:
    if registry_path is None or not registry_path.exists():
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
    legacy_registry: Path | None = None,
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
            "Experiment Accession",
            "Search Route",
        ]
    )

    validated_legacy_pairs = (
        load_validated_legacy_pairs(
            legacy_registry
        )
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

    strong_target_evidence = (
        (ok["Group"] == "Experiment")
        & (ok["Group Confidence"] == "High")
        & (
            ok["Group Evidence"]
            .fillna("")
            .astype(str)
            .str.contains(
                r"Target-specific|Target gene .* present",
                case=False,
                regex=True,
            )
        )
    )

    validated_study_keys = set(
        zip(
            ok.loc[
                strong_target_evidence,
                "Gene",
            ].astype(str).str.casefold(),
            ok.loc[
                strong_target_evidence,
                "Technique",
            ].astype(str),
            ok.loc[
                strong_target_evidence,
                "Database",
            ].astype(str),
            ok.loc[
                strong_target_evidence,
                "Accession",
            ].astype(str),
        )
    )

    def has_strong_study_validation(row):
        key = (
            clean(row["Gene"]).casefold(),
            clean(row["Technique"]),
            clean(row["Database"]),
            clean(row["Accession"]),
        )
        return key in validated_study_keys

    ok["Strong Study Validation"] = ok.apply(
        lambda row: (
            "Yes"
            if has_strong_study_validation(row)
            else "No"
        ),
        axis=1,
    )

    rejected_risky = (
        (ok["Search Route"] == "Legacy symbol fallback")
        & (ok["Risky Symbol"] == "Yes")
        & (ok["Validated Legacy"] != "Yes")
        & (ok["Strong Study Validation"] != "Yes")
    )

    risky_fallback = ok[
        rejected_risky
    ].copy()

    accepted = ok[
        ~rejected_risky
    ].copy()

    accepted = accepted.drop_duplicates(
        subset=[
            "Gene",
            "Technique",
            "Database",
            "Accession",
            "Experiment Accession",
        ]
    )

    summary = (
        accepted.groupby(
            [
                "Gene",
                "Technique",
                "Database",
            ]
        )["Accession"]
        .nunique()
        .reset_index(
            name="Dataset Count"
        )
    )

    sample_metadata = accepted[
        (
            accepted["Experiment Accession"]
            .fillna("")
            .astype(str)
            .str.len()
            > 0
        )
        | (
            accepted["Run Accessions"]
            .fillna("")
            .astype(str)
            .str.len()
            > 0
        )
        | (
            accepted["BioSample Accessions"]
            .fillna("")
            .astype(str)
            .str.len()
            > 0
        )
    ].copy()

    control_experiment = sample_metadata[
        sample_metadata["Group"].isin(
            ["Control", "Experiment"]
        )
    ].copy()

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

        sample_metadata.to_excel(
            writer,
            sheet_name="Sample_Metadata",
            index=False,
        )

        control_experiment.to_excel(
            writer,
            sheet_name="Control_Experiment",
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

    parser.add_argument(
        "--legacy-registry",
        type=Path,
        help=(
            "Optional validated legacy dataset registry used to promote "
            "verified historical gene-accession associations."
        ),
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

        gene_rows = collect_gene(
            gene=gene,
            gene_set=args.gene_set,
        )

        ok_gene_rows = [
            row
            for row in gene_rows
            if row.get("Status") == "OK"
        ]

        if not ok_gene_rows:
            print(
                f"  No hits returned for {gene}; retrying once...",
                flush=True,
            )

            time.sleep(2.0)

            retry_rows = collect_gene(
                gene=gene,
                gene_set=args.gene_set,
            )

            retry_ok_rows = [
                row
                for row in retry_rows
                if row.get("Status") == "OK"
            ]

            if retry_ok_rows:
                print(
                    f"  Retry recovered {len(retry_ok_rows)} rows for {gene}",
                    flush=True,
                )
                gene_rows = retry_rows

        rows.extend(gene_rows)

    dataframe = pd.DataFrame(
        rows
    )

    write_workbook(
        dataframe,
        args.output,
        legacy_registry=args.legacy_registry,
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
