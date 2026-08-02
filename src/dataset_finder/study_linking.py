"""Cross-database accession and publication linking."""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

from dataset_finder.models import DatasetRecord

ACCESSION_PATTERNS: dict[str, re.Pattern[str]] = {
    "geo": re.compile(
        r"\b(?:GSE|GSM|GPL|GDS)\d+\b",
        flags=re.IGNORECASE,
    ),
    "study": re.compile(
        r"\b(?:SRP|ERP|DRP)\d+\b",
        flags=re.IGNORECASE,
    ),
    "experiment": re.compile(
        r"\b(?:SRX|ERX|DRX)\d+\b",
        flags=re.IGNORECASE,
    ),
    "run": re.compile(
        r"\b(?:SRR|ERR|DRR)\d+\b",
        flags=re.IGNORECASE,
    ),
    "bioproject": re.compile(
        r"\b(?:PRJNA|PRJEB|PRJDB)\d+\b",
        flags=re.IGNORECASE,
    ),
    "biosample": re.compile(
        r"\b(?:SAMN|SAME|SAMD)\d+\b",
        flags=re.IGNORECASE,
    ),
    "encode": re.compile(
        r"\bENCSR[A-Z0-9]+\b",
        flags=re.IGNORECASE,
    ),
    "pmid": re.compile(
        r"\bPMID\s*[:#]?\s*(\d{5,10})\b",
        flags=re.IGNORECASE,
    ),
}

DOI_PATTERN = re.compile(
    r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b",
    flags=re.IGNORECASE,
)


@dataclass(frozen=True, slots=True)
class StudyLinks:
    """Normalized identifiers related to one dataset or publication."""

    pubmed_ids: tuple[str, ...] = ()
    dois: tuple[str, ...] = ()
    related_geo_accessions: tuple[str, ...] = ()
    related_study_accessions: tuple[str, ...] = ()
    related_experiment_accessions: tuple[str, ...] = ()
    related_run_accessions: tuple[str, ...] = ()
    related_bioproject_accessions: tuple[str, ...] = ()
    related_biosample_accessions: tuple[str, ...] = ()
    related_encode_accessions: tuple[str, ...] = ()

    @property
    def all_accessions(self) -> tuple[str, ...]:
        """Return every related archive accession."""
        return _unique(
            (
                *self.related_geo_accessions,
                *self.related_study_accessions,
                *self.related_experiment_accessions,
                *self.related_run_accessions,
                *self.related_bioproject_accessions,
                *self.related_biosample_accessions,
                *self.related_encode_accessions,
            )
        )

    @property
    def study_level_accessions(self) -> tuple[str, ...]:
        """Return compact study- and project-level accessions."""
        geo_studies = (
            accession
            for accession in self.related_geo_accessions
            if accession.upper().startswith(
                ("GSE", "GDS")
            )
        )

        return _unique(
            (
                *geo_studies,
                *self.related_study_accessions,
                *self.related_bioproject_accessions,
                *self.related_encode_accessions,
            )
        )


def _flatten(value: Any) -> list[str]:
    """Flatten nested metadata into searchable text."""
    if value is None:
        return []

    if isinstance(value, dict):
        values: list[str] = []

        for key, item in value.items():
            values.extend(_flatten(key))
            values.extend(_flatten(item))

        return values

    if isinstance(value, (list, tuple, set)):
        values = []

        for item in value:
            values.extend(_flatten(item))

        return values

    text = str(value).strip()
    return [text] if text else []


def _unique(values: Iterable[str]) -> tuple[str, ...]:
    """Deduplicate identifiers while preserving source order."""
    result: list[str] = []
    seen: set[str] = set()

    for value in values:
        normalized = str(value).strip()

        if not normalized:
            continue

        identity = normalized.casefold()

        if identity in seen:
            continue

        seen.add(identity)
        result.append(normalized)

    return tuple(result)


def _matches(
    text: str,
    pattern: re.Pattern[str],
) -> tuple[str, ...]:
    """Return normalized regular-expression matches."""
    results: list[str] = []

    for match in pattern.finditer(text):
        value = (
            match.group(1)
            if match.lastindex
            else match.group(0)
        )
        value = value.strip().rstrip(".,;:)")

        if pattern is DOI_PATTERN:
            results.append(value.lower())
        else:
            results.append(value.upper())

    return _unique(results)


def extract_study_links(*values: object) -> StudyLinks:
    """Extract archive accessions, PubMed IDs, and DOIs."""
    text = " ".join(
        fragment
        for value in values
        for fragment in _flatten(value)
    )

    explicit_pmids = _matches(
        text,
        ACCESSION_PATTERNS["pmid"],
    )

    return StudyLinks(
        pubmed_ids=explicit_pmids,
        dois=_matches(text, DOI_PATTERN),
        related_geo_accessions=_matches(
            text,
            ACCESSION_PATTERNS["geo"],
        ),
        related_study_accessions=_matches(
            text,
            ACCESSION_PATTERNS["study"],
        ),
        related_experiment_accessions=_matches(
            text,
            ACCESSION_PATTERNS["experiment"],
        ),
        related_run_accessions=_matches(
            text,
            ACCESSION_PATTERNS["run"],
        ),
        related_bioproject_accessions=_matches(
            text,
            ACCESSION_PATTERNS["bioproject"],
        ),
        related_biosample_accessions=_matches(
            text,
            ACCESSION_PATTERNS["biosample"],
        ),
        related_encode_accessions=_matches(
            text,
            ACCESSION_PATTERNS["encode"],
        ),
    )


def extract_links_from_record(
    record: DatasetRecord,
) -> StudyLinks:
    """Extract links from all normalized and raw record metadata."""
    return extract_study_links(
        record.accession,
        record.project_accession,
        record.experiment_accessions,
        record.sample_accessions,
        record.biosample_accessions,
        record.title,
        record.description,
        record.publication,
        record.evidence_text,
        record.url,
        record.raw_metadata,
    )


def _normalized_title(value: str) -> str:
    """Normalize a title for conservative cross-record grouping."""
    value = value.casefold()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def link_related_records(
    records: Iterable[DatasetRecord],
) -> tuple[DatasetRecord, ...]:
    """Share accessions and publication identifiers across matching studies."""
    from dataclasses import replace

    record_list = list(records)
    groups: dict[str, list[int]] = {}

    for index, record in enumerate(record_list):
        normalized_title = _normalized_title(record.title)

        if len(normalized_title) < 20:
            continue

        groups.setdefault(normalized_title, []).append(index)

    linked_records = list(record_list)

    for indexes in groups.values():
        if len(indexes) < 2:
            continue

        group_records = [
            record_list[index]
            for index in indexes
        ]

        databases = {
            record.database.casefold()
            for record in group_records
            if record.database
        }

        if len(databases) < 2:
            continue

        links = [
            extract_links_from_record(record)
            for record in group_records
        ]

        pubmed_ids = _unique(
            value
            for record, record_links in zip(
                group_records,
                links,
                strict=True,
            )
            for value in (
                *record.pubmed_ids,
                *record_links.pubmed_ids,
            )
        )

        dois = _unique(
            value
            for record, record_links in zip(
                group_records,
                links,
                strict=True,
            )
            for value in (
                *record.dois,
                *record_links.dois,
            )
        )

        related_geo = _unique(
            value
            for record, record_links in zip(
                group_records,
                links,
                strict=True,
            )
            for value in (
                *record.related_geo_accessions,
                *record_links.related_geo_accessions,
            )
        )

        related_studies = _unique(
            value
            for record, record_links in zip(
                group_records,
                links,
                strict=True,
            )
            for value in (
                *record.related_study_accessions,
                *record_links.related_study_accessions,
            )
        )

        related_bioprojects = _unique(
            value
            for record, record_links in zip(
                group_records,
                links,
                strict=True,
            )
            for value in (
                *record.related_bioproject_accessions,
                *record_links.related_bioproject_accessions,
            )
        )

        related_biosamples = _unique(
            value
            for record, record_links in zip(
                group_records,
                links,
                strict=True,
            )
            for value in (
                *record.related_biosample_accessions,
                *record_links.related_biosample_accessions,
            )
        )

        related_accessions = _unique(
            value
            for record, record_links in zip(
                group_records,
                links,
                strict=True,
            )
            for value in (
                record.accession,
                record.project_accession,
                *record.related_accessions,
                *record_links.study_level_accessions,
                *(
                    f"PMID{pmid}"
                    for pmid in (
                        *record.pubmed_ids,
                        *record_links.pubmed_ids,
                    )
                ),
            )
        )

        for index in indexes:
            linked_records[index] = replace(
                linked_records[index],
                pubmed_ids=pubmed_ids,
                dois=dois,
                related_accessions=related_accessions,
                related_geo_accessions=related_geo,
                related_study_accessions=related_studies,
                related_bioproject_accessions=related_bioprojects,
                related_biosample_accessions=related_biosamples,
            )

    return tuple(linked_records)
