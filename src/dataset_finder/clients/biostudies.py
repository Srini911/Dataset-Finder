"""EMBL-EBI BioStudies and ArrayExpress search client."""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import requests

from dataset_finder.models import DatasetRecord

BIOSTUDIES_SEARCH_URL = (
    "https://www.ebi.ac.uk/biostudies/api/v1/search"
)


class BioStudiesClientError(RuntimeError):
    """Raised when a BioStudies request fails."""


class BioStudiesClient:
    """Search public ArrayExpress studies through BioStudies."""

    def __init__(
        self,
        *,
        timeout: float = 45.0,
        session: requests.Session | None = None,
    ) -> None:
        self.timeout = timeout
        self.session = session or requests.Session()
        self.session.headers.update(
            {
                "Accept": "application/json",
                "User-Agent": (
                    "Dataset-Finder/0.3.0 "
                    "(https://github.com/Srini911/Dataset-Finder)"
                ),
            }
        )

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int = 20,
    ) -> list[DatasetRecord]:
        """Search ArrayExpress studies in BioStudies."""
        search_query = f'"{species}" "{query}"'

        try:
            response = self.session.get(
                BIOSTUDIES_SEARCH_URL,
                params={
                    "query": search_query,
                    "collection": "ArrayExpress",
                    "pageSize": max_results,
                },
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
        except (
            requests.RequestException,
            ValueError,
        ) as exc:
            raise BioStudiesClientError(
                f"BioStudies search failed: {exc}"
            ) from exc

        hits = payload.get("hits", [])

        if not isinstance(hits, list):
            return []

        records: list[DatasetRecord] = []

        for hit in hits:
            if not isinstance(hit, dict):
                continue

            record = self._normalize_hit(
                hit,
                fallback_organism=species,
            )

            if record.accession:
                records.append(record)

        return self._deduplicate(records)

    @staticmethod
    def _normalize_hit(
        hit: dict[str, Any],
        *,
        fallback_organism: str,
    ) -> DatasetRecord:
        """Normalize one BioStudies search result."""
        accession = str(
            hit.get("accession", "")
        ).strip()

        title = str(
            hit.get("title", "")
        ).strip()

        content = str(
            hit.get("content", "")
        ).strip()

        author = str(
            hit.get("author", "")
        ).strip()

        release_date = str(
            hit.get("release_date", "")
        ).strip()

        files = hit.get("files")
        sample_count = (
            int(files)
            if isinstance(files, int)
            else None
        )

        evidence = " ".join(
            value
            for value in (
                title,
                content,
                author,
            )
            if value
        )

        return DatasetRecord(
            uid=accession,
            accession=accession,
            title=title,
            organism=fallback_organism,
            study_type="ArrayExpress / BioStudies",
            sample_count=sample_count,
            publication_date=release_date,
            url=(
                "https://www.ebi.ac.uk/biostudies/"
                f"arrayexpress/studies/{quote(accession)}"
            ),
            database="BioStudies",
            project_accession=accession,
            description=content,
            publication=author,
            evidence_text=evidence,
            raw_metadata=hit,
        )

    @staticmethod
    def _deduplicate(
        records: list[DatasetRecord],
    ) -> list[DatasetRecord]:
        """Deduplicate BioStudies records by accession."""
        deduplicated: list[DatasetRecord] = []
        seen: set[str] = set()

        for record in records:
            identity = record.accession.strip().upper()

            if not identity or identity in seen:
                continue

            seen.add(identity)
            deduplicated.append(record)

        return deduplicated
