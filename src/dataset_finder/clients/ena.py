"""European Nucleotide Archive search client."""

from __future__ import annotations

import csv
import io
import urllib.parse

import requests

from dataset_finder.models import DatasetRecord

ENA_TEXT_SEARCH_URL = (
    "https://www.ebi.ac.uk/ena/browser/api/tsv/textsearch"
)
ENA_BROWSER_URL = "https://www.ebi.ac.uk/ena/browser/view"


class ENAClientError(RuntimeError):
    """Raised when an ENA request cannot be completed."""


class ENAClient:
    """Search ENA studies and normalize study-level records."""

    def __init__(
        self,
        *,
        timeout: float = 30.0,
        session: requests.Session | None = None,
    ) -> None:
        self.timeout = timeout
        self.session = session or requests.Session()

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int = 20,
    ) -> list[DatasetRecord]:
        """Search ENA studies using free-text search."""
        species = species.strip()
        query = query.strip()

        if not species:
            raise ValueError("Species cannot be empty.")

        if not query:
            raise ValueError("Query cannot be empty.")

        if max_results < 1:
            raise ValueError(
                "Maximum results must be greater than zero."
            )

        search_text = f"{query} {species}".strip()

        try:
            response = self.session.get(
                ENA_TEXT_SEARCH_URL,
                params={
                    "result": "study",
                    "query": search_text,
                    "limit": max_results,
                },
                timeout=self.timeout,
            )
            response.raise_for_status()
        except requests.RequestException as exc:
            raise ENAClientError(
                f"ENA request failed: {exc}"
            ) from exc

        return self._parse_study_tsv(
            response.text,
            species=species,
            max_results=max_results,
        )

    @classmethod
    def _parse_study_tsv(
        cls,
        text: str,
        *,
        species: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        """Parse ENA study text-search TSV output."""
        if not text.strip():
            return []

        reader = csv.DictReader(
            io.StringIO(text),
            delimiter="\t",
        )

        records: list[DatasetRecord] = []

        for row in reader:
            normalized = {
                str(key or "").strip().casefold(): str(
                    value or ""
                ).strip()
                for key, value in row.items()
            }

            accession = (
                normalized.get("accession", "")
                or normalized.get(
                    "study_accession",
                    "",
                )
            )

            description = (
                normalized.get("description", "")
                or normalized.get("study_description", "")
                or normalized.get("study_title", "")
            )

            if not accession:
                continue

            records.append(
                DatasetRecord(
                    uid=accession,
                    accession=accession,
                    title=description or accession,
                    organism=species,
                    study_type="ENA Study",
                    sample_count=None,
                    publication_date=(
                        normalized.get("first_public", "")
                    ),
                    url=(
                        f"{ENA_BROWSER_URL}/"
                        f"{urllib.parse.quote(accession)}"
                    ),
                    database="ENA",
                    project_accession=accession,
                    description=description,
                    evidence_text=" ".join(
                        value
                        for value in (
                            accession,
                            description,
                            species,
                        )
                        if value
                    ),
                    raw_metadata=dict(row),
                )
            )

            if len(records) >= max_results:
                break

        return records
