"""NCBI BioProject search client."""

from __future__ import annotations

import re
import urllib.parse
from typing import Any

from dataset_finder.clients.ncbi_entrez import (
    NCBIEntrezClient,
)
from dataset_finder.models import DatasetRecord


class NCBIBioProjectClient:
    """Search and normalize NCBI BioProject records."""

    def __init__(
        self,
        *,
        entrez_client: NCBIEntrezClient | None = None,
    ) -> None:
        self.entrez = entrez_client or NCBIEntrezClient()

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int = 20,
    ) -> list[DatasetRecord]:
        """Search BioProject for a gene or regulator."""
        term = (
            f'("{query}"[All Fields]) AND '
            f'"{species}"[Organism]'
        )

        identifiers = self.entrez.search_ids(
            database="bioproject",
            term=term,
            max_results=max_results,
        )
        summaries = self.entrez.summaries(
            database="bioproject",
            identifiers=identifiers,
        )

        return [
            self._normalize_summary(summary)
            for summary in summaries
        ]

    @staticmethod
    def _normalize_summary(
        summary: dict[str, Any],
    ) -> DatasetRecord:
        """Convert a BioProject summary to a dataset record."""
        uid = str(summary.get("uid", ""))

        accession = str(
            summary.get(
                "project_acc",
                summary.get("projectacc", ""),
            )
        ).strip()

        if not accession:
            text = " ".join(
                str(value)
                for value in summary.values()
            )
            match = re.search(
                r"\bPRJ(?:NA|EB|DB)\d+\b",
                text,
                flags=re.IGNORECASE,
            )
            accession = (
                match.group(0).upper()
                if match
                else uid
            )

        title = str(
            summary.get(
                "project_title",
                summary.get("projecttitle", ""),
            )
        ).strip()

        description = str(
            summary.get(
                "project_description",
                summary.get("projectdescription", ""),
            )
        ).strip()

        organism = str(
            summary.get(
                "organism_name",
                summary.get("organismname", ""),
            )
        ).strip()

        publication_date = str(
            summary.get(
                "registration_date",
                summary.get("registrationdate", ""),
            )
        ).strip()

        return DatasetRecord(
            uid=uid,
            accession=accession,
            title=title,
            organism=organism,
            study_type="BioProject",
            sample_count=None,
            publication_date=publication_date,
            url=(
                "https://www.ncbi.nlm.nih.gov/bioproject/"
                f"?term={urllib.parse.quote(accession)}"
            ),
            database="BioProject",
            project_accession=accession,
            description=description,
            evidence_text=" ".join(
                value
                for value in (
                    title,
                    description,
                    organism,
                )
                if value
            ),
            raw_metadata=summary,
        )
