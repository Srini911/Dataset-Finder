"""NCBI Sequence Read Archive search client."""

from __future__ import annotations

import re
import urllib.parse
from typing import Any

from dataset_finder.clients.ncbi_entrez import (
    NCBIEntrezClient,
)
from dataset_finder.models import DatasetRecord


class NCBISRAClient:
    """Search NCBI SRA and normalize study records."""

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
        """Search SRA for one gene or regulator."""
        term = (
            f'("{query}"[All Fields]) AND '
            f'"{species}"[Organism]'
        )

        identifiers = self.entrez.search_ids(
            database="sra",
            term=term,
            max_results=max_results,
        )
        summaries = self.entrez.summaries(
            database="sra",
            identifiers=identifiers,
        )

        records = [
            self._normalize_summary(summary)
            for summary in summaries
        ]

        return self._deduplicate_studies(records)

    @staticmethod
    def _deduplicate_studies(
        records: list[DatasetRecord],
    ) -> list[DatasetRecord]:
        """Remove repeated experiments belonging to the same SRA study."""
        deduplicated: list[DatasetRecord] = []
        seen: set[str] = set()

        for record in records:
            identity = (
                record.project_accession
                or record.accession
                or record.uid
            ).strip().upper()

            if identity in seen:
                continue

            seen.add(identity)
            deduplicated.append(record)

        return deduplicated

    @classmethod
    def _normalize_summary(
        cls,
        summary: dict[str, Any],
    ) -> DatasetRecord:
        """Convert an SRA document summary into a dataset record."""
        uid = str(summary.get("uid", ""))
        title = str(summary.get("title", "")).strip()
        experiment_xml = str(
            summary.get("expxml", "")
        )
        runs_xml = str(summary.get("runs", ""))

        study_accession = cls._first_accession(
            experiment_xml,
            prefixes=("SRP", "ERP", "DRP"),
        )
        experiment_accession = cls._first_accession(
            experiment_xml,
            prefixes=("SRX", "ERX", "DRX"),
        )
        run_accessions = tuple(
            cls._all_accessions(
                runs_xml,
                prefixes=("SRR", "ERR", "DRR"),
            )
        )

        accession = (
            study_accession
            or experiment_accession
            or uid
        )

        organism = cls._extract_attribute(
            experiment_xml,
            "ScientificName",
        )
        sample_count = (
            len(run_accessions)
            if run_accessions
            else None
        )

        return DatasetRecord(
            uid=uid,
            accession=accession,
            title=title,
            organism=organism,
            study_type="Sequence Read Archive",
            sample_count=sample_count,
            publication_date=str(
                summary.get(
                    "createdate",
                    summary.get("updatedate", ""),
                )
            ),
            url=(
                "https://www.ncbi.nlm.nih.gov/sra/"
                f"?term={urllib.parse.quote(accession)}"
            ),
            database="SRA",
            project_accession=study_accession,
            sample_accessions=run_accessions,
            description=cls._strip_xml(experiment_xml),
            evidence_text=cls._strip_xml(
                f"{experiment_xml} {runs_xml}"
            ),
            raw_metadata=summary,
        )

    @staticmethod
    def _extract_attribute(
        text: str,
        attribute: str,
    ) -> str:
        """Extract a quoted XML-like attribute."""
        match = re.search(
            rf'{re.escape(attribute)}="([^"]+)"',
            text,
            flags=re.IGNORECASE,
        )
        return match.group(1).strip() if match else ""

    @staticmethod
    def _strip_xml(text: str) -> str:
        """Convert summary XML fragments into readable text."""
        without_tags = re.sub(r"<[^>]+>", " ", text)
        return re.sub(r"\s+", " ", without_tags).strip()

    @staticmethod
    def _first_accession(
        text: str,
        *,
        prefixes: tuple[str, ...],
    ) -> str:
        """Return the first accession with an accepted prefix."""
        accessions = NCBISRAClient._all_accessions(
            text,
            prefixes=prefixes,
        )
        return accessions[0] if accessions else ""

    @staticmethod
    def _all_accessions(
        text: str,
        *,
        prefixes: tuple[str, ...],
    ) -> list[str]:
        """Extract unique NCBI/ENA-style accessions."""
        prefix_pattern = "|".join(
            re.escape(prefix)
            for prefix in prefixes
        )
        matches = re.findall(
            rf"\b(?:{prefix_pattern})\d+\b",
            text,
            flags=re.IGNORECASE,
        )

        unique: list[str] = []
        seen: set[str] = set()

        for match in matches:
            normalized = match.upper()

            if normalized in seen:
                continue

            seen.add(normalized)
            unique.append(normalized)

        return unique
