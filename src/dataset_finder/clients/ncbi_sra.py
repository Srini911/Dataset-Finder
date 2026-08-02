"""NCBI Sequence Read Archive search client."""

from __future__ import annotations

import re
import urllib.parse
from collections import defaultdict
from dataclasses import replace
from typing import Any

from dataset_finder.clients.ncbi_entrez import NCBIEntrezClient
from dataset_finder.models import DatasetRecord


class NCBISRAClient:
    """Search NCBI SRA and normalize study-level records."""

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

        experiment_records = [
            self._normalize_summary(summary)
            for summary in summaries
        ]

        return self._aggregate_studies(experiment_records)

    @classmethod
    def _normalize_summary(
        cls,
        summary: dict[str, Any],
    ) -> DatasetRecord:
        """Convert one SRA experiment summary into a normalized record."""
        uid = str(summary.get("uid", ""))
        experiment_xml = str(summary.get("expxml", ""))
        runs_xml = str(summary.get("runs", ""))

        study_accession = cls._attribute_from_tag(
            experiment_xml,
            tag="Study",
            attribute="acc",
        )
        study_title = cls._attribute_from_tag(
            experiment_xml,
            tag="Study",
            attribute="name",
        )
        summary_title = cls._text_from_tag(
            experiment_xml,
            tag="Title",
        )
        experiment_accession = cls._attribute_from_tag(
            experiment_xml,
            tag="Experiment",
            attribute="acc",
        )
        sample_accession = cls._attribute_from_tag(
            experiment_xml,
            tag="Sample",
            attribute="acc",
        )
        organism = cls._extract_attribute(
            experiment_xml,
            "ScientificName",
        )
        bioproject = cls._text_from_tag(
            experiment_xml,
            tag="Bioproject",
        )
        biosample = cls._text_from_tag(
            experiment_xml,
            tag="Biosample",
        )
        library_strategy = cls._text_from_tag(
            experiment_xml,
            tag="LIBRARY_STRATEGY",
        )
        library_source = cls._text_from_tag(
            experiment_xml,
            tag="LIBRARY_SOURCE",
        )
        library_selection = cls._text_from_tag(
            experiment_xml,
            tag="LIBRARY_SELECTION",
        )
        platform = cls._extract_attribute(
            experiment_xml,
            "instrument_model",
        )
        layout = cls._library_layout(experiment_xml)

        run_accessions = tuple(
            cls._all_accessions(
                runs_xml,
                prefixes=("SRR", "ERR", "DRR"),
            )
        )

        title = study_title or summary_title
        accession = (
            study_accession
            or experiment_accession
            or uid
        )

        description_parts = [
            summary_title,
            library_strategy,
            library_source,
            library_selection,
            layout,
            platform,
        ]

        return DatasetRecord(
            uid=uid,
            accession=accession,
            title=title,
            organism=organism,
            study_type=library_strategy or "Sequence Read Archive",
            sample_count=len(run_accessions) or None,
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
            project_accession=bioproject or study_accession,
            experiment_accessions=(
                (experiment_accession,)
                if experiment_accession
                else ()
            ),
            sample_accessions=run_accessions,
            biosample_accessions=tuple(
                value
                for value in (
                    sample_accession,
                    biosample,
                )
                if value
            ),
            library_strategy=library_strategy,
            library_source=library_source,
            library_selection=library_selection,
            library_layout=layout,
            platform=platform,
            description=" | ".join(
                value
                for value in description_parts
                if value
            ),
            evidence_text=cls._strip_xml(
                f"{experiment_xml} {runs_xml}"
            ),
            raw_metadata=summary,
        )

    @classmethod
    def _aggregate_studies(
        cls,
        records: list[DatasetRecord],
    ) -> list[DatasetRecord]:
        """Merge SRA experiments and runs belonging to one study."""
        grouped: dict[str, list[DatasetRecord]] = defaultdict(list)

        for record in records:
            identity = (
                record.accession
                or record.project_accession
                or record.uid
            ).strip().upper()

            grouped[identity].append(record)

        aggregated: list[DatasetRecord] = []

        for group in grouped.values():
            first = group[0]

            experiment_accessions = cls._unique(
                value
                for record in group
                for value in record.experiment_accessions
            )
            run_accessions = cls._unique(
                value
                for record in group
                for value in record.sample_accessions
            )
            biosample_accessions = cls._unique(
                value
                for record in group
                for value in record.biosample_accessions
            )

            aggregated.append(
                replace(
                    first,
                    experiment_accessions=experiment_accessions,
                    sample_accessions=run_accessions,
                    biosample_accessions=biosample_accessions,
                    sample_count=len(run_accessions) or None,
                    evidence_text=" ".join(
                        record.evidence_text
                        for record in group
                        if record.evidence_text
                    ),
                    raw_metadata={
                        "experiments": [
                            record.raw_metadata
                            for record in group
                        ]
                    },
                )
            )

        return aggregated

    @staticmethod
    def _unique(values) -> tuple[str, ...]:
        """Deduplicate strings while preserving order."""
        unique: list[str] = []
        seen: set[str] = set()

        for value in values:
            normalized = str(value).strip()

            if not normalized:
                continue

            identity = normalized.upper()

            if identity in seen:
                continue

            seen.add(identity)
            unique.append(normalized)

        return tuple(unique)

    @staticmethod
    def _attribute_from_tag(
        text: str,
        *,
        tag: str,
        attribute: str,
    ) -> str:
        """Extract an attribute from a specific XML-like tag."""
        match = re.search(
            rf"<{re.escape(tag)}\b[^>]*"
            rf'{re.escape(attribute)}="([^"]+)"',
            text,
            flags=re.IGNORECASE,
        )
        return match.group(1).strip() if match else ""

    @staticmethod
    def _text_from_tag(
        text: str,
        *,
        tag: str,
    ) -> str:
        """Extract text enclosed by one XML-like tag."""
        match = re.search(
            rf"<{re.escape(tag)}>(.*?)</{re.escape(tag)}>",
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )
        return (
            re.sub(r"\s+", " ", match.group(1)).strip()
            if match
            else ""
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
    def _library_layout(text: str) -> str:
        """Extract paired- or single-end library layout."""
        if re.search(r"<PAIRED\b", text, flags=re.IGNORECASE):
            return "PAIRED"

        if re.search(r"<SINGLE\b", text, flags=re.IGNORECASE):
            return "SINGLE"

        return ""

    @staticmethod
    def _strip_xml(text: str) -> str:
        """Convert XML fragments into readable text."""
        without_tags = re.sub(r"<[^>]+>", " ", text)
        return re.sub(r"\s+", " ", without_tags).strip()

    @staticmethod
    def _all_accessions(
        text: str,
        *,
        prefixes: tuple[str, ...],
    ) -> list[str]:
        """Extract unique archive accessions."""
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
