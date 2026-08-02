"""NCBI BioSample search client."""

from __future__ import annotations

import re
import urllib.parse
from typing import Any

from dataset_finder.clients.ncbi_entrez import NCBIEntrezClient
from dataset_finder.models import DatasetRecord


class NCBIBioSampleClient:
    """Search and normalize NCBI BioSample records."""

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
        """Search BioSample for biological samples."""
        term = (
            f'("{query}"[All Fields]) AND '
            f'"{species}"[Organism]'
        )

        identifiers = self.entrez.search_ids(
            database="biosample",
            term=term,
            max_results=max_results,
        )

        summaries = self.entrez.summaries(
            database="biosample",
            identifiers=identifiers,
        )

        return [
            self._normalize_summary(summary)
            for summary in summaries
        ]

    @classmethod
    def _normalize_summary(
        cls,
        summary: dict[str, Any],
    ) -> DatasetRecord:
        """Convert one BioSample summary to a dataset record."""
        uid = str(summary.get("uid", "")).strip()

        accession = cls._first_value(
            summary,
            "accession",
            "accessionversion",
            "sample_acc",
            "sampleacc",
        )

        text = " ".join(
            str(value)
            for value in summary.values()
            if value is not None
        )

        if not accession:
            match = re.search(
                r"\bSAM[NED]\d+\b",
                text,
                flags=re.IGNORECASE,
            )
            accession = (
                match.group(0).upper()
                if match
                else uid
            )

        title = cls._first_value(
            summary,
            "title",
            "sample_name",
            "samplename",
        )

        organism = cls._first_value(
            summary,
            "organism",
            "organism_name",
            "organismname",
        )

        publication_date = cls._first_value(
            summary,
            "publication_date",
            "publicationdate",
            "submission_date",
            "submissiondate",
            "createdate",
            "updatedate",
        )

        attributes = cls._extract_attributes(summary)

        description_parts = [
            title,
            organism,
            *(
                f"{key}: {value}"
                for key, value in attributes.items()
                if value
            ),
        ]

        project_accession = cls._extract_accession(
            text,
            prefixes=("PRJNA", "PRJEB", "PRJDB"),
        )

        experiment_accessions = tuple(
            cls._extract_all_accessions(
                text,
                prefixes=("SRX", "ERX", "DRX"),
            )
        )

        run_accessions = tuple(
            cls._extract_all_accessions(
                text,
                prefixes=("SRR", "ERR", "DRR"),
            )
        )

        return DatasetRecord(
            uid=uid,
            accession=accession,
            title=title or accession,
            organism=organism,
            study_type="BioSample",
            sample_count=1,
            publication_date=publication_date,
            url=(
                "https://www.ncbi.nlm.nih.gov/biosample/"
                f"?term={urllib.parse.quote(accession)}"
            ),
            database="BioSample",
            project_accession=project_accession,
            experiment_accessions=experiment_accessions,
            sample_accessions=run_accessions,
            biosample_accessions=(
                (accession,)
                if accession
                else ()
            ),
            description=" | ".join(
                part
                for part in description_parts
                if part
            ),
            tissue=attributes.get("tissue", ""),
            cell_type=attributes.get(
                "cell type",
                attributes.get("cell_type", ""),
            ),
            developmental_stage=attributes.get(
                "developmental stage",
                attributes.get("developmental_stage", ""),
            ),
            sex=attributes.get("sex", ""),
            genotype=attributes.get("genotype", ""),
            strain=attributes.get(
                "strain",
                attributes.get("strain name", ""),
            ),
            treatment=attributes.get(
                "treatment",
                attributes.get("treatment protocol", ""),
            ),
            disease=attributes.get(
                "disease",
                attributes.get("disease state", ""),
            ),
            time_point=attributes.get(
                "time point",
                attributes.get("time_point", ""),
            ),
            evidence_text=" ".join(
                part
                for part in (
                    title,
                    organism,
                    text,
                )
                if part
            ),
            raw_metadata=summary,
        )

    @staticmethod
    def _first_value(
        summary: dict[str, Any],
        *keys: str,
    ) -> str:
        """Return the first non-empty summary value."""
        for key in keys:
            value = summary.get(key)

            if value is not None and str(value).strip():
                return str(value).strip()

        return ""

    @classmethod
    def _extract_attributes(
        cls,
        summary: dict[str, Any],
    ) -> dict[str, str]:
        """Extract BioSample attributes from variable summary formats."""
        attributes: dict[str, str] = {}

        for key, value in summary.items():
            normalized_key = str(key).strip().casefold()

            if isinstance(value, dict):
                for nested_key, nested_value in value.items():
                    cls._store_attribute(
                        attributes,
                        str(nested_key),
                        nested_value,
                    )

            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        name = (
                            item.get("name")
                            or item.get("attribute_name")
                            or item.get("harmonized_name")
                            or item.get("display_name")
                        )
                        item_value = (
                            item.get("value")
                            or item.get("attribute_value")
                        )

                        if name and item_value is not None:
                            cls._store_attribute(
                                attributes,
                                str(name),
                                item_value,
                            )
                    else:
                        cls._parse_attribute_text(
                            attributes,
                            str(item),
                        )

            elif (
                "attribute" in normalized_key
                or normalized_key
                in {
                    "tissue",
                    "sex",
                    "strain",
                    "genotype",
                    "disease",
                    "cell_type",
                    "developmental_stage",
                    "treatment",
                    "time_point",
                }
            ):
                if normalized_key in {
                    "tissue",
                    "sex",
                    "strain",
                    "genotype",
                    "disease",
                    "cell_type",
                    "developmental_stage",
                    "treatment",
                    "time_point",
                }:
                    cls._store_attribute(
                        attributes,
                        normalized_key,
                        value,
                    )
                else:
                    cls._parse_attribute_text(
                        attributes,
                        str(value),
                    )

        return attributes

    @staticmethod
    def _store_attribute(
        attributes: dict[str, str],
        name: str,
        value: object,
    ) -> None:
        """Store one normalized BioSample attribute."""
        key = name.strip().casefold().replace("_", " ")
        normalized_value = str(value).strip()

        if key and normalized_value and key not in attributes:
            attributes[key] = normalized_value

    @classmethod
    def _parse_attribute_text(
        cls,
        attributes: dict[str, str],
        text: str,
    ) -> None:
        """Parse simple key-value BioSample attribute text."""
        for fragment in re.split(r"[;\n|]+", text):
            key, separator, value = fragment.partition("=")

            if not separator:
                key, separator, value = fragment.partition(":")

            if separator:
                cls._store_attribute(
                    attributes,
                    key,
                    value,
                )

    @staticmethod
    def _extract_accession(
        text: str,
        *,
        prefixes: tuple[str, ...],
    ) -> str:
        """Extract the first matching archive accession."""
        values = NCBIBioSampleClient._extract_all_accessions(
            text,
            prefixes=prefixes,
        )
        return values[0] if values else ""

    @staticmethod
    def _extract_all_accessions(
        text: str,
        *,
        prefixes: tuple[str, ...],
    ) -> list[str]:
        """Extract unique accessions while preserving order."""
        prefix_pattern = "|".join(
            re.escape(prefix)
            for prefix in prefixes
        )

        matches = re.findall(
            rf"\b(?:{prefix_pattern})\d+\b",
            text,
            flags=re.IGNORECASE,
        )

        result: list[str] = []
        seen: set[str] = set()

        for match in matches:
            normalized = match.upper()

            if normalized in seen:
                continue

            seen.add(normalized)
            result.append(normalized)

        return result
