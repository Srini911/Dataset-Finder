"""Technique-specific historical dataset discovery."""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, replace
from datetime import UTC, datetime

from dataset_finder.assay_classifier import classify_technique
from dataset_finder.clients.ncbi_entrez import NCBIEntrezClient
from dataset_finder.clients.ncbi_geo import NCBIGEOClient
from dataset_finder.clients.ncbi_sra import NCBISRAClient
from dataset_finder.models import DatasetRecord
from dataset_finder.search_profiles import (
    TECHNIQUE_SEARCH_PROFILES,
    TechniqueSearchProfile,
)


@dataclass(frozen=True, slots=True)
class HistoricalSearchStatus:
    """Result of one gene, technique, and database query."""

    gene_query: str
    technique: str
    database: str
    success: bool
    candidate_count: int
    error: str = ""


@dataclass(frozen=True, slots=True)
class HistoricalSearchResult:
    """Historical records and query-level search statuses."""

    records: tuple[DatasetRecord, ...]
    statuses: tuple[HistoricalSearchStatus, ...]


class HistoricalSearchService:
    """Search SRA and GEO using technique-specific historical queries."""

    def __init__(
        self,
        *,
        entrez_client: NCBIEntrezClient | None = None,
        start_year: int = 2005,
        end_year: int | None = None,
        page_size: int = 100,
        historical_cutoff_year: int = 2015,
    ) -> None:
        self.entrez = entrez_client or NCBIEntrezClient()
        self.start_year = start_year
        self.end_year = (
            end_year
            if end_year is not None
            else datetime.now(UTC).year
        )
        self.page_size = page_size
        self.historical_cutoff_year = historical_cutoff_year

        if self.start_year < 1900:
            raise ValueError(
                "Historical search start year must be 1900 or later."
            )

        if self.end_year < self.start_year:
            raise ValueError(
                "Historical search end year cannot precede start year."
            )

        if self.page_size < 1:
            raise ValueError(
                "Historical search page size must be greater than zero."
            )

    def search(
        self,
        *,
        species: str,
        gene_terms: Iterable[str],
        max_results_per_query: int = 100,
        profiles: tuple[
            TechniqueSearchProfile,
            ...,
        ] = TECHNIQUE_SEARCH_PROFILES,
    ) -> HistoricalSearchResult:
        """Search historical SRA and GEO records for gene-technique pairs."""
        species = species.strip()

        if not species:
            raise ValueError("Species cannot be empty.")

        if max_results_per_query < 1:
            raise ValueError(
                "Maximum results per historical query "
                "must be greater than zero."
            )

        normalized_gene_terms = self._unique_terms(gene_terms)

        if not normalized_gene_terms:
            raise ValueError(
                "At least one gene search term is required."
            )

        records: list[DatasetRecord] = []
        statuses: list[HistoricalSearchStatus] = []

        for gene_term in normalized_gene_terms:
            for profile in profiles:
                query = self._build_query(
                    species=species,
                    gene_term=gene_term,
                    profile=profile,
                )

                for database in ("sra", "gds"):
                    database_name = (
                        "SRA"
                        if database == "sra"
                        else "GEO"
                    )

                    try:
                        discovered = self._search_database(
                            database=database,
                            query=query,
                            gene_term=gene_term,
                            profile=profile,
                            max_results=max_results_per_query,
                        )
                    except Exception as exc:
                        statuses.append(
                            HistoricalSearchStatus(
                                gene_query=gene_term,
                                technique=profile.name,
                                database=database_name,
                                success=False,
                                candidate_count=0,
                                error=str(exc),
                            )
                        )
                        continue

                    records.extend(discovered)
                    statuses.append(
                        HistoricalSearchStatus(
                            gene_query=gene_term,
                            technique=profile.name,
                            database=database_name,
                            success=True,
                            candidate_count=len(discovered),
                        )
                    )

        return HistoricalSearchResult(
            records=tuple(self._deduplicate(records)),
            statuses=tuple(statuses),
        )

    def _search_database(
        self,
        *,
        database: str,
        query: str,
        gene_term: str,
        profile: TechniqueSearchProfile,
        max_results: int,
    ) -> list[DatasetRecord]:
        """Search and normalize one NCBI database."""
        identifiers = self.entrez.search_ids_paged(
            database=database,
            term=query,
            max_results=max_results,
            page_size=self.page_size,
            minimum_date=f"{self.start_year}/01/01",
            maximum_date=f"{self.end_year}/12/31",
            date_type="pdat",
        )

        summaries = self.entrez.summaries(
            database=database,
            identifiers=identifiers,
        )

        if database == "sra":
            normalized = [
                NCBISRAClient._normalize_summary(summary)
                for summary in summaries
            ]
            records = NCBISRAClient._aggregate_studies(
                normalized
            )
        else:
            records = [
                NCBIGEOClient._normalize_summary(summary)
                for summary in summaries
            ]

        return [
            self._annotate_record(
                record=record,
                database=database,
                query=query,
                gene_term=gene_term,
                profile=profile,
            )
            for record in records
        ]

    def _annotate_record(
        self,
        *,
        record: DatasetRecord,
        database: str,
        query: str,
        gene_term: str,
        profile: TechniqueSearchProfile,
    ) -> DatasetRecord:
        """Add search provenance and verified technique evidence."""
        database_name = (
            "SRA"
            if database == "sra"
            else "GEO"
        )

        technique, evidence, evidence_source = (
            self._verified_technique(record)
        )
        technique_match = self._technique_match(
            requested=profile.name,
            verified=technique,
        )

        study_year = self._extract_year(
            record.publication_date
        )

        return replace(
            record,
            database=record.database or database_name,
            technique=technique,
            technique_requested=profile.name,
            technique_match=technique_match,
            technique_search_term="; ".join(profile.terms),
            technique_evidence=evidence,
            technique_evidence_source=evidence_source,
            gene_query_used=gene_term,
            search_query_used=query,
            study_year=study_year,
            historical_study=(
                study_year is not None
                and study_year <= self.historical_cutoff_year
            ),
            raw_metadata={
                **record.raw_metadata,
                "historical_search": {
                    "gene_query": gene_term,
                    "technique_requested": profile.name,
                    "technique_terms": list(profile.terms),
                    "search_query": query,
                    "start_year": self.start_year,
                    "end_year": self.end_year,
                },
            },
        )

    @staticmethod
    def _verified_technique(
        record: DatasetRecord,
    ) -> tuple[str, str, str]:
        """Classify technique with specific assays taking precedence."""
        descriptive_values = (
            record.title,
            record.description,
            record.evidence_text,
        )
        structured_values = (
            record.library_strategy,
            record.study_type,
        )

        descriptive_technique = classify_technique(
            *descriptive_values
        )
        structured_technique = classify_technique(
            *structured_values
        )

        specific_techniques = {
            "CUT_RUN",
            "CUT_TAG",
            "eCLIP",
            "iCLIP",
            "PAR_CLIP",
            "HITS_CLIP",
            "CLIP",
        }

        if descriptive_technique in specific_techniques:
            evidence = " | ".join(
                value
                for value in descriptive_values
                if value
            )

            return (
                descriptive_technique,
                evidence,
                "Specific assay evidence in title or description",
            )

        if structured_technique != "Other_Assays":
            evidence = " | ".join(
                value
                for value in structured_values
                if value
            )

            return (
                structured_technique,
                evidence,
                "Structured repository metadata",
            )

        if descriptive_technique != "Other_Assays":
            evidence = " | ".join(
                value
                for value in descriptive_values
                if value
            )

            return (
                descriptive_technique,
                evidence,
                "Title and descriptive metadata",
            )

        evidence = " | ".join(
            value
            for value in (
                *structured_values,
                *descriptive_values,
            )
            if value
        )

        return (
            "Other_Assays",
            evidence,
            "No supported technique evidence",
        )

    @staticmethod
    def _technique_match(
        *,
        requested: str,
        verified: str,
    ) -> str:
        """Compare the requested profile with the verified assay."""
        requested_map = {
            "RNA-seq": {"RNA_seq", "scRNA_seq", "snRNA_seq"},
            "ChIP-seq": {"ChIP_seq"},
            "CUT&RUN": {"CUT_RUN"},
            "CUT&Tag": {"CUT_TAG"},
            "CLIP-seq": {
                "CLIP",
                "eCLIP",
                "iCLIP",
                "PAR_CLIP",
                "HITS_CLIP",
            },
        }

        expected = requested_map.get(requested, set())

        if verified in expected:
            return "Exact"

        if verified == "Other_Assays":
            return "Unverified"

        return "Mismatch"

    @staticmethod
    def _build_query(
        *,
        species: str,
        gene_term: str,
        profile: TechniqueSearchProfile,
    ) -> str:
        """Build a gene, technique, species, and GEO-series query."""
        technique_query = " OR ".join(
            f'"{term}"[All Fields]'
            for term in profile.terms
        )

        return (
            f'("{gene_term}"[All Fields]) AND '
            f"({technique_query}) AND "
            f'"{species}"[Organism]'
        )

    @staticmethod
    def _extract_year(value: str) -> int | None:
        """Extract a four-digit study year from a date string."""
        match = re.search(
            r"\b(?:19|20)\d{2}\b",
            value or "",
        )

        return int(match.group(0)) if match else None

    @staticmethod
    def _unique_terms(
        values: Iterable[str],
    ) -> tuple[str, ...]:
        """Deduplicate non-empty gene terms while preserving order."""
        terms: list[str] = []
        seen: set[str] = set()

        for value in values:
            normalized = str(value).strip()

            if not normalized:
                continue

            identity = normalized.casefold()

            if identity in seen:
                continue

            seen.add(identity)
            terms.append(normalized)

        return tuple(terms)

    @staticmethod
    def _deduplicate(
        records: list[DatasetRecord],
    ) -> list[DatasetRecord]:
        """Deduplicate equivalent gene-query discoveries."""
        output: list[DatasetRecord] = []
        seen: set[
            tuple[str, str, str, str]
        ] = set()

        for record in records:
            identity = (
                record.database.casefold(),
                record.accession.strip().upper(),
                record.gene_query_used.casefold(),
                record.technique_requested.casefold(),
            )

            if identity in seen:
                continue

            seen.add(identity)
            output.append(record)

        return output
