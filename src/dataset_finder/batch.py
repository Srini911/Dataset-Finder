"""Multi-gene batch search orchestration."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import UTC, datetime

from dataset_finder.assay_classifier import classify_technique
from dataset_finder.clients.flyatlas import (
    FlyAtlasClient,
    FlyAtlasClientError,
    FlyAtlasExpression,
)
from dataset_finder.flybase_resolver import FlyBaseGene, FlyBaseResolver
from dataset_finder.metadata import extract_biological_metadata
from dataset_finder.models import DatasetRecord
from dataset_finder.relevance import assess_relevance
from dataset_finder.search import (
    DatabaseSearchStatus,
    SearchService,
)


@dataclass(frozen=True, slots=True)
class BatchSearchIssue:
    """An error encountered while searching one gene."""

    gene: str
    database: str
    message: str


@dataclass(frozen=True, slots=True)
class BatchDatabaseStatus:
    """Database status for one submitted gene."""

    gene: str
    database: str
    success: bool
    result_count: int
    error: str = ""


@dataclass(frozen=True, slots=True)
class BatchGeneAnnotation:
    """FlyBase and FlyAtlas annotation for one submitted gene."""

    gene: FlyBaseGene
    flyatlas: FlyAtlasExpression


@dataclass(frozen=True, slots=True)
class BatchSearchResult:
    """Results and errors from a multi-gene search."""

    records: tuple[DatasetRecord, ...]
    issues: tuple[BatchSearchIssue, ...]
    database_statuses: tuple[BatchDatabaseStatus, ...]
    gene_annotations: tuple[BatchGeneAnnotation, ...]
    genes: tuple[str, ...]
    gene_set: str
    database: str


class BatchSearchService:
    """Search multiple genes and enrich records with FlyBase metadata."""

    def __init__(
        self,
        *,
        search_service: SearchService | None = None,
        flybase_resolver: FlyBaseResolver | None = None,
        flyatlas_client: FlyAtlasClient | None = None,
    ) -> None:
        self.search_service = search_service or SearchService()
        self.flybase_resolver = flybase_resolver or FlyBaseResolver()
        self.flyatlas_client = flyatlas_client or FlyAtlasClient()

    def search_many(
        self,
        *,
        species: str,
        genes: list[str],
        database: str,
        max_results_per_gene: int,
        gene_set: str = "",
    ) -> BatchSearchResult:
        """Search all genes and return enriched, classified records."""
        if not genes:
            raise ValueError("At least one gene is required.")

        records: list[DatasetRecord] = []
        issues: list[BatchSearchIssue] = []
        database_statuses: list[BatchDatabaseStatus] = []
        gene_annotations: list[BatchGeneAnnotation] = []
        search_date = datetime.now(UTC).date().isoformat()

        for submitted_gene in genes:
            resolved_gene = self.flybase_resolver.resolve(
                submitted_gene
            )

            flyatlas_expression = self._fetch_flyatlas(
                resolved_gene,
            )

            gene_annotations.append(
                BatchGeneAnnotation(
                    gene=resolved_gene,
                    flyatlas=flyatlas_expression,
                )
            )

            query = self._build_search_query(
                submitted_gene,
                resolved_gene,
            )

            candidate_limit = min(
                100,
                max(
                    50,
                    max_results_per_gene * 20,
                ),
            )

            try:
                if hasattr(self.search_service, "search_with_status"):
                    outcome = self.search_service.search_with_status(
                        species=species,
                        query=query,
                        database=database,
                        max_results=candidate_limit,
                    )
                    gene_records = list(outcome.records)

                    for status in outcome.statuses:
                        database_statuses.append(
                            self._batch_status(
                                submitted_gene,
                                status,
                            )
                        )
                else:
                    gene_records = self.search_service.search(
                        species=species,
                        query=query,
                        database=database,
                        max_results=candidate_limit,
                    )
            except Exception as exc:
                issues.append(
                    BatchSearchIssue(
                        gene=submitted_gene,
                        database=database,
                        message=str(exc),
                    )
                )
                continue

            accepted_for_gene = 0

            for record in gene_records:
                if accepted_for_gene >= max_results_per_gene:
                    break

                source_database = (
                    record.database
                    or self._infer_database(record)
                )

                assessment_record = replace(
                    record,
                    database=source_database,
                )

                technique = classify_technique(
                    record.study_type,
                    record.title,
                    record.description,
                    record.evidence_text,
                )

                relevance = assess_relevance(
                    record=assessment_record,
                    submitted_gene=submitted_gene,
                    resolved_gene=resolved_gene,
                )

                if not relevance.accepted:
                    continue

                biological_metadata = extract_biological_metadata(
                    record.title,
                    record.description,
                    record.study_type,
                    record.organism,
                    record.raw_metadata,
                )

                confidence = self._combined_confidence(
                    gene_confidence=self._confidence_label(
                        resolved_gene
                    ),
                    relevance_confidence=relevance.confidence,
                )

                records.append(
                    replace(
                        record,
                        gene=submitted_gene,
                        gene_set=gene_set,
                        official_symbol=(
                            resolved_gene.official_symbol
                            or submitted_gene
                        ),
                        flybase_id=resolved_gene.flybase_id,
                        synonyms=resolved_gene.synonyms,
                        database=source_database,
                        technique=record.technique or technique,
                        tissue=(
                            record.tissue
                            or biological_metadata.tissue
                        ),
                        cell_type=(
                            record.cell_type
                            or biological_metadata.cell_type
                        ),
                        developmental_stage=(
                            record.developmental_stage
                            or biological_metadata.developmental_stage
                        ),
                        sex=record.sex or biological_metadata.sex,
                        genotype=(
                            record.genotype
                            or biological_metadata.genotype
                        ),
                        strain=(
                            record.strain
                            or biological_metadata.strain
                        ),
                        treatment=(
                            record.treatment
                            or biological_metadata.treatment
                        ),
                        control_status=(
                            record.control_status
                            or biological_metadata.control_status
                        ),
                        disease=(
                            record.disease
                            or biological_metadata.disease
                        ),
                        time_point=(
                            record.time_point
                            or biological_metadata.time_point
                        ),
                        perturbation=(
                            record.perturbation
                            or biological_metadata.perturbation
                        ),
                        match_type=relevance.match_type,
                        confidence=confidence,
                        evidence_text=(
                            record.evidence_text
                            or relevance.evidence
                        ),
                        search_date=record.search_date or search_date,
                        flybase_url=resolved_gene.flybase_url,
                        flyatlas_url=resolved_gene.flyatlas_url,
                        flyatlas_brain_male_fpkm=(
                            flyatlas_expression.brain_male_fpkm
                        ),
                        flyatlas_brain_female_fpkm=(
                            flyatlas_expression.brain_female_fpkm
                        ),
                        flyatlas_brain_larval_fpkm=(
                            flyatlas_expression.brain_larval_fpkm
                        ),
                        flyatlas_head_male_fpkm=(
                            flyatlas_expression.head_male_fpkm
                        ),
                        flyatlas_head_female_fpkm=(
                            flyatlas_expression.head_female_fpkm
                        ),
                        flyatlas_top_male_tissue=(
                            flyatlas_expression.top_male_tissue
                        ),
                        flyatlas_top_male_fpkm=(
                            flyatlas_expression.top_male_fpkm
                        ),
                        flyatlas_top_female_tissue=(
                            flyatlas_expression.top_female_tissue
                        ),
                        flyatlas_top_female_fpkm=(
                            flyatlas_expression.top_female_fpkm
                        ),
                        flyatlas_top_larval_tissue=(
                            flyatlas_expression.top_larval_tissue
                        ),
                        flyatlas_top_larval_fpkm=(
                            flyatlas_expression.top_larval_fpkm
                        ),
                    )
                )
                accepted_for_gene += 1

        return BatchSearchResult(
            records=tuple(self._deduplicate(records)),
            issues=tuple(issues),
            database_statuses=tuple(database_statuses),
            gene_annotations=tuple(gene_annotations),
            genes=tuple(genes),
            gene_set=gene_set,
            database=database,
        )

    @staticmethod
    def _batch_status(
        gene: str,
        status: DatabaseSearchStatus,
    ) -> BatchDatabaseStatus:
        """Attach a database search status to its submitted gene."""
        return BatchDatabaseStatus(
            gene=gene,
            database=status.database,
            success=status.success,
            result_count=status.result_count,
            error=status.error,
        )

    def _fetch_flyatlas(
        self,
        resolved_gene: FlyBaseGene,
    ) -> FlyAtlasExpression:
        """Fetch FlyAtlas values without interrupting dataset searches."""
        if not resolved_gene.flybase_id:
            return FlyAtlasExpression(
                flybase_id="",
                symbol="",
            )

        try:
            return self.flyatlas_client.fetch(
                resolved_gene.flybase_id
            )
        except FlyAtlasClientError:
            return FlyAtlasExpression(
                flybase_id=resolved_gene.flybase_id,
                symbol=resolved_gene.official_symbol,
            )

    @staticmethod
    def _build_search_query(
        submitted_gene: str,
        resolved_gene: FlyBaseGene,
    ) -> str:
        """Build a conservative query from symbol and FlyBase ID."""
        terms: list[str] = []

        for value in (
            submitted_gene,
            resolved_gene.official_symbol,
            resolved_gene.flybase_id,
        ):
            value = value.strip()

            if not value:
                continue

            if value.casefold() in {
                existing.casefold()
                for existing in terms
            }:
                continue

            terms.append(value)

        if len(terms) == 1:
            return terms[0]

        return " OR ".join(
            f'"{term}"'
            for term in terms
        )

    @staticmethod
    def _combined_confidence(
        *,
        gene_confidence: str,
        relevance_confidence: str,
    ) -> str:
        """Return the weaker of gene and dataset confidence."""
        rank = {
            "Low": 1,
            "Medium": 2,
            "High": 3,
        }

        return min(
            (gene_confidence, relevance_confidence),
            key=lambda value: rank.get(value, 0),
        )

    @staticmethod
    def _confidence_label(
        resolved_gene: FlyBaseGene,
    ) -> str:
        """Assign a resolution-confidence label."""
        if not resolved_gene.flybase_id:
            return "Low"

        if resolved_gene.match_type == "official_symbol":
            return "High"

        if (
            resolved_gene.match_type == "synonym"
            and not resolved_gene.ambiguous
        ):
            return "High"

        if resolved_gene.match_type == "synonym":
            return "Medium"

        return "Low"

    @staticmethod
    def _infer_database(record: DatasetRecord) -> str:
        """Infer the source database from accession or URL."""
        accession = record.accession.upper()
        url = record.url.casefold()

        if accession.startswith(("GSE", "GSM", "GPL")) or "geo" in url:
            return "GEO"

        if accession.startswith(("ENCSR", "ENCBS", "ENCFF")) or "encodeproject" in url:
            return "ENCODE"

        if accession.startswith(("SRP", "SRR", "SRS", "SRX")) or "/sra" in url:
            return "SRA"

        if accession.startswith(("PRJNA", "PRJEB", "PRJDB")):
            return "BioProject"

        if (
            accession.startswith(("E-", "S-BSST"))
            or "biostudies" in url
        ):
            return "BioStudies"

        return ""

    @staticmethod
    def _deduplicate(
        records: list[DatasetRecord],
    ) -> list[DatasetRecord]:
        """Deduplicate records while preserving gene-specific evidence."""
        deduplicated: list[DatasetRecord] = []
        seen: set[tuple[str, str, str]] = set()

        for record in records:
            identity = (
                record.gene.casefold(),
                record.accession.strip().upper(),
                record.url.strip(),
            )

            if identity in seen:
                continue

            seen.add(identity)
            deduplicated.append(record)

        return deduplicated
