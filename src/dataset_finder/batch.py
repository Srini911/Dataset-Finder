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
from dataset_finder.models import DatasetRecord
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
class BatchSearchResult:
    """Results and errors from a multi-gene search."""

    records: tuple[DatasetRecord, ...]
    issues: tuple[BatchSearchIssue, ...]
    database_statuses: tuple[BatchDatabaseStatus, ...]
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
        search_date = datetime.now(UTC).date().isoformat()

        for submitted_gene in genes:
            resolved_gene = self.flybase_resolver.resolve(
                submitted_gene
            )

            flyatlas_expression = self._fetch_flyatlas(
                resolved_gene,
            )

            query = self._build_search_query(
                submitted_gene,
                resolved_gene,
            )

            try:
                if hasattr(self.search_service, "search_with_status"):
                    outcome = self.search_service.search_with_status(
                        species=species,
                        query=query,
                        database=database,
                        max_results=max_results_per_gene,
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
                        max_results=max_results_per_gene,
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

            for record in gene_records:
                technique = classify_technique(
                    record.study_type,
                    record.title,
                    record.description,
                    record.evidence_text,
                )

                confidence = self._confidence_label(
                    resolved_gene
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
                        database=(
                            record.database
                            or self._infer_database(record)
                        ),
                        technique=record.technique or technique,
                        match_type=resolved_gene.match_type,
                        confidence=confidence,
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

        return BatchSearchResult(
            records=tuple(self._deduplicate(records)),
            issues=tuple(issues),
            database_statuses=tuple(database_statuses),
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
