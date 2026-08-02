"""Tests for multi-gene batch searches."""

from dataset_finder.batch import BatchSearchService
from dataset_finder.models import DatasetRecord


def make_record(accession: str, title: str) -> DatasetRecord:
    """Create a test dataset record."""
    return DatasetRecord(
        uid=accession,
        accession=accession,
        title=title,
        organism="Drosophila melanogaster",
        study_type="Expression profiling by high throughput sequencing",
        sample_count=4,
        publication_date="2026-01-01",
        url=f"https://example.org/{accession}",
    )


class FakeSearchService:
    """Simple deterministic search service."""

    def search(
        self,
        *,
        species: str,
        query: str,
        database: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, database, max_results
        return [
            make_record(
                f"GSE-{query}",
                f"{query} adult brain RNA-seq",
            )
        ]


def test_batch_search_assigns_gene_and_technique() -> None:
    service = BatchSearchService(
        search_service=FakeSearchService(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["bru1", "Hrb98DE"],
        database="geo",
        max_results_per_gene=10,
        gene_set="RBP",
    )

    assert len(result.records) == 2
    assert result.records[0].gene == "bru1"
    assert result.records[0].gene_set == "RBP"
    assert result.records[0].technique == "RNA_seq"
    assert result.records[0].database == "GEO"
    assert not result.issues


class PartiallyFailingSearchService:
    """Search service that fails for one gene."""

    def search(
        self,
        *,
        species: str,
        query: str,
        database: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, database, max_results

        if query == "bad":
            raise RuntimeError("Temporary database failure")

        return [make_record("GSE1", "RNA-seq")]


def test_batch_search_records_gene_errors() -> None:
    service = BatchSearchService(
        search_service=PartiallyFailingSearchService(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["good", "bad"],
        database="geo",
        max_results_per_gene=10,
    )

    assert len(result.records) == 1
    assert len(result.issues) == 1
    assert result.issues[0].gene == "bad"
