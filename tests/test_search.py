"""Tests for dataset search orchestration."""

from __future__ import annotations

import pytest

from dataset_finder.models import DatasetRecord
from dataset_finder.search import SearchService, UnsupportedDatabaseError


class FakeGEOClient:
    """Small test replacement for the live GEO client."""

    def __init__(self) -> None:
        self.received_species: str | None = None
        self.received_query: str | None = None
        self.received_max_results: int | None = None

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        """Return one predictable record."""
        self.received_species = species
        self.received_query = query
        self.received_max_results = max_results

        return [
            DatasetRecord(
                uid="200123456",
                accession="GSE123456",
                title="Example brain RNA-seq dataset",
                organism="Drosophila melanogaster",
                study_type=(
                    "Expression profiling by high throughput sequencing"
                ),
                sample_count=12,
                publication_date="2026/01/10",
                url=(
                    "https://www.ncbi.nlm.nih.gov/geo/query/"
                    "acc.cgi?acc=GSE123456"
                ),
            )
        ]


def test_search_service_calls_geo_client() -> None:
    """The service should delegate GEO searches to the GEO client."""
    client = FakeGEOClient()
    service = SearchService(geo_client=client)

    records = service.search(
        species=" Drosophila melanogaster ",
        query=" brain RNA-seq ",
        database="geo",
        max_results=5,
    )

    assert len(records) == 1
    assert records[0].accession == "GSE123456"
    assert client.received_species == "Drosophila melanogaster"
    assert client.received_query == "brain RNA-seq"
    assert client.received_max_results == 5


def test_all_temporarily_uses_geo_client() -> None:
    """The all option should currently fall back to GEO."""
    client = FakeGEOClient()
    service = SearchService(geo_client=client)

    records = service.search(
        species="Drosophila melanogaster",
        query="brain",
        database="all",
        max_results=3,
    )

    assert len(records) == 1
    assert client.received_max_results == 3


def test_sra_is_not_implemented() -> None:
    """The service should clearly reject SRA until it is implemented."""
    service = SearchService(geo_client=FakeGEOClient())

    with pytest.raises(
        UnsupportedDatabaseError,
        match="SRA search is not implemented",
    ):
        service.search(
            species="Drosophila melanogaster",
            query="brain RNA-seq",
            database="sra",
        )


@pytest.mark.parametrize(
    ("species", "query", "message"),
    [
        ("", "brain", "Species cannot be empty"),
        ("Drosophila melanogaster", "", "Query cannot be empty"),
    ],
)
def test_search_service_validates_text_inputs(
    species: str,
    query: str,
    message: str,
) -> None:
    """The service should reject empty required search values."""
    service = SearchService(geo_client=FakeGEOClient())

    with pytest.raises(ValueError, match=message):
        service.search(
            species=species,
            query=query,
        )


def test_search_service_rejects_invalid_result_limit() -> None:
    """The service should reject nonpositive result limits."""
    service = SearchService(geo_client=FakeGEOClient())

    with pytest.raises(
        ValueError,
        match="Maximum results must be greater than zero",
    ):
        service.search(
            species="Drosophila melanogaster",
            query="brain RNA-seq",
            max_results=0,
        )
