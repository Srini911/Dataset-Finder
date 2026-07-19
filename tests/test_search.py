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


class FakeENCODEClient:
    """Small test replacement for the live ENCODE client."""

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
        """Return one predictable ENCODE record."""
        self.received_species = species
        self.received_query = query
        self.received_max_results = max_results

        return [
            DatasetRecord(
                uid="ENCSR123ABC",
                accession="ENCSR123ABC",
                title="Example CTCF ChIP-seq dataset",
                organism="Homo sapiens",
                study_type="TF ChIP-seq",
                sample_count=2,
                publication_date="2026-01-10",
                url=(
                    "https://www.encodeproject.org/"
                    "experiments/ENCSR123ABC/"
                ),
            )
        ]


def test_all_combines_geo_and_encode_clients() -> None:
    """The all option should combine GEO and ENCODE records."""
    geo_client = FakeGEOClient()
    encode_client = FakeENCODEClient()
    service = SearchService(
        geo_client=geo_client,
        encode_client=encode_client,
    )

    records = service.search(
        species="Homo sapiens",
        query="CTCF",
        database="all",
        max_results=3,
    )

    assert len(records) == 2
    assert records[0].accession == "GSE123456"
    assert records[1].accession == "ENCSR123ABC"

    assert geo_client.received_species == "Homo sapiens"
    assert geo_client.received_query == "CTCF"
    assert geo_client.received_max_results == 3

    assert encode_client.received_species == "Homo sapiens"
    assert encode_client.received_query == "CTCF"
    assert encode_client.received_max_results == 3


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
