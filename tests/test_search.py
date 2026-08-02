"""Tests for dataset search orchestration."""

from __future__ import annotations

import pytest

from dataset_finder.models import DatasetRecord
from dataset_finder.search import SearchService


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


class FakeSRAClient:
    """Return no SRA records during combined-search tests."""

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, query, max_results
        return []


class FakeBioSampleClient:
    """Return no BioSample records during combined-search tests."""

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, query, max_results
        return []


class FakeENAClient:
    """Return no ENA records during combined-search tests."""

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, query, max_results
        return []


class FakeBioStudiesClient:
    """Return no BioStudies records during combined-search tests."""

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, query, max_results
        return []


class FakeBioProjectClient:
    """Return no BioProject records during combined-search tests."""

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, query, max_results
        return []


def test_search_service_calls_geo_client() -> None:
    """The service should delegate GEO searches to the GEO client."""
    client = FakeGEOClient()
    service = SearchService(
        geo_client=client,
        encode_client=FakeENCODEClient(),
        sra_client=FakeSRAClient(),
        bioproject_client=FakeBioProjectClient(),
        biosample_client=FakeBioSampleClient(),
        biostudies_client=FakeBioStudiesClient(),
    )

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


def test_all_combines_supported_clients() -> None:
    """The all option should combine supported client records."""
    geo_client = FakeGEOClient()
    encode_client = FakeENCODEClient()

    service = SearchService(
        geo_client=geo_client,
        encode_client=encode_client,
        sra_client=FakeSRAClient(),
        bioproject_client=FakeBioProjectClient(),
        biosample_client=FakeBioSampleClient(),
        biostudies_client=FakeBioStudiesClient(),
        ena_client=FakeENAClient(),
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


def test_sra_search_is_supported() -> None:
    """The service should route SRA searches to the SRA client."""

    class OneRecordSRAClient:
        def search(
            self,
            *,
            species: str,
            query: str,
            max_results: int,
        ) -> list[DatasetRecord]:
            del species, query, max_results

            return [
                DatasetRecord(
                    uid="101",
                    accession="SRP123456",
                    title="Drosophila brain RNA-seq",
                    organism="Drosophila melanogaster",
                    study_type="Sequence Read Archive",
                    sample_count=2,
                    publication_date="2026/01/01",
                    url=(
                        "https://www.ncbi.nlm.nih.gov/sra/"
                        "?term=SRP123456"
                    ),
                    database="SRA",
                )
            ]

    service = SearchService(
        geo_client=FakeGEOClient(),
        encode_client=FakeENCODEClient(),
        sra_client=OneRecordSRAClient(),
        bioproject_client=FakeBioProjectClient(),
        biostudies_client=FakeBioStudiesClient(),
    )

    records = service.search(
        species="Drosophila melanogaster",
        query="bru1",
        database="sra",
        max_results=5,
    )

    assert len(records) == 1
    assert records[0].accession == "SRP123456"
    assert records[0].database == "SRA"


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
    service = SearchService(
        geo_client=FakeGEOClient(),
        encode_client=FakeENCODEClient(),
        sra_client=FakeSRAClient(),
        bioproject_client=FakeBioProjectClient(),
        biostudies_client=FakeBioStudiesClient(),
    )

    with pytest.raises(ValueError, match=message):
        service.search(
            species=species,
            query=query,
        )


def test_search_service_rejects_invalid_result_limit() -> None:
    """The service should reject nonpositive result limits."""
    service = SearchService(
        geo_client=FakeGEOClient(),
        encode_client=FakeENCODEClient(),
        sra_client=FakeSRAClient(),
        bioproject_client=FakeBioProjectClient(),
        biostudies_client=FakeBioStudiesClient(),
    )

    with pytest.raises(
        ValueError,
        match="Maximum results must be greater than zero",
    ):
        service.search(
            species="Drosophila melanogaster",
            query="brain RNA-seq",
            max_results=0,
        )
