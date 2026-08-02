"""Tests for the NCBI BioProject client."""

from dataset_finder.clients.ncbi_bioproject import (
    NCBIBioProjectClient,
)


class FakeEntrez:
    """Deterministic Entrez client for BioProject tests."""

    def search_ids(
        self,
        *,
        database: str,
        term: str,
        max_results: int,
    ) -> list[str]:
        assert database == "bioproject"
        assert "bru1" in term
        return ["42"]

    def summaries(
        self,
        *,
        database: str,
        identifiers: list[str],
    ) -> list[dict]:
        return [
            {
                "uid": "42",
                "project_acc": "PRJNA123456",
                "project_title": "Fly brain transcriptomics",
                "project_description": "RNA-seq study",
                "organism_name": "Drosophila melanogaster",
                "registration_date": "2026-01-02",
            }
        ]


def test_bioproject_search_normalizes_record() -> None:
    client = NCBIBioProjectClient(
        entrez_client=FakeEntrez(),
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="bru1",
        max_results=5,
    )

    assert len(records) == 1
    assert records[0].accession == "PRJNA123456"
    assert records[0].database == "BioProject"
    assert records[0].organism == "Drosophila melanogaster"
