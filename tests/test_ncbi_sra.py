"""Tests for the NCBI SRA client."""

from dataset_finder.clients.ncbi_sra import NCBISRAClient


class FakeEntrez:
    """Deterministic Entrez client for SRA tests."""

    def search_ids(
        self,
        *,
        database: str,
        term: str,
        max_results: int,
    ) -> list[str]:
        assert database == "sra"
        assert "Drosophila melanogaster" in term
        assert max_results == 5
        return ["101"]

    def summaries(
        self,
        *,
        database: str,
        identifiers: list[str],
    ) -> list[dict]:
        assert database == "sra"
        assert identifiers == ["101"]

        return [
            {
                "uid": "101",
                "title": "bru1 adult brain RNA-seq",
                "createdate": "2026/01/01",
                "expxml": (
                    '<Study acc="SRP123456"/> '
                    '<Experiment acc="SRX123456" '
                    'ScientificName="Drosophila melanogaster"/>'
                ),
                "runs": (
                    '<Run acc="SRR111"/>'
                    '<Run acc="SRR222"/>'
                ),
            }
        ]


def test_sra_search_normalizes_record() -> None:
    client = NCBISRAClient(
        entrez_client=FakeEntrez(),
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="bru1",
        max_results=5,
    )

    assert len(records) == 1
    assert records[0].accession == "SRP123456"
    assert records[0].database == "SRA"
    assert records[0].sample_count == 2
    assert records[0].sample_accessions == (
        "SRR111",
        "SRR222",
    )
