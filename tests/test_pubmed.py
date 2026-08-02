"""Tests for the PubMed client."""

from dataset_finder.clients.pubmed import PubMedClient


class FakeEntrez:
    """Return deterministic PubMed summaries."""

    def search_ids(
        self,
        *,
        database: str,
        term: str,
        max_results: int,
    ) -> list[str]:
        assert database == "pubmed"
        assert "(orb2)" in term
        assert "Drosophila melanogaster" in term
        assert max_results == 5
        return ["38902233"]

    def summaries(
        self,
        *,
        database: str,
        identifiers: list[str],
    ) -> list[dict[str, object]]:
        assert database == "pubmed"
        assert identifiers == ["38902233"]

        return [
            {
                "uid": "38902233",
                "title": (
                    "Orb2 regulation in Drosophila "
                    "GSE263513 PRJNA1097674"
                ),
                "fulljournalname": "Developmental Biology",
                "pubdate": "2024 Jun",
                "authors": [
                    {"name": "Stewart R"},
                    {"name": "Fox D"},
                ],
                "articleids": [
                    {
                        "idtype": "doi",
                        "value": "10.1000/orb2.test",
                    }
                ],
            }
        ]


def test_pubmed_search_returns_publication_record() -> None:
    client = PubMedClient(
        entrez_client=FakeEntrez(),
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="orb2",
        max_results=5,
    )

    assert len(records) == 1

    record = records[0]

    assert record.database == "PubMed"
    assert record.accession == "PMID38902233"
    assert record.pubmed_ids == ("38902233",)
    assert record.dois == ("10.1000/orb2.test",)
    assert record.related_geo_accessions == (
        "GSE263513",
    )
    assert record.related_bioproject_accessions == (
        "PRJNA1097674",
    )
    assert "Stewart R" in record.publication
    assert record.url.endswith("/38902233/")
