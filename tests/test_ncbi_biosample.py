"""Tests for the NCBI BioSample client."""

from dataset_finder.clients.ncbi_biosample import (
    NCBIBioSampleClient,
)


class FakeEntrez:
    """Return deterministic BioSample responses."""

    def search_ids(
        self,
        *,
        database: str,
        term: str,
        max_results: int,
    ) -> list[str]:
        assert database == "biosample"
        assert "Drosophila melanogaster" in term
        assert max_results == 5
        return ["101"]

    def summaries(
        self,
        *,
        database: str,
        identifiers: list[str],
    ) -> list[dict[str, object]]:
        assert database == "biosample"
        assert identifiers == ["101"]

        return [
            {
                "uid": "101",
                "accession": "SAMN40876891",
                "title": "insc-GAL4 orb2 RNAi brain",
                "organism": "Drosophila melanogaster",
                "submissiondate": "2024/04/08",
                "attributes": [
                    {
                        "name": "tissue",
                        "value": "brain",
                    },
                    {
                        "name": "developmental stage",
                        "value": "larva",
                    },
                    {
                        "name": "genotype",
                        "value": "inscGAL4>orb2RNAi",
                    },
                    {
                        "name": "sex",
                        "value": "female",
                    },
                ],
                "links": (
                    "PRJNA1097674 SRX24191199 "
                    "SRR28900001"
                ),
            }
        ]


def test_biosample_search_returns_normalized_record() -> None:
    client = NCBIBioSampleClient(
        entrez_client=FakeEntrez(),
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="orb2",
        max_results=5,
    )

    assert len(records) == 1

    record = records[0]

    assert record.database == "BioSample"
    assert record.accession == "SAMN40876891"
    assert record.organism == "Drosophila melanogaster"
    assert record.sample_count == 1
    assert record.tissue == "brain"
    assert record.developmental_stage == "larva"
    assert record.genotype == "inscGAL4>orb2RNAi"
    assert record.sex == "female"
    assert record.project_accession == "PRJNA1097674"
    assert record.experiment_accessions == (
        "SRX24191199",
    )
    assert record.sample_accessions == (
        "SRR28900001",
    )
    assert record.biosample_accessions == (
        "SAMN40876891",
    )


def test_biosample_falls_back_to_accession_in_text() -> None:
    record = NCBIBioSampleClient._normalize_summary(
        {
            "uid": "202",
            "title": "Drosophila sample SAMN12345678",
            "organism": "Drosophila melanogaster",
        }
    )

    assert record.accession == "SAMN12345678"
