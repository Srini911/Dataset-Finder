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
        return ["101", "102"]

    def summaries(
        self,
        *,
        database: str,
        identifiers: list[str],
    ) -> list[dict]:
        assert database == "sra"
        assert identifiers == ["101", "102"]

        return [
            {
                "uid": "101",
                "createdate": "2022/06/02",
                "expxml": (
                    "<Summary><Title>bru1 rep1</Title>"
                    '<Platform instrument_model="Illumina HiSeq 2500">'
                    "ILLUMINA</Platform></Summary>"
                    '<Experiment acc="SRX15499370"/>'
                    '<Study acc="SRP377648" name="Bruno1 study"/>'
                    '<Organism ScientificName="Drosophila melanogaster"/>'
                    '<Sample acc="SRS13214240"/>'
                    "<Library_descriptor>"
                    "<LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>"
                    "<LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>"
                    "<LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>"
                    "<LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT>"
                    "</Library_descriptor>"
                    "<Bioproject>PRJNA843701</Bioproject>"
                    "<Biosample>SAMN28766850</Biosample>"
                ),
                "runs": '<Run acc="SRR19446263"/>',
            },
            {
                "uid": "102",
                "createdate": "2022/06/02",
                "expxml": (
                    "<Summary><Title>bru1 rep2</Title></Summary>"
                    '<Experiment acc="SRX15499371"/>'
                    '<Study acc="SRP377648" name="Bruno1 study"/>'
                    '<Organism ScientificName="Drosophila melanogaster"/>'
                    '<Sample acc="SRS13214242"/>'
                    "<Library_descriptor>"
                    "<LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>"
                    "<LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>"
                    "<LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>"
                    "<LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT>"
                    "</Library_descriptor>"
                    "<Bioproject>PRJNA843701</Bioproject>"
                    "<Biosample>SAMN28766849</Biosample>"
                ),
                "runs": '<Run acc="SRR19446262"/>',
            },
        ]


def test_sra_search_aggregates_study_metadata() -> None:
    client = NCBISRAClient(
        entrez_client=FakeEntrez(),
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="bru1",
        max_results=5,
    )

    assert len(records) == 1

    record = records[0]

    assert record.accession == "SRP377648"
    assert record.title == "Bruno1 study"
    assert record.project_accession == "PRJNA843701"
    assert record.database == "SRA"
    assert record.sample_count == 2
    assert record.experiment_accessions == (
        "SRX15499370",
        "SRX15499371",
    )
    assert record.sample_accessions == (
        "SRR19446263",
        "SRR19446262",
    )
    assert "SAMN28766850" in record.biosample_accessions
    assert record.library_strategy == "RNA-Seq"
    assert record.library_source == "TRANSCRIPTOMIC"
    assert record.library_selection == "cDNA"
    assert record.library_layout == "PAIRED"
    assert record.platform == "Illumina HiSeq 2500"
