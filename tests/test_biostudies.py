"""Tests for the BioStudies client."""

from dataset_finder.clients.biostudies import BioStudiesClient


class FakeResponse:
    """Minimal fake HTTP response."""

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict:
        return {
            "hits": [
                {
                    "accession": "E-MTAB-493",
                    "type": "study",
                    "title": (
                        "Tissue-wide RNA-seq expression profiles "
                        "of Drosophila melanogaster"
                    ),
                    "author": "Example Author",
                    "files": 6,
                    "release_date": "2011-02-28",
                    "content": (
                        "Drosophila melanogaster RNA-seq "
                        "adult tissues"
                    ),
                }
            ]
        }


class FakeSession:
    """Fake requests session."""

    def __init__(self) -> None:
        self.headers: dict[str, str] = {}
        self.params: dict | None = None

    def get(
        self,
        url: str,
        *,
        params: dict,
        timeout: float,
    ) -> FakeResponse:
        del url, timeout
        self.params = params
        return FakeResponse()


def test_biostudies_search_normalizes_arrayexpress_record() -> None:
    session = FakeSession()
    client = BioStudiesClient(session=session)

    records = client.search(
        species="Drosophila melanogaster",
        query="RNA-seq",
        max_results=3,
    )

    assert len(records) == 1

    record = records[0]

    assert record.accession == "E-MTAB-493"
    assert record.database == "BioStudies"
    assert record.sample_count == 6
    assert record.publication_date == "2011-02-28"
    assert "arrayexpress/studies/E-MTAB-493" in record.url

    assert session.params is not None
    assert session.params["collection"] == "ArrayExpress"
    assert session.params["pageSize"] == 3
