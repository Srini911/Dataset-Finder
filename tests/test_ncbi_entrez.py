from dataset_finder.clients.ncbi_entrez import NCBIEntrezClient


def test_search_ids_passes_date_and_pagination_parameters() -> None:
    class RecordingClient(NCBIEntrezClient):
        def __init__(self) -> None:
            self.parameters = {}

        def request_json(
            self,
            endpoint: str,
            *,
            parameters: dict,
        ) -> dict:
            assert endpoint == "esearch.fcgi"
            self.parameters = parameters

            return {
                "esearchresult": {
                    "idlist": ["10", "11"],
                }
            }

    client = RecordingClient()

    identifiers = client.search_ids(
        database="sra",
        term="orb2",
        max_results=2,
        start=100,
        minimum_date="2005/01/01",
        maximum_date="2026/12/31",
        date_type="pdat",
    )

    assert identifiers == ["10", "11"]
    assert client.parameters["retstart"] == 100
    assert client.parameters["retmax"] == 2
    assert client.parameters["mindate"] == "2005/01/01"
    assert client.parameters["maxdate"] == "2026/12/31"
    assert client.parameters["datetype"] == "pdat"


def test_search_ids_paged_collects_multiple_pages() -> None:
    class PagedClient(NCBIEntrezClient):
        def __init__(self) -> None:
            self.starts: list[int] = []

        def search_ids(
            self,
            *,
            database: str,
            term: str,
            max_results: int,
            start: int = 0,
            minimum_date: str = "",
            maximum_date: str = "",
            date_type: str = "",
        ) -> list[str]:
            del (
                database,
                term,
                minimum_date,
                maximum_date,
                date_type,
            )

            self.starts.append(start)

            pages = {
                0: ["1", "2"],
                2: ["3", "4"],
                4: ["5"],
            }

            return pages.get(start, [])[:max_results]

    client = PagedClient()

    identifiers = client.search_ids_paged(
        database="sra",
        term="orb2",
        max_results=5,
        page_size=2,
        minimum_date="2005/01/01",
        maximum_date="2026/12/31",
        date_type="pdat",
    )

    assert identifiers == ["1", "2", "3", "4", "5"]
    assert client.starts == [0, 2, 4]


def test_search_ids_paged_stops_on_duplicate_page() -> None:
    class DuplicateClient(NCBIEntrezClient):
        def search_ids(
            self,
            *,
            database: str,
            term: str,
            max_results: int,
            start: int = 0,
            minimum_date: str = "",
            maximum_date: str = "",
            date_type: str = "",
        ) -> list[str]:
            del (
                database,
                term,
                max_results,
                start,
                minimum_date,
                maximum_date,
                date_type,
            )

            return ["1", "2"]

    client = DuplicateClient()

    identifiers = client.search_ids_paged(
        database="sra",
        term="orb2",
        max_results=10,
        page_size=2,
    )

    assert identifiers == ["1", "2"]
