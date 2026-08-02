"""Tests for the European Nucleotide Archive client."""

from __future__ import annotations

from typing import Any

import pytest
import requests

from dataset_finder.clients.ena import (
    ENAClient,
    ENAClientError,
)


class FakeResponse:
    """Minimal requests response replacement."""

    def __init__(
        self,
        text: str,
        *,
        error: requests.RequestException | None = None,
    ) -> None:
        self.text = text
        self.error = error

    def raise_for_status(self) -> None:
        if self.error:
            raise self.error


class FakeSession:
    """Return a deterministic ENA response."""

    def __init__(self, response: FakeResponse) -> None:
        self.response = response
        self.calls: list[dict[str, Any]] = []

    def get(
        self,
        url: str,
        *,
        params: dict[str, Any],
        timeout: float,
    ) -> FakeResponse:
        self.calls.append(
            {
                "url": url,
                "params": params,
                "timeout": timeout,
            }
        )
        return self.response


def test_ena_search_returns_normalized_studies() -> None:
    session = FakeSession(
        FakeResponse(
            "accession\tdescription\n"
            "SRP590197\tOrb2 RNA-binding study\n"
            "SRP590257\tOrb2 developmental study\n"
        )
    )

    client = ENAClient(
        session=session,  # type: ignore[arg-type]
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="orb2",
        max_results=5,
    )

    assert len(records) == 2

    first = records[0]

    assert first.database == "ENA"
    assert first.accession == "SRP590197"
    assert first.project_accession == "SRP590197"
    assert first.title == "Orb2 RNA-binding study"
    assert first.organism == "Drosophila melanogaster"
    assert first.study_type == "ENA Study"
    assert first.url.endswith("/SRP590197")

    assert session.calls[0]["params"] == {
        "result": "study",
        "query": "orb2 Drosophila melanogaster",
        "limit": 5,
    }


def test_ena_search_handles_empty_response() -> None:
    client = ENAClient(
        session=FakeSession(
            FakeResponse(""),
        ),  # type: ignore[arg-type]
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="unlikely-query",
        max_results=5,
    )

    assert records == []


def test_ena_search_wraps_request_errors() -> None:
    client = ENAClient(
        session=FakeSession(
            FakeResponse(
                "",
                error=requests.HTTPError(
                    "503 Service Unavailable"
                ),
            )
        ),  # type: ignore[arg-type]
    )

    with pytest.raises(
        ENAClientError,
        match="ENA request failed",
    ):
        client.search(
            species="Drosophila melanogaster",
            query="orb2",
            max_results=5,
        )


@pytest.mark.parametrize(
    ("species", "query", "message"),
    [
        ("", "orb2", "Species cannot be empty"),
        (
            "Drosophila melanogaster",
            "",
            "Query cannot be empty",
        ),
    ],
)
def test_ena_rejects_empty_inputs(
    species: str,
    query: str,
    message: str,
) -> None:
    client = ENAClient()

    with pytest.raises(ValueError, match=message):
        client.search(
            species=species,
            query=query,
        )
