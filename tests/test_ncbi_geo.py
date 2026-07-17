"""Tests for the NCBI GEO client."""

from __future__ import annotations

from typing import Any

import pytest
import requests

from dataset_finder.clients.ncbi_geo import (
    NCBIClientError,
    NCBIGEOClient,
)


class FakeResponse:
    """Minimal requests response replacement."""

    def __init__(
        self,
        payload: dict[str, Any],
        *,
        status_error: requests.RequestException | None = None,
    ) -> None:
        self.payload = payload
        self.status_error = status_error

    def raise_for_status(self) -> None:
        if self.status_error:
            raise self.status_error

    def json(self) -> dict[str, Any]:
        return self.payload


class FakeSession:
    """Return predefined responses for successive requests."""

    def __init__(self, responses: list[FakeResponse]) -> None:
        self.responses = responses
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
        return self.responses.pop(0)


def test_geo_search_returns_normalized_records() -> None:
    """The client should convert NCBI summaries into GEO records."""
    session = FakeSession(
        [
            FakeResponse(
                {
                    "esearchresult": {
                        "idlist": ["200123456"],
                    }
                }
            ),
            FakeResponse(
                {
                    "result": {
                        "200123456": {
                            "uid": "200123456",
                            "accession": "GSE123456",
                            "title": "Drosophila brain RNA sequencing",
                            "taxon": ["Drosophila melanogaster"],
                            "gdstype": ["Expression profiling by high throughput sequencing"],
                            "n_samples": 12,
                            "pdat": "2026/01/10",
                        }
                    }
                }
            ),
        ]
    )

    client = NCBIGEOClient(session=session)  # type: ignore[arg-type]
    records = client.search(
        species="Drosophila melanogaster",
        query="brain RNA-seq",
        max_results=5,
    )

    assert len(records) == 1
    assert records[0].accession == "GSE123456"
    assert records[0].organism == "Drosophila melanogaster"
    assert records[0].sample_count == 12
    assert records[0].url.endswith("?acc=GSE123456")

    assert session.calls[0]["params"]["db"] == "gds"
    assert session.calls[0]["params"]["retmax"] == 5
    assert '"Drosophila melanogaster"[Organism]' in session.calls[0]["params"]["term"]


def test_geo_search_returns_empty_list_when_no_ids_are_found() -> None:
    """The client should handle searches with no results."""
    session = FakeSession(
        [
            FakeResponse(
                {
                    "esearchresult": {
                        "idlist": [],
                    }
                }
            )
        ]
    )

    client = NCBIGEOClient(session=session)  # type: ignore[arg-type]

    assert client.search(
        species="Drosophila melanogaster",
        query="unlikely-query",
    ) == []


def test_geo_search_rejects_empty_species() -> None:
    """Species must contain visible text."""
    client = NCBIGEOClient()

    with pytest.raises(ValueError, match="Species cannot be empty"):
        client.search(species="   ", query="brain RNA-seq")


def test_geo_request_errors_are_wrapped() -> None:
    """Requests errors should become client-specific exceptions."""
    session = FakeSession(
        [
            FakeResponse(
                {},
                status_error=requests.HTTPError("500 Server Error"),
            )
        ]
    )

    client = NCBIGEOClient(session=session)  # type: ignore[arg-type]

    with pytest.raises(NCBIClientError, match="NCBI request failed"):
        client.search(
            species="Drosophila melanogaster",
            query="brain RNA-seq",
        )
