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

    client = NCBIGEOClient(
        session=session,  # type: ignore[arg-type]
        max_attempts=1,
        retry_delay=0,
    )

    with pytest.raises(NCBIClientError, match="NCBI request failed"):
        client.search(
            species="Drosophila melanogaster",
            query="brain RNA-seq",
        )


def test_geo_retries_transient_connection_failure() -> None:
    class FakeResponse:
        status_code = 200

        def raise_for_status(self) -> None:
            return None

        def json(self) -> dict:
            return {
                "esearchresult": {
                    "idlist": [],
                }
            }

    class FlakySession:
        def __init__(self) -> None:
            self.calls = 0

        def get(self, *args, **kwargs):
            del args, kwargs
            self.calls += 1

            if self.calls == 1:
                raise requests.ConnectionError(
                    "Response ended prematurely"
                )

            return FakeResponse()

    session = FlakySession()

    client = NCBIGEOClient(
        session=session,
        max_attempts=3,
        retry_delay=0,
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="orb2",
        max_results=3,
    )

    assert records == []
    assert session.calls == 2


def test_geo_retries_transient_http_status() -> None:
    class FakeResponse:
        def __init__(
            self,
            status_code: int,
            payload: dict | None = None,
        ) -> None:
            self.status_code = status_code
            self._payload = payload or {}

        def raise_for_status(self) -> None:
            if self.status_code >= 400:
                raise requests.HTTPError(
                    f"HTTP {self.status_code}",
                    response=self,
                )

        def json(self) -> dict:
            return self._payload

    class FlakySession:
        def __init__(self) -> None:
            self.calls = 0

        def get(self, *args, **kwargs):
            del args, kwargs
            self.calls += 1

            if self.calls == 1:
                return FakeResponse(503)

            return FakeResponse(
                200,
                {
                    "esearchresult": {
                        "idlist": [],
                    }
                },
            )

    session = FlakySession()

    client = NCBIGEOClient(
        session=session,
        max_attempts=3,
        retry_delay=0,
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="orb2",
        max_results=3,
    )

    assert records == []
    assert session.calls == 2
