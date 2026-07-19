"""Tests for the ENCODE API client."""

from __future__ import annotations

from typing import Any

import pytest
import requests

from dataset_finder.clients.encode import ENCODEClient, ENCODEClientError


class FakeResponse:
    """Minimal fake HTTP response for ENCODE client tests."""

    def __init__(
        self,
        *,
        payload: Any = None,
        status_error: Exception | None = None,
        json_error: Exception | None = None,
    ) -> None:
        self.payload = payload
        self.status_error = status_error
        self.json_error = json_error

    def raise_for_status(self) -> None:
        """Raise the configured HTTP status error."""
        if self.status_error is not None:
            raise self.status_error

    def json(self) -> Any:
        """Return the configured JSON payload."""
        if self.json_error is not None:
            raise self.json_error

        return self.payload


class FakeSession:
    """Small requests-session replacement for offline tests."""

    def __init__(self, response: FakeResponse) -> None:
        self.response = response
        self.headers: dict[str, str] = {}
        self.requested_url: str | None = None
        self.requested_params: list[tuple[str, str]] | None = None
        self.requested_timeout: float | None = None

    def get(
        self,
        url: str,
        *,
        params: list[tuple[str, str]],
        timeout: float,
    ) -> FakeResponse:
        """Return the configured response and record request arguments."""
        self.requested_url = url
        self.requested_params = params
        self.requested_timeout = timeout
        return self.response


def test_encode_search_normalizes_experiment() -> None:
    """ENCODE experiments should become normalized dataset records."""
    response = FakeResponse(
        payload={
            "@graph": [
                {
                    "@id": "/experiments/ENCSR720BJU/",
                    "accession": "ENCSR720BJU",
                    "assay_title": "eCLIP",
                    "biosample_summary": "Homo sapiens K562",
                    "target": {"label": "TARDBP"},
                    "date_released": "2024-01-15",
                    "replicates": [{}, {}],
                }
            ]
        }
    )
    session = FakeSession(response)
    client = ENCODEClient(session=session, timeout=12.5)

    records = client.search(
        species="Homo sapiens",
        query="TARDBP",
        max_results=5,
    )

    assert len(records) == 1

    record = records[0]

    assert record.uid == "ENCSR720BJU"
    assert record.accession == "ENCSR720BJU"
    assert record.title == "TARDBP | eCLIP | Homo sapiens K562"
    assert record.organism == "Homo sapiens"
    assert record.study_type == "eCLIP"
    assert record.sample_count == 2
    assert record.publication_date == "2024-01-15"
    assert record.url == (
        "https://www.encodeproject.org/experiments/ENCSR720BJU/"
    )

    assert session.requested_url == (
        "https://www.encodeproject.org/search/"
    )
    assert session.requested_timeout == 12.5
    assert session.requested_params is not None
    assert ("type", "Experiment") in session.requested_params
    assert ("status", "released") in session.requested_params
    assert ("searchTerm", "TARDBP") in session.requested_params
    assert ("limit", "5") in session.requested_params


def test_encode_search_handles_empty_results() -> None:
    """An empty ENCODE graph should return an empty list."""
    session = FakeSession(FakeResponse(payload={"@graph": []}))
    client = ENCODEClient(session=session)

    records = client.search(
        species="Mus musculus",
        query="CTCF",
        max_results=5,
    )

    assert records == []


def test_encode_search_rejects_malformed_graph() -> None:
    """The client should reject malformed ENCODE search payloads."""
    session = FakeSession(
        FakeResponse(payload={"@graph": {"accession": "invalid"}})
    )
    client = ENCODEClient(session=session)

    with pytest.raises(
        ENCODEClientError,
        match="unexpected search response",
    ):
        client.search(
            species="Homo sapiens",
            query="TARDBP",
        )


def test_encode_search_handles_http_errors() -> None:
    """HTTP request failures should become ENCODE client errors."""
    session = FakeSession(
        FakeResponse(
            status_error=requests.HTTPError("503 Service Unavailable")
        )
    )
    client = ENCODEClient(session=session)

    with pytest.raises(
        ENCODEClientError,
        match="Unable to search ENCODE",
    ):
        client.search(
            species="Homo sapiens",
            query="TARDBP",
        )


def test_encode_search_handles_invalid_json() -> None:
    """Invalid JSON responses should become ENCODE client errors."""
    session = FakeSession(
        FakeResponse(json_error=ValueError("invalid JSON"))
    )
    client = ENCODEClient(session=session)

    with pytest.raises(
        ENCODEClientError,
        match="invalid JSON response",
    ):
        client.search(
            species="Homo sapiens",
            query="TARDBP",
        )


def test_encode_target_path_is_normalized() -> None:
    """A target path should be converted to a readable target name."""
    response = FakeResponse(
        payload={
            "@graph": [
                {
                    "@id": "/experiments/ENCSR123ABC/",
                    "accession": "ENCSR123ABC",
                    "assay_title": "TF ChIP-seq",
                    "biosample_summary": "Homo sapiens HepG2",
                    "target": "/targets/CTCF-human/",
                    "replicates": [],
                }
            ]
        }
    )
    client = ENCODEClient(session=FakeSession(response))

    records = client.search(
        species="Homo sapiens",
        query="CTCF",
        max_results=1,
    )

    assert records[0].title == (
        "CTCF-human | TF ChIP-seq | Homo sapiens HepG2"
    )
