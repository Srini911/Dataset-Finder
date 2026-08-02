"""Tests for the PRIDE Archive client."""

from __future__ import annotations

from typing import Any

import pytest
import requests

from dataset_finder.clients.pride import (
    PRIDEClient,
    PRIDEClientError,
)


class FakeResponse:
    """Minimal JSON response replacement."""

    def __init__(
        self,
        payload: Any,
        *,
        error: requests.RequestException | None = None,
    ) -> None:
        self.payload = payload
        self.error = error

    def raise_for_status(self) -> None:
        if self.error:
            raise self.error

    def json(self) -> Any:
        return self.payload


class FakeSession:
    """Return deterministic PRIDE data."""

    def __init__(self, response: FakeResponse) -> None:
        self.response = response
        self.calls: list[dict[str, Any]] = []
        self.headers: dict[str, str] = {}

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


def test_pride_search_filters_and_normalizes_projects() -> None:
    session = FakeSession(
        FakeResponse(
            [
                {
                    "accession": "PXD123456",
                    "title": "TARDBP proteomics in human neurons",
                    "projectDescription": (
                        "Proteome analysis of TARDBP depletion"
                    ),
                    "sampleProcessingProtocol": (
                        "Human neuronal cells"
                    ),
                    "dataProcessingProtocol": (
                        "Database searching"
                    ),
                    "publicationDate": "2026-01-01",
                    "doi": "10.6019/PXD123456",
                    "organisms": [
                        {
                            "name": "Homo sapiens (human)"
                        }
                    ],
                    "organismParts": [
                        {
                            "name": "brain"
                        }
                    ],
                    "diseases": [
                        {
                            "name": (
                                "amyotrophic lateral sclerosis"
                            )
                        }
                    ],
                    "instruments": [
                        {
                            "name": "Q Exactive"
                        }
                    ],
                    "experimentTypes": [
                        {
                            "name": "Shotgun proteomics"
                        }
                    ],
                    "references": [
                        {
                            "pubmedID": 12345678,
                            "doi": "10.1000/test",
                            "referenceLine": (
                                "Example TARDBP publication"
                            ),
                        }
                    ],
                },
                {
                    "accession": "PXD999999",
                    "title": "Unrelated mouse study",
                    "projectDescription": "Unrelated",
                    "organisms": [
                        {
                            "name": "Mus musculus"
                        }
                    ],
                },
            ]
        )
    )

    client = PRIDEClient(
        session=session,  # type: ignore[arg-type]
    )

    records = client.search(
        species="Homo sapiens",
        query="TARDBP",
        max_results=5,
    )

    assert len(records) == 1

    record = records[0]

    assert record.database == "PRIDE"
    assert record.accession == "PXD123456"
    assert record.organism == "Homo sapiens (human)"
    assert record.platform == "Q Exactive"
    assert record.tissue == "brain"
    assert record.disease == (
        "amyotrophic lateral sclerosis"
    )
    assert record.pubmed_ids == ("12345678",)
    assert "10.6019/pxd123456" in record.dois
    assert "PXD123456" in record.related_accessions


def test_pride_wraps_request_errors() -> None:
    client = PRIDEClient(
        session=FakeSession(
            FakeResponse(
                [],
                error=requests.HTTPError("503"),
            )
        ),  # type: ignore[arg-type]
    )

    with pytest.raises(
        PRIDEClientError,
        match="PRIDE search failed",
    ):
        client.search(
            species="Homo sapiens",
            query="TARDBP",
        )
