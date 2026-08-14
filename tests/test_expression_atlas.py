"""Tests for the Expression Atlas client."""

from __future__ import annotations

from typing import Any

import pytest
import requests

from dataset_finder.clients.expression_atlas import (
    ExpressionAtlasClient,
    ExpressionAtlasClientError,
)


class FakeResponse:
    """Minimal JSON response replacement."""

    def __init__(
        self,
        payload: dict[str, Any],
        *,
        error: requests.RequestException | None = None,
    ) -> None:
        self.payload = payload
        self.error = error

    def raise_for_status(self) -> None:
        if self.error:
            raise self.error

    def json(self) -> dict[str, Any]:
        return self.payload


class FakeSession:
    """Return deterministic Expression Atlas data."""

    def __init__(self, response: FakeResponse) -> None:
        self.response = response
        self.calls: list[dict[str, Any]] = []
        self.headers: dict[str, str] = {}

    def get(
        self,
        url: str,
        *,
        timeout: float,
    ) -> FakeResponse:
        self.calls.append(
            {
                "url": url,
                "timeout": timeout,
            }
        )
        return self.response


def test_expression_atlas_search_filters_locally() -> None:
    session = FakeSession(
        FakeResponse(
            {
                "experiments": [
                    {
                        "experimentAccession": "E-MTAB-ORB2",
                        "experimentDescription": (
                            "Orb2 RNA-seq in Drosophila brain"
                        ),
                        "species": "Drosophila melanogaster",
                        "lastUpdate": "01-08-2026",
                        "rawExperimentType": (
                            "RNASEQ_MRNA_DIFFERENTIAL"
                        ),
                        "technologyType": [
                            "RNA-Seq mRNA"
                        ],
                        "numberOfAssays": 12,
                        "experimentalFactors": [
                            "organism part",
                            "developmental stage",
                        ],
                        "experimentType": "Differential",
                    },
                    {
                        "experimentAccession": "E-MTAB-HUMAN",
                        "experimentDescription": (
                            "Orb2-like human experiment"
                        ),
                        "species": "Homo sapiens",
                        "technologyType": ["RNA-Seq mRNA"],
                    },
                    {
                        "experimentAccession": "E-MTAB-OTHER",
                        "experimentDescription": (
                            "Unrelated fly experiment"
                        ),
                        "species": "Drosophila melanogaster",
                        "technologyType": ["RNA-Seq mRNA"],
                    },
                ]
            }
        )
    )

    client = ExpressionAtlasClient(
        session=session,  # type: ignore[arg-type]
    )

    records = client.search(
        species="Drosophila melanogaster",
        query="orb2",
        max_results=5,
    )

    assert len(records) == 1

    record = records[0]

    assert record.database == "Expression Atlas"
    assert record.accession == "E-MTAB-ORB2"
    assert record.organism == "Drosophila melanogaster"
    assert record.sample_count == 12
    assert "RNA-Seq mRNA" in record.study_type
    assert record.url.endswith(
        "/E-MTAB-ORB2/Results"
    )


def test_expression_atlas_wraps_request_errors() -> None:
    client = ExpressionAtlasClient(
        session=FakeSession(
            FakeResponse(
                {},
                error=requests.HTTPError("503"),
            )
        ),  # type: ignore[arg-type]
    )

    with pytest.raises(
        ExpressionAtlasClientError,
        match="Expression Atlas search failed",
    ):
        client.search(
            species="Drosophila melanogaster",
            query="orb2",
        )
