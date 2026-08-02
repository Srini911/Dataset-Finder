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


def test_parse_geo_sample_soft() -> None:
    from dataset_finder.clients.ncbi_geo import (
        parse_geo_sample_soft,
    )

    soft_text = """^SAMPLE = GSM8193871
!Sample_title = insc-GAL4>wRNAi brains, sample 1
!Sample_source_name_ch1 = brain
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = genotype: inscGAL4>wRNAi
!Sample_treatment_protocol_ch1 = RNAi depletion of white or orb2
!Sample_growth_protocol_ch1 = wandering L3 larvae
!Sample_instrument_model = Illumina NovaSeq 6000
!Sample_library_selection = cDNA
!Sample_library_source = transcriptomic
!Sample_library_strategy = RNA-Seq
!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN40876891
!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX24191199
^SAMPLE = GSM8193872
!Sample_title = insc-GAL4>orb2RNAi brains, sample 2
!Sample_source_name_ch1 = brain
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = genotype: inscGAL4>orb2RNAi
!Sample_treatment_protocol_ch1 = RNAi depletion of white or orb2
!Sample_growth_protocol_ch1 = wandering L3 larvae
!Sample_instrument_model = Illumina NovaSeq 6000
!Sample_library_selection = cDNA
!Sample_library_source = transcriptomic
!Sample_library_strategy = RNA-Seq
!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN40876892
!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX24191200
"""

    metadata = parse_geo_sample_soft(soft_text)

    assert metadata["sample_count"] == 2
    assert metadata["sample_accessions"] == [
        "GSM8193871",
        "GSM8193872",
    ]
    assert metadata["source_names"] == ["brain"]
    assert metadata["characteristics"]["tissue"] == ["brain"]
    assert metadata["characteristics"]["genotype"] == [
        "inscGAL4>wRNAi",
        "inscGAL4>orb2RNAi",
    ]
    assert metadata["library_strategies"] == ["RNA-Seq"]
    assert metadata["library_sources"] == ["transcriptomic"]
    assert metadata["library_selections"] == ["cDNA"]
    assert metadata["instrument_models"] == [
        "Illumina NovaSeq 6000"
    ]
    assert metadata["biosample_accessions"] == [
        "SAMN40876891",
        "SAMN40876892",
    ]
    assert metadata["experiment_accessions"] == [
        "SRX24191199",
        "SRX24191200",
    ]


def test_parse_geo_sample_soft_handles_empty_text() -> None:
    from dataset_finder.clients.ncbi_geo import (
        parse_geo_sample_soft,
    )

    metadata = parse_geo_sample_soft("")

    assert metadata["sample_count"] == 0
    assert metadata["samples"] == []
    assert metadata["characteristics"] == {}


def test_fetch_geo_sample_metadata() -> None:
    class SoftResponse:
        status_code = 200

        def __init__(self, text: str) -> None:
            self.text = text

        def raise_for_status(self) -> None:
            return None

    class SoftSession:
        def __init__(self) -> None:
            self.calls: list[dict[str, Any]] = []

        def get(
            self,
            url: str,
            *,
            params: dict[str, Any],
            timeout: float,
        ) -> SoftResponse:
            self.calls.append(
                {
                    "url": url,
                    "params": params,
                    "timeout": timeout,
                }
            )

            return SoftResponse(
                """^SAMPLE = GSM1
!Sample_title = adult brain sample
!Sample_source_name_ch1 = brain
!Sample_characteristics_ch1 = tissue: brain
!Sample_library_strategy = RNA-Seq
"""
            )

    session = SoftSession()
    client = NCBIGEOClient(
        session=session,  # type: ignore[arg-type]
    )

    metadata = client.fetch_sample_metadata("GSE123")

    assert metadata["sample_count"] == 1
    assert metadata["sample_accessions"] == ["GSM1"]
    assert metadata["source_names"] == ["brain"]
    assert metadata["characteristics"]["tissue"] == ["brain"]
    assert metadata["library_strategies"] == ["RNA-Seq"]

    assert session.calls[0]["params"] == {
        "acc": "GSE123",
        "targ": "gsm",
        "form": "text",
        "view": "full",
    }


def test_fetch_geo_sample_metadata_handles_empty_accession() -> None:
    client = NCBIGEOClient()

    metadata = client.fetch_sample_metadata("")

    assert metadata["sample_count"] == 0
