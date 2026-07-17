"""NCBI GEO DataSets client."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import requests

EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_ACCESSION_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


class NCBIClientError(RuntimeError):
    """Raised when an NCBI request or response cannot be processed."""


@dataclass(frozen=True, slots=True)
class GEORecord:
    """A normalized GEO search result."""

    uid: str
    accession: str
    title: str
    organism: str
    study_type: str
    sample_count: int | None
    publication_date: str
    url: str


class NCBIGEOClient:
    """Search the NCBI GEO DataSets database."""

    def __init__(
        self,
        *,
        timeout: float = 30.0,
        tool: str = "dataset-finder",
        email: str | None = None,
        session: requests.Session | None = None,
    ) -> None:
        self.timeout = timeout
        self.tool = tool
        self.email = email
        self.session = session or requests.Session()

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int = 20,
    ) -> list[GEORecord]:
        """Search GEO DataSets and return normalized records."""
        term = self._build_search_term(species=species, query=query)
        record_ids = self._search_ids(term=term, max_results=max_results)

        if not record_ids:
            return []

        summaries = self._fetch_summaries(record_ids)
        return [self._normalize_summary(summary) for summary in summaries]

    def _build_search_term(self, *, species: str, query: str) -> str:
        species = species.strip()
        query = query.strip()

        if not species:
            raise ValueError("Species cannot be empty.")

        if not query:
            raise ValueError("Query cannot be empty.")

        return f'("{species}"[Organism]) AND ({query}) AND "gse"[Entry Type]'

    def _search_ids(self, *, term: str, max_results: int) -> list[str]:
        payload = self._request_json(
            "esearch.fcgi",
            params={
                "db": "gds",
                "term": term,
                "retmode": "json",
                "retmax": max_results,
                "sort": "date",
            },
        )

        try:
            id_list = payload["esearchresult"]["idlist"]
        except (KeyError, TypeError) as exc:
            raise NCBIClientError("NCBI ESearch returned an unexpected response.") from exc

        return [str(record_id) for record_id in id_list]

    def _fetch_summaries(self, record_ids: list[str]) -> list[dict[str, Any]]:
        payload = self._request_json(
            "esummary.fcgi",
            params={
                "db": "gds",
                "id": ",".join(record_ids),
                "retmode": "json",
            },
        )

        try:
            result = payload["result"]
        except (KeyError, TypeError) as exc:
            raise NCBIClientError("NCBI ESummary returned an unexpected response.") from exc

        summaries: list[dict[str, Any]] = []

        for record_id in record_ids:
            summary = result.get(record_id)

            if isinstance(summary, dict):
                summaries.append(summary)

        return summaries

    def _request_json(
        self,
        endpoint: str,
        *,
        params: dict[str, str | int],
    ) -> dict[str, Any]:
        request_params = dict(params)
        request_params["tool"] = self.tool

        if self.email:
            request_params["email"] = self.email

        try:
            response = self.session.get(
                f"{EUTILS_BASE_URL}/{endpoint}",
                params=request_params,
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
        except requests.RequestException as exc:
            raise NCBIClientError(f"NCBI request failed: {exc}") from exc
        except ValueError as exc:
            raise NCBIClientError("NCBI returned invalid JSON.") from exc

        if not isinstance(payload, dict):
            raise NCBIClientError("NCBI returned an unexpected JSON response.")

        return payload

    @staticmethod
    def _normalize_summary(summary: dict[str, Any]) -> GEORecord:
        accession = str(summary.get("accession", "")).strip()
        organisms = summary.get("taxon", [])
        study_types = summary.get("gdstype", [])

        organism = _join_values(organisms)
        study_type = _join_values(study_types)
        sample_count = _optional_integer(summary.get("n_samples"))

        return GEORecord(
            uid=str(summary.get("uid", "")).strip(),
            accession=accession,
            title=str(summary.get("title", "")).strip(),
            organism=organism,
            study_type=study_type,
            sample_count=sample_count,
            publication_date=str(summary.get("pdat", "")).strip(),
            url=f"{GEO_ACCESSION_URL}?acc={accession}",
        )


def _join_values(value: Any) -> str:
    """Convert an NCBI list-like value into readable text."""
    if isinstance(value, list):
        return ", ".join(str(item).strip() for item in value if str(item).strip())

    return str(value or "").strip()


def _optional_integer(value: Any) -> int | None:
    """Convert an optional value to an integer."""
    if value in (None, ""):
        return None

    try:
        return int(value)
    except (TypeError, ValueError):
        return None
