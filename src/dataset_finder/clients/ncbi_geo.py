"""NCBI GEO DataSets client."""

from __future__ import annotations

import time
from typing import Any

import requests

from dataset_finder.models import DatasetRecord

EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_ACCESSION_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


class NCBIClientError(RuntimeError):
    """Raised when an NCBI request or response cannot be processed."""


# Backward-compatible name for code that imported GEORecord directly.
GEORecord = DatasetRecord


class NCBIGEOClient:
    """Search the NCBI GEO DataSets database."""

    def __init__(
        self,
        *,
        timeout: float = 30.0,
        tool: str = "dataset-finder",
        email: str | None = None,
        session: requests.Session | None = None,
        max_attempts: int = 4,
        retry_delay: float = 0.75,
    ) -> None:
        self.timeout = timeout
        self.tool = tool
        self.email = email
        self.session = session or requests.Session()
        self.max_attempts = max_attempts
        self.retry_delay = retry_delay

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int = 20,
    ) -> list[DatasetRecord]:
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
        """Request NCBI JSON with retries for transient failures."""
        request_params = dict(params)
        request_params["tool"] = self.tool

        if self.email:
            request_params["email"] = self.email

        last_error: Exception | None = None

        for attempt in range(1, self.max_attempts + 1):
            try:
                response = self.session.get(
                    f"{EUTILS_BASE_URL}/{endpoint}",
                    params=request_params,
                    timeout=self.timeout,
                )

                status_code = getattr(
                    response,
                    "status_code",
                    200,
                )

                if status_code in {
                    429,
                    500,
                    502,
                    503,
                    504,
                }:
                    raise requests.HTTPError(
                        f"Transient NCBI HTTP {status_code}",
                        response=response,
                    )

                response.raise_for_status()
                payload = response.json()

                if not isinstance(payload, dict):
                    raise NCBIClientError(
                        "NCBI returned an unexpected JSON response."
                    )

                return payload

            except (
                requests.ConnectionError,
                requests.Timeout,
                requests.HTTPError,
                requests.exceptions.ChunkedEncodingError,
            ) as exc:
                last_error = exc

                if attempt >= self.max_attempts:
                    break

                time.sleep(
                    self.retry_delay * (2 ** (attempt - 1))
                )

            except ValueError as exc:
                raise NCBIClientError(
                    "NCBI returned invalid JSON."
                ) from exc

            except requests.RequestException as exc:
                last_error = exc
                break

        raise NCBIClientError(
            f"NCBI request failed after "
            f"{self.max_attempts} attempts: {last_error}"
        )


    @staticmethod
    def _normalize_summary(summary: dict[str, Any]) -> DatasetRecord:
        accession = str(summary.get("accession", "")).strip()
        organisms = summary.get("taxon", [])
        study_types = summary.get("gdstype", [])

        organism = _join_values(organisms)
        study_type = _join_values(study_types)
        sample_count = _optional_integer(summary.get("n_samples"))

        return DatasetRecord(
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
