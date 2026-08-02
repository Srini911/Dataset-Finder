"""Reusable NCBI Entrez E-utilities client."""

from __future__ import annotations

import os
import time
from typing import Any

import requests

EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class NCBIEntrezError(RuntimeError):
    """Raised when an NCBI Entrez request fails."""


class NCBIEntrezClient:
    """Make retry-aware requests to NCBI E-utilities."""

    def __init__(
        self,
        *,
        email: str | None = None,
        api_key: str | None = None,
        timeout: float = 45.0,
        retries: int = 4,
        session: requests.Session | None = None,
    ) -> None:
        self.email = email or os.environ.get(
            "NCBI_EMAIL",
            "dataset-finder@example.org",
        )
        self.api_key = api_key or os.environ.get(
            "NCBI_API_KEY",
            "",
        )
        self.timeout = timeout
        self.retries = retries
        self.session = session or requests.Session()
        self.session.headers.update(
            {
                "Accept": "application/json",
                "User-Agent": (
                    "Dataset-Finder/0.3.0 "
                    "(https://github.com/Srini911/Dataset-Finder)"
                ),
            }
        )

    def request_json(
        self,
        endpoint: str,
        *,
        parameters: dict[str, Any],
    ) -> dict[str, Any]:
        """Request and decode an NCBI E-utilities JSON response."""
        request_parameters: dict[str, Any] = {
            "tool": "dataset_finder",
            "email": self.email,
            **parameters,
        }

        if self.api_key:
            request_parameters["api_key"] = self.api_key

        last_error: Exception | None = None

        for attempt in range(1, self.retries + 1):
            try:
                response = self.session.get(
                    f"{EUTILS_BASE_URL}/{endpoint}",
                    params=request_parameters,
                    timeout=self.timeout,
                )
                response.raise_for_status()

                payload = response.json()

                if not self.api_key:
                    time.sleep(0.34)

                return payload
            except (
                requests.RequestException,
                ValueError,
            ) as exc:
                last_error = exc

                if attempt < self.retries:
                    time.sleep(min(2 ** (attempt - 1), 8))

        raise NCBIEntrezError(
            f"NCBI request failed after {self.retries} attempts: "
            f"{last_error}"
        )

    def search_ids(
        self,
        *,
        database: str,
        term: str,
        max_results: int,
    ) -> list[str]:
        """Search an Entrez database and return matching UIDs."""
        payload = self.request_json(
            "esearch.fcgi",
            parameters={
                "db": database,
                "term": term,
                "retmode": "json",
                "retmax": max_results,
            },
        )

        result = payload.get("esearchresult", {})
        identifiers = result.get("idlist", [])

        if not isinstance(identifiers, list):
            return []

        return [
            str(identifier)
            for identifier in identifiers
            if identifier
        ]

    def summaries(
        self,
        *,
        database: str,
        identifiers: list[str],
    ) -> list[dict[str, Any]]:
        """Retrieve document summaries for Entrez UIDs."""
        if not identifiers:
            return []

        payload = self.request_json(
            "esummary.fcgi",
            parameters={
                "db": database,
                "id": ",".join(identifiers),
                "retmode": "json",
            },
        )

        result = payload.get("result", {})
        ordered_uids = result.get("uids", identifiers)

        return [
            result[uid]
            for uid in ordered_uids
            if uid in result and isinstance(result[uid], dict)
        ]
