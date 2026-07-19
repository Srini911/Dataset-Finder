"""Client for searching experiments in the ENCODE portal."""

from __future__ import annotations

from typing import Any

import requests

from dataset_finder.models import DatasetRecord


class ENCODEClientError(RuntimeError):
    """Raised when the ENCODE API cannot be queried successfully."""


class ENCODEClient:
    """Search released experiments through the ENCODE REST API."""

    BASE_URL = "https://www.encodeproject.org"

    def __init__(
        self,
        *,
        timeout: float = 30.0,
        session: requests.Session | None = None,
    ) -> None:
        self.timeout = timeout
        self.session = session or requests.Session()
        self.session.headers.update(
            {
                "Accept": "application/json",
                "User-Agent": (
                    "Dataset-Finder/0.1.0 "
                    "(https://github.com/Srini911/Dataset-Finder)"
                ),
            }
        )

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int = 20,
    ) -> list[DatasetRecord]:
        """Search released ENCODE experiments and normalize the results."""
        params = [
            ("type", "Experiment"),
            ("status", "released"),
            (
                "replicates.library.biosample.donor.organism.scientific_name",
                species,
            ),
            ("searchTerm", query),
            ("limit", str(max_results)),
            ("format", "json"),
        ]

        try:
            response = self.session.get(
                f"{self.BASE_URL}/search/",
                params=params,
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
        except requests.RequestException as exc:
            raise ENCODEClientError(
                f"Unable to search ENCODE: {exc}"
            ) from exc
        except ValueError as exc:
            raise ENCODEClientError(
                "ENCODE returned an invalid JSON response."
            ) from exc

        experiments = payload.get("@graph", [])

        if not isinstance(experiments, list):
            raise ENCODEClientError(
                "ENCODE returned an unexpected search response."
            )

        records: list[DatasetRecord] = []

        for experiment in experiments[:max_results]:
            if not isinstance(experiment, dict):
                continue

            records.append(
                self._to_dataset_record(
                    experiment=experiment,
                    requested_species=species,
                )
            )

        return records

    def _to_dataset_record(
        self,
        *,
        experiment: dict[str, Any],
        requested_species: str,
    ) -> DatasetRecord:
        """Convert an ENCODE experiment into a common dataset record."""
        accession = self._text(experiment.get("accession"))
        assay = self._text(
            experiment.get("assay_title")
            or experiment.get("assay_term_name")
        )
        biosample = self._text(experiment.get("biosample_summary"))
        target = self._target_label(experiment.get("target"))

        title_parts = [
            value
            for value in (target, assay, biosample)
            if value
        ]
        title = " | ".join(title_parts)

        if not title:
            title = self._text(experiment.get("description"))

        date_released = self._text(
            experiment.get("date_released")
            or experiment.get("date_created")
        )

        replicates = experiment.get("replicates")
        sample_count = (
            len(replicates)
            if isinstance(replicates, list)
            else None
        )

        experiment_path = self._text(experiment.get("@id"))
        url = (
            f"{self.BASE_URL}{experiment_path}"
            if experiment_path.startswith("/")
            else (
                f"{self.BASE_URL}/experiments/{accession}/"
                if accession
                else self.BASE_URL
            )
        )

        return DatasetRecord(
            uid=accession or experiment_path,
            accession=accession,
            title=title,
            organism=requested_species,
            study_type=assay,
            sample_count=sample_count,
            publication_date=date_released,
            url=url,
        )

    @staticmethod
    def _target_label(target: Any) -> str:
        """Extract a readable target label from an ENCODE target value."""
        if isinstance(target, dict):
            value = target.get("label") or target.get("name")
            return ENCODEClient._text(value)

        if isinstance(target, str):
            parts = [part for part in target.split("/") if part]
            return parts[-1] if parts else target

        return ""

    @staticmethod
    def _text(value: Any) -> str:
        """Return a stripped string or an empty value."""
        return value.strip() if isinstance(value, str) else ""
