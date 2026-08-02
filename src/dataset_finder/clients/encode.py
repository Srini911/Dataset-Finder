"""Client for searching experiments in the ENCODE portal."""

from __future__ import annotations

import re
from typing import Any

import requests

from dataset_finder.models import DatasetRecord


class ENCODEClientError(RuntimeError):
    """Raised when the ENCODE API cannot be queried successfully."""


class ENCODEClient:
    """Search released experiments through the ENCODE REST API."""

    BASE_URL = "https://www.encodeproject.org"
    SUPPORTED_SPECIES = {
        "homo sapiens",
        "mus musculus",
    }

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
                    "Dataset-Finder/0.3.0 "
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
        """Search released ENCODE experiments and normalize results."""
        if species.strip().casefold() not in self.SUPPORTED_SPECIES:
            return []

        simple_terms = self._extract_search_terms(query)

        if not simple_terms:
            return []

        candidates: list[DatasetRecord] = []
        seen: set[str] = set()

        for search_term in simple_terms:
            experiments = self._request_experiments(
                search_term=search_term,
                max_results=max_results,
            )

            for experiment in experiments:
                record = self._to_dataset_record(
                    experiment=experiment,
                    requested_species=species,
                )

                if not self._matches_species(
                    experiment,
                    requested_species=species,
                ):
                    continue

                identity = (
                    record.accession
                    or record.uid
                    or record.url
                ).strip().upper()

                if identity in seen:
                    continue

                seen.add(identity)
                candidates.append(record)

                if len(candidates) >= max_results:
                    return candidates

        return candidates

    def _request_experiments(
        self,
        *,
        search_term: str,
        max_results: int,
    ) -> list[dict[str, Any]]:
        """Request released experiments for one simple search term."""
        params = [
            ("type", "Experiment"),
            ("status", "released"),
            ("searchTerm", search_term),
            ("limit", str(max_results)),
            ("frame", "object"),
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

        return [
            experiment
            for experiment in experiments
            if isinstance(experiment, dict)
        ]

    @staticmethod
    def _extract_search_terms(query: str) -> list[str]:
        """Extract simple terms from a FlyBase-enhanced query."""
        quoted_terms = re.findall(
            r'"([^"]+)"',
            query,
        )

        if quoted_terms:
            raw_terms = quoted_terms
        else:
            raw_terms = re.split(
                r"\s+OR\s+",
                query,
                flags=re.IGNORECASE,
            )

        terms: list[str] = []
        seen: set[str] = set()

        for raw_term in raw_terms:
            term = raw_term.strip().strip('"')

            if not term:
                continue

            identity = term.casefold()

            if identity in seen:
                continue

            seen.add(identity)
            terms.append(term)

        return terms

    @classmethod
    def _matches_species(
        cls,
        experiment: dict[str, Any],
        *,
        requested_species: str,
    ) -> bool:
        """Check returned experiment metadata for the requested species."""
        requested = requested_species.strip().casefold()

        if not requested:
            return True

        metadata_text = cls._metadata_text(experiment).casefold()

        if requested in metadata_text:
            return True

        common_names = {
            "homo sapiens": ("human",),
            "mus musculus": ("mouse", "mice"),
            "drosophila melanogaster": (
                "drosophila",
                "fruit fly",
            ),
            "caenorhabditis elegans": (
                "c. elegans",
                "caenorhabditis",
            ),
        }

        return any(
            alias in metadata_text
            for alias in common_names.get(
                requested,
                (),
            )
        )

    @classmethod
    def _metadata_text(
        cls,
        experiment: dict[str, Any],
    ) -> str:
        """Flatten selected ENCODE metadata into searchable text."""
        values: list[str] = []

        for key in (
            "accession",
            "assay_title",
            "assay_term_name",
            "biosample_summary",
            "description",
            "organism",
            "target",
            "replicates",
        ):
            values.extend(
                cls._flatten_values(
                    experiment.get(key)
                )
            )

        return " ".join(values)

    @classmethod
    def _flatten_values(
        cls,
        value: Any,
    ) -> list[str]:
        """Flatten nested metadata into strings."""
        if isinstance(value, str):
            return [value]

        if isinstance(value, dict):
            flattened: list[str] = []

            for nested_value in value.values():
                flattened.extend(
                    cls._flatten_values(nested_value)
                )

            return flattened

        if isinstance(value, list):
            flattened = []

            for item in value:
                flattened.extend(
                    cls._flatten_values(item)
                )

            return flattened

        return []

    def _to_dataset_record(
        self,
        *,
        experiment: dict[str, Any],
        requested_species: str,
    ) -> DatasetRecord:
        """Convert an ENCODE experiment into a common dataset record."""
        accession = self._text(
            experiment.get("accession")
        )
        assay = self._text(
            experiment.get("assay_title")
            or experiment.get("assay_term_name")
        )
        biosample = self._text(
            experiment.get("biosample_summary")
        )
        target = self._target_label(
            experiment.get("target")
        )

        title_parts = [
            value
            for value in (
                target,
                assay,
                biosample,
            )
            if value
        ]
        title = " | ".join(title_parts)

        if not title:
            title = self._text(
                experiment.get("description")
            )

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

        experiment_path = self._text(
            experiment.get("@id")
        )
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
            database="ENCODE",
            description=self._text(
                experiment.get("description")
            ),
            evidence_text=self._metadata_text(
                experiment
            ),
            raw_metadata=experiment,
        )

    @staticmethod
    def _target_label(target: Any) -> str:
        """Extract a readable target label."""
        if isinstance(target, dict):
            value = (
                target.get("label")
                or target.get("name")
            )
            return ENCODEClient._text(value)

        if isinstance(target, str):
            parts = [
                part
                for part in target.split("/")
                if part
            ]
            return (
                parts[-1]
                if parts
                else target
            )

        return ""

    @staticmethod
    def _text(value: Any) -> str:
        """Return a stripped string or an empty value."""
        return (
            value.strip()
            if isinstance(value, str)
            else ""
        )
