"""EMBL-EBI Expression Atlas search client."""

from __future__ import annotations

import re
from typing import Any
from urllib.parse import quote

import requests

from dataset_finder.models import DatasetRecord

EXPRESSION_ATLAS_EXPERIMENTS_URL = (
    "https://www.ebi.ac.uk/gxa/json/experiments"
)

EXPRESSION_ATLAS_EXPERIMENT_URL = (
    "https://www.ebi.ac.uk/gxa/experiments"
)


class ExpressionAtlasClientError(RuntimeError):
    """Raised when an Expression Atlas request fails."""


class ExpressionAtlasClient:
    """Search Expression Atlas experiment metadata."""

    def __init__(
        self,
        *,
        timeout: float = 45.0,
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
        """Search Expression Atlas and locally filter experiments."""
        species = species.strip()
        query = query.strip()

        if not species:
            raise ValueError("Species cannot be empty.")

        if not query:
            raise ValueError("Query cannot be empty.")

        if max_results < 1:
            raise ValueError(
                "Maximum results must be greater than zero."
            )

        try:
            response = self.session.get(
                EXPRESSION_ATLAS_EXPERIMENTS_URL,
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
        except (
            requests.RequestException,
            ValueError,
        ) as exc:
            raise ExpressionAtlasClientError(
                f"Expression Atlas search failed: {exc}"
            ) from exc

        experiments = payload.get("experiments", [])

        if not isinstance(experiments, list):
            return []

        records: list[DatasetRecord] = []

        for experiment in experiments:
            if not isinstance(experiment, dict):
                continue

            if not self._matches(
                experiment,
                species=species,
                query=query,
            ):
                continue

            record = self._normalize_experiment(experiment)

            if record.accession:
                records.append(record)

            if len(records) >= max_results:
                break

        return self._deduplicate(records)

    @classmethod
    def _matches(
        cls,
        experiment: dict[str, Any],
        *,
        species: str,
        query: str,
    ) -> bool:
        """Require both organism and query evidence."""
        experiment_species = str(
            experiment.get("species", "")
        ).strip()

        if experiment_species.casefold() != species.casefold():
            return False

        searchable_text = " ".join(
            cls._flatten_text(
                {
                    "accession": experiment.get(
                        "experimentAccession",
                        "",
                    ),
                    "description": experiment.get(
                        "experimentDescription",
                        "",
                    ),
                    "species": experiment_species,
                    "technology": experiment.get(
                        "technologyType",
                        [],
                    ),
                    "factors": experiment.get(
                        "experimentalFactors",
                        [],
                    ),
                    "projects": experiment.get(
                        "experimentProjects",
                        [],
                    ),
                    "type": experiment.get(
                        "experimentType",
                        "",
                    ),
                    "raw_type": experiment.get(
                        "rawExperimentType",
                        "",
                    ),
                }
            )
        )

        return cls._contains_query(
            searchable_text,
            query,
        )

    @staticmethod
    def _contains_query(
        text: str,
        query: str,
    ) -> bool:
        """Match a direct query or Boolean query components."""
        text = text.casefold()

        quoted_terms = re.findall(
            r'"([^"]+)"',
            query,
        )

        if quoted_terms:
            return any(
                term.casefold() in text
                for term in quoted_terms
            )

        tokens = [
            token
            for token in re.findall(
                r"[A-Za-z0-9_.-]+",
                query,
            )
            if token.casefold()
            not in {
                "and",
                "or",
                "not",
                "all",
                "fields",
                "organism",
            }
        ]

        return any(
            token.casefold() in text
            for token in tokens
            if len(token) >= 2
        )

    @classmethod
    def _normalize_experiment(
        cls,
        experiment: dict[str, Any],
    ) -> DatasetRecord:
        """Convert one Expression Atlas experiment."""
        accession = str(
            experiment.get(
                "experimentAccession",
                "",
            )
        ).strip()

        description = str(
            experiment.get(
                "experimentDescription",
                "",
            )
        ).strip()

        species = str(
            experiment.get("species", "")
        ).strip()

        technology_types = cls._string_values(
            experiment.get("technologyType")
        )

        experimental_factors = cls._string_values(
            experiment.get("experimentalFactors")
        )

        projects = cls._string_values(
            experiment.get("experimentProjects")
        )

        experiment_type = str(
            experiment.get("experimentType", "")
        ).strip()

        raw_experiment_type = str(
            experiment.get("rawExperimentType", "")
        ).strip()

        study_type = "; ".join(
            value
            for value in (
                *technology_types,
                experiment_type,
            )
            if value
        )

        assay_count = cls._optional_integer(
            experiment.get("numberOfAssays")
        )

        publication_date = str(
            experiment.get(
                "lastUpdate",
                experiment.get("loadDate", ""),
            )
        ).strip()

        evidence_text = " ".join(
            value
            for value in (
                accession,
                description,
                species,
                study_type,
                raw_experiment_type,
                "; ".join(experimental_factors),
                "; ".join(projects),
            )
            if value
        )

        return DatasetRecord(
            uid=accession,
            accession=accession,
            title=description or accession,
            organism=species,
            study_type=(
                study_type
                or raw_experiment_type
                or "Expression Atlas"
            ),
            sample_count=assay_count,
            publication_date=publication_date,
            url=(
                f"{EXPRESSION_ATLAS_EXPERIMENT_URL}/"
                f"{quote(accession)}/Results"
            ),
            database="Expression Atlas",
            project_accession=accession,
            description=description,
            evidence_text=evidence_text,
            raw_metadata=experiment,
        )

    @classmethod
    def _flatten_text(
        cls,
        value: Any,
    ) -> list[str]:
        """Flatten nested Atlas values into strings."""
        if value is None:
            return []

        if isinstance(value, dict):
            fragments: list[str] = []

            for key, item in value.items():
                fragments.extend(
                    cls._flatten_text(key)
                )
                fragments.extend(
                    cls._flatten_text(item)
                )

            return fragments

        if isinstance(value, (list, tuple, set)):
            fragments = []

            for item in value:
                fragments.extend(
                    cls._flatten_text(item)
                )

            return fragments

        text = str(value).strip()
        return [text] if text else []

    @staticmethod
    def _string_values(
        value: Any,
    ) -> tuple[str, ...]:
        """Normalize optional list-like metadata."""
        if isinstance(value, list):
            values = value
        elif value in (None, ""):
            values = []
        else:
            values = [value]

        result: list[str] = []
        seen: set[str] = set()

        for item in values:
            normalized = str(item).strip()

            if not normalized:
                continue

            identity = normalized.casefold()

            if identity in seen:
                continue

            seen.add(identity)
            result.append(normalized)

        return tuple(result)

    @staticmethod
    def _optional_integer(
        value: Any,
    ) -> int | None:
        """Convert an optional value to an integer."""
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _deduplicate(
        records: list[DatasetRecord],
    ) -> list[DatasetRecord]:
        """Deduplicate records by experiment accession."""
        result: list[DatasetRecord] = []
        seen: set[str] = set()

        for record in records:
            identity = record.accession.upper()

            if not identity or identity in seen:
                continue

            seen.add(identity)
            result.append(record)

        return result
