"""PRIDE Archive proteomics project search client."""

from __future__ import annotations

import re
from typing import Any
from urllib.parse import quote

import requests

from dataset_finder.models import DatasetRecord
from dataset_finder.study_linking import (
    extract_study_links,
)

PRIDE_PROJECTS_URL = (
    "https://www.ebi.ac.uk/pride/ws/archive/v2/projects"
)

PRIDE_PROJECT_URL = (
    "https://www.ebi.ac.uk/pride/archive/projects"
)


class PRIDEClientError(RuntimeError):
    """Raised when a PRIDE request fails."""


class PRIDEClient:
    """Search PRIDE Archive proteomics projects."""

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
        """Search PRIDE and strictly filter returned projects."""
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

        candidate_limit = min(
            200,
            max(50, max_results * 10),
        )

        try:
            response = self.session.get(
                PRIDE_PROJECTS_URL,
                params={
                    "keyword": query,
                    "page": 0,
                    "pageSize": candidate_limit,
                },
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
        except (
            requests.RequestException,
            ValueError,
        ) as exc:
            raise PRIDEClientError(
                f"PRIDE search failed: {exc}"
            ) from exc

        projects = self._project_list(payload)

        records: list[DatasetRecord] = []

        for project in projects:
            if not isinstance(project, dict):
                continue

            if not self._matches(
                project,
                species=species,
                query=query,
            ):
                continue

            record = self._normalize_project(project)

            if record.accession:
                records.append(record)

            if len(records) >= max_results:
                break

        return self._deduplicate(records)

    @staticmethod
    def _project_list(
        payload: Any,
    ) -> list[Any]:
        """Support current and wrapped PRIDE responses."""
        if isinstance(payload, list):
            return payload

        if not isinstance(payload, dict):
            return []

        embedded = payload.get("_embedded", {})

        if isinstance(embedded, dict):
            for key in (
                "projects",
                "projectList",
                "compactprojects",
            ):
                value = embedded.get(key)

                if isinstance(value, list):
                    return value

        for key in (
            "projects",
            "content",
            "results",
        ):
            value = payload.get(key)

            if isinstance(value, list):
                return value

        return []

    @classmethod
    def _matches(
        cls,
        project: dict[str, Any],
        *,
        species: str,
        query: str,
    ) -> bool:
        """Require query evidence and compatible organism metadata."""
        organism_names = cls._named_values(
            project.get("organisms")
        )

        searchable_text = " ".join(
            cls._flatten_text(
                {
                    "accession": project.get(
                        "accession",
                        "",
                    ),
                    "title": project.get("title", ""),
                    "description": project.get(
                        "projectDescription",
                        "",
                    ),
                    "sample_processing": project.get(
                        "sampleProcessingProtocol",
                        "",
                    ),
                    "data_processing": project.get(
                        "dataProcessingProtocol",
                        "",
                    ),
                    "tags": project.get(
                        "projectTags",
                        [],
                    ),
                    "keywords": project.get(
                        "keywords",
                        [],
                    ),
                    "organisms": organism_names,
                    "organism_parts": project.get(
                        "organismParts",
                        [],
                    ),
                    "diseases": project.get(
                        "diseases",
                        [],
                    ),
                    "references": project.get(
                        "references",
                        [],
                    ),
                    "attributes": project.get(
                        "additionalAttributes",
                        [],
                    ),
                }
            )
        )

        if organism_names:
            species_matches = any(
                cls._organism_matches(
                    organism,
                    species,
                )
                for organism in organism_names
            )
        else:
            species_matches = (
                species.casefold()
                in searchable_text.casefold()
            )

        return (
            species_matches
            and cls._contains_query(
                searchable_text,
                query,
            )
        )

    @staticmethod
    def _organism_matches(
        organism: str,
        species: str,
    ) -> bool:
        """Compare organism names while ignoring parenthetical aliases."""
        normalized_organism = re.sub(
            r"\s*\([^)]*\)\s*",
            "",
            organism,
        ).strip().casefold()

        normalized_species = re.sub(
            r"\s*\([^)]*\)\s*",
            "",
            species,
        ).strip().casefold()

        return (
            normalized_organism == normalized_species
            or normalized_species in normalized_organism
            or normalized_organism in normalized_species
        )

    @staticmethod
    def _contains_query(
        text: str,
        query: str,
    ) -> bool:
        """Match direct or expanded Boolean query terms."""
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
    def _normalize_project(
        cls,
        project: dict[str, Any],
    ) -> DatasetRecord:
        """Convert one PRIDE project to a dataset record."""
        accession = str(
            project.get("accession", "")
        ).strip()

        title = str(
            project.get("title", "")
        ).strip()

        description = str(
            project.get("projectDescription", "")
        ).strip()

        sample_processing = str(
            project.get(
                "sampleProcessingProtocol",
                "",
            )
        ).strip()

        data_processing = str(
            project.get(
                "dataProcessingProtocol",
                "",
            )
        ).strip()

        organism_names = cls._named_values(
            project.get("organisms")
        )

        organism = "; ".join(organism_names)

        organism_parts = cls._named_values(
            project.get("organismParts")
        )

        diseases = cls._named_values(
            project.get("diseases")
        )

        instruments = cls._named_values(
            project.get("instruments")
        )

        experiment_types = cls._named_values(
            project.get("experimentTypes")
        )

        references = project.get("references", [])

        links = extract_study_links(
            project,
        )

        reference_pmids = cls._reference_values(
            references,
            "pubmedID",
        )

        reference_dois = tuple(
            value.casefold()
            for value in cls._reference_values(
                references,
                "doi",
            )
        )

        project_doi = str(
            project.get("doi", "")
        ).strip().casefold()

        pubmed_ids = cls._unique(
            (
                *reference_pmids,
                *links.pubmed_ids,
            )
        )

        dois = cls._unique(
            (
                project_doi,
                *reference_dois,
                *links.dois,
            )
        )

        publication_date = str(
            project.get(
                "publicationDate",
                project.get("submissionDate", ""),
            )
        ).strip()

        publication = cls._reference_lines(
            references
        )

        evidence_text = " ".join(
            value
            for value in (
                accession,
                title,
                description,
                sample_processing,
                data_processing,
                organism,
                "; ".join(organism_parts),
                "; ".join(diseases),
                "; ".join(instruments),
                "; ".join(experiment_types),
                publication,
            )
            if value
        )

        related_accessions = cls._unique(
            (
                accession,
                *links.study_level_accessions,
                *(
                    f"PMID{pmid}"
                    for pmid in pubmed_ids
                ),
            )
        )

        return DatasetRecord(
            uid=accession,
            accession=accession,
            title=title or accession,
            organism=organism,
            study_type=(
                "; ".join(experiment_types)
                or "Proteomics"
            ),
            sample_count=None,
            publication_date=publication_date,
            url=(
                f"{PRIDE_PROJECT_URL}/"
                f"{quote(accession)}"
            ),
            database="PRIDE",
            project_accession=accession,
            platform="; ".join(instruments),
            description=description,
            tissue="; ".join(organism_parts),
            disease="; ".join(diseases),
            publication=publication,
            pubmed_ids=pubmed_ids,
            dois=dois,
            related_accessions=related_accessions,
            related_geo_accessions=(
                links.related_geo_accessions
            ),
            related_study_accessions=(
                links.related_study_accessions
            ),
            related_bioproject_accessions=(
                links.related_bioproject_accessions
            ),
            related_biosample_accessions=(
                links.related_biosample_accessions
            ),
            evidence_text=evidence_text,
            raw_metadata=project,
        )

    @classmethod
    def _flatten_text(
        cls,
        value: Any,
    ) -> list[str]:
        """Flatten nested PRIDE values."""
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
    def _named_values(
        value: Any,
    ) -> tuple[str, ...]:
        """Extract unique name values from CV objects."""
        if not isinstance(value, list):
            return ()

        result: list[str] = []
        seen: set[str] = set()

        for item in value:
            if isinstance(item, dict):
                name = str(
                    item.get("name", "")
                ).strip()
            else:
                name = str(item).strip()

            if not name:
                continue

            identity = name.casefold()

            if identity in seen:
                continue

            seen.add(identity)
            result.append(name)

        return tuple(result)

    @staticmethod
    def _reference_values(
        references: Any,
        key: str,
    ) -> tuple[str, ...]:
        """Extract one field from PRIDE references."""
        if not isinstance(references, list):
            return ()

        values: list[str] = []

        for reference in references:
            if not isinstance(reference, dict):
                continue

            value = reference.get(key)

            if value not in (None, ""):
                values.append(str(value).strip())

        return PRIDEClient._unique(values)

    @staticmethod
    def _reference_lines(
        references: Any,
    ) -> str:
        """Return readable publication reference lines."""
        if not isinstance(references, list):
            return ""

        values: list[str] = []

        for reference in references:
            if not isinstance(reference, dict):
                continue

            line = str(
                reference.get(
                    "referenceLine",
                    "",
                )
            ).strip()

            if line:
                values.append(line)

        return " | ".join(
            PRIDEClient._unique(values)
        )

    @staticmethod
    def _unique(
        values: Any,
    ) -> tuple[str, ...]:
        """Deduplicate strings while preserving order."""
        result: list[str] = []
        seen: set[str] = set()

        for value in values:
            normalized = str(value).strip()

            if not normalized:
                continue

            identity = normalized.casefold()

            if identity in seen:
                continue

            seen.add(identity)
            result.append(normalized)

        return tuple(result)

    @staticmethod
    def _deduplicate(
        records: list[DatasetRecord],
    ) -> list[DatasetRecord]:
        """Deduplicate PRIDE projects."""
        result: list[DatasetRecord] = []
        seen: set[str] = set()

        for record in records:
            identity = record.accession.upper()

            if not identity or identity in seen:
                continue

            seen.add(identity)
            result.append(record)

        return result
