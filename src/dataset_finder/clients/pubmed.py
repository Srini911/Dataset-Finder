"""NCBI PubMed search client."""

from __future__ import annotations

import urllib.parse
from typing import Any

from dataset_finder.clients.ncbi_entrez import NCBIEntrezClient
from dataset_finder.models import DatasetRecord
from dataset_finder.study_linking import extract_study_links


class PubMedClient:
    """Search PubMed and normalize publication records."""

    def __init__(
        self,
        *,
        entrez_client: NCBIEntrezClient | None = None,
    ) -> None:
        self.entrez = entrez_client or NCBIEntrezClient()

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int = 20,
    ) -> list[DatasetRecord]:
        """Search PubMed for publications related to a query."""
        term = (
            f"({query}) AND "
            f'("{species}"[Title/Abstract] OR '
            f'"{species}"[MeSH Terms])'
        )

        identifiers = self.entrez.search_ids(
            database="pubmed",
            term=term,
            max_results=max_results,
        )

        summaries = self.entrez.summaries(
            database="pubmed",
            identifiers=identifiers,
        )

        return [
            self._normalize_summary(summary)
            for summary in summaries
        ]

    @classmethod
    def _normalize_summary(
        cls,
        summary: dict[str, Any],
    ) -> DatasetRecord:
        """Convert a PubMed document summary to a dataset record."""
        pmid = str(summary.get("uid", "")).strip()
        title = str(summary.get("title", "")).strip()
        journal = str(
            summary.get(
                "fulljournalname",
                summary.get("source", ""),
            )
        ).strip()
        publication_date = str(
            summary.get(
                "pubdate",
                summary.get("sortpubdate", ""),
            )
        ).strip()

        authors = cls._authors(summary.get("authors"))
        article_ids = cls._article_ids(
            summary.get("articleids")
        )

        doi = (
            article_ids.get("doi", "")
            or article_ids.get("elocationid", "")
        )

        metadata_links = extract_study_links(
            title,
            summary,
        )

        publication = " | ".join(
            value
            for value in (
                authors,
                journal,
            )
            if value
        )

        return DatasetRecord(
            uid=pmid,
            accession=f"PMID{pmid}",
            title=title or f"PubMed record {pmid}",
            organism="",
            study_type="Publication",
            sample_count=None,
            publication_date=publication_date,
            url=(
                "https://pubmed.ncbi.nlm.nih.gov/"
                f"{urllib.parse.quote(pmid)}/"
            ),
            database="PubMed",
            description=journal,
            publication=publication,
            pubmed_ids=(pmid,) if pmid else (),
            dois=(doi.lower(),) if doi else metadata_links.dois,
            related_accessions=metadata_links.all_accessions,
            related_geo_accessions=(
                metadata_links.related_geo_accessions
            ),
            related_study_accessions=(
                metadata_links.related_study_accessions
            ),
            related_bioproject_accessions=(
                metadata_links.related_bioproject_accessions
            ),
            related_biosample_accessions=(
                metadata_links.related_biosample_accessions
            ),
            evidence_text=" ".join(
                value
                for value in (
                    title,
                    authors,
                    journal,
                    doi,
                )
                if value
            ),
            raw_metadata=summary,
        )

    @staticmethod
    def _authors(value: Any) -> str:
        """Convert PubMed author objects into readable text."""
        if not isinstance(value, list):
            return ""

        names: list[str] = []

        for author in value:
            if isinstance(author, dict):
                name = str(
                    author.get(
                        "name",
                        author.get("authtype", ""),
                    )
                ).strip()
            else:
                name = str(author).strip()

            if name:
                names.append(name)

        return "; ".join(names)

    @staticmethod
    def _article_ids(value: Any) -> dict[str, str]:
        """Normalize PubMed article identifier objects."""
        if not isinstance(value, list):
            return {}

        identifiers: dict[str, str] = {}

        for item in value:
            if not isinstance(item, dict):
                continue

            id_type = str(
                item.get("idtype", "")
            ).strip().casefold()
            identifier = str(
                item.get("value", "")
            ).strip()

            if id_type and identifier:
                identifiers[id_type] = identifier

        return identifiers
