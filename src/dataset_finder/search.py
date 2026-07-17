"""Dataset search orchestration."""

from __future__ import annotations

from dataset_finder.clients.ncbi_geo import NCBIGEOClient
from dataset_finder.models import DatasetRecord


class UnsupportedDatabaseError(ValueError):
    """Raised when a requested database is not yet supported."""


class SearchService:
    """Coordinate searches across supported public databases."""

    def __init__(
        self,
        *,
        geo_client: NCBIGEOClient | None = None,
    ) -> None:
        self.geo_client = geo_client or NCBIGEOClient()

    def search(
        self,
        *,
        species: str,
        query: str,
        database: str = "geo",
        max_results: int = 20,
    ) -> list[DatasetRecord]:
        """Search the requested database and return normalized records."""
        species = species.strip()
        query = query.strip()
        database = database.strip().lower()

        if not species:
            raise ValueError("Species cannot be empty.")

        if not query:
            raise ValueError("Query cannot be empty.")

        if max_results < 1:
            raise ValueError("Maximum results must be greater than zero.")

        if database == "sra":
            raise UnsupportedDatabaseError("SRA search is not implemented yet.")

        if database not in {"geo", "all"}:
            raise UnsupportedDatabaseError(
                f"Unsupported database: {database}"
            )

        return self.geo_client.search(
            species=species,
            query=query,
            max_results=max_results,
        )
