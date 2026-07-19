"""Dataset search orchestration."""

from __future__ import annotations

from dataset_finder.clients.encode import ENCODEClient
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
        encode_client: ENCODEClient | None = None,
    ) -> None:
        self.geo_client = geo_client or NCBIGEOClient()
        self.encode_client = encode_client or ENCODEClient()

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

        if database == "geo":
            return self.geo_client.search(
                species=species,
                query=query,
                max_results=max_results,
            )

        if database == "encode":
            return self.encode_client.search(
                species=species,
                query=query,
                max_results=max_results,
            )

        if database == "sra":
            raise UnsupportedDatabaseError(
                "SRA search is not implemented yet."
            )

        if database == "all":
            geo_records = self.geo_client.search(
                species=species,
                query=query,
                max_results=max_results,
            )
            encode_records = self.encode_client.search(
                species=species,
                query=query,
                max_results=max_results,
            )
            return self._interleave_records(
                geo_records,
                encode_records,
                max_results=max_results,
            )

        raise UnsupportedDatabaseError(
            f"Unsupported database: {database}"
        )

    @staticmethod
    def _interleave_records(
        geo_records: list[DatasetRecord],
        encode_records: list[DatasetRecord],
        *,
        max_results: int,
    ) -> list[DatasetRecord]:
        """Interleave and deduplicate records from multiple databases."""
        merged: list[DatasetRecord] = []
        seen: set[tuple[str, str]] = set()

        longest = max(len(geo_records), len(encode_records))

        for index in range(longest):
            for records in (geo_records, encode_records):
                if index >= len(records):
                    continue

                record = records[index]
                identity = (
                    record.accession.strip().upper(),
                    record.url.strip(),
                )

                if identity in seen:
                    continue

                seen.add(identity)
                merged.append(record)

                if len(merged) >= max_results:
                    return merged

        return merged
