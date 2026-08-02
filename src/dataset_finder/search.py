"""Dataset search orchestration."""

from __future__ import annotations

from collections.abc import Iterable

from dataset_finder.clients.encode import ENCODEClient
from dataset_finder.clients.ncbi_bioproject import (
    NCBIBioProjectClient,
)
from dataset_finder.clients.ncbi_geo import NCBIGEOClient
from dataset_finder.clients.ncbi_sra import NCBISRAClient
from dataset_finder.models import DatasetRecord


class UnsupportedDatabaseError(ValueError):
    """Raised when a requested database is unsupported."""


class SearchService:
    """Coordinate searches across public databases."""

    def __init__(
        self,
        *,
        geo_client: NCBIGEOClient | None = None,
        encode_client: ENCODEClient | None = None,
        sra_client: NCBISRAClient | None = None,
        bioproject_client: NCBIBioProjectClient | None = None,
    ) -> None:
        self.geo_client = geo_client or NCBIGEOClient()
        self.encode_client = encode_client or ENCODEClient()
        self.sra_client = sra_client or NCBISRAClient()
        self.bioproject_client = (
            bioproject_client
            or NCBIBioProjectClient()
        )

    def search(
        self,
        *,
        species: str,
        query: str,
        database: str = "geo",
        max_results: int = 20,
    ) -> list[DatasetRecord]:
        """Search the requested database and normalize results."""
        species = species.strip()
        query = query.strip()
        database = database.strip().lower()

        if not species:
            raise ValueError("Species cannot be empty.")

        if not query:
            raise ValueError("Query cannot be empty.")

        if max_results < 1:
            raise ValueError(
                "Maximum results must be greater than zero."
            )

        clients = {
            "geo": self.geo_client,
            "encode": self.encode_client,
            "sra": self.sra_client,
            "bioproject": self.bioproject_client,
        }

        if database in clients:
            return clients[database].search(
                species=species,
                query=query,
                max_results=max_results,
            )

        if database == "all":
            result_groups: list[list[DatasetRecord]] = []

            for client in clients.values():
                try:
                    records = client.search(
                        species=species,
                        query=query,
                        max_results=max_results,
                    )
                except Exception:
                    continue

                result_groups.append(records)

            return self._interleave_records(
                result_groups,
                max_results=max_results,
            )

        raise UnsupportedDatabaseError(
            f"Unsupported database: {database}"
        )

    @staticmethod
    def _interleave_records(
        *record_groups: Iterable[list[DatasetRecord]],
        max_results: int,
    ) -> list[DatasetRecord]:
        """Interleave and deduplicate multiple result groups."""
        if (
            len(record_groups) == 1
            and not isinstance(record_groups[0], list)
        ):
            groups = list(record_groups[0])
        elif (
            len(record_groups) == 1
            and record_groups
            and record_groups[0]
            and isinstance(record_groups[0][0], list)
        ):
            groups = list(record_groups[0])
        else:
            groups = list(record_groups)

        merged: list[DatasetRecord] = []
        seen: set[tuple[str, str]] = set()

        longest = max(
            (len(records) for records in groups),
            default=0,
        )

        for index in range(longest):
            for records in groups:
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
