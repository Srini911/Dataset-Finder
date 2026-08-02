"""Dataset search orchestration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from dataset_finder.clients.biostudies import BioStudiesClient
from dataset_finder.clients.ena import ENAClient
from dataset_finder.clients.encode import ENCODEClient
from dataset_finder.clients.ncbi_bioproject import NCBIBioProjectClient
from dataset_finder.clients.ncbi_biosample import NCBIBioSampleClient
from dataset_finder.clients.ncbi_geo import NCBIGEOClient
from dataset_finder.clients.ncbi_sra import NCBISRAClient
from dataset_finder.models import DatasetRecord


class SearchClient(Protocol):
    """Protocol implemented by database search clients."""

    def search(
        self,
        *,
        species: str,
        query: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        """Search one public database."""


@dataclass(frozen=True, slots=True)
class DatabaseSearchStatus:
    """Outcome of one database search."""

    database: str
    success: bool
    result_count: int
    error: str = ""


@dataclass(frozen=True, slots=True)
class SearchOutcome:
    """Records and source-level statuses from a search."""

    records: tuple[DatasetRecord, ...]
    statuses: tuple[DatabaseSearchStatus, ...]


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
        biosample_client: NCBIBioSampleClient | None = None,
        biostudies_client: BioStudiesClient | None = None,
        ena_client: ENAClient | None = None,
    ) -> None:
        self.geo_client = geo_client or NCBIGEOClient()
        self.encode_client = encode_client or ENCODEClient()
        self.sra_client = sra_client or NCBISRAClient()
        self.bioproject_client = bioproject_client or NCBIBioProjectClient()
        self.biosample_client = biosample_client or NCBIBioSampleClient()
        self.biostudies_client = biostudies_client or BioStudiesClient()
        self.ena_client = ena_client or ENAClient()

    def _clients(self) -> dict[str, SearchClient]:
        return {
            "geo": self.geo_client,
            "encode": self.encode_client,
            "sra": self.sra_client,
            "bioproject": self.bioproject_client,
            "biosample": self.biosample_client,
            "biostudies": self.biostudies_client,
            "ena": self.ena_client,
        }

    @staticmethod
    def _validate(
        *,
        species: str,
        query: str,
        max_results: int,
    ) -> tuple[str, str]:
        species = species.strip()
        query = query.strip()

        if not species:
            raise ValueError("Species cannot be empty.")

        if not query:
            raise ValueError("Query cannot be empty.")

        if max_results < 1:
            raise ValueError("Maximum results must be greater than zero.")

        return species, query

    def search_with_status(
        self,
        *,
        species: str,
        query: str,
        database: str = "geo",
        max_results: int = 20,
    ) -> SearchOutcome:
        """Search databases and preserve every database outcome."""
        species, query = self._validate(
            species=species,
            query=query,
            max_results=max_results,
        )
        database = database.strip().lower()
        clients = self._clients()

        if database == "all":
            selected_clients = clients
        elif database in clients:
            selected_clients = {
                database: clients[database],
            }
        else:
            raise UnsupportedDatabaseError(
                f"Unsupported database: {database}"
            )

        result_groups: list[list[DatasetRecord]] = []
        statuses: list[DatabaseSearchStatus] = []

        for database_name, client in selected_clients.items():
            try:
                records = client.search(
                    species=species,
                    query=query,
                    max_results=max_results,
                )
            except Exception as exc:
                statuses.append(
                    DatabaseSearchStatus(
                        database=database_name,
                        success=False,
                        result_count=0,
                        error=str(exc),
                    )
                )
                continue

            result_groups.append(records)
            statuses.append(
                DatabaseSearchStatus(
                    database=database_name,
                    success=True,
                    result_count=len(records),
                )
            )

        records = self._interleave_records(
            result_groups,
            max_results=max_results,
        )

        return SearchOutcome(
            records=tuple(records),
            statuses=tuple(statuses),
        )

    def search(
        self,
        *,
        species: str,
        query: str,
        database: str = "geo",
        max_results: int = 20,
    ) -> list[DatasetRecord]:
        """Search databases and return normalized records."""
        outcome = self.search_with_status(
            species=species,
            query=query,
            database=database,
            max_results=max_results,
        )

        if database.strip().lower() != "all":
            failed = [
                status
                for status in outcome.statuses
                if not status.success
            ]

            if failed:
                raise RuntimeError(failed[0].error)

        return list(outcome.records)

    @staticmethod
    def _interleave_records(
        record_groups: list[list[DatasetRecord]],
        *,
        max_results: int,
    ) -> list[DatasetRecord]:
        """Interleave and deduplicate multiple result groups."""
        merged: list[DatasetRecord] = []
        seen: set[tuple[str, str]] = set()

        longest = max(
            (len(records) for records in record_groups),
            default=0,
        )

        for index in range(longest):
            for records in record_groups:
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
