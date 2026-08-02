"""Resolve curated Drosophila symbols using packaged FlyBase mappings."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from importlib.resources import files


@dataclass(frozen=True, slots=True)
class FlyBaseGene:
    """Resolved FlyBase gene information."""

    submitted_symbol: str
    official_symbol: str
    flybase_id: str
    current_fullname: str
    synonyms: tuple[str, ...]
    secondary_flybase_ids: tuple[str, ...]
    annotation_id: str
    match_type: str
    ambiguous: bool

    @property
    def flybase_url(self) -> str:
        """Return the FlyBase gene-report URL."""
        if not self.flybase_id:
            return ""

        return (
            "https://flybase.org/reports/"
            f"{self.flybase_id}.html"
        )

    @property
    def flyatlas_url(self) -> str:
        """Return the FlyAtlas 2 gene-results URL."""
        if self.flybase_id:
            return (
                "https://motif.mvls.gla.ac.uk/FlyAtlas2/"
                "index.html?search=gene&gene="
                f"{self.flybase_id}&idtype=fbgn"
            )

        if self.official_symbol:
            return (
                "https://motif.mvls.gla.ac.uk/FlyAtlas2/"
                "index.html?search=gene&gene="
                f"{self.official_symbol}&idtype=symbol"
            )

        return ""

    @property
    def flyatlas_download_url(self) -> str:
        """Return the FlyAtlas 2 direct gene-table download URL."""
        if not self.flybase_id:
            return ""

        return (
            "https://motif.mvls.gla.ac.uk/FA2Direct/"
            "index.html?fbgn="
            f"{self.flybase_id}&tableOut=gene"
        )

    @property
    def search_terms(self) -> tuple[str, ...]:
        """Return conservative search terms for external databases."""
        values = [
            self.official_symbol,
            self.submitted_symbol,
            self.flybase_id,
        ]

        if len(self.submitted_symbol) >= 3:
            values.extend(self.synonyms)

        unique: list[str] = []
        seen: set[str] = set()

        for value in values:
            value = value.strip()

            if not value:
                continue

            identity = value.casefold()

            if identity in seen:
                continue

            seen.add(identity)
            unique.append(value)

        return tuple(unique)


class FlyBaseResolver:
    """Resolve symbols from the packaged compact FlyBase index."""

    def __init__(self) -> None:
        (
            self._exact_records,
            self._casefold_records,
        ) = self._load_records()

    @staticmethod
    def _split(value: str) -> tuple[str, ...]:
        return tuple(
            item.strip()
            for item in value.split("|")
            if item.strip()
        )

    def _load_records(
        self,
    ) -> tuple[
        dict[str, FlyBaseGene],
        dict[str, list[FlyBaseGene]],
    ]:
        resource = (
            files("dataset_finder")
            .joinpath("data")
            .joinpath("flybase")
            .joinpath("drosophila_gene_index.tsv")
        )

        exact_records: dict[str, FlyBaseGene] = {}
        casefold_records: dict[str, list[FlyBaseGene]] = {}

        with resource.open(
            "r",
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.DictReader(
                handle,
                delimiter="\t",
            )

            for row in reader:
                submitted_symbol = row[
                    "submitted_symbol"
                ].strip()

                record = FlyBaseGene(
                    submitted_symbol=submitted_symbol,
                    official_symbol=row[
                        "official_symbol"
                    ].strip(),
                    flybase_id=row[
                        "flybase_id"
                    ].strip(),
                    current_fullname=row[
                        "current_fullname"
                    ].strip(),
                    synonyms=self._split(
                        row["symbol_synonyms"]
                    ),
                    secondary_flybase_ids=self._split(
                        row["secondary_flybase_ids"].replace(
                            ",",
                            "|",
                        )
                    ),
                    annotation_id=row[
                        "annotation_id"
                    ].strip(),
                    match_type=row[
                        "match_type"
                    ].strip(),
                    ambiguous=(
                        row["ambiguous"]
                        .strip()
                        .casefold()
                        == "yes"
                    ),
                )

                exact_records[submitted_symbol] = record
                casefold_records.setdefault(
                    submitted_symbol.casefold(),
                    [],
                ).append(record)

        return exact_records, casefold_records

    def resolve(self, symbol: str) -> FlyBaseGene:
        """Resolve one submitted symbol."""
        submitted = symbol.strip()
        exact_record = self._exact_records.get(
            submitted
        )

        if exact_record is not None:
            return exact_record

        folded_records = self._casefold_records.get(
            submitted.casefold(),
            [],
        )

        if len(folded_records) == 1:
            return folded_records[0]

        return FlyBaseGene(
            submitted_symbol=submitted,
            official_symbol="",
            flybase_id="",
            current_fullname="",
            synonyms=(),
            secondary_flybase_ids=(),
            annotation_id="",
            match_type="unresolved",
            ambiguous=False,
        )

    def resolve_many(
        self,
        symbols: list[str],
    ) -> dict[str, FlyBaseGene]:
        """Resolve multiple submitted symbols."""
        return {
            symbol: self.resolve(symbol)
            for symbol in symbols
        }
