"""FlyAtlas 2 expression-table client."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass

import requests

FLYATLAS_DOWNLOAD_URL = (
    "https://motif.mvls.gla.ac.uk/FA2Direct/index.html"
)


class FlyAtlasClientError(RuntimeError):
    """Raised when a FlyAtlas request or response cannot be processed."""


@dataclass(frozen=True, slots=True)
class FlyAtlasExpression:
    """Selected FlyAtlas tissue-expression measurements."""

    flybase_id: str
    symbol: str
    brain_male_fpkm: float | None = None
    brain_female_fpkm: float | None = None
    brain_larval_fpkm: float | None = None
    head_male_fpkm: float | None = None
    head_female_fpkm: float | None = None
    top_male_tissue: str = ""
    top_male_fpkm: float | None = None
    top_female_tissue: str = ""
    top_female_fpkm: float | None = None
    top_larval_tissue: str = ""
    top_larval_fpkm: float | None = None


class FlyAtlasClient:
    """Download and parse FlyAtlas 2 gene-expression tables."""

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
                "Accept": "text/plain",
                "User-Agent": (
                    "Dataset-Finder/0.3.0 "
                    "(https://github.com/Srini911/Dataset-Finder)"
                ),
            }
        )
        self._cache: dict[str, FlyAtlasExpression] = {}

    def fetch(self, flybase_id: str) -> FlyAtlasExpression:
        """Fetch expression information for one FlyBase gene ID."""
        normalized_id = flybase_id.strip()

        if not normalized_id:
            return FlyAtlasExpression(
                flybase_id="",
                symbol="",
            )

        if normalized_id in self._cache:
            return self._cache[normalized_id]

        try:
            response = self.session.get(
                FLYATLAS_DOWNLOAD_URL,
                params={
                    "fbgn": normalized_id,
                    "tableOut": "gene",
                },
                timeout=self.timeout,
            )
            response.raise_for_status()
        except requests.RequestException as exc:
            raise FlyAtlasClientError(
                f"FlyAtlas request failed for {normalized_id}: {exc}"
            ) from exc

        expression = self.parse_table(response.text)

        if not expression.flybase_id:
            expression = FlyAtlasExpression(
                flybase_id=normalized_id,
                symbol=expression.symbol,
                brain_male_fpkm=expression.brain_male_fpkm,
                brain_female_fpkm=expression.brain_female_fpkm,
                brain_larval_fpkm=expression.brain_larval_fpkm,
                head_male_fpkm=expression.head_male_fpkm,
                head_female_fpkm=expression.head_female_fpkm,
                top_male_tissue=expression.top_male_tissue,
                top_male_fpkm=expression.top_male_fpkm,
                top_female_tissue=expression.top_female_tissue,
                top_female_fpkm=expression.top_female_fpkm,
                top_larval_tissue=expression.top_larval_tissue,
                top_larval_fpkm=expression.top_larval_fpkm,
            )

        self._cache[normalized_id] = expression
        return expression

    @classmethod
    def parse_table(cls, text: str) -> FlyAtlasExpression:
        """Parse one FlyAtlas direct-download gene table."""
        lines = [
            line
            for line in text.splitlines()
            if line.strip()
        ]

        metadata: dict[str, str] = {}
        tissue_start: int | None = None

        for index, line in enumerate(lines):
            columns = line.split("\t")

            if columns[0].strip() == "Tissue":
                tissue_start = index + 1
                break

            if len(columns) >= 2:
                metadata[columns[0].strip()] = columns[1].strip()

        if tissue_start is None:
            raise FlyAtlasClientError(
                "FlyAtlas response did not contain a tissue table."
            )

        tissue_rows = list(
            csv.reader(
                io.StringIO(
                    "\n".join(lines[tissue_start:])
                ),
                delimiter="\t",
            )
        )

        brain_values = cls._find_tissue(
            tissue_rows,
            "Brain / CNS",
        )
        head_values = cls._find_tissue(
            tissue_rows,
            "Head",
        )

        top_male = cls._top_tissue(
            tissue_rows,
            value_index=1,
        )
        top_female = cls._top_tissue(
            tissue_rows,
            value_index=4,
        )
        top_larval = cls._top_tissue(
            tissue_rows,
            value_index=9,
        )

        return FlyAtlasExpression(
            flybase_id=metadata.get("FlyBase ID", ""),
            symbol=metadata.get("Symbol", ""),
            brain_male_fpkm=cls._number(brain_values, 1),
            brain_female_fpkm=cls._number(brain_values, 4),
            brain_larval_fpkm=cls._number(brain_values, 9),
            head_male_fpkm=cls._number(head_values, 1),
            head_female_fpkm=cls._number(head_values, 4),
            top_male_tissue=top_male[0],
            top_male_fpkm=top_male[1],
            top_female_tissue=top_female[0],
            top_female_fpkm=top_female[1],
            top_larval_tissue=top_larval[0],
            top_larval_fpkm=top_larval[1],
        )

    @staticmethod
    def _find_tissue(
        rows: list[list[str]],
        tissue_name: str,
    ) -> list[str]:
        for row in rows:
            if row and row[0].strip() == tissue_name:
                return row

        return []

    @staticmethod
    def _number(
        row: list[str],
        index: int,
    ) -> float | None:
        if index >= len(row):
            return None

        value = row[index].strip()

        if not value or value == "-":
            return None

        try:
            return float(value)
        except ValueError:
            return None

    @classmethod
    def _top_tissue(
        cls,
        rows: list[list[str]],
        *,
        value_index: int,
    ) -> tuple[str, float | None]:
        best_tissue = ""
        best_value: float | None = None

        for row in rows:
            if not row:
                continue

            tissue = row[0].strip()

            if not tissue or tissue == "Whole body":
                continue

            value = cls._number(row, value_index)

            if value is None:
                continue

            if best_value is None or value > best_value:
                best_tissue = tissue
                best_value = value

        return best_tissue, best_value
