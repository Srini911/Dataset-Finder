"""Shared data models for Dataset Finder."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True, slots=True)
class DatasetRecord:
    """A normalized public functional-genomics dataset record."""

    uid: str
    accession: str
    title: str
    organism: str
    study_type: str
    sample_count: int | None
    publication_date: str
    url: str

    database: str = ""
    gene: str = ""
    gene_set: str = ""
    official_symbol: str = ""
    flybase_id: str = ""
    synonyms: tuple[str, ...] = field(default_factory=tuple)
    technique: str = ""
    technique_subtype: str = ""
    project_accession: str = ""
    sample_accessions: tuple[str, ...] = field(default_factory=tuple)
    description: str = ""
    tissue: str = ""
    developmental_stage: str = ""
    sex: str = ""
    genotype: str = ""
    perturbation: str = ""
    control: str = ""
    publication: str = ""
    evidence_text: str = ""
    match_type: str = ""
    confidence: str = ""
    search_date: str = ""
    flyatlas_url: str = ""
    flybase_url: str = ""
    raw_metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a serialization-friendly dictionary."""
        result = asdict(self)
        result["synonyms"] = "; ".join(self.synonyms)
        result["sample_accessions"] = "; ".join(self.sample_accessions)
        result["raw_metadata"] = str(self.raw_metadata)
        return result
