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
    experiment_accessions: tuple[str, ...] = field(default_factory=tuple)
    sample_accessions: tuple[str, ...] = field(default_factory=tuple)
    biosample_accessions: tuple[str, ...] = field(default_factory=tuple)
    library_strategy: str = ""
    library_source: str = ""
    library_selection: str = ""
    library_layout: str = ""
    platform: str = ""
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
    flyatlas_brain_male_fpkm: float | None = None
    flyatlas_brain_female_fpkm: float | None = None
    flyatlas_brain_larval_fpkm: float | None = None
    flyatlas_head_male_fpkm: float | None = None
    flyatlas_head_female_fpkm: float | None = None
    flyatlas_top_male_tissue: str = ""
    flyatlas_top_male_fpkm: float | None = None
    flyatlas_top_female_tissue: str = ""
    flyatlas_top_female_fpkm: float | None = None
    flyatlas_top_larval_tissue: str = ""
    flyatlas_top_larval_fpkm: float | None = None
    raw_metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a serialization-friendly dictionary."""
        result = asdict(self)
        result["synonyms"] = "; ".join(self.synonyms)
        result["experiment_accessions"] = "; ".join(
            self.experiment_accessions
        )
        result["sample_accessions"] = "; ".join(
            self.sample_accessions
        )
        result["biosample_accessions"] = "; ".join(
            self.biosample_accessions
        )
        result["raw_metadata"] = str(self.raw_metadata)
        return result
