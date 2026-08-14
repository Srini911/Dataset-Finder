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
    technique_requested: str = ""
    technique_match: str = ""
    technique_search_term: str = ""
    technique_evidence: str = ""
    technique_evidence_source: str = ""
    gene_query_used: str = ""
    search_query_used: str = ""
    study_year: int | None = None
    historical_study: bool = False
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
    cell_type: str = ""
    developmental_stage: str = ""
    sex: str = ""
    genotype: str = ""
    strain: str = ""
    treatment: str = ""
    control_status: str = ""
    disease: str = ""
    time_point: str = ""
    perturbation: str = ""
    control: str = ""
    publication: str = ""
    pubmed_ids: tuple[str, ...] = field(default_factory=tuple)
    dois: tuple[str, ...] = field(default_factory=tuple)
    related_accessions: tuple[str, ...] = field(default_factory=tuple)
    related_geo_accessions: tuple[str, ...] = field(default_factory=tuple)
    related_study_accessions: tuple[str, ...] = field(default_factory=tuple)
    related_bioproject_accessions: tuple[str, ...] = field(
        default_factory=tuple
    )
    related_biosample_accessions: tuple[str, ...] = field(
        default_factory=tuple
    )
    evidence_text: str = ""
    match_type: str = ""
    confidence: str = ""
    ranking_score: int = 0
    ranking_reasons: tuple[str, ...] = field(default_factory=tuple)
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
        result["pubmed_ids"] = "; ".join(self.pubmed_ids)
        result["dois"] = "; ".join(self.dois)
        result["related_accessions"] = "; ".join(
            self.related_accessions
        )
        result["related_geo_accessions"] = "; ".join(
            self.related_geo_accessions
        )
        result["related_study_accessions"] = "; ".join(
            self.related_study_accessions
        )
        result["related_bioproject_accessions"] = "; ".join(
            self.related_bioproject_accessions
        )
        result["related_biosample_accessions"] = "; ".join(
            self.related_biosample_accessions
        )
        result["raw_metadata"] = str(self.raw_metadata)
        return result
