"""Tests for historical technique-specific discovery."""

from __future__ import annotations

from typing import Any

import pytest

from dataset_finder.historical_search import (
    HistoricalSearchService,
)
from dataset_finder.models import DatasetRecord
from dataset_finder.search_profiles import (
    TechniqueSearchProfile,
)

RNA_PROFILE = TechniqueSearchProfile(
    name="RNA-seq",
    terms=(
        "RNA-seq",
        "transcriptome sequencing",
    ),
)


class FakeEntrezClient:
    """Return deterministic historical NCBI records."""

    def __init__(self) -> None:
        self.search_calls: list[dict[str, Any]] = []

    def search_ids_paged(
        self,
        *,
        database: str,
        term: str,
        max_results: int,
        page_size: int,
        minimum_date: str,
        maximum_date: str,
        date_type: str,
    ) -> list[str]:
        self.search_calls.append(
            {
                "database": database,
                "term": term,
                "max_results": max_results,
                "page_size": page_size,
                "minimum_date": minimum_date,
                "maximum_date": maximum_date,
                "date_type": date_type,
            }
        )

        return (
            ["101"]
            if database == "sra"
            else ["201"]
        )

    def summaries(
        self,
        *,
        database: str,
        identifiers: list[str],
    ) -> list[dict[str, Any]]:
        del identifiers

        if database == "sra":
            return [
                {
                    "uid": "101",
                    "createdate": "2012/03/14",
                    "expxml": (
                        '<Study acc="SRP000101" '
                        'name="orb2 larval brain transcriptome"/>'
                        '<Experiment acc="SRX000101"/>'
                        '<Sample acc="SRS000101" '
                        'ScientificName="Drosophila melanogaster"/>'
                        "<Bioproject>PRJNA101</Bioproject>"
                        "<Biosample>SAMN000101</Biosample>"
                        "<LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>"
                        "<LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>"
                        "<LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>"
                    ),
                    "runs": '<Run acc="SRR000101"/>',
                }
            ]

        return [
            {
                "uid": "201",
                "accession": "GSE000201",
                "title": "orb2 RNA-seq in larval brain",
                "taxon": ["Drosophila melanogaster"],
                "gdstype": [
                    "Expression profiling by "
                    "high throughput sequencing"
                ],
                "n_samples": 4,
                "pdat": "2010/06/01",
            }
        ]


def test_historical_search_queries_sra_and_geo() -> None:
    entrez = FakeEntrezClient()
    service = HistoricalSearchService(
        entrez_client=entrez,  # type: ignore[arg-type]
        start_year=2005,
        end_year=2026,
        page_size=50,
    )

    result = service.search(
        species="Drosophila melanogaster",
        gene_terms=["orb2"],
        max_results_per_query=100,
        profiles=(RNA_PROFILE,),
    )

    assert len(result.records) == 2
    assert {
        record.database
        for record in result.records
    } == {"SRA", "GEO"}

    assert {
        call["database"]
        for call in entrez.search_calls
    } == {"sra", "gds"}

    for call in entrez.search_calls:
        assert call["minimum_date"] == "2005/01/01"
        assert call["maximum_date"] == "2026/12/31"
        assert call["date_type"] == "pdat"
        assert '"orb2"[All Fields]' in call["term"]
        assert '"RNA-seq"[All Fields]' in call["term"]


def test_historical_search_preserves_real_accessions() -> None:
    service = HistoricalSearchService(
        entrez_client=FakeEntrezClient(),  # type: ignore[arg-type]
        start_year=2005,
        end_year=2026,
    )

    result = service.search(
        species="Drosophila melanogaster",
        gene_terms=["orb2"],
        profiles=(RNA_PROFILE,),
    )

    sra_record = next(
        record
        for record in result.records
        if record.database == "SRA"
    )

    assert sra_record.uid == "101"
    assert sra_record.accession == "SRP000101"
    assert sra_record.experiment_accessions == (
        "SRX000101",
    )
    assert sra_record.sample_accessions == (
        "SRR000101",
    )
    assert sra_record.project_accession == "PRJNA101"


def test_historical_search_records_provenance_and_year() -> None:
    service = HistoricalSearchService(
        entrez_client=FakeEntrezClient(),  # type: ignore[arg-type]
        start_year=2005,
        end_year=2026,
        historical_cutoff_year=2015,
    )

    result = service.search(
        species="Drosophila melanogaster",
        gene_terms=["orb2"],
        profiles=(RNA_PROFILE,),
    )

    for record in result.records:
        assert record.gene_query_used == "orb2"
        assert record.technique_requested == "RNA-seq"
        assert "RNA-seq" in record.technique_search_term
        assert record.search_query_used
        assert record.historical_study is True
        assert record.study_year in {2010, 2012}


def test_structured_metadata_verifies_rna_seq() -> None:
    service = HistoricalSearchService(
        entrez_client=FakeEntrezClient(),  # type: ignore[arg-type]
        start_year=2005,
        end_year=2026,
    )

    result = service.search(
        species="Drosophila melanogaster",
        gene_terms=["orb2"],
        profiles=(RNA_PROFILE,),
    )

    for record in result.records:
        assert record.technique == "RNA_seq"
        assert record.technique_evidence
        assert record.technique_evidence_source == (
            "Structured repository metadata"
        )


def test_historical_search_rejects_invalid_year_range() -> None:
    with pytest.raises(
        ValueError,
        match="cannot precede",
    ):
        HistoricalSearchService(
            start_year=2026,
            end_year=2005,
        )


def test_specific_cut_run_title_overrides_generic_chip_strategy() -> None:
    record = DatasetRecord(
        uid="1",
        accession="SRP1",
        title="Anti-GFP Cut & Run in Drosophila wing tissue",
        organism="Drosophila melanogaster",
        study_type="ChIP-Seq",
        sample_count=2,
        publication_date="2023/01/01",
        url="https://example.org/SRP1",
        library_strategy="ChIP-Seq",
    )

    technique, evidence, source = (
        HistoricalSearchService._verified_technique(record)
    )

    assert technique == "CUT_RUN"
    assert "Cut & Run" in evidence
    assert source == (
        "Specific assay evidence in title or description"
    )


def test_requested_and_verified_technique_match() -> None:
    assert HistoricalSearchService._technique_match(
        requested="CUT&RUN",
        verified="CUT_RUN",
    ) == "Exact"

    assert HistoricalSearchService._technique_match(
        requested="RNA-seq",
        verified="ChIP_seq",
    ) == "Mismatch"

    assert HistoricalSearchService._technique_match(
        requested="CLIP-seq",
        verified="eCLIP",
    ) == "Exact"

    assert HistoricalSearchService._technique_match(
        requested="ChIP-seq",
        verified="Other_Assays",
    ) == "Unverified"
