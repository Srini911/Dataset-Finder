"""Tests for gene-to-dataset relevance validation."""

from dataset_finder.flybase_resolver import FlyBaseGene, FlyBaseResolver
from dataset_finder.models import DatasetRecord
from dataset_finder.relevance import assess_relevance


def make_record(
    text: str,
    *,
    database: str = "GEO",
) -> DatasetRecord:
    return DatasetRecord(
        uid="1",
        accession="GSE1",
        title=text,
        organism="Drosophila melanogaster",
        study_type="RNA-seq",
        sample_count=4,
        publication_date="2026-01-01",
        url="https://example.org/GSE1",
        database=database,
    )


def test_official_symbol_match_is_high_confidence() -> None:
    gene = FlyBaseResolver().resolve("bru1")

    result = assess_relevance(
        record=make_record("bru1 brain RNA-seq"),
        submitted_gene="bru1",
        resolved_gene=gene,
    )

    assert result.accepted
    assert result.match_type == "Official symbol"
    assert result.confidence == "High"


def test_flybase_identifier_match_is_high_confidence() -> None:
    gene = FlyBaseResolver().resolve("bru1")

    result = assess_relevance(
        record=make_record("Study targeting FBgn0000114"),
        submitted_gene="bru1",
        resolved_gene=gene,
    )

    assert result.accepted
    assert result.match_type == "FlyBase identifier"


def test_short_symbol_does_not_match_species_abbreviation() -> None:
    gene = FlyBaseResolver().resolve("D")

    result = assess_relevance(
        record=make_record(
            "Oxidative stress in D. melanogaster",
        ),
        submitted_gene="D",
        resolved_gene=gene,
    )

    assert not result.accepted


def test_short_symbol_requires_gene_context() -> None:
    gene = FlyBaseResolver().resolve("D")

    result = assess_relevance(
        record=make_record(
            "RNA-seq analysis of D mutant embryos",
        ),
        submitted_gene="D",
        resolved_gene=gene,
    )

    assert result.accepted
    assert result.match_type in {
        "Official symbol",
        "Submitted symbol",
    }


def test_unverified_database_match_is_rejected() -> None:
    gene = FlyBaseResolver().resolve("caz")

    result = assess_relevance(
        record=make_record(
            "Adult fly transcriptome",
            database="BioStudies",
        ),
        submitted_gene="caz",
        resolved_gene=gene,
    )

    assert not result.accepted


def test_long_synonym_can_support_match() -> None:
    gene = FlyBaseResolver().resolve("caz")

    result = assess_relevance(
        record=make_record("Drosophila dFUS neuronal study"),
        submitted_gene="caz",
        resolved_gene=gene,
    )

    assert result.accepted
    assert result.match_type == "FlyBase synonym"


def test_compact_gene_name_matches_bruno1_spelling() -> None:
    gene = FlyBaseResolver().resolve("bru1")

    result = assess_relevance(
        record=make_record(
            "Bruno1-mediated repression of oskar translation"
        ),
        submitted_gene="bru1",
        resolved_gene=gene,
    )

    assert result.accepted
    assert result.match_type in {
        "FlyBase gene name",
        "FlyBase synonym",
    }


def test_orb_does_not_match_orb2() -> None:
    gene = FlyBaseResolver().resolve("orb")

    result = assess_relevance(
        record=make_record(
            "The ORB2 RNA-binding protein regulates maternal transcripts"
        ),
        submitted_gene="orb",
        resolved_gene=gene,
    )

    assert not result.accepted


def test_orb2_matches_orb2() -> None:
    gene = FlyBaseResolver().resolve("orb2")

    result = assess_relevance(
        record=make_record(
            "The ORB2 RNA-binding protein regulates maternal transcripts"
        ),
        submitted_gene="orb2",
        resolved_gene=gene,
    )

    assert result.accepted
    assert result.match_type in {
        "Official symbol",
        "Submitted symbol",
        "FlyBase synonym",
    }


def test_distinctive_historical_query_provenance_can_support_match() -> None:
    resolved = FlyBaseGene(
        submitted_symbol="CG7804",
        official_symbol="cocoon",
        flybase_id="FBgn0036496",
        current_fullname="cocoon",
        synonyms=("CG7804",),
        secondary_flybase_ids=(),
        annotation_id="",
        match_type="synonym",
        ambiguous=False,
    )

    record = DatasetRecord(
        uid="4211124",
        accession="SRP110269",
        title=(
            "Fast evolution of gained essential function by a young gene "
            "through gained interaction with other essential genes [RNA-Seq]"
        ),
        organism="Drosophila melanogaster",
        study_type="RNA-Seq",
        sample_count=1,
        publication_date="2017/01/01",
        url="https://example.org/SRP110269",
        database="SRA",
        gene_query_used="CG7804",
    )

    result = assess_relevance(
        record=record,
        submitted_gene="CG7804",
        resolved_gene=resolved,
    )

    assert result.accepted is True
    assert result.match_type == "Historical query provenance"
    assert result.evidence == "CG7804"


def test_short_ambiguous_query_provenance_is_not_enough() -> None:
    resolved = FlyBaseGene(
        submitted_symbol="sm",
        official_symbol="sm",
        flybase_id="FBgn0003435",
        current_fullname="smooth",
        synonyms=("sm",),
        secondary_flybase_ids=(),
        annotation_id="",
        match_type="official_symbol",
        ambiguous=False,
    )

    record = DatasetRecord(
        uid="30736048",
        accession="SRP475297",
        title="SM2698B systemic infection",
        organism="Drosophila melanogaster",
        study_type="RNA-Seq",
        sample_count=1,
        publication_date="2023/01/01",
        url="https://example.org/SRP475297",
        database="SRA",
        gene_query_used="sm",
    )

    result = assess_relevance(
        record=record,
        submitted_gene="sm",
        resolved_gene=resolved,
    )

    assert result.accepted is False
