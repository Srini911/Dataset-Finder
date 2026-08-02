"""Tests for gene-to-dataset relevance validation."""

from dataset_finder.flybase_resolver import FlyBaseResolver
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
