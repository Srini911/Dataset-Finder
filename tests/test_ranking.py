"""Tests for biological dataset ranking."""

from dataset_finder.models import DatasetRecord
from dataset_finder.ranking import (
    rank_records,
    score_dataset_record,
)


def make_record(
    *,
    accession: str,
    title: str,
    match_type: str = "Official symbol",
    confidence: str = "High",
    technique_match: str = "Exact",
    tissue: str = "",
    sex: str = "",
    historical_study: bool = False,
    database: str = "GEO",
) -> DatasetRecord:
    return DatasetRecord(
        uid=accession,
        accession=accession,
        title=title,
        organism="Drosophila melanogaster",
        study_type="RNA-seq",
        sample_count=2,
        publication_date="2014/01/01",
        url=f"https://example.org/{accession}",
        database=database,
        match_type=match_type,
        confidence=confidence,
        technique_match=technique_match,
        tissue=tissue,
        sex=sex,
        historical_study=historical_study,
        study_year=2014,
    )


def test_brain_dataset_scores_above_non_neural_dataset() -> None:
    brain = make_record(
        accession="GSE_BRAIN",
        title="RNA-seq of adult Drosophila brain neurons",
        tissue="brain",
    )

    muscle = make_record(
        accession="GSE_MUSCLE",
        title="RNA-seq of adult flight muscle",
        tissue="muscle",
    )

    brain_score, _ = score_dataset_record(brain)
    muscle_score, _ = score_dataset_record(muscle)

    assert brain_score > muscle_score


def test_technique_mismatch_receives_large_penalty() -> None:
    exact = make_record(
        accession="GSE_EXACT",
        title="brain RNA-seq",
        technique_match="Exact",
    )

    mismatch = make_record(
        accession="GSE_MISMATCH",
        title="brain ChIP-seq",
        technique_match="Mismatch",
    )

    exact_score, _ = score_dataset_record(exact)
    mismatch_score, _ = score_dataset_record(mismatch)

    assert exact_score > mismatch_score
    assert exact_score - mismatch_score >= 100


def test_historical_neural_dataset_can_rank_above_new_generic_record() -> None:
    historical = make_record(
        accession="GSE_OLD",
        title="Drosophila brain RNA-binding protein RNA-seq",
        tissue="brain",
        historical_study=True,
    )

    generic = make_record(
        accession="GSE_NEW",
        title="Drosophila whole organism RNA-seq",
    )

    ranked = rank_records(
        [
            generic,
            historical,
        ]
    )

    assert ranked[0].accession == "GSE_OLD"
    assert ranked[0].ranking_score > ranked[1].ranking_score
    assert ranked[0].ranking_reasons


def test_sample_level_record_is_penalized() -> None:
    series = make_record(
        accession="GSE12345",
        title="brain RNA-seq",
    )

    sample = make_record(
        accession="GSM12345",
        title="brain RNA-seq",
    )

    series_score, _ = score_dataset_record(series)
    sample_score, _ = score_dataset_record(sample)

    assert series_score > sample_score
