"""Biological ranking of discovered datasets."""

from __future__ import annotations

import re
from dataclasses import replace

from dataset_finder.models import DatasetRecord

STRONG_NEURAL_PATTERNS = (
    r"\bbrain\b",
    r"\bcentral nervous system\b",
    r"\bcns\b",
    r"\bmushroom bod(?:y|ies)\b",
    r"\bkenyon cells?\b",
    r"\boptic lobes?\b",
    r"\bantennal lobes?\b",
    r"\bventral nerve cord\b",
    r"\bneurons?\b",
    r"\bneuronal\b",
    r"\bglia\b",
    r"\bglial\b",
    r"\bneural stem cells?\b",
    r"\bneuroblast(?:s)?\b",
)

PERTURBATION_PATTERNS = (
    r"\brnai\b",
    r"\bknockdown\b",
    r"\bknock[- ]?out\b",
    r"\bdepletion\b",
    r"\bmutant\b",
    r"\boverexpression\b",
    r"\boverexpress",
    r"\bcrispr\b",
    r"\bimmunoprecipitation\b",
    r"\bip\b",
    r"\bchip\b",
)

MATCH_TYPE_SCORES = {
    "FlyBase identifier": 100,
    "Official symbol": 90,
    "Submitted symbol": 85,
    "FlyBase gene name": 80,
    "FlyBase synonym": 70,
}

CONFIDENCE_SCORES = {
    "High": 25,
    "Medium": 10,
    "Low": 0,
}

TECHNIQUE_MATCH_SCORES = {
    "Exact": 60,
    "Compatible": 40,
    "Unverified": -10,
    "Mismatch": -100,
}


def _record_text(record: DatasetRecord) -> str:
    values = (
        record.title,
        record.description,
        record.evidence_text,
        record.tissue,
        record.cell_type,
        record.developmental_stage,
        record.genotype,
        record.treatment,
        record.perturbation,
        record.study_type,
        record.library_strategy,
    )

    return " ".join(
        str(value).strip()
        for value in values
        if value and str(value).strip()
    )


def _matches_any(
    text: str,
    patterns: tuple[str, ...],
) -> bool:
    return any(
        re.search(pattern, text, flags=re.IGNORECASE)
        for pattern in patterns
    )


def score_dataset_record(
    record: DatasetRecord,
) -> tuple[int, tuple[str, ...]]:
    """Score one accepted record for Ray Lab usefulness."""
    score = 0
    reasons: list[str] = []

    match_score = MATCH_TYPE_SCORES.get(
        record.match_type,
        0,
    )

    if match_score:
        score += match_score
        reasons.append(
            f"gene evidence: {record.match_type} (+{match_score})"
        )

    confidence_score = CONFIDENCE_SCORES.get(
        record.confidence,
        0,
    )

    if confidence_score:
        score += confidence_score
        reasons.append(
            f"{record.confidence.lower()} confidence "
            f"(+{confidence_score})"
        )

    technique_score = TECHNIQUE_MATCH_SCORES.get(
        record.technique_match,
        0,
    )

    if technique_score:
        score += technique_score
        reasons.append(
            f"technique {record.technique_match.lower()} "
            f"({technique_score:+d})"
        )

    text = _record_text(record)

    if _matches_any(text, STRONG_NEURAL_PATTERNS):
        score += 55
        reasons.append("brain/neural evidence (+55)")

    if record.sex.strip():
        score += 15
        reasons.append("sex metadata (+15)")

    if _matches_any(text, PERTURBATION_PATTERNS):
        score += 25
        reasons.append("experimental perturbation/target evidence (+25)")

    if record.historical_study:
        score += 15
        reasons.append("historical study (+15)")

    if record.pubmed_ids:
        score += 10
        reasons.append("linked publication (+10)")

    if record.database.casefold() in {"geo", "sra"}:
        score += 10
        reasons.append("primary sequencing repository (+10)")

    accession = record.accession.upper()

    if accession.startswith(("GSM", "SRR", "ERR", "DRR")):
        score -= 15
        reasons.append("sample/run-level record (-15)")

    return score, tuple(reasons)


def rank_records(
    records: list[DatasetRecord],
) -> list[DatasetRecord]:
    """Rank accepted records from most to least useful."""
    scored: list[DatasetRecord] = []

    for record in records:
        score, reasons = score_dataset_record(record)

        scored.append(
            replace(
                record,
                ranking_score=score,
                ranking_reasons=reasons,
            )
        )

    return sorted(
        scored,
        key=lambda record: (
            -record.ranking_score,
            -int(record.historical_study),
            -(record.study_year or 0),
            record.database.casefold(),
            record.accession.casefold(),
        ),
    )
