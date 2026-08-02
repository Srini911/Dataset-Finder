"""Validate whether dataset metadata supports a gene match."""

from __future__ import annotations

import re
from dataclasses import dataclass

from dataset_finder.flybase_resolver import FlyBaseGene
from dataset_finder.models import DatasetRecord


@dataclass(frozen=True, slots=True)
class RelevanceAssessment:
    """Evidence supporting one gene-to-dataset match."""

    accepted: bool
    match_type: str
    confidence: str
    evidence: str


def _contains_exact_term(
    text: str,
    term: str,
    *,
    case_sensitive: bool = False,
) -> bool:
    """Match a complete term without adjacent letters or digits."""
    if not term:
        return False

    pattern = (
        rf"(?<![A-Za-z0-9])"
        rf"{re.escape(term)}"
        rf"(?![A-Za-z0-9])"
    )

    flags = 0 if case_sensitive else re.IGNORECASE

    return bool(re.search(pattern, text, flags=flags))


def _normalized_compact(value: str) -> str:
    """Normalize names for forms such as Bruno 1 versus Bruno1."""
    return re.sub(
        r"[^A-Za-z0-9]+",
        "",
        value,
    ).casefold()


def _contains_compact_name(
    text: str,
    name: str,
) -> bool:
    """Match a multi-character name while ignoring spaces and punctuation."""
    compact_name = _normalized_compact(name)

    if len(compact_name) < 3:
        return False

    compact_text = _normalized_compact(text)
    return compact_name in compact_text


def _short_symbol_context_match(
    text: str,
    symbol: str,
) -> bool:
    """Require explicit gene-like context for one- or two-character symbols."""
    if not symbol or len(symbol) > 2:
        return False

    escaped = re.escape(symbol)

    patterns = (
        rf"\bgene\s+{escaped}\b",
        rf"\b{escaped}\s+gene\b",
        rf"\b{escaped}[-\s]?RNAi\b",
        rf"\b{escaped}\s+mutant\b",
        rf"\b{escaped}\s+knockdown\b",
        rf"\b{escaped}\s+knockout\b",
        rf"\b{escaped}\s+overexpression\b",
        rf"\b{escaped}\s+depletion\b",
        rf"\b{escaped}\s+loss[-\s]of[-\s]function\b",
    )

    return any(
        re.search(
            pattern,
            text,
            flags=re.IGNORECASE,
        )
        for pattern in patterns
    )


def record_search_text(record: DatasetRecord) -> str:
    """Combine searchable record metadata."""
    values = [
        record.title,
        record.description,
        record.evidence_text,
        record.study_type,
        record.accession,
        record.project_accession,
        record.publication,
        record.tissue,
        record.genotype,
        record.perturbation,
    ]

    return " ".join(
        str(value).strip()
        for value in values
        if value and str(value).strip()
    )


def assess_relevance(
    *,
    record: DatasetRecord,
    submitted_gene: str,
    resolved_gene: FlyBaseGene,
) -> RelevanceAssessment:
    """Assess whether metadata genuinely supports a gene match."""
    text = record_search_text(record)

    if (
        resolved_gene.flybase_id
        and _contains_exact_term(
            text,
            resolved_gene.flybase_id,
        )
    ):
        return RelevanceAssessment(
            accepted=True,
            match_type="FlyBase identifier",
            confidence="High",
            evidence=resolved_gene.flybase_id,
        )

    official_symbol = resolved_gene.official_symbol

    if official_symbol:
        if len(official_symbol) <= 2:
            official_match = _short_symbol_context_match(
                text,
                official_symbol,
            )
        else:
            official_match = _contains_exact_term(
                text,
                official_symbol,
            )

        if official_match:
            return RelevanceAssessment(
                accepted=True,
                match_type="Official symbol",
                confidence="High",
                evidence=official_symbol,
            )

    if submitted_gene:
        if len(submitted_gene) <= 2:
            submitted_match = _short_symbol_context_match(
                text,
                submitted_gene,
            )
        else:
            submitted_match = _contains_exact_term(
                text,
                submitted_gene,
            )

        if submitted_match:
            return RelevanceAssessment(
                accepted=True,
                match_type="Submitted symbol",
                confidence="High",
                evidence=submitted_gene,
            )

    current_fullname = getattr(
        resolved_gene,
        "current_fullname",
        "",
    )

    if (
        current_fullname
        and _contains_compact_name(
            text,
            current_fullname,
        )
    ):
        return RelevanceAssessment(
            accepted=True,
            match_type="FlyBase gene name",
            confidence="High",
            evidence=current_fullname,
        )

    for synonym in resolved_gene.synonyms:
        if len(_normalized_compact(synonym)) < 3:
            continue

        if (
            _contains_exact_term(text, synonym)
            or _contains_compact_name(text, synonym)
        ):
            return RelevanceAssessment(
                accepted=True,
                match_type="FlyBase synonym",
                confidence="Medium",
                evidence=synonym,
            )

    return RelevanceAssessment(
        accepted=False,
        match_type="No metadata evidence",
        confidence="Low",
        evidence="",
    )
