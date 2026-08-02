"""Biological metadata extraction and normalization."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class BiologicalMetadata:
    """Normalized biological metadata extracted from dataset records."""

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


TISSUE_PATTERNS = {
    "brain": (
        r"\bbrain\b",
        r"\bcentral nervous system\b",
        r"\bcns\b",
    ),
    "head": (r"\bhead\b",),
    "eye": (r"\beye\b", r"\bretina\b"),
    "muscle": (
        r"\bmuscle\b",
        r"\bflight muscle\b",
        r"\bindirect flight musculature\b",
        r"\bifm\b",
    ),
    "gut": (r"\bgut\b",),
    "midgut": (r"\bmidgut\b",),
    "hindgut": (r"\bhindgut\b",),
    "fat body": (r"\bfat body\b",),
    "ovary": (r"\bovary\b", r"\bovarian\b"),
    "testis": (r"\btestis\b", r"\btestes\b"),
    "wing disc": (
        r"\bwing disc\b",
        r"\bwing imaginal disc\b",
    ),
    "imaginal disc": (r"\bimaginal disc\b",),
    "salivary gland": (r"\bsalivary gland\b",),
    "heart": (r"\bheart\b",),
    "hemocyte": (r"\bhemocyte\b", r"\bhaemocyte\b"),
}


CELL_TYPE_PATTERNS = {
    "neuron": (
        r"\bneurons?\b",
        r"\bneuronal\b",
    ),
    "glia": (r"\bglia\b", r"\bglial\b"),
    "astrocyte-like glia": (
        r"\bastrocyte[- ]like\b",
        r"\bastrocyte[- ]like glia\b",
    ),
    "hemocyte": (r"\bhemocyte\b", r"\bhaemocyte\b"),
    "enterocyte": (r"\benterocyte\b",),
    "intestinal stem cell": (
        r"\bintestinal stem cell\b",
        r"\bisc\b",
    ),
    "enteroblast": (r"\benteroblast\b",),
    "photoreceptor": (r"\bphotoreceptor\b",),
    "adipocyte": (r"\badipocyte\b",),
    "germ cell": (r"\bgerm cell\b",),
}


DEVELOPMENTAL_STAGE_PATTERNS = {
    "embryo": (r"\bembryo\b", r"\bembryonic\b"),
    "larva": (
        r"\blarva\b",
        r"\blarval\b",
        r"\bfirst instar\b",
        r"\bsecond instar\b",
        r"\bthird instar\b",
    ),
    "pupa": (r"\bpupa\b", r"\bpupal\b"),
    "adult": (r"\badult\b",),
}


SEX_PATTERNS = {
    "male": (r"\bmale\b", r"\bmales\b"),
    "female": (r"\bfemale\b", r"\bfemales\b"),
    "mixed": (
        r"\bmixed sex\b",
        r"\bmixed-sex\b",
        r"\bboth sexes\b",
    ),
}


STRAIN_PATTERNS = {
    "w1118": (r"\bw1118\b", r"\bw\[1118\]\b"),
    "Canton-S": (r"\bcanton[- ]s\b",),
    "Oregon-R": (r"\boregon[- ]r\b",),
    "DGRP-551": (r"\bdgrp[- ]551\b",),
}


PERTURBATION_PATTERNS = {
    "RNAi": (r"\brnai\b", r"\brna interference\b"),
    "knockdown": (r"\bknockdown\b", r"\bknock-down\b"),
    "knockout": (r"\bknockout\b", r"\bknock-out\b"),
    "overexpression": (
        r"\boverexpression\b",
        r"\bover-expression\b",
    ),
    "CRISPR": (r"\bcrispr\b",),
    "mutant": (r"\bmutant\b", r"\bmutation\b"),
}


CONTROL_PATTERNS = {
    "control": (
        r"\bcontrol\b",
        r"\buntreated\b",
        r"\bvehicle-treated\b",
        r"\bwild type\b",
        r"\bwild-type\b",
    ),
    "treated": (
        r"\btreated\b",
        r"\btreatment\b",
        r"\bexposed\b",
        r"\bstimulated\b",
    ),
}


TIME_POINT_PATTERN = re.compile(
    r"\b\d+(?:\.\d+)?\s*"
    r"(?:h|hr|hrs|hour|hours|d|day|days|"
    r"min|mins|minute|minutes|week|weeks)\b",
    flags=re.IGNORECASE,
)


def _flatten_text(value: Any) -> list[str]:
    """Convert nested values into searchable text fragments."""
    if value is None:
        return []

    if isinstance(value, dict):
        fragments: list[str] = []

        for key, item in value.items():
            fragments.extend(_flatten_text(key))
            fragments.extend(_flatten_text(item))

        return fragments

    if isinstance(value, (list, tuple, set)):
        fragments = []

        for item in value:
            fragments.extend(_flatten_text(item))

        return fragments

    text = str(value).strip()
    return [text] if text else []


def _combined_text(values: tuple[object, ...]) -> str:
    """Join arbitrary metadata values into normalized searchable text."""
    fragments: list[str] = []

    for value in values:
        fragments.extend(_flatten_text(value))

    return " ".join(fragments).casefold()


def _first_match(
    text: str,
    patterns: dict[str, tuple[str, ...]],
) -> str:
    """Return the first normalized label supported by the text."""
    for label, expressions in patterns.items():
        if any(
            re.search(
                expression,
                text,
                flags=re.IGNORECASE,
            )
            for expression in expressions
        ):
            return label

    return ""


def extract_biological_metadata(
    *texts: object,
) -> BiologicalMetadata:
    """Extract normalized biological metadata from arbitrary values."""
    text = _combined_text(texts)

    time_point_match = TIME_POINT_PATTERN.search(text)

    return BiologicalMetadata(
        tissue=_first_match(text, TISSUE_PATTERNS),
        cell_type=_first_match(text, CELL_TYPE_PATTERNS),
        developmental_stage=_first_match(
            text,
            DEVELOPMENTAL_STAGE_PATTERNS,
        ),
        sex=_first_match(text, SEX_PATTERNS),
        strain=_first_match(text, STRAIN_PATTERNS),
        control_status=_first_match(
            text,
            CONTROL_PATTERNS,
        ),
        time_point=(
            time_point_match.group(0)
            if time_point_match
            else ""
        ),
        perturbation=_first_match(
            text,
            PERTURBATION_PATTERNS,
        ),
    )
