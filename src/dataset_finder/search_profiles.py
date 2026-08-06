"""Search profiles for technique-specific dataset discovery."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class TechniqueSearchProfile:
    """Technique label and repository search terminology."""

    name: str
    terms: tuple[str, ...]


TECHNIQUE_SEARCH_PROFILES = (
    TechniqueSearchProfile(
        name="RNA-seq",
        terms=(
            "RNA-seq",
            "RNA seq",
            "RNAseq",
            "transcriptome sequencing",
            "transcriptomic sequencing",
            "expression profiling by high throughput sequencing",
        ),
    ),
    TechniqueSearchProfile(
        name="ChIP-seq",
        terms=(
            "ChIP-seq",
            "ChIP seq",
            "ChIPseq",
            "chromatin immunoprecipitation sequencing",
            "genome binding occupancy profiling",
        ),
    ),
    TechniqueSearchProfile(
        name="CUT&RUN",
        terms=(
            "CUT&RUN",
            "CUT and RUN",
            "CUT-AND-RUN",
            "CUT N RUN",
            "CUT-RUN",
            "cleavage under targets and release using nuclease",
        ),
    ),
    TechniqueSearchProfile(
        name="CUT&Tag",
        terms=(
            "CUT&Tag",
            "CUT and Tag",
            "CUT-TAG",
            "CUT N TAG",
            "CUT-AND-TAG",
            "cleavage under targets and tagmentation",
        ),
    ),
    TechniqueSearchProfile(
        name="CLIP-seq",
        terms=(
            "CLIP",
            "CLIP-seq",
            "CLIP seq",
            "HITS-CLIP",
            "iCLIP",
            "PAR-CLIP",
            "eCLIP",
            "crosslinking immunoprecipitation",
        ),
    ),
)


def technique_profile(name: str) -> TechniqueSearchProfile:
    """Return one technique profile by normalized name."""
    normalized = name.strip().casefold()

    for profile in TECHNIQUE_SEARCH_PROFILES:
        if profile.name.casefold() == normalized:
            return profile

    raise ValueError(
        f"Unsupported technique profile: {name}"
    )
