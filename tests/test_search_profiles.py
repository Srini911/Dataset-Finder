"""Tests for technique-specific search profiles."""

import pytest

from dataset_finder.search_profiles import (
    TECHNIQUE_SEARCH_PROFILES,
    technique_profile,
)


def test_expected_technique_profiles_are_available() -> None:
    names = {
        profile.name
        for profile in TECHNIQUE_SEARCH_PROFILES
    }

    assert names == {
        "RNA-seq",
        "ChIP-seq",
        "CUT&RUN",
        "CUT&Tag",
        "CLIP-seq",
    }


def test_rna_seq_profile_contains_historical_synonyms() -> None:
    profile = technique_profile("RNA-seq")

    assert "RNAseq" in profile.terms
    assert (
        "expression profiling by high throughput sequencing"
        in profile.terms
    )


def test_clip_profile_includes_eclip() -> None:
    profile = technique_profile("clip-seq")

    assert "eCLIP" in profile.terms
    assert "HITS-CLIP" in profile.terms


def test_unknown_profile_is_rejected() -> None:
    with pytest.raises(
        ValueError,
        match="Unsupported technique profile",
    ):
        technique_profile("unknown")
