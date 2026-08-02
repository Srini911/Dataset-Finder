"""Tests for built-in curated gene sets."""

import pytest

from dataset_finder.builtin_gene_sets import (
    BuiltinGeneSetError,
    available_gene_sets,
    load_builtin_gene_set,
)


def test_available_gene_sets() -> None:
    assert available_gene_sets() == ("rbp", "tf")


def test_load_rbp_gene_set() -> None:
    genes = load_builtin_gene_set("rbp")

    assert "bru1" in genes
    assert "Hrb98DE" in genes
    assert "caz" in genes
    assert len(genes) > 100


def test_load_tf_gene_set() -> None:
    genes = load_builtin_gene_set("tf")

    assert "abd-A" in genes
    assert "D" in genes
    assert "zfh2" in genes
    assert len(genes) > 200


def test_gene_set_name_is_case_insensitive() -> None:
    assert load_builtin_gene_set("RBP") == load_builtin_gene_set("rbp")


def test_unknown_gene_set_is_rejected() -> None:
    with pytest.raises(BuiltinGeneSetError, match="Unknown built-in gene set"):
        load_builtin_gene_set("unknown")
