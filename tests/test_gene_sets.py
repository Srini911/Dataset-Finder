"""Tests for multi-gene input handling."""

from pathlib import Path

import pytest

from dataset_finder.gene_sets import (
    GeneInputError,
    collect_genes,
    deduplicate_genes,
    load_gene_file,
    split_gene_values,
)


def test_split_gene_values_supports_common_separators() -> None:
    genes = split_gene_values(
        [
            "bru1,Hrb98DE",
            "mod;sm",
            "snf\ntra2",
        ]
    )

    assert genes == [
        "bru1",
        "Hrb98DE",
        "mod",
        "sm",
        "snf",
        "tra2",
    ]


def test_deduplicate_genes_preserves_order() -> None:
    genes = deduplicate_genes(
        ["bru1", "Hrb98DE", "BRU1", "mod", "hrb98de"]
    )

    assert genes == ["bru1", "Hrb98DE", "mod"]


def test_load_gene_file(tmp_path: Path) -> None:
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text(
        "bru1\nHrb98DE\nmod\nbru1\n",
        encoding="utf-8",
    )

    assert load_gene_file(gene_file) == [
        "bru1",
        "Hrb98DE",
        "mod",
    ]


def test_load_gene_file_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(GeneInputError, match="does not exist"):
        load_gene_file(tmp_path / "missing.txt")


def test_load_gene_file_rejects_empty_file(tmp_path: Path) -> None:
    gene_file = tmp_path / "empty.txt"
    gene_file.write_text("\n\n", encoding="utf-8")

    with pytest.raises(GeneInputError, match="contains no gene symbols"):
        load_gene_file(gene_file)


def test_collect_genes_combines_sources() -> None:
    genes = collect_genes(
        query="bru1",
        genes=["Hrb98DE,mod", "sm"],
    )

    assert genes == [
        "bru1",
        "Hrb98DE",
        "mod",
        "sm",
    ]


def test_collect_genes_requires_input() -> None:
    with pytest.raises(GeneInputError, match="Provide at least one gene"):
        collect_genes()
