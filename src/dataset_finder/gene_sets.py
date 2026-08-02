"""Utilities for loading and normalizing multiple gene symbols."""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path


class GeneInputError(ValueError):
    """Raised when gene input cannot be loaded or validated."""


def normalize_gene(value: str) -> str:
    """Normalize one gene symbol without changing biological capitalization."""
    return value.strip()


def split_gene_values(values: Iterable[str]) -> list[str]:
    """Split whitespace-, comma-, semicolon-, or newline-separated gene values."""
    genes: list[str] = []

    for value in values:
        if value is None:
            continue

        parts = re.split(r"[,;\n\r\t]+", str(value))

        for part in parts:
            gene = normalize_gene(part)

            if gene:
                genes.append(gene)

    return genes


def deduplicate_genes(genes: Iterable[str]) -> list[str]:
    """Remove duplicate gene symbols while preserving input order."""
    unique_genes: list[str] = []
    seen: set[str] = set()

    for gene in genes:
        normalized = normalize_gene(gene)

        if not normalized:
            continue

        identity = normalized.casefold()

        if identity in seen:
            continue

        seen.add(identity)
        unique_genes.append(normalized)

    return unique_genes


def load_gene_file(path: str | Path) -> list[str]:
    """Load gene symbols from a plain-text file."""
    gene_path = Path(path)

    if not gene_path.exists():
        raise GeneInputError(f"Gene file does not exist: {gene_path}")

    if not gene_path.is_file():
        raise GeneInputError(f"Gene path is not a file: {gene_path}")

    try:
        content = gene_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise GeneInputError(
            f"Unable to read gene file {gene_path}: {exc}"
        ) from exc

    genes = split_gene_values([content])
    genes = deduplicate_genes(genes)

    if not genes:
        raise GeneInputError(f"Gene file contains no gene symbols: {gene_path}")

    return genes


def collect_genes(
    *,
    query: str | None = None,
    genes: Iterable[str] | None = None,
    gene_file: str | Path | None = None,
) -> list[str]:
    """Collect genes from CLI query, explicit values, and an optional file."""
    collected: list[str] = []

    if query:
        collected.extend(split_gene_values([query]))

    if genes:
        collected.extend(split_gene_values(genes))

    if gene_file:
        collected.extend(load_gene_file(gene_file))

    collected = deduplicate_genes(collected)

    if not collected:
        raise GeneInputError(
            "Provide at least one gene using --query, --genes, or --gene-file."
        )

    return collected
