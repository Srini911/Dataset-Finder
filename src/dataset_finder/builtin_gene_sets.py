"""Built-in curated gene sets distributed with Dataset Finder."""

from __future__ import annotations

from importlib.resources import files


class BuiltinGeneSetError(ValueError):
    """Raised when a built-in gene set cannot be loaded."""


GENE_SET_FILES = {
    "rbp": "drosophila_rbps.txt",
    "tf": "drosophila_tfs.txt",
}


def available_gene_sets() -> tuple[str, ...]:
    """Return available built-in gene-set names."""
    return tuple(GENE_SET_FILES)


def load_builtin_gene_set(name: str) -> list[str]:
    """Load a curated built-in gene list."""
    normalized_name = name.strip().casefold()

    if normalized_name not in GENE_SET_FILES:
        choices = ", ".join(available_gene_sets())
        raise BuiltinGeneSetError(
            f"Unknown built-in gene set: {name}. Available sets: {choices}."
        )

    resource = (
        files("dataset_finder")
        .joinpath("data")
        .joinpath("gene_sets")
        .joinpath(GENE_SET_FILES[normalized_name])
    )

    try:
        content = resource.read_text(encoding="utf-8")
    except OSError as exc:
        raise BuiltinGeneSetError(
            f"Unable to load built-in gene set {normalized_name}: {exc}"
        ) from exc

    genes = [
        line.strip()
        for line in content.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

    if not genes:
        raise BuiltinGeneSetError(
            f"Built-in gene set is empty: {normalized_name}"
        )

    return genes
