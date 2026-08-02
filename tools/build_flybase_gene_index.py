"""Build a compact FlyBase gene index for Dataset Finder gene sets."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path

from dataset_finder.builtin_gene_sets import load_builtin_gene_set

ROOT = Path(__file__).resolve().parents[1]

ANNOTATION_PATH = (
    ROOT
    / "data"
    / "flybase"
    / "fbgn_annotation_ID_fb_2026_02.tsv.gz"
)

SYNONYM_PATH = (
    ROOT
    / "data"
    / "flybase"
    / "fb_synonym_fb_2026_02.tsv.gz"
)

OUTPUT_PATH = (
    ROOT
    / "src"
    / "dataset_finder"
    / "data"
    / "flybase"
    / "drosophila_gene_index.tsv"
)


# Curated historical symbols whose aliases are shared by multiple genes.
# Values are current FlyBase symbols.
CURATED_SYMBOL_OVERRIDES = {
    "h": "hry",
}


def split_synonyms(value: str) -> list[str]:
    """Split FlyBase synonym fields conservatively."""
    if not value:
        return []

    normalized = value.replace(";", "|")

    return [
        item.strip()
        for item in normalized.split("|")
        if item.strip()
    ]


def load_annotation_rows() -> tuple[
    dict[str, dict[str, str]],
    dict[str, set[str]],
]:
    """Load exact and case-insensitive FlyBase annotation mappings."""
    exact_rows: dict[str, dict[str, str]] = {}
    casefold_symbols: dict[str, set[str]] = defaultdict(set)

    with gzip.open(
        ANNOTATION_PATH,
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row in reader:
            if not row or row[0].startswith("#"):
                continue

            row += [""] * (6 - len(row))

            symbol = row[0].strip()
            organism = row[1].strip()

            if organism != "Dmel" or not symbol:
                continue

            exact_rows[symbol] = {
                "official_symbol": symbol,
                "flybase_id": row[2].strip(),
                "secondary_flybase_ids": row[3].strip(),
                "annotation_id": row[4].strip(),
                "secondary_annotation_ids": row[5].strip(),
            }
            casefold_symbols[symbol.casefold()].add(symbol)

    return exact_rows, casefold_symbols


def load_synonym_rows() -> tuple[
    dict[str, dict[str, object]],
    dict[str, set[str]],
    dict[str, set[str]],
]:
    """Load exact symbols, case-insensitive symbols, and aliases."""
    by_exact_symbol: dict[str, dict[str, object]] = {}
    by_casefold_symbol: dict[str, set[str]] = defaultdict(set)
    alias_to_symbols: dict[str, set[str]] = defaultdict(set)

    with gzip.open(
        SYNONYM_PATH,
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row in reader:
            if not row or row[0].startswith("#"):
                continue

            row += [""] * (6 - len(row))

            primary_fbid = row[0].strip()
            organism = row[1].strip()
            current_symbol = row[2].strip()
            current_fullname = row[3].strip()
            fullname_synonyms = split_synonyms(row[4].strip())
            symbol_synonyms = split_synonyms(row[5].strip())

            if (
                organism != "Dmel"
                or not primary_fbid.startswith("FBgn")
                or not current_symbol
            ):
                continue

            entry = {
                "flybase_id": primary_fbid,
                "official_symbol": current_symbol,
                "current_fullname": current_fullname,
                "symbol_synonyms": symbol_synonyms,
                "fullname_synonyms": fullname_synonyms,
            }

            by_exact_symbol[current_symbol] = entry
            by_casefold_symbol[
                current_symbol.casefold()
            ].add(current_symbol)

            for alias in (
                current_symbol,
                *symbol_synonyms,
            ):
                if alias:
                    alias_to_symbols[
                        alias.casefold()
                    ].add(current_symbol)

    return (
        by_exact_symbol,
        by_casefold_symbol,
        alias_to_symbols,
    )


def main() -> None:
    """Build the compact curated FlyBase index."""
    if not ANNOTATION_PATH.exists():
        raise FileNotFoundError(ANNOTATION_PATH)

    if not SYNONYM_PATH.exists():
        raise FileNotFoundError(SYNONYM_PATH)

    requested_genes: list[str] = []
    seen_requested: set[str] = set()

    for gene_set in ("rbp", "tf"):
        for gene in load_builtin_gene_set(gene_set):
            identity = gene.casefold()

            if identity in seen_requested:
                continue

            seen_requested.add(identity)
            requested_genes.append(gene)

    annotation_exact, annotation_casefold = (
        load_annotation_rows()
    )
    (
        synonym_exact,
        synonym_casefold,
        alias_to_symbols,
    ) = load_synonym_rows()

    output_rows: list[dict[str, str]] = []
    unresolved: list[str] = []
    ambiguous: list[str] = []

    for submitted_symbol in requested_genes:
        folded = submitted_symbol.casefold()
        candidate_symbols: set[str] = set()

        # Exact case-sensitive official symbols always win.
        if submitted_symbol in synonym_exact:
            candidate_symbols.add(submitted_symbol)

        if submitted_symbol in annotation_exact:
            candidate_symbols.add(submitted_symbol)

        # Collect case-insensitive and alias matches only for diagnostics
        # and fallback resolution.
        candidate_symbols.update(
            synonym_casefold.get(folded, set())
        )
        candidate_symbols.update(
            annotation_casefold.get(folded, set())
        )
        candidate_symbols.update(
            alias_to_symbols.get(folded, set())
        )

        if not candidate_symbols:
            unresolved.append(submitted_symbol)
            output_rows.append(
                {
                    "submitted_symbol": submitted_symbol,
                    "official_symbol": "",
                    "flybase_id": "",
                    "current_fullname": "",
                    "symbol_synonyms": "",
                    "secondary_flybase_ids": "",
                    "annotation_id": "",
                    "match_type": "unresolved",
                    "ambiguous": "No",
                }
            )
            continue

        if len(candidate_symbols) > 1:
            ambiguous.append(submitted_symbol)

        override_symbol = CURATED_SYMBOL_OVERRIDES.get(
            submitted_symbol
        )

        if (
            override_symbol
            and override_symbol in synonym_exact
        ):
            selected_symbol = override_symbol
        elif submitted_symbol in synonym_exact:
            selected_symbol = submitted_symbol
        elif submitted_symbol in annotation_exact:
            selected_symbol = submitted_symbol
        else:
            exact_case_candidates = [
                symbol
                for symbol in candidate_symbols
                if symbol == submitted_symbol
            ]

            if exact_case_candidates:
                selected_symbol = exact_case_candidates[0]
            else:
                selected_symbol = sorted(
                    candidate_symbols,
                    key=lambda value: (
                        value.casefold(),
                        value,
                    ),
                )[0]

        synonym_entry = synonym_exact.get(
            selected_symbol,
            {},
        )
        annotation_entry = annotation_exact.get(
            selected_symbol,
            {},
        )

        official_symbol = str(
            synonym_entry.get(
                "official_symbol",
                annotation_entry.get(
                    "official_symbol",
                    submitted_symbol,
                ),
            )
        )

        flybase_id = str(
            synonym_entry.get(
                "flybase_id",
                annotation_entry.get(
                    "flybase_id",
                    "",
                ),
            )
        )

        symbol_synonyms = synonym_entry.get(
            "symbol_synonyms",
            [],
        )

        if not isinstance(symbol_synonyms, list):
            symbol_synonyms = []

        match_type = (
            "official_symbol"
            if official_symbol == submitted_symbol
            else "synonym"
        )

        output_rows.append(
            {
                "submitted_symbol": submitted_symbol,
                "official_symbol": official_symbol,
                "flybase_id": flybase_id,
                "current_fullname": str(
                    synonym_entry.get(
                        "current_fullname",
                        "",
                    )
                ),
                "symbol_synonyms": "|".join(
                    sorted(
                        {
                            value
                            for value in symbol_synonyms
                            if value
                        },
                        key=lambda value: (
                            value.casefold(),
                            value,
                        ),
                    )
                ),
                "secondary_flybase_ids": str(
                    annotation_entry.get(
                        "secondary_flybase_ids",
                        "",
                    )
                ),
                "annotation_id": str(
                    annotation_entry.get(
                        "annotation_id",
                        "",
                    )
                ),
                "match_type": match_type,
                "ambiguous": (
                    "Yes"
                    if len(candidate_symbols) > 1
                    else "No"
                ),
            }
        )

    OUTPUT_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "submitted_symbol",
        "official_symbol",
        "flybase_id",
        "current_fullname",
        "symbol_synonyms",
        "secondary_flybase_ids",
        "annotation_id",
        "match_type",
        "ambiguous",
    ]

    with OUTPUT_PATH.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Requested genes: {len(requested_genes)}")
    print(f"Index rows: {len(output_rows)}")
    print(
        "Resolved: "
        f"{len(output_rows) - len(unresolved)}"
    )
    print(f"Unresolved: {len(unresolved)}")
    print(f"Ambiguous: {len(ambiguous)}")

    if unresolved:
        print("Unresolved genes:")
        print(", ".join(unresolved))

    if ambiguous:
        print("Ambiguous genes:")
        print(", ".join(ambiguous))

    print(f"Saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
