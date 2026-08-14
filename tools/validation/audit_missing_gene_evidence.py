"""Audit remaining legacy misses for explicit gene evidence."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from dataset_finder.flybase_resolver import FlyBaseResolver


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def distinctive(term: str) -> bool:
    compact = re.sub(r"[^A-Za-z0-9]+", "", term)

    if re.fullmatch(r"FBgn\d+", term, flags=re.IGNORECASE):
        return True

    if re.fullmatch(r"CG\d+", term, flags=re.IGNORECASE):
        return True

    return len(compact) >= 4


def contains_term(text: str, term: str) -> bool:
    if not term:
        return False

    pattern = (
        rf"(?<![A-Za-z0-9])"
        rf"{re.escape(term)}"
        rf"(?![A-Za-z0-9])"
    )

    return bool(
        re.search(
            pattern,
            text,
            flags=re.IGNORECASE,
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--missing-report",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--resolved-old",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    missing = pd.read_excel(
        args.missing_report,
        sheet_name="True_Missing",
    )

    resolved = pd.read_excel(
        args.resolved_old,
        sheet_name="Resolved_Old_Results",
    )

    resolver = FlyBaseResolver()

    resolved["gene_key"] = (
        resolved["Gene"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
    )

    resolved["study_key"] = (
        resolved["study_accession"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

    rows: list[dict[str, object]] = []

    for _, target in missing.iterrows():
        gene = clean(target["Gene"])
        study = clean(
            target["Expected Study Accession"]
        ).upper()

        gene_info = resolver.resolve(gene)

        terms = [
            gene,
            gene_info.official_symbol,
            gene_info.flybase_id,
            *gene_info.secondary_flybase_ids,
            gene_info.annotation_id,
            gene_info.current_fullname,
            *gene_info.synonyms,
        ]

        unique_terms: list[str] = []
        seen: set[str] = set()

        for term in terms:
            term = clean(term)

            if not term:
                continue

            key = term.casefold()

            if key in seen:
                continue

            seen.add(key)
            unique_terms.append(term)

        old_rows = resolved[
            (
                resolved["gene_key"]
                == gene.casefold()
            )
            & (
                resolved["study_key"]
                == study
            )
        ]

        old_text_parts: list[str] = []

        for _, old_row in old_rows.iterrows():
            for column in (
                "study_title",
                "raw_title",
                "library_strategy",
                "library_source",
                "library_selection",
                "all_detected_accessions",
            ):
                if column in old_row.index:
                    value = clean(
                        old_row[column]
                    )

                    if value:
                        old_text_parts.append(value)

        text = " | ".join(old_text_parts)

        matched_terms = [
            term
            for term in unique_terms
            if contains_term(text, term)
        ]

        strong_matches = [
            term
            for term in matched_terms
            if distinctive(term)
        ]

        if strong_matches:
            verdict = "EXPLICIT_STRONG_GENE_EVIDENCE"
        elif matched_terms:
            verdict = "ONLY_SHORT_OR_AMBIGUOUS_EVIDENCE"
        elif old_rows.empty:
            verdict = "NO_MATCHING_OLD_RESOLVED_ROW"
        else:
            verdict = "NO_EXPLICIT_GENE_EVIDENCE"

        rows.append(
            {
                "Gene Set": target["Gene Set"],
                "Gene": gene,
                "Technique": target["Technique"],
                "Expected Study Accession": study,
                "Official Symbol": (
                    gene_info.official_symbol
                ),
                "FlyBase ID": gene_info.flybase_id,
                "Annotation ID": (
                    gene_info.annotation_id
                ),
                "Matched Gene Terms": (
                    "; ".join(matched_terms)
                ),
                "Strong Matched Terms": (
                    "; ".join(strong_matches)
                ),
                "Evidence Verdict": verdict,
                "Study Title": clean(
                    target["Study Title"]
                ),
                "Old Metadata Text": text,
            }
        )

    audit = pd.DataFrame(rows)

    with pd.ExcelWriter(
        args.output,
        engine="xlsxwriter",
    ) as writer:
        audit.to_excel(
            writer,
            sheet_name="Missing_Gene_Evidence",
            index=False,
        )

        summary = (
            audit.groupby(
                [
                    "Gene Set",
                    "Evidence Verdict",
                ]
            )
            .size()
            .reset_index(name="Count")
        )

        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

    print("\n===== REMAINING MISS GENE-EVIDENCE AUDIT =====\n")

    print(
        audit["Evidence Verdict"]
        .value_counts()
        .to_string()
    )

    print("\nBY GENE SET\n")

    print(
        pd.crosstab(
            audit["Gene Set"],
            audit["Evidence Verdict"],
        ).to_string()
    )

    print(
        "\nSTRONG GENE EVIDENCE — THESE ARE REAL "
        "PIPELINE RECOVERY TARGETS\n"
    )

    strong = audit[
        audit["Evidence Verdict"]
        == "EXPLICIT_STRONG_GENE_EVIDENCE"
    ]

    if strong.empty:
        print("None")
    else:
        print(
            strong[
                [
                    "Gene Set",
                    "Gene",
                    "Technique",
                    "Expected Study Accession",
                    "Strong Matched Terms",
                    "Study Title",
                ]
            ].to_string(index=False)
        )

    print(
        "\nNO EXPLICIT GENE EVIDENCE — "
        "DO NOT FORCE THESE INTO THE PIPELINE\n"
    )

    weak = audit[
        audit["Evidence Verdict"]
        == "NO_EXPLICIT_GENE_EVIDENCE"
    ]

    if weak.empty:
        print("None")
    else:
        print(
            weak[
                [
                    "Gene Set",
                    "Gene",
                    "Technique",
                    "Expected Study Accession",
                    "Study Title",
                ]
            ].to_string(index=False)
        )

    print("\nReport:", args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
