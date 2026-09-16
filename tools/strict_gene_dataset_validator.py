from pathlib import Path
import re
import pandas as pd

ROOT = Path("results_2026_09_13_metadata")
INDEX = Path("src/dataset_finder/data/flybase/drosophila_gene_index.tsv")

CONFIG = {
    "RBP": ROOT / "Drosophila_RBP_Metadata_FINAL.xlsx",
    "TF": ROOT / "Drosophila_TF_Metadata_FINAL.xlsx",
}

OUT = ROOT / "STRICT_GENE_DATASET_VALIDATION_AUDIT.xlsx"

ASSAY_TECHNIQUES = {"ChIP_seq", "CUT_RUN", "CUT_TAG", "CLIP"}

GENE_CONTEXT = re.compile(
    r"\b("
    r"rnai|rna[- ]?i|knockdown|knock[- ]?down|knockout|knock[- ]?out|"
    r"mutant|mutation|deplet(?:e|ed|ion)|overexpress(?:ion|ed)?|"
    r"transgenic|crispr|cas9|sirna|dsrna|shrna|"
    r"anti[- ]?|antibody|chip|clip|cut&?run|cut&?tag|"
    r"flag|gfp|ha[- ]?tag|tagged|ip\b|immunoprecip"
    r")",
    re.I,
)

BAD_SHORT_CONTEXT = {
    "ac": [
        r"h3k27[- ]?ac",
        r"h3k9[- ]?ac",
        r"h4k\d+[- ]?ac",
        r"acetyl",
    ],
    "h": [
        r"\b\d+(?:\.\d+)?\s*h\b",
        r"\bhour(?:s)?\b",
    ],
}


def clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def norm(x):
    return re.sub(r"\s+", " ", clean(x)).strip()


def boundary_pattern(term):
    term = clean(term)
    if not term:
        return None

    return re.compile(
        r"(?<![A-Za-z0-9])" +
        re.escape(term) +
        r"(?![A-Za-z0-9])",
        re.I,
    )


def split_aliases(x):
    return [
        a.strip()
        for a in clean(x).split("|")
        if a.strip()
    ]


def build_gene_info(index):
    info = {}

    for _, r in index.iterrows():
        gene = clean(r["submitted_symbol"])

        aliases = []

        for value in [
            r.get("official_symbol"),
            r.get("flybase_id"),
            r.get("current_fullname"),
            r.get("annotation_id"),
        ]:
            value = clean(value)
            if value:
                aliases.append(value)

        aliases.extend(split_aliases(r.get("symbol_synonyms")))

        aliases = list(dict.fromkeys(aliases))

        strong = []
        weak = []

        for a in aliases:
            alnum = re.sub(r"[^A-Za-z0-9]", "", a)

            if (
                a.startswith("FBgn")
                or re.fullmatch(r"CG\d+", a, re.I)
                or len(alnum) >= 4
            ):
                strong.append(a)
            else:
                weak.append(a)

        info[gene] = {
            "official_symbol": clean(r.get("official_symbol")) or gene,
            "flybase_id": clean(r.get("flybase_id")),
            "full_name": clean(r.get("current_fullname")),
            "strong_aliases": list(dict.fromkeys(strong)),
            "weak_aliases": list(dict.fromkeys(weak)),
        }

    return info


def evidence_text(row):
    fields = [
        "Experiment Title",
        "Study Title",
        "Metadata Evidence",
        "Group Evidence",
    ]

    return " | ".join(
        norm(row.get(c))
        for c in fields
        if norm(row.get(c))
    )


def find_matches(text, aliases):
    hits = []

    for alias in aliases:
        p = boundary_pattern(alias)

        if p and p.search(text):
            hits.append(alias)

    return list(dict.fromkeys(hits))


def bad_short_collision(gene, text):
    for pattern in BAD_SHORT_CONTEXT.get(gene, []):
        if re.search(pattern, text, re.I):
            return True

    return False


def classify_sample(row, gene, gi):
    text = evidence_text(row)

    if not text:
        return "REVIEW", "No sample-level evidence text", ""

    strong_hits = find_matches(text, gi["strong_aliases"])
    weak_hits = find_matches(text, gi["weak_aliases"])

    if bad_short_collision(gene, text) and not strong_hits:
        return (
            "REJECTED",
            "Short-symbol collision with unrelated terminology",
            "; ".join(weak_hits),
        )

    technique = clean(row.get("Technique"))

    if strong_hits:
        if technique in ASSAY_TECHNIQUES:
            if GENE_CONTEXT.search(text):
                return (
                    "VALIDATED",
                    "Strong identifier/alias with assay-target context",
                    "; ".join(strong_hits),
                )

            return (
                "REVIEW",
                "Strong gene evidence present but assay-target role is unclear",
                "; ".join(strong_hits),
            )

        if technique == "RNA_seq":
            if GENE_CONTEXT.search(text):
                return (
                    "VALIDATED",
                    "Strong gene evidence with perturbation/genotype context",
                    "; ".join(strong_hits),
                )

            return (
                "REVIEW",
                "Gene is present but RNA-seq perturbation role is unclear",
                "; ".join(strong_hits),
            )

        return (
            "REVIEW",
            "Strong gene evidence but technique context is unresolved",
            "; ".join(strong_hits),
        )

    if weak_hits:
        return (
            "REVIEW",
            "Only short/weak gene alias evidence",
            "; ".join(weak_hits),
        )

    return (
        "REJECTED",
        "No gene-specific evidence in sample metadata",
        "",
    )


def validate(kind, path, gene_info):
    samples = pd.read_excel(path, sheet_name="Sample_Metadata")

    rows = []

    for _, r in samples.iterrows():
        gene = clean(r.get("Gene"))

        if gene not in gene_info:
            continue

        gi = gene_info[gene]

        status, reason, matched = classify_sample(
            r,
            gene,
            gi,
        )

        rows.append({
            "Panel": kind,
            "Gene": gene,
            "Official Symbol": gi["official_symbol"],
            "FlyBase ID": gi["flybase_id"],
            "Full Gene Name": gi["full_name"],
            "Technique": clean(r.get("Technique")),
            "Study Accession": clean(r.get("Study Accession")),
            "Experiment Accession": clean(r.get("Experiment Accession")),
            "Run Accessions": clean(r.get("Run Accessions")),
            "Experiment Title": clean(r.get("Experiment Title")),
            "Study Title": clean(r.get("Study Title")),
            "Original Group": clean(r.get("Group")),
            "Original Group Evidence": clean(r.get("Group Evidence")),
            "Validation Status": status,
            "Validation Reason": reason,
            "Matched Gene Evidence": matched,
        })

    sample_audit = pd.DataFrame(rows)

    study_rows = []

    keys = [
        "Panel",
        "Gene",
        "Official Symbol",
        "FlyBase ID",
        "Full Gene Name",
        "Technique",
        "Study Accession",
    ]

    for key, g in sample_audit.groupby(keys, dropna=False):
        statuses = set(g["Validation Status"])

        if "VALIDATED" in statuses:
            final = "VALIDATED"
        elif "REVIEW" in statuses:
            final = "REVIEW"
        else:
            final = "REJECTED"

        validated_n = int(
            g["Validation Status"].eq("VALIDATED").sum()
        )
        review_n = int(
            g["Validation Status"].eq("REVIEW").sum()
        )
        rejected_n = int(
            g["Validation Status"].eq("REJECTED").sum()
        )

        study_rows.append({
            **dict(zip(keys, key)),
            "Study Validation Status": final,
            "Validated Samples": validated_n,
            "Review Samples": review_n,
            "Rejected Samples": rejected_n,
            "Total Samples": len(g),
        })

    study_audit = pd.DataFrame(study_rows)

    return sample_audit, study_audit


def main():
    index = pd.read_csv(INDEX, sep="\t")
    gene_info = build_gene_info(index)

    all_samples = []
    all_studies = []

    for kind, path in CONFIG.items():
        sample_audit, study_audit = validate(
            kind,
            path,
            gene_info,
        )

        all_samples.append(sample_audit)
        all_studies.append(study_audit)

        print(f"\n===== {kind} =====")
        print("Sample-level:")
        print(
            sample_audit["Validation Status"]
            .value_counts()
            .to_string()
        )

        print("\nGene-study-level:")
        print(
            study_audit["Study Validation Status"]
            .value_counts()
            .to_string()
        )

        print("\nGenes with >=1 validated study:")
        print(
            study_audit.loc[
                study_audit["Study Validation Status"].eq("VALIDATED"),
                "Gene",
            ].nunique()
        )

    samples = pd.concat(all_samples, ignore_index=True)
    studies = pd.concat(all_studies, ignore_index=True)

    with pd.ExcelWriter(OUT, engine="openpyxl") as writer:
        studies.to_excel(
            writer,
            sheet_name="Gene_Study_Audit",
            index=False,
        )

        samples.to_excel(
            writer,
            sheet_name="Sample_Audit",
            index=False,
        )

        for ws in writer.book.worksheets:
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions

    print("\nCREATED:", OUT)
    print("\nIMPORTANT: No final workbook was overwritten.")


if __name__ == "__main__":
    main()
