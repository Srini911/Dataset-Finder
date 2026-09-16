from pathlib import Path
import pandas as pd

ROOT = Path("results_2026_09_13_metadata")

AUDIT = ROOT / "STRICT_GENE_DATASET_VALIDATION_AUDIT.xlsx"
INDEX = Path("src/dataset_finder/data/flybase/drosophila_gene_index.tsv")

PANELS = {
    "RBP": {
        "genes": Path("src/dataset_finder/data/gene_sets/drosophila_rbps.txt"),
        "metadata": ROOT / "Drosophila_RBP_Metadata_FINAL.xlsx",
        "comparisons": ROOT / "Drosophila_RBP_Control_Experiment_Sex_FINAL.xlsx",
        "output": ROOT / "Drosophila_RBP_TRUTH_VALIDATED_FINAL.xlsx",
    },
    "TF": {
        "genes": Path("src/dataset_finder/data/gene_sets/drosophila_tfs.txt"),
        "metadata": ROOT / "Drosophila_TF_Metadata_FINAL.xlsx",
        "comparisons": ROOT / "Drosophila_TF_Control_Experiment_Sex_FINAL.xlsx",
        "output": ROOT / "Drosophila_TF_TRUTH_VALIDATED_FINAL.xlsx",
    },
}

TECH_MAP = {
    "CUT_RUN": "CUT&RUN",
    "CUT_TAG": "CUT&Tag",
    "ChIP_seq": "ChIP-seq",
    "CLIP": "CLIP",
    "RNA_seq": "RNA-seq",
}

TECHNIQUES = [
    "CUT&RUN",
    "CUT&Tag",
    "ChIP-seq",
    "CLIP",
    "RNA-seq",
]

FINAL_COLUMNS = [
    "Gene Name",
    "Gene Symbol",
    "FlyBase ID",
    "Full Gene Name",
    "Technique",
    "Study Name",
    "GEO/SRA Number",
    "Experiment Accession Code",
    "Control Accession Code",
    "Sex Label",
    "Validation Status",
]

REVIEW_COLUMNS = [
    "Gene",
    "Official Symbol",
    "FlyBase ID",
    "Full Gene Name",
    "Technique",
    "Study Accession",
    "Study Validation Status",
    "Validated Samples",
    "Review Samples",
    "Rejected Samples",
    "Total Samples",
]

def clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()

def read_genes(path):
    return [
        x.strip()
        for x in path.read_text().splitlines()
        if x.strip()
    ]

def load_index():
    df = pd.read_csv(INDEX, sep="\t")
    df["submitted_symbol"] = df["submitted_symbol"].astype(str).str.strip()
    return df.set_index("submitted_symbol", drop=False)

def annotation(gene, idx):
    if gene not in idx.index:
        return gene, "", ""

    r = idx.loc[gene]

    if isinstance(r, pd.DataFrame):
        r = r.iloc[0]

    return (
        clean(r.get("official_symbol")) or gene,
        clean(r.get("flybase_id")),
        clean(r.get("current_fullname")),
    )

def main():
    audit = pd.read_excel(AUDIT, sheet_name="Gene_Study_Audit")

    audit["Gene"] = audit["Gene"].fillna("").astype(str).str.strip()
    audit["Study Accession"] = (
        audit["Study Accession"].fillna("").astype(str).str.strip()
    )

    idx = load_index()

    for panel, cfg in PANELS.items():
        genes = read_genes(cfg["genes"])

        panel_audit = audit[audit["Panel"].eq(panel)].copy()

        validated = panel_audit[
            panel_audit["Study Validation Status"].eq("VALIDATED")
        ].copy()

        review = panel_audit[
            panel_audit["Study Validation Status"].eq("REVIEW")
        ].copy()

        rejected = panel_audit[
            panel_audit["Study Validation Status"].eq("REJECTED")
        ].copy()

        metadata = pd.read_excel(
            cfg["metadata"],
            sheet_name="All_Datasets",
        )

        metadata["Gene"] = (
            metadata["Gene"].fillna("").astype(str).str.strip()
        )

        metadata["Study Accession"] = (
            metadata["Study Accession"]
            .fillna("")
            .astype(str)
            .str.strip()
        )

        comparisons = pd.read_excel(
            cfg["comparisons"],
            sheet_name="Summary",
        )

        for c in [
            "Gene Symbol",
            "Technique",
            "GEO/SRA Number",
        ]:
            comparisons[c] = (
                comparisons[c].fillna("").astype(str).str.strip()
            )

        valid_keys = set(
            zip(
                validated["Gene"],
                validated["Technique"],
                validated["Study Accession"],
            )
        )

        rows = []

        for gene in genes:
            official, flybase, fullname = annotation(gene, idx)

            gene_valid = validated[
                validated["Gene"].eq(gene)
            ].copy()

            if gene_valid.empty:
                rows.append({
                    "Gene Name": gene,
                    "Gene Symbol": official,
                    "FlyBase ID": flybase,
                    "Full Gene Name": fullname,
                    "Technique": "",
                    "Study Name": "",
                    "GEO/SRA Number": "",
                    "Experiment Accession Code": "",
                    "Control Accession Code": "",
                    "Sex Label": "Unspecified",
                    "Validation Status": "No validated dataset found",
                })
                continue

            for _, vr in gene_valid.iterrows():
                raw_tech = clean(vr["Technique"])
                tech = TECH_MAP.get(raw_tech, raw_tech)
                study = clean(vr["Study Accession"])

                m = metadata[
                    metadata["Gene"].eq(gene)
                    & metadata["Technique"].eq(raw_tech)
                    & metadata["Study Accession"].eq(study)
                ].copy()

                study_name = ""

                if not m.empty:
                    names = [
                        clean(x)
                        for x in m["Study Title"].tolist()
                        if clean(x)
                    ]

                    if names:
                        study_name = names[0]

                comp = comparisons[
                    comparisons["Gene Symbol"].eq(official)
                    & comparisons["Technique"].eq(tech)
                    & comparisons["GEO/SRA Number"].eq(study)
                ].copy()

                if comp.empty:
                    rows.append({
                        "Gene Name": gene,
                        "Gene Symbol": official,
                        "FlyBase ID": flybase,
                        "Full Gene Name": fullname,
                        "Technique": tech,
                        "Study Name": study_name,
                        "GEO/SRA Number": study,
                        "Experiment Accession Code": "",
                        "Control Accession Code": "",
                        "Sex Label": "Unspecified",
                        "Validation Status": "Validated dataset",
                    })

                else:
                    for _, cr in comp.iterrows():
                        rows.append({
                            "Gene Name": gene,
                            "Gene Symbol": official,
                            "FlyBase ID": flybase,
                            "Full Gene Name": fullname,
                            "Technique": tech,
                            "Study Name": study_name or clean(
                                cr.get("Study Name")
                            ),
                            "GEO/SRA Number": study,
                            "Experiment Accession Code": clean(
                                cr.get("Experiment Accession Code")
                            ),
                            "Control Accession Code": clean(
                                cr.get("Control Accession Code")
                            ),
                            "Sex Label": clean(
                                cr.get("Sex Label")
                            ) or "Unspecified",
                            "Validation Status": (
                                "Validated comparison"
                            ),
                        })

        final = pd.DataFrame(rows, columns=FINAL_COLUMNS)

        final = final.drop_duplicates().reset_index(drop=True)

        coverage_rows = []

        for gene in genes:
            official, flybase, fullname = annotation(gene, idx)

            gv = validated[validated["Gene"].eq(gene)]
            gr = review[review["Gene"].eq(gene)]
            gj = rejected[rejected["Gene"].eq(gene)]

            if not gv.empty:
                status = "Validated dataset"
            elif not gr.empty:
                status = "Needs review"
            else:
                status = "No validated dataset found"

            coverage_rows.append({
                "Gene Name": gene,
                "Gene Symbol": official,
                "FlyBase ID": flybase,
                "Full Gene Name": fullname,
                "Coverage Status": status,
                "Validated Studies": len(gv),
                "Review Studies": len(gr),
                "Rejected Studies": len(gj),
            })

        coverage = pd.DataFrame(coverage_rows)

        review_out = review[
            [c for c in REVIEW_COLUMNS if c in review.columns]
        ].copy()

        rejected_out = rejected[
            [c for c in REVIEW_COLUMNS if c in rejected.columns]
        ].copy()

        with pd.ExcelWriter(
            cfg["output"],
            engine="openpyxl",
        ) as writer:

            coverage.to_excel(
                writer,
                sheet_name="Gene_Coverage",
                index=False,
            )

            final.to_excel(
                writer,
                sheet_name="All_Genes",
                index=False,
            )

            for tech in TECHNIQUES:
                x = final[final["Technique"].eq(tech)].copy()

                x.to_excel(
                    writer,
                    sheet_name=tech,
                    index=False,
                )

            review_out.to_excel(
                writer,
                sheet_name="Needs_Review",
                index=False,
            )

            rejected_out.to_excel(
                writer,
                sheet_name="Rejected_Associations",
                index=False,
            )

            for ws in writer.book.worksheets:
                ws.freeze_panes = "A2"
                ws.auto_filter.ref = ws.dimensions

        print(f"\n===== {panel} =====")
        print("Panel genes:", len(genes))
        print(
            "Genes represented:",
            final["Gene Name"].nunique(),
        )
        print(
            "Genes with validated datasets:",
            validated["Gene"].nunique(),
        )
        print(
            "Genes needing review only:",
            len(
                set(review["Gene"])
                - set(validated["Gene"])
            ),
        )
        print(
            "Genes with no validated/review study:",
            len(
                set(genes)
                - set(validated["Gene"])
                - set(review["Gene"])
            ),
        )
        print("Final rows:", len(final))
        print("Created:", cfg["output"])

if __name__ == "__main__":
    main()
