import re
from pathlib import Path

import pandas as pd


ROOT = Path("results_2026_09_13_metadata")

FILES = {
    "RBP": ROOT / "Drosophila_RBP_Metadata_FINAL.xlsx",
    "TF": ROOT / "Drosophila_TF_Metadata_FINAL.xlsx",
}

OUT = {
    "RBP": ROOT / "Drosophila_RBP_Control_Experiment_Sex_FINAL.xlsx",
    "TF": ROOT / "Drosophila_TF_Control_Experiment_Sex_FINAL.xlsx",
}

TECHS = ["CUT_RUN", "CUT_TAG", "ChIP_seq", "CLIP", "RNA_seq"]

LABEL = {
    "CUT_RUN": "CUT&RUN",
    "CUT_TAG": "CUT&Tag",
    "ChIP_seq": "ChIP-seq",
    "CLIP": "CLIP",
    "RNA_seq": "RNA-seq",
}

COLS = [
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
]


def text(x):
    return "" if pd.isna(x) else str(x).strip()


def title(row):
    return text(row.get("Experiment Title", ""))


def run_list(frame):
    values = []

    for value in frame["Run Accessions"]:
        values.extend(
            x.upper()
            for x in re.findall(
                r"(?:SRR|ERR|DRR)\d+",
                text(value),
                re.I,
            )
        )

    if not values:
        for value in frame["Experiment Accession"]:
            values.extend(
                x.upper()
                for x in re.findall(
                    r"(?:SRX|ERX|DRX)\d+",
                    text(value),
                    re.I,
                )
            )

    return ", ".join(dict.fromkeys(values))


def normalize_sex(value):
    value = text(value).lower()

    if value == "female":
        return "Female"

    if value == "male":
        return "Male"

    if value == "mixed":
        return "Mixed"

    return "Unspecified"


def frame_sex(frame):
    vals = {
        normalize_sex(x)
        for x in frame["Sex"]
        if normalize_sex(x) != "Unspecified"
    }

    if vals == {"Female"}:
        return "Female"

    if vals == {"Male"}:
        return "Male"

    if vals == {"Mixed"}:
        return "Mixed"

    if vals == {"Female", "Male"}:
        return "Mixed"

    return "Unspecified"


def extract_timepoint(value):
    value = text(value).lower()

    patterns = [
        r"\b(\d+(?:\.\d+)?\s*-\s*\d+(?:\.\d+)?\s*h)\b",
        r"\b(\d+(?:\.\d+)?\s*h)\b",
        r"\b(\d+\s*-\s*\d+\s*min)\b",
    ]

    for pattern in patterns:
        m = re.search(pattern, value, re.I)
        if m:
            return re.sub(r"\s+", "", m.group(1).lower())

    return ""


def is_rna_sample(row):
    t = title(row).lower()
    strategy = text(row.get("Library Strategy", "")).lower()

    bad = [
        r"\bmerip\b",
        r"\brip[- ]?seq\b",
        r"\bripseq\b",
        r"\bclip\b",
        r"\bicl[i]?p\b",
        r"\binput\b",
        r"\bigg\b",
        r"\bchip\b",
    ]

    if any(re.search(p, t, re.I) for p in bad):
        return False

    return "rna-seq" in t or "rna-seq" in strategy or "rna seq" in t


def is_chip_experiment(row):
    t = title(row).lower()

    if re.search(r"\binput\b|\bigg\b", t, re.I):
        return False

    return bool(
        re.search(
            r"\bchip\b|chip-seq|chip-nexus",
            t,
            re.I,
        )
    )


def target_supported(row, gene, technique):
    t = title(row).lower()
    gene_l = text(gene).lower()

    if technique == "RNA_seq":
        return True

    if technique == "ChIP_seq":
        if not is_chip_experiment(row):
            return False

        # Avoid unsafe matching of one-character gene symbols.
        if len(gene_l) <= 1:
            return False

        official = text(
            row.get("Official Symbol", "")
        ).lower()

        candidates = {
            x for x in [gene_l, official]
            if len(x) > 1
        }

        return any(
            re.search(
                rf"(?<![a-z0-9]){re.escape(x)}(?![a-z0-9])",
                t,
                re.I,
            )
            for x in candidates
        )

    return True

def compatible_sex(exp, controls):
    exp_sex = frame_sex(exp)

    if exp_sex not in {"Female", "Male"}:
        return controls

    matched = controls[
        controls["Sex"].apply(normalize_sex) == exp_sex
    ]

    return matched if not matched.empty else controls.iloc[0:0]


def split_by_time(exp, controls):
    exp = exp.copy()
    controls = controls.copy()

    exp["_time"] = exp["Experiment Title"].apply(extract_timepoint)
    controls["_time"] = controls["Experiment Title"].apply(extract_timepoint)

    exp_times = [x for x in exp["_time"].unique() if x]
    ctrl_times = {x for x in controls["_time"].unique() if x}

    if not exp_times:
        return [(exp, controls)]

    pairs = []

    for tp in exp_times:
        e = exp[exp["_time"] == tp]
        c = controls[controls["_time"] == tp]

        if not e.empty and not c.empty:
            pairs.append((e, c))

    if pairs:
        return pairs

    if ctrl_times:
        return []

    return [(exp, controls)]


def prepare_rna(block):
    block = block[
        block.apply(is_rna_sample, axis=1)
    ].copy()

    exp = block[block["Group"] == "Experiment"].copy()
    ctrl = block[block["Group"] == "Control"].copy()

    if exp.empty or ctrl.empty:
        return []

    ctrl = compatible_sex(exp, ctrl)

    if ctrl.empty:
        return []

    return split_by_time(exp, ctrl)


def prepare_chip(block, gene):
    exp = block[block["Group"] == "Experiment"].copy()
    ctrl = block[block["Group"] == "Control"].copy()

    if exp.empty or ctrl.empty:
        return []

    exp = exp[
        exp.apply(
            lambda r: target_supported(
                r,
                gene,
                "ChIP_seq",
            ),
            axis=1,
        )
    ].copy()

    if exp.empty:
        return []

    ctrl = compatible_sex(exp, ctrl)

    if ctrl.empty:
        return []

    return split_by_time(exp, ctrl)


def prepare_other(block):
    exp = block[block["Group"] == "Experiment"].copy()
    ctrl = block[block["Group"] == "Control"].copy()

    if exp.empty or ctrl.empty:
        return []

    ctrl = compatible_sex(exp, ctrl)

    if ctrl.empty:
        return []

    return [(exp, ctrl)]


def make_row(block, exp, ctrl, technique):
    first = block.iloc[0]

    sex = frame_sex(exp)

    if sex == "Unspecified":
        sex = frame_sex(ctrl)

    study_name = text(first.get("Study Title", ""))

    if not study_name:
        study_name = text(first.get("Title", ""))

    return {
        "Gene Name": text(first.get("Gene", "")),
        "Gene Symbol": text(first.get("Official Symbol", "")),
        "FlyBase ID": text(first.get("FlyBase ID", "")),
        "Full Gene Name": text(first.get("Full Name", "")),
        "Technique": LABEL[technique],
        "Study Name": study_name,
        "GEO/SRA Number": text(
            first.get("Study Accession", "")
        ),
        "Experiment Accession Code": run_list(exp),
        "Control Accession Code": run_list(ctrl),
        "Sex Label": sex,
    }


def build(kind):
    print(f"\nReading {FILES[kind]}")

    df = pd.read_excel(
        FILES[kind],
        sheet_name="Sample_Metadata",
    )

    df = df[
        df["Group"].isin(["Control", "Experiment"])
    ].copy()

    rows = []

    grouped = df.groupby(
        ["Gene", "Technique", "Study Accession"],
        dropna=False,
        sort=True,
    )

    for (gene, technique, study), block in grouped:
        if technique not in TECHS:
            continue

        excluded = {
            ("TF", "EcR", "ChIP_seq", "SRP226806"),
            ("TF", "h", "ChIP_seq", "SRP656693"),
            ("TF", "z", "RNA_seq", "SRP484798"),
            ("TF", "z", "RNA_seq", "SRP527914"),
        }

        if (
            kind,
            text(gene),
            text(technique),
            text(study),
        ) in excluded:
            continue

        groups = set(block["Group"])

        if not {"Control", "Experiment"} <= groups:
            continue

        if technique == "RNA_seq":
            pairs = prepare_rna(block)

        elif technique == "ChIP_seq":
            pairs = prepare_chip(block, gene)

        else:
            pairs = prepare_other(block)

        for exp, ctrl in pairs:
            if exp.empty or ctrl.empty:
                continue

            row = make_row(
                block,
                exp,
                ctrl,
                technique,
            )

            if (
                row["Experiment Accession Code"]
                and row["Control Accession Code"]
            ):
                rows.append(row)

    result = pd.DataFrame(rows, columns=COLS)

    if not result.empty:
        result = result.drop_duplicates().sort_values(
            [
                "Technique",
                "Gene Name",
                "GEO/SRA Number",
                "Experiment Accession Code",
            ],
            kind="stable",
        )

    with pd.ExcelWriter(
        OUT[kind],
        engine="openpyxl",
    ) as writer:
        result.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

        for technique in TECHS:
            label = LABEL[technique]

            sub = result[
                result["Technique"] == label
            ].copy()

            sub.to_excel(
                writer,
                sheet_name=label,
                index=False,
            )

            ws = writer.book[label]
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions

        ws = writer.book["Summary"]
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

        for ws in writer.book.worksheets:
            widths = {
                "A": 18,
                "B": 18,
                "C": 18,
                "D": 36,
                "E": 14,
                "F": 65,
                "G": 18,
                "H": 55,
                "I": 55,
                "J": 15,
            }

            for col, width in widths.items():
                ws.column_dimensions[col].width = width

    print("CREATED:", OUT[kind])
    print("Rows:", len(result))
    print(
        "Genes:",
        result["Gene Name"].nunique()
        if not result.empty
        else 0,
    )

    if not result.empty:
        print("\nTechniques:")
        print(result["Technique"].value_counts().to_string())

        print("\nSex:")
        print(result["Sex Label"].value_counts().to_string())


for kind in ["RBP", "TF"]:
    build(kind)
