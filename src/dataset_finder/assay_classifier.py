"""Classify functional-genomics records by experimental technique."""

from __future__ import annotations

import re
from collections.abc import Iterable

TECHNIQUE_PATTERNS: dict[str, tuple[str, ...]] = {
    "CUT_RUN": (
        r"\bcut\s*&\s*run\b",
        r"\bcut[-\s]?and[-\s]?run\b",
        r"\bcutnrun\b",
        r"cleavage under targets and release using nuclease",
    ),
    "CUT_TAG": (
        r"\bcut\s*&\s*tag\b",
        r"\bcut[-\s]?and[-\s]?tag\b",
        r"\bcutntag\b",
        r"cleavage under targets and tagmentation",
    ),
    "eCLIP": (
        r"\beclip\b",
        r"enhanced clip",
        r"enhanced crosslinking immunoprecipitation",
    ),
    "iCLIP": (
        r"\biclip\b",
        r"individual nucleotide resolution clip",
    ),
    "PAR_CLIP": (
        r"\bpar[-\s]?clip\b",
        r"photoactivatable ribonucleoside",
    ),
    "HITS_CLIP": (
        r"\bhits[-\s]?clip\b",
        r"high-throughput sequencing of rna isolated",
    ),
    "CLIP": (
        r"\bclip[-\s]?seq\b",
        r"\bclip sequencing\b",
        r"crosslinking immunoprecipitation",
    ),
    "ChIP_seq": (
        r"\bchip[-\s]?seq\b",
        r"\bchip sequencing\b",
        r"chromatin immunoprecipitation sequencing",
        r"genome binding/occupancy profiling by high throughput sequencing",
    ),
    "ATAC_seq": (
        r"\batac[-\s]?seq\b",
        r"assay for transposase-accessible chromatin",
    ),
    "scRNA_seq": (
        r"\bscrna[-\s]?seq\b",
        r"\bsingle[-\s]?cell rna",
        r"\bsingle cell transcriptom",
        r"\b10x genomics\b",
    ),
    "snRNA_seq": (
        r"\bsnrna[-\s]?seq\b",
        r"\bsingle[-\s]?nucleus rna",
        r"\bsingle nucleus transcriptom",
    ),
    "Spatial": (
        r"\bspatial transcriptom",
        r"\bvisium\b",
        r"\bmerfish\b",
        r"\bslide[-\s]?seq\b",
    ),
    "RNA_seq": (
        r"\brna[-\s]?seq\b",
        r"\btranscriptome sequencing\b",
        r"\btranscriptomic sequencing\b",
        r"\bmrna sequencing\b",
        r"expression profiling by high throughput sequencing",
    ),
    "Microarray": (
        r"\bmicroarray\b",
        r"\bexpression array\b",
        r"\barray profiling\b",
        r"expression profiling by array",
        r"genome binding/occupancy profiling by array",
    ),
    "Proteomics": (
        r"\bproteom",
        r"\bmass spectrom",
        r"\blc[-\s]?ms\b",
    ),
}


TECHNIQUE_ORDER = tuple(TECHNIQUE_PATTERNS) + ("Other_Assays",)


def normalize_search_text(values: Iterable[object]) -> str:
    """Combine arbitrary metadata values into normalized searchable text."""
    return " ".join(
        str(value).strip()
        for value in values
        if value is not None and str(value).strip()
    ).casefold()


def classify_technique(*values: object) -> str:
    """Return the most specific matching experimental technique."""
    text = normalize_search_text(values)

    for technique, patterns in TECHNIQUE_PATTERNS.items():
        if any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns):
            return technique

    return "Other_Assays"
