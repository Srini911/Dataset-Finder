"""Tests for experimental-technique classification."""

import pytest

from dataset_finder.assay_classifier import classify_technique


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("single-cell RNA-seq of adult fly brain", "scRNA_seq"),
        ("single nucleus RNA sequencing", "snRNA_seq"),
        ("enhanced eCLIP experiment", "eCLIP"),
        ("PAR-CLIP binding analysis", "PAR_CLIP"),
        ("HITS-CLIP dataset", "HITS_CLIP"),
        ("CLIP-seq experiment", "CLIP"),
        ("CUT&RUN profiling", "CUT_RUN"),
        ("CUT and Tag assay", "CUT_TAG"),
        ("ChIP-seq experiment", "ChIP_seq"),
        ("ATAC-seq chromatin accessibility", "ATAC_seq"),
        ("bulk RNA-seq transcriptome", "RNA_seq"),
        ("spatial transcriptomics", "Spatial"),
        ("gene expression microarray", "Microarray"),
        ("LC-MS proteomics", "Proteomics"),
        ("unclassified functional assay", "Other_Assays"),
    ],
)
def test_classify_technique(text: str, expected: str) -> None:
    assert classify_technique(text) == expected


def test_specific_clip_type_precedes_generic_clip() -> None:
    assert classify_technique("eCLIP crosslinking immunoprecipitation") == "eCLIP"


def test_single_cell_precedes_generic_rna_seq() -> None:
    assert classify_technique("single-cell RNA-seq") == "scRNA_seq"


@pytest.mark.parametrize(
    ("study_type", "expected"),
    [
        (
            "Expression profiling by high throughput sequencing",
            "RNA_seq",
        ),
        (
            "Genome binding/occupancy profiling by high throughput sequencing",
            "ChIP_seq",
        ),
        (
            "Expression profiling by array",
            "Microarray",
        ),
    ],
)
def test_classify_geo_controlled_study_types(
    study_type: str,
    expected: str,
) -> None:
    assert classify_technique(study_type) == expected
