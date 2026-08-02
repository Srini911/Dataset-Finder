"""Tests for biological metadata extraction."""

from dataset_finder.metadata import (
    BiologicalMetadata,
    extract_biological_metadata,
)


def test_biological_metadata_defaults_are_empty() -> None:
    metadata = BiologicalMetadata()

    assert metadata.tissue == ""
    assert metadata.cell_type == ""
    assert metadata.developmental_stage == ""
    assert metadata.sex == ""
    assert metadata.genotype == ""
    assert metadata.strain == ""
    assert metadata.treatment == ""
    assert metadata.control_status == ""
    assert metadata.disease == ""
    assert metadata.time_point == ""
    assert metadata.perturbation == ""


def test_extracts_core_drosophila_metadata() -> None:
    metadata = extract_biological_metadata(
        "RNA-seq of adult male Drosophila brain neurons",
        "w1118 control collected after 24 hours",
    )

    assert metadata.tissue == "brain"
    assert metadata.cell_type == "neuron"
    assert metadata.developmental_stage == "adult"
    assert metadata.sex == "male"
    assert metadata.strain == "w1118"
    assert metadata.control_status == "control"
    assert metadata.time_point == "24 hours"


def test_extracts_perturbation_and_tissue() -> None:
    metadata = extract_biological_metadata(
        "Orb2 RNAi knockdown in larval brain",
    )

    assert metadata.tissue == "brain"
    assert metadata.developmental_stage == "larva"
    assert metadata.perturbation == "RNAi"


def test_extracts_nested_metadata_values() -> None:
    metadata = extract_biological_metadata(
        {
            "source_name": "adult female ovary",
            "characteristics": [
                "strain: Canton-S",
                "treatment: overexpression",
            ],
        }
    )

    assert metadata.tissue == "ovary"
    assert metadata.developmental_stage == "adult"
    assert metadata.sex == "female"
    assert metadata.strain == "Canton-S"
    assert metadata.perturbation == "overexpression"


def test_unknown_metadata_remains_empty() -> None:
    metadata = extract_biological_metadata(
        "functional genomics experiment",
    )

    assert metadata == BiologicalMetadata()
