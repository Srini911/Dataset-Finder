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


def test_extracts_genotype_treatment_and_disease() -> None:
    metadata = extract_biological_metadata(
        "Adult male brain RNA-seq from homozygous mutant flies",
        "oxidative stress model of Parkinson disease",
    )

    assert metadata.genotype == "homozygous"
    assert metadata.treatment == "oxidative stress"
    assert metadata.disease == "Parkinson disease"


def test_extracts_wild_type_heat_shock() -> None:
    metadata = extract_biological_metadata(
        "Wild-type Drosophila exposed to heat shock",
    )

    assert metadata.genotype == "wild type"
    assert metadata.treatment == "heat shock"


def test_extracts_transgenic_infection_model() -> None:
    metadata = extract_biological_metadata(
        "Transgenic GAL4 flies after bacterial infection",
    )

    assert metadata.genotype == "transgenic"
    assert metadata.treatment == "infection"


def test_extracts_alzheimer_disease() -> None:
    metadata = extract_biological_metadata(
        "Drosophila model of Alzheimer's disease",
    )

    assert metadata.disease == "Alzheimer disease"


def test_specific_zygosity_precedes_generic_mutant() -> None:
    metadata = extract_biological_metadata(
        "heterozygous mutant Drosophila adults",
    )

    assert metadata.genotype == "heterozygous"


def test_extracts_plural_brain_and_depletion() -> None:
    metadata = extract_biological_metadata(
        "Effect of Orb2 depletion on mRNA expression "
        "in Drosophila larval brains"
    )

    assert metadata.tissue == "brain"
    assert metadata.developmental_stage == "larva"
    assert metadata.perturbation == "knockdown"


def test_maternal_to_zygotic_transition_maps_to_embryo() -> None:
    metadata = extract_biological_metadata(
        "Orb2 regulation during the maternal-to-zygotic transition"
    )

    assert metadata.developmental_stage == "embryo"


def test_gal4_driver_genotype_is_transgenic() -> None:
    metadata = extract_biological_metadata(
        {
            "characteristics": {
                "genotype": [
                    "inscGAL4>wRNAi",
                    "inscGAL4>orb2RNAi",
                ]
            }
        }
    )

    assert metadata.genotype == "transgenic"
