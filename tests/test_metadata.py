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
    assert metadata.perturbation == "RNAi; knockdown"


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


def test_male_and_female_samples_are_mixed() -> None:
    metadata = extract_biological_metadata(
        {
            "samples": [
                {"sex": "male"},
                {"sex": "female"},
            ]
        }
    )

    assert metadata.sex == "mixed"


def test_control_and_treated_samples_are_mixed() -> None:
    metadata = extract_biological_metadata(
        {
            "sample_titles": [
                "untreated control replicate",
                "drug-treated replicate",
            ]
        }
    )

    assert metadata.control_status == "mixed"


def test_wild_type_and_mutant_genotypes_are_mixed() -> None:
    metadata = extract_biological_metadata(
        {
            "characteristics": {
                "genotype": [
                    "wild type",
                    "orb2 mutant",
                ]
            }
        }
    )

    assert metadata.genotype == "mixed"


def test_single_genotype_is_not_marked_mixed() -> None:
    metadata = extract_biological_metadata(
        {
            "characteristics": {
                "genotype": [
                    "inscGAL4>orb2RNAi",
                    "inscGAL4>orb2RNAi",
                ]
            }
        }
    )

    assert metadata.genotype == "transgenic"


def test_single_sex_is_not_marked_mixed() -> None:
    metadata = extract_biological_metadata(
        "adult female brain samples",
    )

    assert metadata.sex == "female"


def test_zygosity_with_mutant_word_is_not_mixed() -> None:
    homozygous = extract_biological_metadata(
        "homozygous mutant flies",
    )
    heterozygous = extract_biological_metadata(
        "heterozygous mutant flies",
    )

    assert homozygous.genotype == "homozygous"
    assert heterozygous.genotype == "heterozygous"


def test_extracts_multiple_time_points() -> None:
    metadata = extract_biological_metadata(
        "Samples collected at 3 days, 6 days, and 9 days",
    )

    assert metadata.time_point == "3 days; 6 days; 9 days"


def test_extracts_multiple_strains() -> None:
    metadata = extract_biological_metadata(
        "Comparison of w1118, Canton-S, and DGRP-551 flies",
    )

    assert metadata.strain == "w1118; Canton-S; DGRP-551"


def test_extracts_multiple_perturbations() -> None:
    metadata = extract_biological_metadata(
        "RNAi knockdown and CRISPR knockout experiments",
    )

    assert metadata.perturbation == (
        "RNAi; knockdown; knockout; CRISPR"
    )


def test_specific_gut_tissue_suppresses_generic_gut() -> None:
    metadata = extract_biological_metadata(
        "RNA-seq from adult midgut and hindgut",
    )

    assert metadata.tissue == "midgut; hindgut"


def test_wing_disc_suppresses_generic_imaginal_disc() -> None:
    metadata = extract_biological_metadata(
        "Drosophila wing imaginal disc samples",
    )

    assert metadata.tissue == "wing disc"


def test_multiple_developmental_stages_are_retained() -> None:
    metadata = extract_biological_metadata(
        "Embryonic, larval, pupal, and adult samples",
    )

    assert metadata.developmental_stage == (
        "embryo; larva; pupa; adult"
    )


def test_specific_glial_cell_type_suppresses_generic_glia() -> None:
    metadata = extract_biological_metadata(
        "astrocyte-like glia from adult brains",
    )

    assert metadata.cell_type == "astrocyte-like glia"
