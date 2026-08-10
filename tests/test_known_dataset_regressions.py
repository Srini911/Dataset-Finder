"""Known historically validated datasets used as regression targets."""

KNOWN_DATASETS = {
    ("CG7804", "RNA_seq"): {
        "SRP110269",
    },
    ("CG7804", "ChIP_seq"): {
        "SRP110266",
    },
    ("Hrb98DE", "RNA_seq"): {
        "SRP186005",
        "SRP001537",
    },
    ("bru1", "RNA_seq"): {
        "SRP377648",
        "SRP050336",
    },
    ("snf", "RNA_seq"): {
        "SRP055034",
    },
    ("snf", "ChIP_seq"): {
        "SRP131779",
    },
}


def test_known_dataset_registry_has_unique_targets() -> None:
    for accessions in KNOWN_DATASETS.values():
        assert accessions
        assert len(accessions) == len(set(accessions))


def test_known_dataset_registry_uses_study_level_accessions() -> None:
    for accessions in KNOWN_DATASETS.values():
        for accession in accessions:
            assert accession.startswith(
                (
                    "SRP",
                    "ERP",
                    "DRP",
                    "GSE",
                )
            )
