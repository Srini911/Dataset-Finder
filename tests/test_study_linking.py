"""Tests for cross-database study linking."""

from dataset_finder.models import DatasetRecord
from dataset_finder.study_linking import (
    extract_links_from_record,
    extract_study_links,
)


def test_extracts_archive_and_publication_links() -> None:
    links = extract_study_links(
        """
        GEO GSE263513, BioProject PRJNA1097674,
        BioSample SAMN40876891, SRP590197,
        SRX24191199, SRR28900001,
        PMID: 38902233 and doi:10.1016/j.test.2026.01.001
        """
    )

    assert links.pubmed_ids == ("38902233",)
    assert links.dois == (
        "10.1016/j.test.2026.01.001",
    )
    assert links.related_geo_accessions == (
        "GSE263513",
    )
    assert links.related_study_accessions == (
        "SRP590197",
    )
    assert links.related_experiment_accessions == (
        "SRX24191199",
    )
    assert links.related_run_accessions == (
        "SRR28900001",
    )
    assert links.related_bioproject_accessions == (
        "PRJNA1097674",
    )
    assert links.related_biosample_accessions == (
        "SAMN40876891",
    )


def test_extracts_links_from_dataset_record() -> None:
    record = DatasetRecord(
        uid="1",
        accession="GSE263513",
        title="Study linked to PRJNA1097674",
        organism="Drosophila melanogaster",
        study_type="RNA-seq",
        sample_count=2,
        publication_date="2024",
        url="https://example.org",
        raw_metadata={
            "relations": [
                "SAMN40876891",
                "SRX24191199",
            ]
        },
    )

    links = extract_links_from_record(record)

    assert "GSE263513" in links.all_accessions
    assert "PRJNA1097674" in links.all_accessions
    assert "SAMN40876891" in links.all_accessions
    assert "SRX24191199" in links.all_accessions


def test_links_records_with_matching_titles_across_databases() -> None:
    geo = DatasetRecord(
        uid="1",
        accession="GSE299109",
        title="Orb2 regulation during development",
        organism="Drosophila melanogaster",
        study_type="RNA-seq",
        sample_count=4,
        publication_date="2025",
        url="https://example.org/GSE299109",
        database="GEO",
    )

    publication = DatasetRecord(
        uid="41206465",
        accession="PMID41206465",
        title="Orb2 regulation during development",
        organism="",
        study_type="Publication",
        sample_count=None,
        publication_date="2026",
        url="https://pubmed.ncbi.nlm.nih.gov/41206465/",
        database="PubMed",
        pubmed_ids=("41206465",),
        dois=("10.1000/orb2",),
    )

    from dataset_finder.study_linking import link_related_records

    linked = link_related_records(
        [geo, publication]
    )

    for record in linked:
        assert record.pubmed_ids == ("41206465",)
        assert record.dois == ("10.1000/orb2",)
        assert "GSE299109" in record.related_accessions
        assert "PMID41206465" in record.related_accessions
        assert record.related_geo_accessions == (
            "GSE299109",
        )


def test_does_not_link_same_title_from_only_one_database() -> None:
    first = DatasetRecord(
        uid="1",
        accession="GSE1",
        title="Repeated study title",
        organism="Drosophila melanogaster",
        study_type="RNA-seq",
        sample_count=2,
        publication_date="2025",
        url="https://example.org/GSE1",
        database="GEO",
    )

    second = DatasetRecord(
        uid="2",
        accession="GSE2",
        title="Repeated study title",
        organism="Drosophila melanogaster",
        study_type="RNA-seq",
        sample_count=2,
        publication_date="2025",
        url="https://example.org/GSE2",
        database="GEO",
    )

    from dataset_finder.study_linking import link_related_records

    linked = link_related_records([first, second])

    assert linked[0].related_accessions == ()
    assert linked[1].related_accessions == ()


def test_study_level_accessions_exclude_samples_and_experiments() -> None:
    links = extract_study_links(
        """
        GSE263513 GSM8193871 SRP590197 SRX24191199
        SRR28900001 PRJNA1097674 SAMN40876891
        """
    )

    assert links.study_level_accessions == (
        "GSE263513",
        "SRP590197",
        "PRJNA1097674",
    )


def test_cross_record_links_remain_compact() -> None:
    geo = DatasetRecord(
        uid="1",
        accession="GSE299109",
        title="A shared Orb2 study title",
        organism="Drosophila melanogaster",
        study_type="RNA-seq",
        sample_count=2,
        publication_date="2025",
        url="https://example.org/GSE299109",
        database="GEO",
        sample_accessions=("GSM1",),
        experiment_accessions=("SRX1",),
        biosample_accessions=("SAMN1",),
    )

    publication = DatasetRecord(
        uid="41206465",
        accession="PMID41206465",
        title="A shared Orb2 study title",
        organism="",
        study_type="Publication",
        sample_count=None,
        publication_date="2026",
        url="https://pubmed.ncbi.nlm.nih.gov/41206465/",
        database="PubMed",
        pubmed_ids=("41206465",),
    )

    from dataset_finder.study_linking import link_related_records

    linked = link_related_records([geo, publication])

    for record in linked:
        assert record.related_accessions == (
            "GSE299109",
            "PMID41206465",
        )
