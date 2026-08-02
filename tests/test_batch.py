"""Tests for multi-gene batch searches."""

from dataset_finder.batch import BatchSearchService
from dataset_finder.models import DatasetRecord


def make_record(accession: str, title: str) -> DatasetRecord:
    """Create a test dataset record."""
    return DatasetRecord(
        uid=accession,
        accession=accession,
        title=title,
        organism="Drosophila melanogaster",
        study_type="Expression profiling by high throughput sequencing",
        sample_count=4,
        publication_date="2026-01-01",
        url=f"https://example.org/{accession}",
    )


class FakeSearchService:
    """Simple deterministic search service."""

    def search(
        self,
        *,
        species: str,
        query: str,
        database: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, database, max_results
        return [
            make_record(
                f"GSE-{query}",
                f"{query} adult brain RNA-seq",
            )
        ]


def test_batch_search_assigns_gene_and_technique() -> None:
    service = BatchSearchService(
        search_service=FakeSearchService(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["bru1", "Hrb98DE"],
        database="geo",
        max_results_per_gene=10,
        gene_set="RBP",
    )

    assert len(result.records) == 2
    assert result.records[0].gene == "bru1"
    assert result.records[0].gene_set == "RBP"
    assert result.records[0].technique == "RNA_seq"
    assert result.records[0].database == "GEO"
    assert not result.issues


class PartiallyFailingSearchService:
    """Search service that fails for one gene."""

    def search(
        self,
        *,
        species: str,
        query: str,
        database: str,
        max_results: int,
    ) -> list[DatasetRecord]:
        del species, database, max_results

        if query == "bad":
            raise RuntimeError("Temporary database failure")

        return [
            make_record(
                "GSE1",
                f"{query} RNA-seq",
            )
        ]


def test_batch_search_records_gene_errors() -> None:
    service = BatchSearchService(
        search_service=PartiallyFailingSearchService(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["good", "bad"],
        database="geo",
        max_results_per_gene=10,
    )

    assert len(result.records) == 1
    assert len(result.issues) == 1
    assert result.issues[0].gene == "bad"


class FakeFlyBaseGene:
    """Minimal resolved FlyBase record used by tests."""

    submitted_symbol = "h"
    official_symbol = "hry"
    flybase_id = "FBgn0001168"
    current_fullname = "hairy"
    synonyms = ("hairy", "h")
    match_type = "synonym"
    ambiguous = True
    flybase_url = (
        "https://flybase.org/reports/"
        "FBgn0001168.html"
    )
    flyatlas_url = (
        "https://motif.mvls.gla.ac.uk/FlyAtlas2/"
        "index.html?search=gene&gene="
        "FBgn0001168&idtype=fbgn"
    )


class FakeFlyBaseResolver:
    """Return deterministic FlyBase metadata."""

    def resolve(self, symbol: str) -> FakeFlyBaseGene:
        del symbol
        return FakeFlyBaseGene()


def test_batch_search_adds_flybase_metadata() -> None:
    service = BatchSearchService(
        search_service=FakeSearchService(),
        flybase_resolver=FakeFlyBaseResolver(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=5,
        gene_set="TF",
    )

    record = result.records[0]

    assert record.gene == "h"
    assert record.official_symbol == "hry"
    assert record.flybase_id == "FBgn0001168"
    assert record.synonyms == ("hairy", "h")
    assert record.match_type == "FlyBase identifier"
    assert record.confidence == "Medium"
    assert record.flybase_url.endswith(
        "FBgn0001168.html"
    )
    assert "FBgn0001168" in record.flyatlas_url
    assert "idtype=fbgn" in record.flyatlas_url


def test_batch_search_query_includes_resolved_flybase_terms() -> None:
    class CapturingSearchService:
        def __init__(self) -> None:
            self.query = ""

        def search(
            self,
            *,
            species: str,
            query: str,
            database: str,
            max_results: int,
        ) -> list[DatasetRecord]:
            del species, database, max_results
            self.query = query
            return [
                    make_record(
                        "GSE1",
                        "hairy gene RNA-seq",
                    )
                ]

    search_service = CapturingSearchService()

    service = BatchSearchService(
        search_service=search_service,
        flybase_resolver=FakeFlyBaseResolver(),
    )

    service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=5,
    )

    assert '"h"' in search_service.query
    assert '"hry"' in search_service.query
    assert '"FBgn0001168"' in search_service.query


class FakeFlyAtlasClient:
    """Return deterministic FlyAtlas expression values."""

    def fetch(self, flybase_id: str):
        from dataset_finder.clients.flyatlas import FlyAtlasExpression

        return FlyAtlasExpression(
            flybase_id=flybase_id,
            symbol="hry",
            brain_male_fpkm=12.5,
            brain_female_fpkm=10.25,
            brain_larval_fpkm=7.5,
            head_male_fpkm=15.0,
            head_female_fpkm=14.0,
            top_male_tissue="Testis",
            top_male_fpkm=100.0,
            top_female_tissue="Ovary",
            top_female_fpkm=120.0,
            top_larval_tissue="Hindgut",
            top_larval_fpkm=40.0,
        )


def test_batch_search_adds_flyatlas_expression() -> None:
    service = BatchSearchService(
        search_service=FakeSearchService(),
        flybase_resolver=FakeFlyBaseResolver(),
        flyatlas_client=FakeFlyAtlasClient(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=5,
    )

    record = result.records[0]

    assert record.flyatlas_brain_male_fpkm == 12.5
    assert record.flyatlas_brain_female_fpkm == 10.25
    assert record.flyatlas_top_male_tissue == "Testis"
    assert record.flyatlas_top_female_tissue == "Ovary"
    assert record.flyatlas_top_larval_tissue == "Hindgut"


def test_batch_overfetches_candidates_before_relevance_filtering() -> None:
    class CandidateSearchService:
        def __init__(self) -> None:
            self.max_results = 0

        def search(
            self,
            *,
            species: str,
            query: str,
            database: str,
            max_results: int,
        ) -> list[DatasetRecord]:
            del species, query, database
            self.max_results = max_results

            return [
                make_record(
                    f"GSE{index}",
                    (
                        "unrelated transcriptome"
                        if index < 4
                        else "hairy gene RNA-seq"
                    ),
                )
                for index in range(5)
            ]

    search_service = CandidateSearchService()

    service = BatchSearchService(
        search_service=search_service,
        flybase_resolver=FakeFlyBaseResolver(),
        flyatlas_client=FakeFlyAtlasClient(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=1,
    )

    assert search_service.max_results == 50
    assert len(result.records) == 1
    assert result.records[0].gene == "h"


def test_batch_limits_final_accepted_results_per_gene() -> None:
    class ManyAcceptedSearchService:
        def search(
            self,
            *,
            species: str,
            query: str,
            database: str,
            max_results: int,
        ) -> list[DatasetRecord]:
            del species, query, database, max_results

            return [
                make_record(
                    f"GSE{index}",
                    "hairy gene RNA-seq",
                )
                for index in range(10)
            ]

    service = BatchSearchService(
        search_service=ManyAcceptedSearchService(),
        flybase_resolver=FakeFlyBaseResolver(),
        flyatlas_client=FakeFlyAtlasClient(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=3,
    )

    assert len(result.records) == 3


def test_batch_search_extracts_biological_metadata() -> None:
    class MetadataSearchService:
        def search(
            self,
            *,
            species: str,
            query: str,
            database: str,
            max_results: int,
        ) -> list[DatasetRecord]:
            del species, query, database, max_results

            return [
                DatasetRecord(
                    uid="GSE_META",
                    accession="GSE_META",
                    title=(
                        "RNA-seq of adult male Drosophila "
                        "brain neurons after hairy RNAi"
                    ),
                    organism="Drosophila melanogaster",
                    study_type=(
                        "Expression profiling by "
                        "high throughput sequencing"
                    ),
                    sample_count=6,
                    publication_date="2026/01/01",
                    url="https://example.org/GSE_META",
                    description=(
                        "w1118 control and treated samples "
                        "collected after 24 hours"
                    ),
                )
            ]

    service = BatchSearchService(
        search_service=MetadataSearchService(),
        flybase_resolver=FakeFlyBaseResolver(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=5,
    )

    assert len(result.records) == 1

    record = result.records[0]

    assert record.tissue == "brain"
    assert record.cell_type == "neuron"
    assert record.developmental_stage == "adult"
    assert record.sex == "male"
    assert record.strain == "w1118"
    assert record.control_status == "control"
    assert record.time_point == "24 hours"
    assert record.perturbation == "RNAi"


def test_existing_record_metadata_is_preserved() -> None:
    class MetadataSearchService:
        def search(
            self,
            *,
            species: str,
            query: str,
            database: str,
            max_results: int,
        ) -> list[DatasetRecord]:
            del species, query, database, max_results

            return [
                DatasetRecord(
                    uid="GSE_EXISTING",
                    accession="GSE_EXISTING",
                    title="hairy gene study in adult male brain",
                    organism="Drosophila melanogaster",
                    study_type="RNA-seq",
                    sample_count=4,
                    publication_date="2026/01/01",
                    url="https://example.org/GSE_EXISTING",
                    tissue="head",
                    sex="female",
                    strain="Canton-S",
                )
            ]

    service = BatchSearchService(
        search_service=MetadataSearchService(),
        flybase_resolver=FakeFlyBaseResolver(),
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=5,
    )

    record = result.records[0]

    assert record.tissue == "head"
    assert record.sex == "female"
    assert record.strain == "Canton-S"


class FakeGEOSampleMetadataClient:
    """Return deterministic GEO sample-level metadata."""

    def __init__(self) -> None:
        self.accessions: list[str] = []

    def fetch_sample_metadata(
        self,
        accession: str,
    ) -> dict[str, object]:
        self.accessions.append(accession)

        return {
            "sample_accessions": [
                "GSM8193871",
                "GSM8193872",
            ],
            "biosample_accessions": [
                "SAMN40876891",
                "SAMN40876892",
            ],
            "experiment_accessions": [
                "SRX24191199",
                "SRX24191200",
            ],
            "source_names": ["brain"],
            "characteristics": {
                "tissue": ["brain"],
                "genotype": [
                    "inscGAL4>wRNAi",
                    "inscGAL4>orb2RNAi",
                ],
            },
            "treatment_protocols": [
                "RNAi depletion of white or orb2",
            ],
            "growth_protocols": [
                "wandering L3 larvae",
            ],
            "library_strategies": ["RNA-Seq"],
            "library_sources": ["transcriptomic"],
            "library_selections": ["cDNA"],
            "instrument_models": [
                "Illumina NovaSeq 6000",
            ],
        }


def test_batch_enriches_accepted_geo_record_with_sample_metadata() -> None:
    class GEOSearchService:
        def search(
            self,
            *,
            species: str,
            query: str,
            database: str,
            max_results: int,
        ) -> list[DatasetRecord]:
            del species, query, database, max_results

            return [
                DatasetRecord(
                    uid="GSE263513",
                    accession="GSE263513",
                    title=(
                        "Effect of hairy RNAi depletion "
                        "in Drosophila larval brains"
                    ),
                    organism="Drosophila melanogaster",
                    study_type=(
                        "Expression profiling by "
                        "high throughput sequencing"
                    ),
                    sample_count=2,
                    publication_date="2024/04/09",
                    url="https://example.org/GSE263513",
                    database="GEO",
                )
            ]

    geo_client = FakeGEOSampleMetadataClient()

    service = BatchSearchService(
        search_service=GEOSearchService(),
        flybase_resolver=FakeFlyBaseResolver(),
        geo_client=geo_client,
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=5,
    )

    assert len(result.records) == 1
    assert geo_client.accessions == ["GSE263513"]

    record = result.records[0]

    assert record.tissue == "brain"
    assert record.developmental_stage == "larva"
    assert record.genotype == "transgenic"
    assert record.perturbation == "RNAi"
    assert record.sample_accessions == (
        "GSM8193871",
        "GSM8193872",
    )
    assert record.biosample_accessions == (
        "SAMN40876891",
        "SAMN40876892",
    )
    assert record.experiment_accessions == (
        "SRX24191199",
        "SRX24191200",
    )
    assert record.library_strategy == "RNA-Seq"
    assert record.library_source == "transcriptomic"
    assert record.library_selection == "cDNA"
    assert record.platform == "Illumina NovaSeq 6000"


def test_rejected_geo_record_does_not_fetch_sample_metadata() -> None:
    class IrrelevantGEOSearchService:
        def search(
            self,
            *,
            species: str,
            query: str,
            database: str,
            max_results: int,
        ) -> list[DatasetRecord]:
            del species, query, database, max_results

            return [
                DatasetRecord(
                    uid="GSE_OTHER",
                    accession="GSE_OTHER",
                    title="Unrelated adult fly transcriptome",
                    organism="Drosophila melanogaster",
                    study_type="RNA-seq",
                    sample_count=4,
                    publication_date="2026/01/01",
                    url="https://example.org/GSE_OTHER",
                    database="GEO",
                )
            ]

    geo_client = FakeGEOSampleMetadataClient()

    service = BatchSearchService(
        search_service=IrrelevantGEOSearchService(),
        flybase_resolver=FakeFlyBaseResolver(),
        geo_client=geo_client,
    )

    result = service.search_many(
        species="Drosophila melanogaster",
        genes=["h"],
        database="geo",
        max_results_per_gene=5,
    )

    assert result.records == ()
    assert geo_client.accessions == []
