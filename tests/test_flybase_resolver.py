"""Tests for packaged FlyBase gene resolution."""

from dataset_finder.flybase_resolver import FlyBaseResolver


def test_resolve_known_rbp() -> None:
    resolver = FlyBaseResolver()
    gene = resolver.resolve("bru1")

    assert gene.official_symbol
    assert gene.flybase_id.startswith("FBgn")
    assert gene.flybase_url.startswith(
        "https://flybase.org/reports/"
    )


def test_resolve_known_tf() -> None:
    resolver = FlyBaseResolver()
    gene = resolver.resolve("abd-A")

    assert gene.official_symbol == "abd-A"
    assert gene.flybase_id == "FBgn0000014"
    assert gene.match_type == "official_symbol"


def test_resolution_is_case_insensitive() -> None:
    resolver = FlyBaseResolver()

    assert (
        resolver.resolve("HRB98DE").flybase_id
        == resolver.resolve("Hrb98DE").flybase_id
    )


def test_unresolved_symbol_returns_safe_record() -> None:
    resolver = FlyBaseResolver()
    gene = resolver.resolve("not_a_real_gene")

    assert gene.submitted_symbol == "not_a_real_gene"
    assert gene.official_symbol == ""
    assert gene.flybase_id == ""
    assert gene.match_type == "unresolved"


def test_search_terms_include_flybase_id() -> None:
    resolver = FlyBaseResolver()
    gene = resolver.resolve("caz")

    assert "caz" in gene.search_terms
    assert gene.flybase_id in gene.search_terms


def test_exact_official_symbols_win_over_ambiguous_synonyms() -> None:
    resolver = FlyBaseResolver()

    for symbol in ("D", "ss", "z"):
        gene = resolver.resolve(symbol)

        assert gene.official_symbol == symbol
        assert gene.match_type == "official_symbol"


def test_historical_h_symbol_resolves_to_hairy() -> None:
    resolver = FlyBaseResolver()
    gene = resolver.resolve("h")

    assert gene.official_symbol == "hry"
    assert gene.flybase_id == "FBgn0001168"
    assert gene.current_fullname == "hairy"
    assert gene.match_type == "synonym"


def test_flyatlas_urls_use_flybase_identifier() -> None:
    resolver = FlyBaseResolver()
    gene = resolver.resolve("bru1")

    assert gene.flybase_id == "FBgn0000114"
    assert "FBgn0000114" in gene.flyatlas_url
    assert "idtype=fbgn" in gene.flyatlas_url
    assert "FBgn0000114" in gene.flyatlas_download_url
    assert "tableOut=gene" in gene.flyatlas_download_url
