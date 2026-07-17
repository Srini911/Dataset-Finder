"""Tests for the Dataset Finder command-line interface."""

from __future__ import annotations

import pytest

from dataset_finder.cli import build_parser, main
from dataset_finder.clients.ncbi_geo import GEORecord


def test_parser_program_name() -> None:
    """The parser should expose the correct program name."""
    parser = build_parser()

    assert parser.prog == "dataset-finder"


def test_search_arguments_are_parsed() -> None:
    """The search command should parse its required arguments."""
    parser = build_parser()

    args = parser.parse_args(
        [
            "search",
            "--species",
            "Drosophila melanogaster",
            "--query",
            "brain RNA-seq",
        ]
    )

    assert args.command == "search"
    assert args.species == "Drosophila melanogaster"
    assert args.query == "brain RNA-seq"
    assert args.database == "geo"
    assert args.max_results == 20


def test_search_rejects_nonpositive_max_results() -> None:
    """The search command should reject invalid result limits."""
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "search",
                "--species",
                "Drosophila melanogaster",
                "--query",
                "brain RNA-seq",
                "--max-results",
                "0",
            ]
        )


def test_geo_search_command_output(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The GEO search command should print normalized records."""
    fake_records = [
        GEORecord(
            uid="200123456",
            accession="GSE123456",
            title="Drosophila brain RNA sequencing",
            organism="Drosophila melanogaster",
            study_type="Expression profiling by high throughput sequencing",
            sample_count=12,
            publication_date="2026/01/10",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123456",
        )
    ]

    monkeypatch.setattr(
        "dataset_finder.search.NCBIGEOClient.search",
        lambda self, **kwargs: fake_records,
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "dataset-finder",
            "search",
            "--species",
            "Drosophila melanogaster",
            "--query",
            "brain RNA-seq",
            "--database",
            "geo",
            "--max-results",
            "5",
        ],
    )

    exit_code = main()
    output = capsys.readouterr().out

    assert exit_code == 0
    assert "Dataset Finder GEO search" in output
    assert "Results found: 1" in output
    assert "GSE123456" in output
    assert "Drosophila brain RNA sequencing" in output
    assert "Samples: 12" in output


def test_geo_search_handles_no_results(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The CLI should report an empty GEO search clearly."""
    monkeypatch.setattr(
        "dataset_finder.search.NCBIGEOClient.search",
        lambda self, **kwargs: [],
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "dataset-finder",
            "search",
            "--species",
            "Drosophila melanogaster",
            "--query",
            "unlikely query",
        ],
    )

    exit_code = main()
    output = capsys.readouterr().out

    assert exit_code == 0
    assert "Results found: 0" in output
    assert "No matching GEO Series records were found." in output


def test_sra_search_is_not_implemented(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The CLI should reject SRA searches until support is added."""
    monkeypatch.setattr(
        "sys.argv",
        [
            "dataset-finder",
            "search",
            "--species",
            "Drosophila melanogaster",
            "--query",
            "brain RNA-seq",
            "--database",
            "sra",
        ],
    )

    exit_code = main()
    captured = capsys.readouterr()

    assert exit_code == 2
    assert "SRA search is not implemented yet." in captured.err
