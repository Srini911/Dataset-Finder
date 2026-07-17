"""Tests for the Dataset Finder command-line interface."""

from __future__ import annotations

import pytest

from dataset_finder.cli import build_parser, main


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
    assert args.database == "all"
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


def test_search_command_output(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The search command should print a validated request summary."""
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
            "10",
        ],
    )

    exit_code = main()
    output = capsys.readouterr().out

    assert exit_code == 0
    assert "Dataset Finder search" in output
    assert "Species: Drosophila melanogaster" in output
    assert "Query: brain RNA-seq" in output
    assert "Database: geo" in output
    assert "Maximum results: 10" in output
