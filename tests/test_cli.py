"""Tests for the Dataset Finder command-line interface."""

from dataset_finder.cli import build_parser


def test_parser_program_name() -> None:
    """The command-line parser should use the expected program name."""
    parser = build_parser()
    assert parser.prog == "dataset-finder"
