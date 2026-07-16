"""Command-line interface for Dataset Finder."""

from __future__ import annotations

import argparse

from dataset_finder import __version__


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="dataset-finder",
        description=(
            "Discover and organize public functional-genomics datasets "
            "across species."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def main() -> int:
    """Run the Dataset Finder command-line interface."""
    parser = build_parser()
    parser.parse_args()
    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
