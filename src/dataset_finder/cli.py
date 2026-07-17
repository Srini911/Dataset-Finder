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

    subparsers = parser.add_subparsers(dest="command")

    search_parser = subparsers.add_parser(
        "search",
        help="Search for public functional-genomics datasets.",
        description="Prepare a structured dataset search request.",
    )
    search_parser.add_argument(
        "--species",
        required=True,
        help='Scientific species name, for example "Drosophila melanogaster".',
    )
    search_parser.add_argument(
        "--query",
        required=True,
        help='Search terms, for example "brain RNA-seq".',
    )
    search_parser.add_argument(
        "--database",
        choices=("geo", "sra", "all"),
        default="all",
        help="Database to search. Default: all.",
    )
    search_parser.add_argument(
        "--max-results",
        type=positive_integer,
        default=20,
        help="Maximum number of results. Default: 20.",
    )

    return parser


def positive_integer(value: str) -> int:
    """Parse a positive integer argument."""
    parsed_value = int(value)

    if parsed_value < 1:
        raise argparse.ArgumentTypeError("value must be greater than zero")

    return parsed_value


def run_search(args: argparse.Namespace) -> int:
    """Run a dataset search request."""
    species = args.species.strip()
    query = args.query.strip()

    if not species:
        raise ValueError("Species cannot be empty.")

    if not query:
        raise ValueError("Query cannot be empty.")

    print("Dataset Finder search")
    print(f"Species: {species}")
    print(f"Query: {query}")
    print(f"Database: {args.database}")
    print(f"Maximum results: {args.max_results}")
    print("Status: search request validated; database integration is next.")

    return 0


def main() -> int:
    """Run the Dataset Finder command-line interface."""
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "search":
        return run_search(args)

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
