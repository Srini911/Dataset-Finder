"""Command-line interface for Dataset Finder."""

from __future__ import annotations

import argparse
import sys

from dataset_finder import __version__
from dataset_finder.clients.ncbi_geo import NCBIClientError
from dataset_finder.models import DatasetRecord
from dataset_finder.search import SearchService, UnsupportedDatabaseError


def positive_integer(value: str) -> int:
    """Parse a positive integer argument."""
    try:
        parsed_value = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "value must be an integer"
        ) from exc

    if parsed_value < 1:
        raise argparse.ArgumentTypeError(
            "value must be greater than zero"
        )

    return parsed_value


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
        description="Search supported public biological databases.",
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
        default="geo",
        help="Database to search. Default: geo.",
    )
    search_parser.add_argument(
        "--max-results",
        type=positive_integer,
        default=20,
        help="Maximum number of results. Default: 20.",
    )

    return parser


def print_records(records: list[DatasetRecord]) -> None:
    """Print normalized dataset records."""
    for index, record in enumerate(records, start=1):
        print(f"{index}. {record.accession or record.uid}")
        print(f"   Title: {record.title or 'Not available'}")
        print(f"   Organism: {record.organism or 'Not available'}")
        print(f"   Study type: {record.study_type or 'Not available'}")
        print(
            "   Samples: "
            f"{record.sample_count if record.sample_count is not None else 'Not available'}"
        )
        print(
            "   Publication date: "
            f"{record.publication_date or 'Not available'}"
        )
        print(f"   URL: {record.url}")


def run_search(args: argparse.Namespace) -> int:
    """Run a dataset search."""
    species = args.species.strip()
    query = args.query.strip()

    if args.database == "all":
        print(
            "Combined database search is not implemented yet; searching GEO only.",
            file=sys.stderr,
        )

    service = SearchService()

    try:
        records = service.search(
            species=species,
            query=query,
            database=args.database,
            max_results=args.max_results,
        )
    except UnsupportedDatabaseError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    except NCBIClientError as exc:
        print(f"Dataset Finder error: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"Dataset Finder error: {exc}", file=sys.stderr)
        return 2

    print("Dataset Finder GEO search")
    print(f"Species: {species}")
    print(f"Query: {query}")
    print(f"Results found: {len(records)}")

    if not records:
        print("No matching GEO Series records were found.")
        return 0

    print()
    print_records(records)
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
