"""Command-line interface for Dataset Finder."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dataset_finder import __version__
from dataset_finder.batch import BatchSearchService
from dataset_finder.exporters import export_csv, export_excel, export_json
from dataset_finder.gene_sets import GeneInputError, collect_genes
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
            "Discover public RBP and TF functional-genomics datasets "
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
        description="Search one or many genes across supported databases.",
    )
    search_parser.add_argument(
        "--species",
        required=True,
        help='Scientific species name, for example "Drosophila melanogaster".',
    )
    search_parser.add_argument(
        "--query",
        help="One gene or legacy free-text query.",
    )
    search_parser.add_argument(
        "--genes",
        nargs="+",
        help="Multiple genes separated by spaces or commas.",
    )
    search_parser.add_argument(
        "--gene-file",
        type=Path,
        help="Plain-text file containing one or more gene symbols.",
    )
    search_parser.add_argument(
        "--gene-set-name",
        default="",
        help="Label stored in exports, for example RBP or TF.",
    )
    search_parser.add_argument(
        "--database",
        choices=("geo", "encode", "all"),
        default="geo",
        help="Database to search. Default: geo.",
    )
    search_parser.add_argument(
        "--max-results",
        type=positive_integer,
        default=20,
        help="Maximum results per gene. Default: 20.",
    )
    search_parser.add_argument(
        "--format",
        choices=("table", "csv", "json", "xlsx"),
        default="table",
        help="Output format. Default: table.",
    )
    search_parser.add_argument(
        "--output",
        type=Path,
        help="Output file path.",
    )

    return parser


def print_records(records: list[DatasetRecord]) -> None:
    """Print normalized dataset records."""
    for index, record in enumerate(records, start=1):
        print(f"{index}. {record.accession or record.uid}")

        if record.gene:
            print(f"   Gene: {record.gene}")

        print(f"   Title: {record.title or 'Not available'}")
        print(f"   Organism: {record.organism or 'Not available'}")
        print(f"   Study type: {record.study_type or 'Not available'}")

        if record.technique:
            print(f"   Technique: {record.technique}")

        print(
            "   Samples: "
            f"{record.sample_count if record.sample_count is not None else 'Not available'}"
        )
        print(
            "   Publication date: "
            f"{record.publication_date or 'Not available'}"
        )
        print(f"   URL: {record.url}")


def default_output_path(output_format: str) -> Path:
    """Return the default output filename for an export format."""
    extension = "xlsx" if output_format == "xlsx" else output_format
    return Path(f"dataset_finder_results.{extension}")


def export_records(
    records: list[DatasetRecord],
    *,
    output_format: str,
    output_path: Path | None,
) -> Path:
    """Export simple record lists to CSV or JSON."""
    path = output_path or default_output_path(output_format)

    if output_format == "csv":
        return export_csv(records, path)

    if output_format == "json":
        return export_json(records, path)

    raise ValueError(f"Unsupported export format: {output_format}")


def run_search(args: argparse.Namespace) -> int:
    """Run a single- or multi-gene dataset search."""
    species = args.species.strip()
    database = args.database.strip().lower()

    try:
        genes = collect_genes(
            query=args.query,
            genes=args.genes,
            gene_file=args.gene_file,
        )
    except GeneInputError as exc:
        print(f"Dataset Finder error: {exc}", file=sys.stderr)
        return 2

    if len(genes) == 1 and args.format != "xlsx":
        service = SearchService()

        try:
            records = service.search(
                species=species,
                query=genes[0],
                database=database,
                max_results=args.max_results,
            )
        except UnsupportedDatabaseError as exc:
            print(str(exc), file=sys.stderr)
            return 2
        except ValueError as exc:
            print(f"Dataset Finder error: {exc}", file=sys.stderr)
            return 2
        except Exception as exc:
            print(f"Dataset Finder error: {exc}", file=sys.stderr)
            return 1

        print(f"Dataset Finder {database.upper()} search")
        print(f"Species: {species}")
        print(f"Query: {genes[0]}")
        print(f"Results found: {len(records)}")

        if args.format in {"csv", "json"}:
            try:
                exported_path = export_records(
                    records,
                    output_format=args.format,
                    output_path=args.output,
                )
            except OSError as exc:
                print(
                    f"Dataset Finder export error: {exc}",
                    file=sys.stderr,
                )
                return 1

            print(f"Exported results: {exported_path}")
            return 0

        if not records:
            if database == "geo":
                print("No matching GEO Series records were found.")
            else:
                print(f"No matching {database.upper()} records were found.")
            return 0

        print()
        print_records(records)
        return 0

    batch_service = BatchSearchService()

    try:
        result = batch_service.search_many(
            species=species,
            genes=genes,
            database=database,
            max_results_per_gene=args.max_results,
            gene_set=args.gene_set_name.strip(),
        )
    except ValueError as exc:
        print(f"Dataset Finder error: {exc}", file=sys.stderr)
        return 2

    print("Dataset Finder multi-gene search")
    print(f"Species: {species}")
    print(f"Genes searched: {len(genes)}")
    print(f"Results found: {len(result.records)}")
    print(f"Search errors: {len(result.issues)}")

    if args.format == "xlsx":
        output_path = args.output or default_output_path("xlsx")

        try:
            exported_path = export_excel(
                result,
                output_path,
            )
        except OSError as exc:
            print(
                f"Dataset Finder export error: {exc}",
                file=sys.stderr,
            )
            return 1

        print(f"Exported workbook: {exported_path}")
        return 0

    records = list(result.records)

    if args.format in {"csv", "json"}:
        exported_path = export_records(
            records,
            output_format=args.format,
            output_path=args.output,
        )
        print(f"Exported results: {exported_path}")
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
