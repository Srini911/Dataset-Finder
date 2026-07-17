# Dataset Finder

> Discover, organize, and curate public functional-genomics datasets across species.

Dataset Finder is an open-source bioinformatics toolkit for discovering and organizing publicly available datasets from major biological databases.

It is designed to work across species, tissues, diseases, genes, and experimental assays without hard-coding the workflow to a single organism.

## Project status

Dataset Finder is currently under active development.

The first releases will focus on:

- reproducible command-line searches
- structured metadata collection
- assay classification
- quality scoring
- Excel, CSV, and JSON outputs
- methods-ready reports

## Planned database support

Initial support is planned for:

- GEO
- SRA
- BioProject
- BioSample
- PubMed

Future integrations may include:

- ENA
- BioStudies
- Expression Atlas
- ProteomeXchange

## Planned assay support

Dataset Finder will classify experimental techniques from public metadata, including:

- bulk RNA-seq
- single-cell RNA-seq
- single-nucleus RNA-seq
- ATAC-seq
- ChIP-seq
- CUT&RUN
- CUT&Tag
- CLIP-based assays
- spatial transcriptomics
- other functional-genomics assays

## Guiding principle

**Search once. Discover everywhere.**

The software is intended to remain:

- species independent
- tissue independent
- assay independent
- reproducible
- publication ready
- community driven

## Installation

Dataset Finder currently requires Python 3.11 or newer.

```bash
git clone https://github.com/Srini911/Dataset-Finder.git
cd Dataset-Finder

python3 -m venv .venv
source .venv/bin/activate

pip install -e ".[dev]"
```

## Current command-line interface

Check the installed version:

```bash
dataset-finder --version
```

Display command help:

```bash
dataset-finder --help
```

## Planned search usage

The following examples show the intended future interface and are not implemented yet.

```bash
dataset-finder search \
  --species "Mus musculus" \
  --tissue brain \
  --gene Tardbp \
  --assays RNA-seq scRNA-seq
```

```bash
dataset-finder search \
  --species "Homo sapiens" \
  --tissue hippocampus \
  --disease "Alzheimer disease" \
  --assays any
```

## Planned outputs

Each completed search is intended to produce:

- dataset inventory
- accession identifiers
- direct database links
- sample and assay metadata
- publication references
- evidence and confidence classifications
- Excel workbook
- CSV tables
- JSON records
- reproducible search report
- methods-ready summary

## Repository structure

```text
config/       Species, assay, tissue, and database configuration
docs/         User and developer documentation
examples/     Reproducible example searches
outputs/      Locally generated search results
scripts/      Development and maintenance utilities
src/          Dataset Finder Python package
tests/        Automated tests
```

## Development checks

Run linting:

```bash
ruff check src tests
```

Run tests:

```bash
pytest
```

## Citation and acknowledgment

If Dataset Finder contributes to your research, publication, thesis, report, software, or teaching, please cite the software or acknowledge its creator.

**Creator and lead developer:** Srinivas Amla

Suggested acknowledgment:

> Public functional-genomics datasets were identified and organized using Dataset Finder, developed by Srinivas Amla.

Formal citation metadata is provided in `CITATION.cff`.

## Contributing

Community contributions are welcome. Please read `CONTRIBUTING.md` before submitting a pull request.

## License

Dataset Finder is distributed under the MIT License. See `LICENSE`.

## Author

**Srinivas Amla**  
University of Massachusetts Boston

Research interests include bioinformatics, computational biology, functional genomics, machine learning, and public biological-data integration.
