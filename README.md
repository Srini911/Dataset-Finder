# Dataset Finder

> Discover, organize, and curate public functional-genomics datasets across multiple biological databases.

Dataset Finder is an open-source Python bioinformatics toolkit for discovering publicly available functional-genomics datasets from major biological repositories.

The project provides a unified command-line interface for searching and normalizing dataset metadata across species, genes, tissues, diseases, regulators, and experimental assays.

## Current Capabilities

- Search GEO datasets
- Search ENCODE experiments
- Search GEO and ENCODE through a common command-line interface
- Normalize results into a shared dataset-record model
- Validate code with pytest and Ruff

## Project Status

Dataset Finder is under active development.

### Implemented

- GEO dataset search
- ENCODE experiment search
- Combined database search
- Python command-line interface
- Normalized metadata model
- Automated test suite

### In Progress

- Documentation improvements
- Export utilities
- Ranking and relevance scoring
- Expanded metadata normalization

### Planned

- SRA integration
- BioProject integration
- BioSample integration
- PubMed integration
- ENA integration
- Expression Atlas integration

## Supported Databases

| Database | Status |
|---|---|
| GEO | Implemented |
| ENCODE | Implemented |
| SRA | Planned |
| BioProject | Planned |
| BioSample | Planned |
| PubMed | Planned |
| ENA | Future |
| Expression Atlas | Future |
| ProteomeXchange | Future |

## Supported Assay Types

Dataset Finder is intended to support functional-genomics studies including:

- bulk RNA-seq
- single-cell RNA-seq
- single-nucleus RNA-seq
- ATAC-seq
- ChIP-seq
- CUT&RUN
- CUT&Tag
- eCLIP and other CLIP-based assays
- spatial transcriptomics
- perturbation-based transcriptomics

## Guiding Principle

**Search once. Discover everywhere.**

## Installation

Dataset Finder requires Python 3.11 or newer.

```bash
git clone https://github.com/Srini911/Dataset-Finder.git
cd Dataset-Finder

python3 -m venv .venv
source .venv/bin/activate

pip install -e ".[dev]"
```

## Command-Line Interface

Check the installed version:

```bash
dataset-finder --version
```

Display general help:

```bash
dataset-finder --help
```

Display search help:

```bash
dataset-finder search --help
```

## Quick Start

### Search GEO

```bash
dataset-finder search \
  --database geo \
  --species "Mus musculus" \
  --query "Tardbp" \
  --max-results 5
```

### Search ENCODE

```bash
dataset-finder search \
  --database encode \
  --species "Homo sapiens" \
  --query "TARDBP" \
  --max-results 5
```

### Search All Implemented Databases

```bash
dataset-finder search \
  --database all \
  --species "Homo sapiens" \
  --query "TARDBP" \
  --max-results 10
```

The current combined search queries GEO and ENCODE and returns normalized dataset records.

## Search Result Metadata

Search results may include:

- accession identifier
- dataset title
- organism
- assay or study type
- sample or replicate count
- release or publication date when available
- direct database URL

## Planned Outputs

Future releases are intended to provide:

- dataset inventories
- accession identifiers
- direct database links
- sample and assay metadata
- publication references
- relevance scores
- Excel workbooks
- CSV tables
- JSON records
- reproducible search reports
- methods-ready summaries

## Architecture

```text
Command-Line Interface
        |
        v
SearchService
        |
        +-- GEO Client
        |
        +-- ENCODE Client
        |
        +-- Future Database Clients
```

All database clients normalize results into a common dataset-record model.

## Repository Structure

```text
src/dataset_finder/
├── clients/       Database clients
├── exporters/     Export functionality
├── utils/         Shared utilities
├── cli.py         Command-line interface
├── models.py      Shared data models
├── ranking.py     Ranking utilities
└── search.py      Search orchestration

tests/             Automated tests
```

## Development

Run linting:

```bash
ruff check src tests
```

Run tests:

```bash
pytest
```

Run all checks:

```bash
ruff check src tests && pytest
```

## Roadmap

### Version 0.2

- GEO search
- ENCODE search
- Combined database search
- ENCODE client tests
- Improved documentation

### Version 0.3

- SRA integration
- CSV export
- JSON export
- Excel export
- Relevance ranking

### Version 0.4

- BioProject integration
- PubMed integration
- ENA integration
- Advanced assay filters
- Tissue and disease filters

### Version 1.0

- Unified multi-database search
- Harmonized metadata
- Publication-ready reports
- Extensible client architecture
- Stable public release

## Citation and Acknowledgment

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
