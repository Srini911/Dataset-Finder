cat > README.md <<'EOF'
# Dataset Finder

> Search and export public functional-genomics datasets from GEO and ENCODE through a unified command-line interface.

Dataset Finder is an open-source Python bioinformatics toolkit for discovering, organizing, and exporting publicly available functional-genomics datasets from major biological repositories.

It provides a consistent command-line interface for searching supported databases and normalizes returned metadata into a common dataset-record format for downstream analysis.

## Features

- Search GEO datasets
- Search ENCODE experiments
- Search GEO and ENCODE through one interface
- Discover datasets across multiple species
- Normalize metadata into a shared record structure
- Export search results as CSV
- Export search results as JSON
- Display results directly in the terminal
- Automated testing with pytest
- Code-quality validation with Ruff
- Continuous Integration with GitHub Actions

## Project Status

Dataset Finder is under active development.

### Implemented

- GEO dataset search
- ENCODE experiment search
- Combined GEO and ENCODE search
- Python command-line interface
- Normalized dataset metadata model
- Terminal table output
- CSV export
- JSON export
- Automated test suite
- GitHub Actions continuous integration

### In Progress

- Ranking and relevance scoring
- Expanded metadata normalization
- Enhanced biological search filters
- Documentation improvements
- Additional database integrations

### Planned

- SRA integration
- BioProject integration
- BioSample integration
- PubMed integration
- ENA integration
- Expression Atlas integration
- Excel export
- Reproducible search reports

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

## Target Assay Coverage

Dataset Finder is designed to support the discovery of functional-genomics studies involving:

- Bulk RNA sequencing
- Single-cell RNA sequencing
- Single-nucleus RNA sequencing
- ATAC-seq
- ChIP-seq
- CUT&RUN
- CUT&Tag
- eCLIP and other CLIP-based assays
- Spatial transcriptomics
- Perturbation-based transcriptomics

## Guiding Principle

**Search once. Discover everywhere.**

## Requirements

Dataset Finder requires:

- Python 3.11 or newer
- Internet access for querying public databases
- Git for cloning the repository

## Installation

Clone the repository:

```bash
git clone https://github.com/Srini911/Dataset-Finder.git
cd Dataset-Finder
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install Dataset Finder and its development dependencies:

```bash
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

Confirm that the command is available:

```bash
dataset-finder --version
```

## Command-Line Interface

Display general help:

```bash
dataset-finder --help
```

Display search-specific help:

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
  --query "CTCF" \
  --max-results 10
```

The combined search queries GEO and ENCODE and returns normalized dataset records through a common interface.

## Export Search Results

Dataset Finder supports terminal, CSV, and JSON output formats.

### Export Results as CSV

```bash
dataset-finder search \
  --species "Homo sapiens" \
  --query "CTCF" \
  --format csv \
  --output results.csv
```

Example terminal output:

```text
Dataset Finder GEO search
Species: Homo sapiens
Query: CTCF
Results found: 20
Exported results: /path/to/results.csv
```

### Export Results as JSON

```bash
dataset-finder search \
  --species "Drosophila melanogaster" \
  --query "Fru" \
  --format json \
  --output drosophila_fru_results.json
```

When `--output` is omitted, Dataset Finder creates a default file in the current working directory.

### Display Results in the Terminal

```bash
dataset-finder search \
  --database geo \
  --species "Drosophila melanogaster" \
  --query "fruitless" \
  --format table \
  --max-results 5
```

## Available Output Formats

| Format | CLI value | Status |
|---|---|---|
| Terminal table | `table` | Implemented |
| CSV | `csv` | Implemented |
| JSON | `json` | Implemented |
| Excel workbook | Not yet available | Planned |
| Reproducible report | Not yet available | Planned |

## Search Result Metadata

Normalized search results may include:

| Field | Description |
|---|---|
| `uid` | Database-specific unique identifier |
| `accession` | Public dataset or experiment accession |
| `title` | Dataset title |
| `organism` | Study organism or organisms |
| `study_type` | Assay or study classification |
| `sample_count` | Number of samples or replicates when available |
| `publication_date` | Release or publication date when available |
| `url` | Direct link to the source database record |

Metadata availability may differ between databases because GEO and ENCODE expose different record structures.

## Example CSV Structure

```text
uid,accession,title,organism,study_type,sample_count,publication_date,url
200304737,GSE304737,Example dataset,Drosophila melanogaster,Genome binding/occupancy profiling by high throughput sequencing,7,2026/02/09,https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304737
```

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
        |
        v
Normalized Dataset Records
        |
        +-- Terminal Output
        +-- CSV Export
        +-- JSON Export
```

Each supported database client converts source-specific metadata into the common `DatasetRecord` model.

## Repository Structure

```text
Dataset-Finder/
├── .github/
│   └── workflows/          GitHub Actions workflows
├── src/
│   └── dataset_finder/
│       ├── clients/        Database-specific clients
│       ├── exporters/      CSV and JSON exporters
│       ├── utils/          Shared utilities
│       ├── cli.py          Command-line interface
│       ├── models.py       Shared dataset models
│       ├── ranking.py      Ranking utilities
│       └── search.py       Search orchestration
├── tests/                  Automated tests
├── CHANGELOG.md
├── CITATION.cff
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
├── LICENSE
├── README.md
└── SECURITY.md
```

## Development

Install the project with development dependencies:

```bash
pip install -e ".[dev]"
```

Run Ruff:

```bash
ruff check src tests
```

Run the test suite:

```bash
pytest
```

Run all local validation checks:

```bash
ruff check src tests
pytest
git diff --check
```

## Roadmap

### Completed Foundation

- GEO dataset search
- ENCODE experiment search
- Combined database search
- Unified command-line interface
- Normalized metadata model
- CSV export
- JSON export
- Automated testing
- Continuous Integration

### Next Milestone

- SRA integration
- Relevance ranking
- Improved metadata normalization
- Assay filters
- Tissue filters
- Disease filters
- Excel export

### Future Development

- BioProject integration
- BioSample integration
- PubMed integration
- ENA integration
- Expression Atlas integration
- Publication-linked dataset discovery
- Reproducible search reports
- Methods-ready dataset summaries
- Stable extensible database-client framework

## Limitations

Dataset Finder is under active development.

Current limitations include:

- Only GEO and ENCODE are currently implemented.
- Search relevance depends on metadata exposed by the source database.
- Metadata completeness varies between repositories.
- Search results should be reviewed before being used in publications or downstream analyses.
- Database APIs may change independently of Dataset Finder.

## Citation and Acknowledgment

When Dataset Finder contributes to research, a publication, thesis, report, software project, or teaching material, please cite the software or acknowledge its creator.

**Creator and lead developer:** Srinivas Amla

Suggested acknowledgment:

> Public functional-genomics datasets were identified and organized using Dataset Finder, developed by Srinivas Amla.

Formal citation metadata is available in [`CITATION.cff`](CITATION.cff).

## Contributing

Community contributions are welcome.

Before contributing, read [`CONTRIBUTING.md`](CONTRIBUTING.md) for development and pull-request guidance.

Please also review:

- [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md)
- [`SECURITY.md`](SECURITY.md)

## License

Dataset Finder is distributed under the MIT License. See [`LICENSE`](LICENSE).

## Author

**Srinivas Amla**
University of Massachusetts Boston

Research interests include bioinformatics, computational biology, functional genomics, machine learning, and public biological-data integration.
EOF
