# Dataset Finder

Dataset Finder is an open-source Python toolkit for discovering, integrating, annotating, and exporting public functional genomics datasets from major biological repositories.

The project provides a unified command-line interface for searching biological databases, normalizing heterogeneous metadata, resolving gene identifiers, enriching gene annotations, classifying experimental techniques, evaluating gene-to-dataset relevance, and exporting standardized results for downstream analysis.

Dataset Finder was initially developed for large-scale screening of RNA-binding proteins, transcription factors, and other regulatory genes in *Drosophila melanogaster*. The software also supports general single-query searches across other species where the selected repository provides coverage.

## Features

Dataset Finder currently provides:

- Single-gene dataset searches
- Multi-gene batch searches
- Built-in Drosophila RNA-binding protein gene set
- Built-in Drosophila transcription factor gene set
- Custom gene lists supplied on the command line
- Plain-text gene-file input
- GEO search
- SRA search
- BioProject search
- BioStudies and ArrayExpress search
- ENCODE search for supported organisms
- FlyBase gene-symbol and identifier resolution
- FlyBase synonym and historical-name matching
- FlyAtlas tissue-expression enrichment
- Repository metadata normalization
- Gene-to-dataset relevance assessment
- Experimental technique classification
- Match-type and confidence reporting
- Study-level SRA aggregation
- Terminal, CSV, JSON, and Excel output
- Database-level success and failure reporting
- Automatic retries for transient NCBI failures
- Automated tests with pytest
- Static analysis with Ruff
- Continuous integration with GitHub Actions

## Built-in Gene Sets

The package includes two curated *Drosophila melanogaster* gene collections:

| Gene set | CLI value | Number of genes |
|----------|-----------|-----------------|
| RNA-binding proteins | `rbp` | 129 |
| Transcription factors | `tf` | 292 |

## Supported Data Sources

| Source | Status | Role |
|--------|--------|------|
| GEO | Implemented | Dataset discovery |
| SRA | Implemented | Sequencing-study discovery |
| BioProject | Implemented | Project-level metadata |
| BioStudies / ArrayExpress | Implemented | Study discovery |
| ENCODE | Implemented | Functional genomics experiments for supported organisms |
| FlyBase | Implemented | Gene annotation and identifier resolution |
| FlyAtlas | Implemented | Tissue-expression enrichment |
| BioSample | Metadata enrichment | Parsed from SRA records when available |
| PubMed | Planned | Publication discovery |
| ENA | Planned | Sequencing-data discovery |
| Expression Atlas | Planned | Expression-study integration |
| ProteomeXchange | Planned | Proteomics-data integration |

Repository coverage differs by organism. The current ENCODE portal primarily provides human and mouse experiment coverage, while the Drosophila workflow relies more heavily on GEO, SRA, BioProject, BioStudies, FlyBase, and FlyAtlas.

## Supported Experimental Techniques

Dataset Finder classifies records into the following categories:

- RNA-seq
- Single-cell RNA-seq
- Single-nucleus RNA-seq
- ChIP-seq
- CUT&RUN
- CUT&Tag
- ATAC-seq
- CLIP-seq
- eCLIP
- iCLIP
- PAR-CLIP
- HITS-CLIP
- Spatial transcriptomics
- Microarray
- Proteomics
- Other assays

The classifier also recognizes controlled GEO descriptions such as:

- `Expression profiling by high throughput sequencing`
- `Genome binding/occupancy profiling by high throughput sequencing`
- `Expression profiling by array`

## Requirements

- Python 3.11 or newer
- Git
- Internet access

The project is developed and tested primarily on macOS and Linux.

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

Install Dataset Finder and the development dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Confirm that the command is available:

```bash
dataset-finder --version
dataset-finder search --help
```

## Command-Line Interface

The main command is:

```bash
dataset-finder search
```

Available search inputs include:

- `--query`
- `--genes`
- `--gene-file`
- `--gene-set`
- `--gene-set-name`

Available databases include:

- `geo`
- `encode`
- `sra`
- `bioproject`
- `biostudies`
- `all`

Available output formats include:

- `table`
- `csv`
- `json`
- `xlsx`

Display all search options:

```bash
dataset-finder search --help
```

## Single-Query Search

Use `--query` for one direct repository search.

```bash
dataset-finder search \
  --species "Drosophila melanogaster" \
  --query orb2 \
  --database geo \
  --max-results 5
```

## Multi-Gene Search

Use `--genes` to search multiple genes in one run:

```bash
dataset-finder search \
  --species "Drosophila melanogaster" \
  --genes bru1 caz Hrb98DE orb orb2 \
  --database all \
  --max-results 3
```

In batch mode, `--max-results` limits the number of accepted records per gene. Dataset Finder may retrieve a larger internal candidate pool before relevance filtering.

## Built-in RBP Search

```bash
dataset-finder search \
  --species "Drosophila melanogaster" \
  --gene-set rbp \
  --database all \
  --max-results 3 \
  --format xlsx \
  --output Drosophila_RBP_Dataset_Screening.xlsx
```

## Built-in TF Search

```bash
dataset-finder search \
  --species "Drosophila melanogaster" \
  --gene-set tf \
  --database all \
  --max-results 3 \
  --format xlsx \
  --output Drosophila_TF_Dataset_Screening.xlsx
```

## Gene-File Input

A gene file may contain one symbol per line:

```text
bru1
Hrb98DE
orb
orb2
caz
```

Search the file:

```bash
dataset-finder search \
  --species "Drosophila melanogaster" \
  --gene-file genes.txt \
  --gene-set-name CUSTOM_SET \
  --database all \
  --max-results 3 \
  --format xlsx \
  --output custom_gene_screening.xlsx
```

## Output Formats

### Terminal

```bash
dataset-finder search \
  --species "Drosophila melanogaster" \
  --genes orb orb2 \
  --database geo \
  --max-results 3 \
  --format table
```

### CSV

```bash
dataset-finder search \
  --species "Homo sapiens" \
  --query CTCF \
  --database geo \
  --format csv \
  --output CTCF_results.csv
```

### JSON

```bash
dataset-finder search \
  --species "Mus musculus" \
  --query Tardbp \
  --database geo \
  --format json \
  --output Tardbp_results.json
```

### Excel

```bash
dataset-finder search \
  --species "Drosophila melanogaster" \
  --genes bru1 caz orb2 \
  --gene-set-name RBP_TEST \
  --database all \
  --max-results 5 \
  --format xlsx \
  --output RBP_Test_Workbook.xlsx
```

## Search Workflow

For each submitted gene, Dataset Finder performs the following steps:

1. Normalize the submitted gene input.
2. Resolve the gene through FlyBase.
3. Retrieve the official symbol, FlyBase identifier, full name, and known synonyms.
4. Build a repository search query using accepted identifiers.
5. Search one or more selected repositories.
6. Normalize repository-specific metadata into the shared `DatasetRecord` model.
7. Evaluate whether the returned record contains sufficient gene evidence.
8. Classify the experimental technique.
9. Retrieve FlyAtlas tissue-expression values when available.
10. Export the accepted records and search-status information.

```text
Gene, gene list, or built-in gene set
                  |
                  v
          FlyBase resolution
                  |
                  v
     Repository-specific searches
   GEO | SRA | BioProject | BioStudies | ENCODE
                  |
                  v
         Metadata normalization
                  |
                  v
       Gene relevance assessment
                  |
                  v
       Technique classification
                  |
                  v
        FlyAtlas enrichment
                  |
                  v
      Terminal | CSV | JSON | Excel
```

## FlyBase Gene Resolution

For *Drosophila melanogaster*, Dataset Finder resolves:

- Submitted symbol
- Current official symbol
- FlyBase gene identifier
- Current full gene name
- Symbol synonyms
- Historical symbols
- Secondary FlyBase identifiers
- Annotation identifier
- Resolution type
- Alias ambiguity status
- FlyBase record URL

Exact official symbols are prioritized over ambiguous synonyms.

Example:

```text
h -> hry -> FBgn0001168 -> hairy
```

## FlyAtlas Enrichment

FlyAtlas enrichment may include:

- Adult male brain FPKM
- Adult female brain FPKM
- Larval brain FPKM
- Adult male head FPKM
- Adult female head FPKM
- Highest-expression adult male tissue
- Highest-expression adult female tissue
- Highest-expression larval tissue
- FlyAtlas record URL

Gene annotation and FlyAtlas information are retained even when no dataset passes strict relevance filtering.

## Gene Relevance Assessment

Repository search matches are not automatically treated as valid biological matches.

Dataset Finder evaluates available metadata for evidence including:

- FlyBase identifier
- Official symbol
- Submitted symbol
- Current full gene name
- Accepted FlyBase synonyms
- Dataset title
- Study description
- Repository metadata

Short symbols require stronger context to reduce accidental matches.

Examples:

- `D` does not match the species abbreviation in `D. melanogaster`
- `orb` does not match an `ORB2`-only record
- Generic BioStudies results are excluded when no gene evidence is present

A gene may legitimately return zero accepted datasets.

## SRA Normalization

SRA results are aggregated at the study level to avoid repeated rows for experiments belonging to the same study.

Normalized SRA metadata may include:

- SRP study accession
- Study title
- BioProject accession
- SRX experiment accessions
- SRR run accessions
- BioSample accessions
- Library strategy
- Library source
- Library selection
- Library layout
- Sequencing platform
- Run count
- Release date
- Direct SRA URL

## Excel Workbook Structure

Excel exports contain the following worksheets:

| Worksheet | Description |
|-----------|-------------|
| `README` | Workbook description and search context |
| `Gene_Summary` | Per-gene counts by database and technique |
| `Gene_Annotations` | FlyBase resolution and FlyAtlas expression |
| `All_Datasets` | All accepted normalized dataset records |
| `Database_Status` | Per-gene and per-database search status |
| `Errors` | Search failures and error messages |
| Technique worksheets | Accepted datasets grouped by assay |

Technique worksheets include:

- `CUT_RUN`
- `CUT_TAG`
- `eCLIP`
- `iCLIP`
- `PAR_CLIP`
- `HITS_CLIP`
- `CLIP`
- `ChIP_seq`
- `ATAC_seq`
- `scRNA_seq`
- `snRNA_seq`
- `Spatial`
- `RNA_seq`
- `Microarray`
- `Proteomics`
- `Other_Assays`

## Database Status Reporting

The `Database_Status` worksheet records:

- Submitted gene
- Database name
- Success or failure
- Number of returned records
- Error message

A temporary repository failure is distinguishable from a successful search with zero results.

## Network Reliability

The GEO client retries transient failures such as:

- Interrupted connections
- Premature response termination
- Request timeouts
- HTTP 429
- HTTP 500
- HTTP 502
- HTTP 503
- HTTP 504

## Repository Structure

```text
Dataset-Finder/
├── .github/
│   └── workflows/
├── src/
│   └── dataset_finder/
│       ├── clients/
│       ├── data/
│       ├── exporters/
│       ├── assay_classifier.py
│       ├── batch.py
│       ├── builtin_gene_sets.py
│       ├── cli.py
│       ├── flybase_resolver.py
│       ├── gene_sets.py
│       ├── models.py
│       ├── relevance.py
│       └── search.py
├── tests/
├── tools/
├── AUTHORS.md
├── CHANGELOG.md
├── CITATION.cff
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
├── LICENSE
├── README.md
├── SECURITY.md
└── pyproject.toml
```

## Development

```bash
python -m pip install -e ".[dev]"
ruff check src tests tools
pytest -q
git diff --check
```

## Continuous Integration

GitHub Actions validates commits and pull requests using Ruff, automated tests, and supported Python versions.

## Limitations

- Repository metadata quality varies between databases.
- Some datasets may not expose enough metadata for strict gene validation.
- A zero-result gene does not necessarily mean no public data exists.
- ENCODE organism coverage differs from GEO and SRA coverage.
- External APIs may temporarily fail or change independently of Dataset Finder.
- Search results should be reviewed before publication or downstream analysis.

## Roadmap

Planned work includes:

- PubMed integration
- ENA integration
- Expression Atlas integration
- ProteomeXchange integration
- Additional metadata normalization
- Tissue and cell-type filtering
- Disease-related filtering
- Dataset ranking improvements
- Citation export
- Reproducible search reports

## Contributing

Before contributing, review:

- `CONTRIBUTING.md`
- `CODE_OF_CONDUCT.md`
- `SECURITY.md`

## Citation

If Dataset Finder contributes to a publication, thesis, report, software project, or teaching material, cite the software using the metadata provided in `CITATION.cff`.

Suggested acknowledgement:

> Public functional genomics datasets were identified and organized using Dataset Finder, developed by Srinivas Amla.

## License

Dataset Finder is distributed under the MIT License. See `LICENSE` for the complete license text.

## Author

**Srinivas Amla**

University of Massachusetts Boston
