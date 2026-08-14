# Changelog

All notable changes to Dataset Finder will be documented in this file.

The format is based on Keep a Changelog, and the project follows semantic
versioning.

## [Unreleased]

### Added

- Expanded public-database support for SRA, BioProject, BioSample, BioStudies/ArrayExpress, ENA, PubMed, Expression Atlas, and PRIDE
- Packaged Drosophila RNA-binding protein and transcription-factor gene sets
- FlyBase gene resolution and FlyAtlas annotation
- Biological metadata extraction for tissue, sex, genotype, treatment, disease, and related study information
- Historical technique-specific GEO and SRA searches
- Dataset ranking, relevance scoring, and known-dataset regression checks
- Focused Drosophila RBP/TF screening workflow with canonical identifier and historical-alias searches
- Validation utilities for historical dataset recovery and ambiguous gene-symbol handling

### Changed

- Improved SRA study aggregation and metadata normalization
- Improved multi-gene recall and gene-to-dataset relevance filtering
- Added neural-dataset highlighting and richer Excel exports
- Reorganized validation utilities under `tools/validation`
- Renamed the focused screening utility to `screen_drosophila_gene_sets.py`
- Made NCBI email and validated legacy-registry configuration portable

## [0.3.0] - 2026-08-02

### Added

- Multi-gene dataset screening
- Built-in Drosophila gene-set support
- CSV and JSON export
- Expanded repository integration and project documentation

### Changed

- Improved combined repository searches and result handling
- Updated continuous-integration and project metadata

### Verified

- Automated tests and Ruff static-analysis checks

## [0.2.0] - 2026-07-19

### Added

- ENCODE REST API client
- ENCODE experiment search support
- Combined GEO and ENCODE searching
- ENCODE client unit tests
- Expanded installation, usage, architecture, roadmap, and citation documentation
- GitHub Actions continuous-integration workflow
- Formal citation metadata
- Community contribution, security, authorship, and conduct documentation

### Changed

- Updated the package version from 0.1.0 to 0.2.0
- Updated project metadata to include ENCODE
- Limited CLI database choices to implemented user-facing integrations
- Interleaved combined GEO and ENCODE results so one source does not
  automatically hide the other
- Clarified implemented and planned capabilities

### Verified

- Ruff static-analysis checks
- Automated unit tests on Python 3.11, 3.12, and 3.13 through GitHub Actions

## [0.1.0]

### Added

- Initial Dataset Finder Python package
- GEO dataset search
- Shared dataset-record model
- Search-service orchestration
- Command-line interface
- Initial testing and project configuration

[Unreleased]: https://github.com/Srini911/Dataset-Finder/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/Srini911/Dataset-Finder/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/Srini911/Dataset-Finder/releases/tag/v0.2.0
