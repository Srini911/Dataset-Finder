# Changelog

All notable changes to Dataset Finder will be documented in this file.

The format is based on Keep a Changelog, and the project follows semantic
versioning.

## [Unreleased]

### Planned

- Structured export utilities for CSV, JSON, and Excel
- Dataset relevance ranking and quality scoring
- Additional public biological-database integrations
- Structured filters for tissues, diseases, assays, and regulators

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

[Unreleased]: https://github.com/Srini911/Dataset-Finder/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Srini911/Dataset-Finder/releases/tag/v0.2.0
