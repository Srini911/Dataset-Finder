# Contributing to Dataset Finder

Thank you for considering a contribution to Dataset Finder.

## Development setup

Fork or clone the repository and create a virtual environment:

```bash
git clone https://github.com/Srini911/Dataset-Finder.git
cd Dataset-Finder
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

## Branch workflow

Create a focused branch from the latest `main` branch:

```bash
git switch main
git pull --ff-only origin main
git switch -c feature/brief-description
```

Use `fix/brief-description` for bug fixes and
`docs/brief-description` for documentation-only changes.

## Code standards

Contributions should:

- support Python 3.11 or newer;
- use type annotations for public functions and methods;
- keep database-specific behavior inside the appropriate client module;
- return normalized `DatasetRecord` objects from search clients;
- include clear error messages for network and validation failures;
- avoid committing generated outputs, credentials, caches, or virtual environments.

## Required checks

Run these commands before submitting a pull request:

```bash
ruff check src tests
pytest
```

New behavior should include appropriate automated tests.

## Pull requests

A pull request should:

- explain the problem and the proposed solution;
- describe important implementation decisions;
- include tests for new or changed behavior;
- update the README or changelog when user-facing behavior changes;
- remain focused on one logical change.

## Bug reports and feature requests

Use GitHub Issues and include:

- the Dataset Finder version;
- the Python version;
- the operating system;
- the exact command that was run;
- the complete error message;
- the expected and observed behavior.

Do not include passwords, API tokens, private data, or unpublished
sensitive research information.
