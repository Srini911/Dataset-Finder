"""Shared data models for Dataset Finder."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class DatasetRecord:
    """A normalized public-dataset search result."""

    uid: str
    accession: str
    title: str
    organism: str
    study_type: str
    sample_count: int | None
    publication_date: str
    url: str
