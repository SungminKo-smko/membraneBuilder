"""Shared fixtures and path configuration for membrane builder tests."""

import sys
from pathlib import Path

# Ensure the src package is importable
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pytest


@pytest.fixture
def sample_pdb_path() -> Path:
    """Return the path to the sample protein PDB fixture file."""
    path = Path(__file__).parent / "fixtures" / "sample_protein.pdb"
    assert path.exists(), f"Sample PDB fixture not found at {path}"
    return path
