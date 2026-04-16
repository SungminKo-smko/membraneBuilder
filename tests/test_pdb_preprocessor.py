"""tests/test_pdb_preprocessor.py — Tests for pdb_preprocessor.extract_chains.

Uses the example PDB file at examples/7CFN.pdb (resolved relative to this
test file so pytest can be run from any working directory).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Make sure the package is importable when running without installation
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC_DIR = _REPO_ROOT / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from membrane_builder_mcp.pdb_preprocessor import extract_chains  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_EXAMPLES_DIR = _REPO_ROOT / "examples"
_SAMPLE_PDB = _EXAMPLES_DIR / "7CFN.pdb"


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _count_record(pdb_file: Path, record_type: str) -> int:
    """Count lines whose first 6 characters match *record_type* (strip-matched)."""
    count = 0
    with pdb_file.open() as fh:
        for line in fh:
            if line[:6].strip() == record_type:
                count += 1
    return count


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_extract_chain_r(tmp_path: Path) -> None:
    """Extract chain R only — verify basic success fields."""
    out = tmp_path / "chain_R.pdb"

    result = await extract_chains(
        pdb_path=str(_SAMPLE_PDB),
        chains=["R"],
        output_path=str(out),
    )

    assert result["success"] is True
    assert result["chains_found"] == ["R"]
    assert result["atom_count"] > 0
    assert out.exists()


@pytest.mark.asyncio
async def test_extract_with_ligands(tmp_path: Path) -> None:
    """Extract chain R together with ligands FX0, CLR, PLM."""
    out = tmp_path / "chain_R_ligands.pdb"

    result = await extract_chains(
        pdb_path=str(_SAMPLE_PDB),
        chains=["R"],
        keep_ligands=["FX0", "CLR", "PLM"],
        remove_waters=True,
        output_path=str(out),
    )

    assert result["success"] is True

    # All three requested ligands must appear in the output
    found = set(result["ligands_found"])
    for lig in ("FX0", "CLR", "PLM"):
        assert lig in found, f"Expected ligand {lig!r} not found; got {found}"

    assert result["hetatm_count"] > 0

    # HOH must not appear in the output
    with out.open() as fh:
        content = fh.read()
    assert "HOH" not in content, "HOH water records should have been removed"


@pytest.mark.asyncio
async def test_remove_waters(tmp_path: Path) -> None:
    """HOH lines must be absent when remove_waters=True (the default)."""
    out = tmp_path / "no_water.pdb"

    await extract_chains(
        pdb_path=str(_SAMPLE_PDB),
        chains=["R"],
        remove_waters=True,
        output_path=str(out),
    )

    assert out.exists()
    with out.open() as fh:
        for line in fh:
            # No HETATM line should contain HOH in the residue name field (cols 17-19)
            if line[:6].strip() == "HETATM":
                resname = line[17:20].strip()
                assert resname != "HOH", f"HOH found in output: {line.rstrip()!r}"


@pytest.mark.asyncio
async def test_output_path(tmp_path: Path) -> None:
    """A custom output_path must be honoured."""
    custom_out = tmp_path / "subdir" / "my_output.pdb"

    result = await extract_chains(
        pdb_path=str(_SAMPLE_PDB),
        chains=["R"],
        output_path=str(custom_out),
    )

    assert result["success"] is True
    assert Path(result["output_path"]).resolve() == custom_out.resolve()
    assert custom_out.exists()


@pytest.mark.asyncio
async def test_missing_file() -> None:
    """extract_chains must raise FileNotFoundError for a non-existent PDB."""
    with pytest.raises(FileNotFoundError):
        await extract_chains(
            pdb_path="/nonexistent/path/no_such_file.pdb",
            chains=["A"],
        )
