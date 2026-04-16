"""tests/test_openmm_runner.py — Tests for openmm_runner.run_openmm_simulation.

Tests that require OpenMM are skipped automatically when the library is not
installed in the current Python environment.
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

from membrane_builder_mcp.openmm_runner import (  # noqa: E402
    _OPENMM_AVAILABLE,
    run_openmm_simulation,
)

# ---------------------------------------------------------------------------
# Paths to example GROMACS files
# ---------------------------------------------------------------------------
_GROMACS_DIR = _REPO_ROOT / "examples" / "output_7cfn" / "gromacs"
_GRO_FILE = _GROMACS_DIR / "membrane_7CFN.gro"
_TOP_FILE = _GROMACS_DIR / "membrane_7CFN.top"


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_openmm_import_flag() -> None:
    """_OPENMM_AVAILABLE must be a plain bool regardless of whether OpenMM is installed."""
    assert isinstance(_OPENMM_AVAILABLE, bool)


@pytest.mark.asyncio
async def test_run_without_openmm() -> None:
    """When OpenMM is absent, run_openmm_simulation must raise ImportError with install instructions."""
    if _OPENMM_AVAILABLE:
        pytest.skip("OpenMM is installed — this test only applies when it is absent.")

    with pytest.raises(ImportError, match="conda install"):
        await run_openmm_simulation(
            gro_path=str(_GRO_FILE),
            top_path=str(_TOP_FILE),
        )


@pytest.mark.asyncio
async def test_missing_gro_file(tmp_path: Path) -> None:
    """FileNotFoundError must be raised for a non-existent .gro file."""
    if not _OPENMM_AVAILABLE:
        # The public function raises ImportError before the file check; call
        # the internal _run_simulation directly to exercise the file guard.
        from membrane_builder_mcp.openmm_runner import _run_simulation

        with pytest.raises((ImportError, FileNotFoundError)):
            _run_simulation(
                gro_path="/nonexistent/path/missing.gro",
                top_path=str(_TOP_FILE),
                output_dir=tmp_path,
                minimize=False,
                min_max_iterations=0,
                nvt_steps=0,
                npt_steps=0,
                production_steps=0,
                temperature=310.0,
                pressure=1.0,
                timestep=0.002,
                nonbonded_cutoff=1.2,
                report_interval=10,
                checkpoint_interval=50,
                platform="CPU",
                precision="mixed",
                gromacs_include_dir=None,
            )
    else:
        from membrane_builder_mcp.openmm_runner import _run_simulation

        with pytest.raises(FileNotFoundError, match=r"\.gro"):
            _run_simulation(
                gro_path="/nonexistent/path/missing.gro",
                top_path=str(_TOP_FILE),
                output_dir=tmp_path,
                minimize=False,
                min_max_iterations=0,
                nvt_steps=0,
                npt_steps=0,
                production_steps=0,
                temperature=310.0,
                pressure=1.0,
                timestep=0.002,
                nonbonded_cutoff=1.2,
                report_interval=10,
                checkpoint_interval=50,
                platform="CPU",
                precision="mixed",
                gromacs_include_dir=None,
            )


@pytest.mark.asyncio
async def test_missing_top_file(tmp_path: Path) -> None:
    """FileNotFoundError must be raised for a non-existent .top file."""
    from membrane_builder_mcp.openmm_runner import _run_simulation

    if not _OPENMM_AVAILABLE:
        with pytest.raises((ImportError, FileNotFoundError)):
            _run_simulation(
                gro_path=str(_GRO_FILE),
                top_path="/nonexistent/path/missing.top",
                output_dir=tmp_path,
                minimize=False,
                min_max_iterations=0,
                nvt_steps=0,
                npt_steps=0,
                production_steps=0,
                temperature=310.0,
                pressure=1.0,
                timestep=0.002,
                nonbonded_cutoff=1.2,
                report_interval=10,
                checkpoint_interval=50,
                platform="CPU",
                precision="mixed",
                gromacs_include_dir=None,
            )
    else:
        with pytest.raises(FileNotFoundError, match=r"\.top"):
            _run_simulation(
                gro_path=str(_GRO_FILE),
                top_path="/nonexistent/path/missing.top",
                output_dir=tmp_path,
                minimize=False,
                min_max_iterations=0,
                nvt_steps=0,
                npt_steps=0,
                production_steps=0,
                temperature=310.0,
                pressure=1.0,
                timestep=0.002,
                nonbonded_cutoff=1.2,
                report_interval=10,
                checkpoint_interval=50,
                platform="CPU",
                precision="mixed",
                gromacs_include_dir=None,
            )


@pytest.mark.skipif(not _OPENMM_AVAILABLE, reason="OpenMM not installed")
def test_gromacs_file_loading() -> None:
    """Verify OpenMM can load the example GROMACS files and create a system.

    The full system (~234k atoms) is too large for MD on CPU in test time,
    so we only verify file loading + system creation here.
    """
    if not _GRO_FILE.is_file():
        pytest.skip(f"Example .gro file not found: {_GRO_FILE}")
    if not _TOP_FILE.is_file():
        pytest.skip(f"Example .top file not found: {_TOP_FILE}")

    from openmm.app import GromacsGroFile, GromacsTopFile, PME, HBonds
    from openmm.unit import nanometer

    gro = GromacsGroFile(str(_GRO_FILE))
    top = GromacsTopFile(str(_TOP_FILE), periodicBoxVectors=gro.getPeriodicBoxVectors())

    # Verify positions loaded
    assert len(gro.positions) > 200000, f"Expected ~234k atoms, got {len(gro.positions)}"

    # Verify system creation with PME (the most critical step)
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.2 * nanometer,
        constraints=HBonds,
    )
    assert system.getNumParticles() == len(gro.positions)
