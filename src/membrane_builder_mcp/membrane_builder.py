"""Core module for generating Packmol input files and building lipid bilayer membranes."""

from __future__ import annotations

import asyncio
import logging
import math
import os
import shutil
import subprocess
from pathlib import Path

from .defaults import (
    DEFAULT_BILAYER_THICKNESS,
    DEFAULT_LIPID_AREA,
    DEFAULT_LIPID_COMPOSITION,
    DEFAULT_TOLERANCE,
    DEFAULT_WATER_DENSITY,
    DEFAULT_WATER_PADDING,
    DEFAULT_XY_PADDING,
    LIPID_DB,
)
from .pdb_utils import estimate_cross_section, parse_pdb

logger = logging.getLogger(__name__)

LIPIDS_DIR = Path(__file__).parent / "lipids"


def get_lipid_pdb_path(lipid_name: str, leaflet: str) -> Path:
    """Return the path to a lipid PDB template file.

    Args:
        lipid_name: Lipid type key as it appears in LIPID_DB (e.g. "POPC", "cholesterol").
        leaflet: Which leaflet -- "upper" or "lower".

    Returns:
        Absolute path to the lipid PDB file.

    Raises:
        ValueError: If the lipid name is unknown or the leaflet is invalid.
        FileNotFoundError: If the template PDB file does not exist.
    """
    if lipid_name not in LIPID_DB:
        raise ValueError(
            f"Unknown lipid '{lipid_name}'. Available: {list(LIPID_DB.keys())}"
        )
    if leaflet not in ("upper", "lower"):
        raise ValueError(f"Leaflet must be 'upper' or 'lower', got '{leaflet}'")

    filename = f"{lipid_name}_{leaflet}.pdb"
    path = LIPIDS_DIR / filename

    if not path.exists():
        raise FileNotFoundError(f"Lipid template not found: {path}")

    return path


def calculate_lipid_counts(
    membrane_x: float,
    membrane_y: float,
    cross_section: float,
    lipid_composition: dict[str, float],
    lipid_area: float,
) -> dict[str, int]:
    """Calculate the number of each lipid type per leaflet.

    Args:
        membrane_x: Membrane X dimension in Angstroms.
        membrane_y: Membrane Y dimension in Angstroms.
        cross_section: Protein cross-sectional area in Angstroms squared.
        lipid_composition: Dict mapping lipid name to its mole fraction (must sum to ~1.0).
        lipid_area: Area per lipid in Angstroms squared.

    Returns:
        Dict mapping lipid name to count per leaflet (e.g. {"POPC": 150, "POPE": 50}).
    """
    available_area = membrane_x * membrane_y - cross_section
    if available_area <= 0:
        logger.warning(
            "Available area (%.1f A^2) is non-positive; protein may be too large "
            "for the specified membrane dimensions.",
            available_area,
        )
        return {name: 0 for name in lipid_composition}

    total_lipids = available_area / lipid_area

    counts: dict[str, int] = {}
    for name, fraction in lipid_composition.items():
        counts[name] = max(int(math.floor(total_lipids * fraction)), 0)

    return counts


def generate_packmol_input(
    protein_pdb_path: str,
    output_path: str,
    lipid_counts: dict[str, int],
    membrane_x: float,
    membrane_y: float,
    z_center: float,
    bilayer_thickness: float,
    water_padding: float,
    water_density: float,
    tolerance: float,
    n_water_upper: int,
    n_water_lower: int,
) -> str:
    """Generate Packmol input file content as a string.

    Args:
        protein_pdb_path: Path to the protein PDB file.
        output_path: Path for the Packmol output PDB file.
        lipid_counts: Dict mapping lipid name to count per leaflet.
        membrane_x: Membrane X dimension in Angstroms.
        membrane_y: Membrane Y dimension in Angstroms.
        z_center: Z coordinate of the bilayer center.
        bilayer_thickness: Thickness of the lipid bilayer in Angstroms.
        water_padding: Water layer thickness above/below bilayer in Angstroms.
        water_density: Water density in molecules per cubic Angstrom.
        tolerance: Minimum distance between molecules in Angstroms.
        n_water_upper: Number of water molecules above the upper leaflet.
        n_water_lower: Number of water molecules below the lower leaflet.

    Returns:
        Packmol input file content as a string.
    """
    half_thickness = bilayer_thickness / 2.0
    xmin = -membrane_x / 2.0
    xmax = membrane_x / 2.0
    ymin = -membrane_y / 2.0
    ymax = membrane_y / 2.0

    z_upper_top = z_center + half_thickness
    z_lower_bottom = z_center - half_thickness
    z_water_top = z_upper_top + water_padding
    z_water_bottom = z_lower_bottom - water_padding

    water_pdb_path = str(LIPIDS_DIR / "water.pdb")

    lines: list[str] = []
    lines.append(f"tolerance {tolerance:.1f}")
    lines.append("filetype pdb")
    lines.append(f"output {output_path}")
    lines.append("seed -1")
    lines.append("")

    # Protein (fixed)
    lines.append("# Protein (fixed)")
    lines.append(f"structure {protein_pdb_path}")
    lines.append("  number 1")
    cx, cy, cz = 0.0, 0.0, z_center  # protein will be centered at origin in XY
    lines.append(f"  fixed {cx:.3f} {cy:.3f} {cz:.3f} 0. 0. 0.")
    lines.append("  center")
    lines.append("end structure")
    lines.append("")

    # Lipids -- upper leaflet
    for lipid_name, count in lipid_counts.items():
        if count <= 0:
            continue
        upper_pdb = str(get_lipid_pdb_path(lipid_name, "upper"))
        lines.append(f"# {lipid_name} upper leaflet")
        lines.append(f"structure {upper_pdb}")
        lines.append(f"  number {count}")
        lines.append(
            f"  inside box {xmin:.3f} {ymin:.3f} {z_center:.3f} "
            f"{xmax:.3f} {ymax:.3f} {z_upper_top:.3f}"
        )
        lines.append("end structure")
        lines.append("")

    # Lipids -- lower leaflet
    for lipid_name, count in lipid_counts.items():
        if count <= 0:
            continue
        lower_pdb = str(get_lipid_pdb_path(lipid_name, "lower"))
        lines.append(f"# {lipid_name} lower leaflet")
        lines.append(f"structure {lower_pdb}")
        lines.append(f"  number {count}")
        lines.append(
            f"  inside box {xmin:.3f} {ymin:.3f} {z_lower_bottom:.3f} "
            f"{xmax:.3f} {ymax:.3f} {z_center:.3f}"
        )
        lines.append("end structure")
        lines.append("")

    # Water upper
    if n_water_upper > 0:
        lines.append("# Water upper (above upper leaflet)")
        lines.append(f"structure {water_pdb_path}")
        lines.append(f"  number {n_water_upper}")
        lines.append(
            f"  inside box {xmin:.3f} {ymin:.3f} {z_upper_top:.3f} "
            f"{xmax:.3f} {ymax:.3f} {z_water_top:.3f}"
        )
        lines.append("end structure")
        lines.append("")

    # Water lower
    if n_water_lower > 0:
        lines.append("# Water lower (below lower leaflet)")
        lines.append(f"structure {water_pdb_path}")
        lines.append(f"  number {n_water_lower}")
        lines.append(
            f"  inside box {xmin:.3f} {ymin:.3f} {z_water_bottom:.3f} "
            f"{xmax:.3f} {ymax:.3f} {z_lower_bottom:.3f}"
        )
        lines.append("end structure")
        lines.append("")

    return "\n".join(lines)


async def build_membrane(
    protein_pdb_path: str,
    output_dir: str = "./output",
    lipid_composition: dict[str, float] | None = None,
    membrane_x_size: float | None = None,
    membrane_y_size: float | None = None,
    bilayer_thickness: float = DEFAULT_BILAYER_THICKNESS,
    water_padding: float = DEFAULT_WATER_PADDING,
    lipid_area: float = DEFAULT_LIPID_AREA,
    water_density: float = DEFAULT_WATER_DENSITY,
    tolerance: float = DEFAULT_TOLERANCE,
) -> dict:
    """Build a lipid bilayer membrane around a protein using Packmol.

    This is the main orchestration function that:
    1. Parses the protein PDB
    2. Calculates membrane dimensions and lipid counts
    3. Generates the Packmol input file
    4. Runs Packmol (if available)

    Args:
        protein_pdb_path: Path to the protein PDB file.
        output_dir: Directory for output files.
        lipid_composition: Dict mapping lipid names to mole fractions.
            Defaults to DEFAULT_LIPID_COMPOSITION (pure POPC).
        membrane_x_size: X dimension of the membrane in Angstroms.
            If None, auto-calculated from protein bounding box.
        membrane_y_size: Y dimension of the membrane in Angstroms.
            If None, auto-calculated from protein bounding box.
        bilayer_thickness: Bilayer thickness in Angstroms.
        water_padding: Water layer padding in Angstroms.
        lipid_area: Area per lipid in Angstroms squared.
        water_density: Water molecule density in molecules per cubic Angstrom.
        tolerance: Minimum distance between molecules in Angstroms.

    Returns:
        Dict with build results including lipid counts, water counts,
        box dimensions, and packmol execution status.
    """
    if lipid_composition is None:
        lipid_composition = dict(DEFAULT_LIPID_COMPOSITION)

    # Step a: Parse protein PDB
    pdb_data = parse_pdb(protein_pdb_path)
    atoms = pdb_data["atoms"]
    center = pdb_data["center"]
    bbox = pdb_data["bbox"]
    z_center = center[2]

    cross_section = estimate_cross_section(atoms, z_center, bilayer_thickness)

    # Step b: Calculate membrane dimensions
    if membrane_x_size is None:
        membrane_x_size = bbox["size"][0] + DEFAULT_XY_PADDING * 2
    if membrane_y_size is None:
        membrane_y_size = bbox["size"][1] + DEFAULT_XY_PADDING * 2

    # Step c: Calculate lipid counts per leaflet
    lipid_counts = calculate_lipid_counts(
        membrane_x=membrane_x_size,
        membrane_y=membrane_y_size,
        cross_section=cross_section,
        lipid_composition=lipid_composition,
        lipid_area=lipid_area,
    )

    # Step d: Calculate water molecule counts
    water_volume_upper = membrane_x_size * membrane_y_size * water_padding
    water_volume_lower = membrane_x_size * membrane_y_size * water_padding
    n_water_upper = int(math.floor(water_volume_upper * water_density))
    n_water_lower = int(math.floor(water_volume_lower * water_density))

    # Step e: Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Step f & g: Generate and write packmol input
    output_pdb_path = os.path.join(output_dir, "membrane.pdb")
    inp_content = generate_packmol_input(
        protein_pdb_path=os.path.abspath(protein_pdb_path),
        output_path=os.path.abspath(output_pdb_path),
        lipid_counts=lipid_counts,
        membrane_x=membrane_x_size,
        membrane_y=membrane_y_size,
        z_center=z_center,
        bilayer_thickness=bilayer_thickness,
        water_padding=water_padding,
        water_density=water_density,
        tolerance=tolerance,
        n_water_upper=n_water_upper,
        n_water_lower=n_water_lower,
    )

    inp_file_path = os.path.join(output_dir, "membrane.inp")
    with open(inp_file_path, "w") as f:
        f.write(inp_content)

    logger.info("Packmol input written to %s", inp_file_path)

    # Step h: Check if packmol is available
    packmol_available = shutil.which("packmol") is not None

    packmol_success: bool | None = None
    packmol_log: str = ""

    if not packmol_available:
        logger.warning(
            "Packmol is not installed or not found in PATH. "
            "Input file was generated but Packmol was not run."
        )
        packmol_log = "Packmol is not installed. Please install it and run manually."
    else:
        # Step i: Run packmol via asyncio.to_thread
        def _run_packmol() -> subprocess.CompletedProcess:
            return subprocess.run(
                ["packmol"],
                input=inp_content,
                capture_output=True,
                text=True,
                timeout=600,
            )

        try:
            result = await asyncio.to_thread(_run_packmol)
            stdout = result.stdout or ""
            stderr = result.stderr or ""

            # Step j: Check convergence
            success_markers = ("Success", "ENDED WITHOUT SYMMETRY")
            packmol_success = any(marker in stdout for marker in success_markers)

            # Last 20 lines of stdout
            stdout_lines = stdout.strip().splitlines()
            packmol_log = "\n".join(stdout_lines[-20:])

            if not packmol_success:
                logger.warning("Packmol did not converge. stderr: %s", stderr)
                if stderr:
                    packmol_log += "\n--- stderr ---\n" + stderr

        except subprocess.TimeoutExpired:
            packmol_success = False
            packmol_log = "Packmol timed out after 600 seconds."
            logger.error("Packmol timed out.")
        except Exception as exc:
            packmol_success = False
            packmol_log = f"Packmol execution failed: {exc}"
            logger.error("Packmol execution failed: %s", exc)

    # Step k: Build result dict
    half_thickness = bilayer_thickness / 2.0
    box_dimensions = {
        "x": membrane_x_size,
        "y": membrane_y_size,
        "z_min": z_center - half_thickness - water_padding,
        "z_max": z_center + half_thickness + water_padding,
    }

    return {
        "output_path": os.path.abspath(output_pdb_path),
        "input_path": os.path.abspath(inp_file_path),
        "n_lipids_upper": dict(lipid_counts),
        "n_lipids_lower": dict(lipid_counts),
        "n_water_upper": n_water_upper,
        "n_water_lower": n_water_lower,
        "box_dimensions": box_dimensions,
        "protein_center": center,
        "cross_section_area": cross_section,
        "packmol_success": packmol_success,
        "packmol_log": packmol_log,
    }
