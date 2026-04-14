"""PDB file parsing and protein structure analysis utilities for Packmol membrane builder."""

import numpy as np


def parse_pdb(pdb_path: str) -> dict:
    """Parse a PDB file and extract atom information, center of mass, and bounding box.

    Args:
        pdb_path: Path to the PDB file.

    Returns:
        Dictionary with keys:
            - 'atoms': list of dicts with keys (name, resname, chain, resnum, x, y, z)
            - 'center': tuple (cx, cy, cz) simple average of all atom coordinates
            - 'bbox': dict with 'min', 'max', and 'size' tuples
    """
    atoms = []

    with open(pdb_path, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue

            atom = {
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resnum": int(line[22:26].strip()),
                "x": float(line[30:38].strip()),
                "y": float(line[38:46].strip()),
                "z": float(line[46:54].strip()),
            }
            atoms.append(atom)

    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms])

    center = tuple(coords.mean(axis=0).tolist())

    coord_min = tuple(coords.min(axis=0).tolist())
    coord_max = tuple(coords.max(axis=0).tolist())
    size = tuple((coords.max(axis=0) - coords.min(axis=0)).tolist())

    bbox = {
        "min": coord_min,
        "max": coord_max,
        "size": size,
    }

    return {
        "atoms": atoms,
        "center": center,
        "bbox": bbox,
    }


def estimate_cross_section(atoms: list, z_center: float, thickness: float) -> float:
    """Estimate the cross-sectional area of atoms in the XY plane at a given Z slice.

    Filters atoms within z_center +/- thickness/2, projects onto the XY plane,
    and counts occupied cells on a 2A resolution grid.

    Args:
        atoms: List of atom dicts (each with 'x', 'y', 'z' keys).
        z_center: Center of the Z slice.
        thickness: Thickness of the Z slice.

    Returns:
        Cross-sectional area in Angstroms squared.
    """
    half = thickness / 2.0
    z_lo = z_center - half
    z_hi = z_center + half

    filtered = [a for a in atoms if z_lo <= a["z"] <= z_hi]

    if not filtered:
        return 0.0

    coords = np.array([[a["x"], a["y"]] for a in filtered])

    grid_resolution = 2.0  # Angstroms

    # Snap each XY coordinate to a grid cell index
    grid_indices = np.floor(coords / grid_resolution).astype(int)

    # Count unique occupied cells
    unique_cells = set(map(tuple, grid_indices))
    cell_area = grid_resolution ** 2

    return len(unique_cells) * cell_area


def suggest_membrane_params(pdb_path: str, xy_padding: float = 15.0) -> dict:
    """Suggest membrane building parameters based on protein structure analysis.

    Args:
        pdb_path: Path to the PDB file.
        xy_padding: Padding to add on each side in the XY plane (Angstroms).

    Returns:
        Dictionary with suggested parameters:
            - 'membrane_x_size': suggested X dimension of the membrane
            - 'membrane_y_size': suggested Y dimension of the membrane
            - 'n_lipids_per_leaflet': suggested number of lipids per leaflet
            - 'protein_center': center of the protein
            - 'protein_bbox': bounding box of the protein
            - 'cross_section_area': estimated cross-sectional area at the membrane plane
    """
    pdb_data = parse_pdb(pdb_path)
    atoms = pdb_data["atoms"]
    center = pdb_data["center"]
    bbox = pdb_data["bbox"]

    membrane_x = bbox["size"][0] + 2.0 * xy_padding
    membrane_y = bbox["size"][1] + 2.0 * xy_padding

    cross_section = estimate_cross_section(atoms, z_center=center[2], thickness=10.0)

    membrane_area = membrane_x * membrane_y
    area_per_lipid = 65.0  # approximate area per lipid in Angstroms squared
    n_lipids = int((membrane_area - cross_section) / area_per_lipid)

    return {
        "membrane_x_size": round(membrane_x, 1),
        "membrane_y_size": round(membrane_y, 1),
        "n_lipids_per_leaflet": max(n_lipids, 0),
        "protein_center": center,
        "protein_bbox": bbox,
        "cross_section_area": round(cross_section, 1),
    }
