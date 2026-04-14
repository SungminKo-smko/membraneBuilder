"""MCP server for the Packmol membrane builder."""

import json
from dataclasses import asdict

from mcp.server.fastmcp import FastMCP

from .defaults import LIPID_DB
from .membrane_builder import build_membrane
from .pdb_utils import suggest_membrane_params

mcp = FastMCP("packmol-membrane-builder")


@mcp.tool()
async def build_membrane_tool(
    protein_pdb_path: str,
    output_dir: str = "./output",
    lipid_composition: dict[str, float] | None = None,
    membrane_x_size: float | None = None,
    membrane_y_size: float | None = None,
    bilayer_thickness: float = 40.0,
    water_padding: float = 15.0,
    lipid_area: float = 65.0,
    water_density: float = 0.0334,
    tolerance: float = 2.0,
) -> str:
    """Build a lipid bilayer membrane around a protein using Packmol.

    Args:
        protein_pdb_path: Path to the protein PDB file.
        output_dir: Output directory for generated files.
        lipid_composition: Lipid composition ratios (sum=1.0), e.g. {"POPC": 0.7, "POPE": 0.3}.
        membrane_x_size: Membrane X size in angstroms. Auto-calculated if not specified.
        membrane_y_size: Membrane Y size in angstroms. Auto-calculated if not specified.
        bilayer_thickness: Bilayer thickness in angstroms.
        water_padding: Water layer padding in angstroms.
        lipid_area: Area per lipid in angstroms squared.
        water_density: Water density in molecules per cubic angstrom.
        tolerance: Packmol tolerance in angstroms.
    """
    try:
        result = await build_membrane(
            protein_pdb_path=protein_pdb_path,
            output_dir=output_dir,
            lipid_composition=lipid_composition,
            membrane_x_size=membrane_x_size,
            membrane_y_size=membrane_y_size,
            bilayer_thickness=bilayer_thickness,
            water_padding=water_padding,
            lipid_area=lipid_area,
            water_density=water_density,
            tolerance=tolerance,
        )
        return json.dumps(result, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e)}, indent=2)


@mcp.tool()
async def list_available_lipids() -> str:
    """List all available lipid types with their properties."""
    try:
        lipids = {
            name: asdict(info)
            for name, info in LIPID_DB.items()
        }
        return json.dumps(lipids, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e)}, indent=2)


@mcp.tool()
async def analyze_protein_structure(protein_pdb_path: str) -> str:
    """Analyze a protein PDB file and suggest membrane building parameters.

    Args:
        protein_pdb_path: Path to the protein PDB file.
    """
    try:
        params = suggest_membrane_params(protein_pdb_path)
        return json.dumps(params, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e)}, indent=2)


if __name__ == "__main__":
    mcp.run(transport="stdio")
