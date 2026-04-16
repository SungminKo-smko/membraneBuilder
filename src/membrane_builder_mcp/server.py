"""MCP server for membrane building using packmol-memgen and GROMACS conversion."""

import json

from mcp.server.fastmcp import FastMCP

from membrane_builder_mcp import membrane_builder, converter, pdb_preprocessor, openmm_runner

mcp = FastMCP("membrane-builder")


@mcp.tool()
async def build_membrane(
    protein_pdb_path: str,
    output_dir: str = "./output",
    lipids: str = "POPC",
    ratio: str = "1",
    distxy_fix: float | None = None,
    dist: float = 15.0,
    dist_wat: float = 17.5,
    salt: bool = True,
    salt_concentration: float = 0.15,
    parametrize: bool = True,
    ffprot: str = "ff19SB",
    ffwat: str = "tip3p",
    fflip: str = "lipid21",
    preoriented: bool = False,
    keep_ligands: bool = True,
    keep_files: bool = True,
    ligand_params: list[str] | None = None,
    gaff2: bool = True,
    convert_to_gromacs: bool = True,
) -> str:
    """Build a membrane around a protein using packmol-memgen, then optionally convert to GROMACS format.

    Args:
        protein_pdb_path: Path to the input protein PDB file (required).
        output_dir: Directory where output files will be written.
        lipids: Lipid type(s) to use, e.g. "POPC" or "POPC:POPE".
        ratio: Ratio of lipid types, e.g. "1" or "3:1".
        distxy_fix: Fixed XY distance (Angstroms) for the membrane patch. Auto-detected if None.
        dist: Minimum distance (Angstroms) between protein and box edge in Z.
        dist_wat: Water layer thickness (Angstroms) above and below the membrane.
        salt: Whether to add salt ions to neutralize the system.
        salt_concentration: Salt concentration in mol/L (default 0.15 M).
        parametrize: Whether to run LEaP parametrization to produce AMBER topology.
        ffprot: AMBER force field for the protein (e.g. ff19SB, ff14SB).
        ffwat: Water model (e.g. tip3p, tip4pew).
        fflip: Lipid force field (e.g. lipid21, lipid17).
        preoriented: Set True if the protein is already oriented along the Z-axis.
        keep_ligands: Whether to retain non-standard residues/ligands.
        keep_files: Whether to keep intermediate build files.
        convert_to_gromacs: If True, automatically convert AMBER output to GROMACS format.

    Returns:
        JSON string with all output file paths and success status.
    """
    try:
        result = await membrane_builder.build_membrane(
            protein_pdb_path=protein_pdb_path,
            output_dir=output_dir,
            lipids=lipids,
            ratio=ratio,
            distxy_fix=distxy_fix,
            dist=dist,
            dist_wat=dist_wat,
            salt=salt,
            salt_concentration=salt_concentration,
            parametrize=parametrize,
            ffprot=ffprot,
            ffwat=ffwat,
            fflip=fflip,
            preoriented=preoriented,
            keep_ligands=keep_ligands,
            keep_files=keep_files,
            ligand_params=ligand_params,
            gaff2=gaff2,
        )

        if convert_to_gromacs and result.get("success") and parametrize:
            amber_prmtop = result.get("prmtop")
            amber_inpcrd = result.get("inpcrd")
            if amber_prmtop and amber_inpcrd:
                try:
                    gromacs_result = await converter.convert_amber_to_gromacs(
                        amber_prmtop=amber_prmtop,
                        amber_inpcrd=amber_inpcrd,
                        output_dir=output_dir,
                        output_prefix="system",
                    )
                    result["gromacs"] = gromacs_result
                except Exception as gromacs_err:
                    result["gromacs"] = {
                        "success": False,
                        "error": str(gromacs_err),
                    }
            else:
                result["gromacs"] = {
                    "success": False,
                    "error": "AMBER prmtop/inpcrd not found in build_membrane output; skipping GROMACS conversion.",
                }

        return json.dumps(result, indent=2)

    except Exception as exc:
        return json.dumps({"success": False, "error": str(exc)}, indent=2)


@mcp.tool()
async def list_available_lipids() -> str:
    """Return a list of lipids supported by packmol-memgen.

    Returns:
        JSON string containing a list of lipid dicts with name, charge, full_name, and comment fields.
    """
    try:
        lipids = await membrane_builder.list_available_lipids()
        return json.dumps(lipids, indent=2)
    except Exception as exc:
        return json.dumps({"success": False, "error": str(exc)}, indent=2)


@mcp.tool()
async def analyze_protein(protein_pdb_path: str) -> str:
    """Analyze a protein PDB structure and return structural statistics and recommended build parameters.

    Args:
        protein_pdb_path: Path to the protein PDB file (required).

    Returns:
        JSON string with atom count, residue count, chain IDs, bounding box, and recommended parameters.
    """
    try:
        result = await membrane_builder.analyze_protein(protein_pdb_path=protein_pdb_path)
        return json.dumps(result, indent=2)
    except Exception as exc:
        return json.dumps({"success": False, "error": str(exc)}, indent=2)


@mcp.tool()
async def convert_to_gromacs(
    amber_prmtop: str,
    amber_inpcrd: str,
    output_dir: str | None = None,
    output_prefix: str = "system",
) -> str:
    """Convert AMBER topology and coordinate files to GROMACS format using ParmEd.

    Args:
        amber_prmtop: Path to the AMBER parameter/topology (.prmtop) file (required).
        amber_inpcrd: Path to the AMBER coordinate (.inpcrd / .rst7) file (required).
        output_dir: Directory to write GROMACS output files. Defaults to the same directory as the input files.
        output_prefix: Prefix for the output files (default "system").

    Returns:
        JSON string with paths to the generated GROMACS files and success status.
    """
    try:
        result = await converter.convert_amber_to_gromacs(
            amber_prmtop=amber_prmtop,
            amber_inpcrd=amber_inpcrd,
            output_dir=output_dir,
            output_prefix=output_prefix,
        )
        return json.dumps(result, indent=2)
    except Exception as exc:
        return json.dumps({"success": False, "error": str(exc)}, indent=2)


@mcp.tool()
async def extract_protein_chains(
    pdb_path: str,
    chains: list[str],
    keep_ligands: list[str] | None = None,
    remove_waters: bool = True,
    output_path: str | None = None,
) -> str:
    """Extract specific chains and ligands from a PDB file.

    Use this to prepare a protein structure by keeping only the chains and
    ligands of interest.  For example, extract GPCR19 (chain R) with its
    ligands from a multi-chain cryo-EM complex before membrane building.

    Args:
        pdb_path: Path to the source PDB file (required).
        chains: List of chain IDs to keep, e.g. ["R"]. Case-sensitive.
        keep_ligands: Whitelist of HETATM residue names to keep, e.g. ["FX0", "CLR", "PLM"].
            None means keep all non-water HETATMs in the selected chains.
        remove_waters: Drop HOH records. Defaults to True.
        output_path: Destination file. None auto-generates <stem>_extracted.pdb.

    Returns:
        JSON string with output_path, atom_count, hetatm_count, chains_found, ligands_found.
    """
    try:
        result = await pdb_preprocessor.extract_chains(
            pdb_path=pdb_path,
            chains=chains,
            keep_ligands=keep_ligands,
            remove_waters=remove_waters,
            output_path=output_path,
        )
        return json.dumps(result, indent=2)
    except Exception as exc:
        return json.dumps({"success": False, "error": str(exc)}, indent=2)


@mcp.tool()
async def run_md_simulation(
    gro_path: str,
    top_path: str,
    output_dir: str = "./md_output",
    minimize: bool = True,
    min_max_iterations: int = 0,
    nvt_steps: int = 50000,
    npt_steps: int = 500000,
    production_steps: int = 5000000,
    temperature: float = 310.0,
    pressure: float = 1.0,
    timestep: float = 0.002,
    nonbonded_cutoff: float = 1.2,
    report_interval: int = 5000,
    checkpoint_interval: int = 25000,
    platform: str = "CUDA",
    precision: str = "mixed",
    gromacs_include_dir: str | None = None,
) -> str:
    """Run an OpenMM MD simulation using GROMACS input files (.gro + .top).

    Executes a full membrane protein MD pipeline:
    1. Energy minimization
    2. NVT equilibration (with position restraints)
    3. NPT equilibration (semiisotropic barostat, position restraints)
    4. Production MD (restraints released)

    Args:
        gro_path: Path to the GROMACS .gro coordinate file (required).
        top_path: Path to the GROMACS .top topology file (required).
        output_dir: Directory for output files.
        minimize: Whether to run energy minimization first.
        min_max_iterations: Max minimization iterations.
        nvt_steps: NVT equilibration steps (default 50000 = 100ps at 2fs).
        npt_steps: NPT equilibration steps (default 500000 = 1ns at 2fs).
        production_steps: Production MD steps (default 5000000 = 10ns at 2fs).
        temperature: Temperature in Kelvin (default 310 K).
        pressure: Pressure in bar (default 1.0 bar).
        timestep: Integration timestep in ps (default 0.002 = 2fs).
        nonbonded_cutoff: PME cutoff in nm (default 1.2).
        report_interval: Steps between trajectory/energy writes.
        checkpoint_interval: Steps between checkpoint saves.
        platform: OpenMM platform — CUDA, OpenCL, or CPU. Auto-fallback.
        precision: GPU precision — mixed, single, or double.
        gromacs_include_dir: Path to GROMACS include dir for force field files.

    Returns:
        JSON string with success, stage, output_files, minimized_energy, final_temperature.
    """
    try:
        result = await openmm_runner.run_openmm_simulation(
            gro_path=gro_path,
            top_path=top_path,
            output_dir=output_dir,
            minimize=minimize,
            min_max_iterations=min_max_iterations,
            nvt_steps=nvt_steps,
            npt_steps=npt_steps,
            production_steps=production_steps,
            temperature=temperature,
            pressure=pressure,
            timestep=timestep,
            nonbonded_cutoff=nonbonded_cutoff,
            report_interval=report_interval,
            checkpoint_interval=checkpoint_interval,
            platform=platform,
            precision=precision,
            gromacs_include_dir=gromacs_include_dir,
        )
        return json.dumps(result, indent=2)
    except Exception as exc:
        return json.dumps({"success": False, "error": str(exc)}, indent=2)


if __name__ == "__main__":
    mcp.run(transport="stdio")
