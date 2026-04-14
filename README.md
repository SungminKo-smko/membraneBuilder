# membrane-builder-mcp

An MCP (Model Context Protocol) server that builds lipid bilayer membrane systems around proteins using AmberTools25 `packmol-memgen`. Provide a protein PDB file and the server invokes `packmol-memgen` to embed the protein in a solvated, ionized membrane, then optionally converts the output to GROMACS format via ParmEd.

## Features

- Wraps the AmberTools25 `packmol-memgen` CLI as an MCP server
- Supports 3000+ lipid types from the Lipid21 force field library
- Mixed bilayer composition with configurable leaflet ratios
- Automatic solvation and ion placement
- Optional AMBER to GROMACS format conversion (ParmEd)
- Protein structure analysis tool

## Project Structure

```
src/membrane_builder_mcp/
├── server.py           - MCP server exposing 4 tools
├── membrane_builder.py - packmol-memgen CLI wrapper
└── converter.py        - AMBER to GROMACS conversion (ParmEd)
```

## Prerequisites

- Python 3.10 or later
- Miniconda or Anaconda
- AmberTools25 installed in a conda environment

### Installing AmberTools25

```bash
conda create --name AmberTools25 python=3.12
conda activate AmberTools25
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install dacase::ambertools-dac=25
```

## Installation

```bash
git clone https://github.com/SungminKo-smko/membraneBuilder.git
cd membraneBuilder
pip install -e .
```

## Environment Variables

| Variable | Default | Description |
|---|---|---|
| `CONDA_BASE` | `~/miniconda3` | Path to Miniconda installation |
| `MEMBRANE_CONDA_ENV` | `AmberTools25` | Name of the AmberTools conda environment |

## MCP Configuration

Add the following to your Claude Desktop or Claude Code MCP configuration:

```json
{
  "mcpServers": {
    "membrane-builder": {
      "command": "python",
      "args": ["-m", "membrane_builder_mcp.server"],
      "env": {
        "CONDA_BASE": "/path/to/miniconda3"
      }
    }
  }
}
```

Replace `/path/to/miniconda3` with the absolute path to your Miniconda or Anaconda installation.

## Tools

### build_membrane

The primary tool. Embeds a protein in a lipid bilayer, adds water and ions, and returns AMBER topology and coordinate files. Optionally converts output to GROMACS format.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `protein_pdb_path` | `str` | required | Path to the input protein PDB file |
| `lipids` | `str` | `"POPC"` | Lipid name(s), colon-separated for mixed bilayers (e.g. `"DOPE:DOPG"`) |
| `ratio` | `str` | `"1"` | Molar ratio(s) matching the lipids list (e.g. `"3:1"`) |
| `dist` | `float` | `10.0` | Minimum distance in angstroms from protein to bilayer edge |
| `dist_wat` | `float` | `17.5` | Water layer thickness in angstroms |
| `salt` | `float` | `0.15` | Salt concentration in mol/L |
| `parametrize` | `bool` | `true` | Run tleap to generate AMBER topology files |
| `ffprot` | `str` | `"ff19SB"` | Protein force field |
| `ffwat` | `str` | `"tip3p"` | Water model |
| `fflip` | `str` | `"lipid21"` | Lipid force field |
| `preoriented` | `bool` | `false` | Skip OPM orientation step if protein is already oriented |
| `keep_ligands` | `bool` | `false` | Retain non-standard residues/ligands in the system |
| `convert_to_gromacs` | `bool` | `false` | Convert AMBER output to GROMACS `.gro`/`.top` via ParmEd |

Output files generated:

- `bilayer_*.pdb` - full membrane system PDB
- `*_lipid.prmtop` / `*_lipid.inpcrd` - AMBER topology and coordinates
- `system.gro` / `system.top` - GROMACS format (when `convert_to_gromacs=true`)

### list_available_lipids

Returns the list of lipid types supported by `packmol-memgen` (3000+ entries from the Lipid21 library). Takes no parameters.

### analyze_protein

Analyzes a protein PDB file and reports atom count, residue count, chain identifiers, and bounding box dimensions. Useful for reviewing the structure before building a membrane.

| Parameter | Type | Description |
|---|---|---|
| `protein_pdb_path` | `str` | Path to the protein PDB file |

### convert_to_gromacs

Standalone AMBER to GROMACS conversion using ParmEd. Converts an existing `.prmtop`/`.inpcrd` pair to `.gro`/`.top`.

| Parameter | Type | Description |
|---|---|---|
| `prmtop_path` | `str` | Path to the AMBER topology file |
| `inpcrd_path` | `str` | Path to the AMBER coordinate file |
| `output_dir` | `str` | Directory for GROMACS output files |

## Example Usage

Building a pure POPC bilayer (default settings):

```
build_membrane(protein_pdb_path="/data/protein.pdb")
```

Building a mixed DOPE:DOPG bilayer at 3:1 ratio with GROMACS output:

```
build_membrane(
    protein_pdb_path="/data/protein.pdb",
    lipids="DOPE:DOPG",
    ratio="3:1",
    convert_to_gromacs=true
)
```

Equivalent direct CLI invocation:

```bash
# Pure POPC
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1 --parametrize

# Mixed bilayer
packmol-memgen --pdb protein.pdb --lipids DOPE:DOPG --ratio 3:1 --parametrize
```

## License

MIT
