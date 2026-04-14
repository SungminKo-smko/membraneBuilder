# membrane-builder-mcp

An MCP (Model Context Protocol) server that builds lipid bilayer membranes around proteins using [Packmol](https://m3g.github.io/packmol/). Provide a protein PDB file and the server automatically calculates membrane dimensions, generates Packmol input, and produces a solvated membrane system.

## Features

- Automatic membrane sizing based on protein bounding box
- Configurable lipid composition with support for mixed bilayers
- Built-in lipid library: POPC, POPE, DPPC, and cholesterol
- Protein structure analysis with suggested membrane parameters
- Water box generation with adjustable padding and density
- Runs as a stdio-based MCP server compatible with Claude Desktop and Claude Code

## Prerequisites

- Python 3.10 or later
- [Packmol](https://m3g.github.io/packmol/download.shtml) installed and available on your `PATH`
- [uv](https://docs.astral.sh/uv/) package manager

## Installation

```bash
git clone https://github.com/<your-org>/membrane-builder-mcp.git
cd membrane-builder-mcp
uv sync
```

To include development dependencies (pytest):

```bash
uv sync --extra dev
```

## Packmol Installation

Packmol must be installed separately. Download and build instructions are available at:

https://m3g.github.io/packmol/download.shtml

Verify your installation:

```bash
packmol < /dev/null
```

## Usage

### Running the MCP Server

```bash
python -m membrane_builder_mcp.server
```

The server communicates over stdio using the MCP protocol.

### MCP Configuration

Add the following to your Claude Desktop or Claude Code MCP configuration:

```json
{
  "mcpServers": {
    "membrane-builder": {
      "command": "uv",
      "args": ["--directory", "/path/to/membrane-builder-mcp", "run", "python", "-m", "membrane_builder_mcp.server"]
    }
  }
}
```

Replace `/path/to/membrane-builder-mcp` with the absolute path to this repository.

## Tools

### build_membrane

Build a lipid bilayer membrane around a protein using Packmol.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `protein_pdb_path` | `str` | (required) | Path to the protein PDB file |
| `output_dir` | `str` | `"./output"` | Output directory for generated files |
| `lipid_composition` | `dict[str, float]` | `{"POPC": 1.0}` | Lipid composition ratios (values must sum to 1.0) |
| `membrane_x_size` | `float` | auto | Membrane X dimension in angstroms |
| `membrane_y_size` | `float` | auto | Membrane Y dimension in angstroms |
| `bilayer_thickness` | `float` | `40.0` | Bilayer thickness in angstroms |
| `water_padding` | `float` | `15.0` | Water layer padding in angstroms |
| `lipid_area` | `float` | `65.0` | Area per lipid in angstroms squared |
| `water_density` | `float` | `0.0334` | Water density in molecules per cubic angstrom |
| `tolerance` | `float` | `2.0` | Packmol distance tolerance in angstroms |

### list_available_lipids

List all available lipid types with their properties (description, headgroup, atom count, and area per lipid). Takes no parameters.

### analyze_protein_structure

Analyze a protein PDB file and suggest membrane building parameters.

| Parameter | Type | Description |
|---|---|---|
| `protein_pdb_path` | `str` | Path to the protein PDB file |

## Example Usage Scenarios

**Simple membrane with defaults** -- provide only the protein PDB file and the server auto-calculates dimensions and uses a pure POPC bilayer:

```
build_membrane(protein_pdb_path="/data/1abc.pdb")
```

**Mixed lipid composition** -- build a membrane with 70% POPC and 30% POPE:

```
build_membrane(
    protein_pdb_path="/data/1abc.pdb",
    lipid_composition={"POPC": 0.7, "POPE": 0.3}
)
```

**Analyze first, then build** -- inspect the protein structure to review suggested parameters before building:

```
analyze_protein_structure(protein_pdb_path="/data/1abc.pdb")
```

**Custom dimensions** -- specify explicit membrane size and increased water padding:

```
build_membrane(
    protein_pdb_path="/data/1abc.pdb",
    membrane_x_size=120.0,
    membrane_y_size=120.0,
    water_padding=20.0
)
```

## Running Tests

```bash
uv run pytest tests/ -v
```

## License

MIT
