"""Tests for the membrane_builder module."""

import pytest
import pytest_asyncio

from membrane_builder_mcp.membrane_builder import (
    build_membrane,
    calculate_lipid_counts,
    generate_packmol_input,
    get_lipid_pdb_path,
)


class TestGetLipidPdbPath:
    """Tests for get_lipid_pdb_path function."""

    def test_returns_existing_path_for_valid_lipid(self):
        path = get_lipid_pdb_path("POPC", "upper")
        assert path.exists()
        assert path.suffix == ".pdb"

    def test_upper_and_lower_leaflets(self):
        for leaflet in ("upper", "lower"):
            path = get_lipid_pdb_path("POPC", leaflet)
            assert path.exists()

    def test_all_known_lipids(self):
        known_lipids = ["POPC", "POPE", "DPPC", "cholesterol"]
        for lipid in known_lipids:
            for leaflet in ("upper", "lower"):
                path = get_lipid_pdb_path(lipid, leaflet)
                assert path.exists(), f"Missing PDB for {lipid}_{leaflet}"

    def test_raises_for_unknown_lipid(self):
        with pytest.raises(ValueError, match="Unknown lipid"):
            get_lipid_pdb_path("FAKE_LIPID", "upper")

    def test_raises_for_invalid_leaflet(self):
        with pytest.raises(ValueError, match="Leaflet must be"):
            get_lipid_pdb_path("POPC", "middle")


class TestCalculateLipidCounts:
    """Tests for calculate_lipid_counts function."""

    def test_single_lipid_composition(self):
        counts = calculate_lipid_counts(
            membrane_x=100.0,
            membrane_y=100.0,
            cross_section=0.0,
            lipid_composition={"POPC": 1.0},
            lipid_area=65.0,
        )
        assert "POPC" in counts
        assert counts["POPC"] > 0
        # 100*100 / 65 = ~153
        expected = int(10000.0 / 65.0)
        assert counts["POPC"] == expected

    def test_mixed_composition(self):
        counts = calculate_lipid_counts(
            membrane_x=100.0,
            membrane_y=100.0,
            cross_section=0.0,
            lipid_composition={"POPC": 0.7, "POPE": 0.3},
            lipid_area=65.0,
        )
        assert "POPC" in counts
        assert "POPE" in counts
        assert counts["POPC"] > counts["POPE"]
        # Check proportions are roughly correct
        total = counts["POPC"] + counts["POPE"]
        assert total > 0
        popc_frac = counts["POPC"] / total
        assert 0.6 < popc_frac < 0.8

    def test_cross_section_reduces_count(self):
        counts_no_cs = calculate_lipid_counts(
            membrane_x=100.0,
            membrane_y=100.0,
            cross_section=0.0,
            lipid_composition={"POPC": 1.0},
            lipid_area=65.0,
        )
        counts_with_cs = calculate_lipid_counts(
            membrane_x=100.0,
            membrane_y=100.0,
            cross_section=2000.0,
            lipid_composition={"POPC": 1.0},
            lipid_area=65.0,
        )
        assert counts_with_cs["POPC"] < counts_no_cs["POPC"]

    def test_returns_zero_when_protein_fills_membrane(self):
        counts = calculate_lipid_counts(
            membrane_x=10.0,
            membrane_y=10.0,
            cross_section=200.0,  # exceeds 10*10=100
            lipid_composition={"POPC": 1.0},
            lipid_area=65.0,
        )
        assert counts["POPC"] == 0


class TestGeneratePackmolInput:
    """Tests for generate_packmol_input function."""

    def test_returns_string(self, sample_pdb_path):
        content = generate_packmol_input(
            protein_pdb_path=str(sample_pdb_path),
            output_path="/tmp/test_output.pdb",
            lipid_counts={"POPC": 100},
            membrane_x=80.0,
            membrane_y=80.0,
            z_center=0.0,
            bilayer_thickness=40.0,
            water_padding_upper=15.0,
            water_padding_lower=15.0,
            water_density=0.0334,
            tolerance=2.0,
            n_water_upper=500,
            n_water_lower=500,
        )
        assert isinstance(content, str)
        assert len(content) > 0

    def test_contains_expected_keywords(self, sample_pdb_path):
        content = generate_packmol_input(
            protein_pdb_path=str(sample_pdb_path),
            output_path="/tmp/test_output.pdb",
            lipid_counts={"POPC": 100},
            membrane_x=80.0,
            membrane_y=80.0,
            z_center=0.0,
            bilayer_thickness=40.0,
            water_padding_upper=15.0,
            water_padding_lower=15.0,
            water_density=0.0334,
            tolerance=2.0,
            n_water_upper=500,
            n_water_lower=500,
        )
        assert "tolerance" in content
        assert "filetype pdb" in content
        assert "output" in content
        assert "structure" in content
        assert "end structure" in content

    def test_contains_protein_structure_block(self, sample_pdb_path):
        content = generate_packmol_input(
            protein_pdb_path=str(sample_pdb_path),
            output_path="/tmp/test_output.pdb",
            lipid_counts={"POPC": 50},
            membrane_x=80.0,
            membrane_y=80.0,
            z_center=0.0,
            bilayer_thickness=40.0,
            water_padding_upper=15.0,
            water_padding_lower=15.0,
            water_density=0.0334,
            tolerance=2.0,
            n_water_upper=500,
            n_water_lower=500,
        )
        assert str(sample_pdb_path) in content
        assert "fixed" in content
        assert "center" in content

    def test_contains_water_blocks(self, sample_pdb_path):
        content = generate_packmol_input(
            protein_pdb_path=str(sample_pdb_path),
            output_path="/tmp/test_output.pdb",
            lipid_counts={"POPC": 50},
            membrane_x=80.0,
            membrane_y=80.0,
            z_center=0.0,
            bilayer_thickness=40.0,
            water_padding_upper=15.0,
            water_padding_lower=15.0,
            water_density=0.0334,
            tolerance=2.0,
            n_water_upper=500,
            n_water_lower=500,
        )
        assert "Water upper" in content
        assert "Water lower" in content
        assert "water.pdb" in content

    def test_no_water_blocks_when_zero(self, sample_pdb_path):
        content = generate_packmol_input(
            protein_pdb_path=str(sample_pdb_path),
            output_path="/tmp/test_output.pdb",
            lipid_counts={"POPC": 50},
            membrane_x=80.0,
            membrane_y=80.0,
            z_center=0.0,
            bilayer_thickness=40.0,
            water_padding_upper=15.0,
            water_padding_lower=15.0,
            water_density=0.0334,
            tolerance=2.0,
            n_water_upper=0,
            n_water_lower=0,
        )
        assert "Water upper" not in content
        assert "Water lower" not in content

    def test_multiple_lipid_types(self, sample_pdb_path):
        content = generate_packmol_input(
            protein_pdb_path=str(sample_pdb_path),
            output_path="/tmp/test_output.pdb",
            lipid_counts={"POPC": 70, "POPE": 30},
            membrane_x=80.0,
            membrane_y=80.0,
            z_center=0.0,
            bilayer_thickness=40.0,
            water_padding_upper=15.0,
            water_padding_lower=15.0,
            water_density=0.0334,
            tolerance=2.0,
            n_water_upper=500,
            n_water_lower=500,
        )
        assert "POPC upper leaflet" in content
        assert "POPC lower leaflet" in content
        assert "POPE upper leaflet" in content
        assert "POPE lower leaflet" in content


@pytest.mark.asyncio
class TestBuildMembrane:
    """Tests for build_membrane async function."""

    async def test_runs_with_sample_protein(self, sample_pdb_path, tmp_path):
        result = await build_membrane(
            protein_pdb_path=str(sample_pdb_path),
            output_dir=str(tmp_path / "output"),
        )
        assert isinstance(result, dict)

    async def test_creates_output_directory_with_inp_file(self, sample_pdb_path, tmp_path):
        output_dir = tmp_path / "output"
        result = await build_membrane(
            protein_pdb_path=str(sample_pdb_path),
            output_dir=str(output_dir),
        )
        assert output_dir.exists()
        inp_files = list(output_dir.glob("*.inp"))
        assert len(inp_files) == 1

    async def test_returns_expected_keys(self, sample_pdb_path, tmp_path):
        result = await build_membrane(
            protein_pdb_path=str(sample_pdb_path),
            output_dir=str(tmp_path / "output"),
        )
        expected_keys = {
            "output_path",
            "input_path",
            "n_lipids_upper",
            "n_lipids_lower",
            "n_water_upper",
            "n_water_lower",
            "box_dimensions",
            "protein_center",
            "cross_section_area",
            "packmol_success",
            "packmol_log",
        }
        assert set(result.keys()) == expected_keys

    async def test_lipid_counts_are_dicts(self, sample_pdb_path, tmp_path):
        result = await build_membrane(
            protein_pdb_path=str(sample_pdb_path),
            output_dir=str(tmp_path / "output"),
        )
        assert isinstance(result["n_lipids_upper"], dict)
        assert isinstance(result["n_lipids_lower"], dict)

    async def test_water_counts_are_non_negative(self, sample_pdb_path, tmp_path):
        result = await build_membrane(
            protein_pdb_path=str(sample_pdb_path),
            output_dir=str(tmp_path / "output"),
        )
        assert result["n_water_upper"] >= 0
        assert result["n_water_lower"] >= 0

    async def test_custom_lipid_composition(self, sample_pdb_path, tmp_path):
        result = await build_membrane(
            protein_pdb_path=str(sample_pdb_path),
            output_dir=str(tmp_path / "output"),
            lipid_composition={"POPC": 0.7, "POPE": 0.3},
        )
        assert "POPC" in result["n_lipids_upper"]
        assert "POPE" in result["n_lipids_upper"]

    async def test_inp_file_content_is_valid(self, sample_pdb_path, tmp_path):
        output_dir = tmp_path / "output"
        result = await build_membrane(
            protein_pdb_path=str(sample_pdb_path),
            output_dir=str(output_dir),
        )
        inp_path = result["input_path"]
        with open(inp_path, "r") as f:
            content = f.read()
        assert "tolerance" in content
        assert "filetype pdb" in content
        assert "end structure" in content
