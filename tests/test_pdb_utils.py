"""Tests for the pdb_utils module."""

import pytest

from membrane_builder_mcp.pdb_utils import (
    estimate_cross_section,
    parse_pdb,
    suggest_membrane_params,
)


class TestParsePdb:
    """Tests for parse_pdb function."""

    def test_returns_expected_keys(self, sample_pdb_path):
        result = parse_pdb(str(sample_pdb_path))
        assert "atoms" in result
        assert "center" in result
        assert "bbox" in result

    def test_atoms_is_list_of_dicts(self, sample_pdb_path):
        result = parse_pdb(str(sample_pdb_path))
        atoms = result["atoms"]
        assert isinstance(atoms, list)
        assert len(atoms) > 0
        for atom in atoms:
            assert isinstance(atom, dict)
            for key in ("name", "resname", "chain", "resnum", "x", "y", "z"):
                assert key in atom

    def test_atom_count_matches_fixture(self, sample_pdb_path):
        result = parse_pdb(str(sample_pdb_path))
        # The fixture has 30 ATOM records
        assert len(result["atoms"]) == 30

    def test_center_is_float_tuple(self, sample_pdb_path):
        result = parse_pdb(str(sample_pdb_path))
        center = result["center"]
        assert isinstance(center, tuple)
        assert len(center) == 3
        for val in center:
            assert isinstance(val, float)

    def test_bbox_has_min_max_size(self, sample_pdb_path):
        result = parse_pdb(str(sample_pdb_path))
        bbox = result["bbox"]
        assert "min" in bbox
        assert "max" in bbox
        assert "size" in bbox
        # min should be less than or equal to max in each dimension
        for i in range(3):
            assert bbox["min"][i] <= bbox["max"][i]
            assert bbox["size"][i] >= 0

    def test_residue_names_from_fixture(self, sample_pdb_path):
        result = parse_pdb(str(sample_pdb_path))
        resnames = {a["resname"] for a in result["atoms"]}
        assert "ALA" in resnames
        assert "GLY" in resnames
        assert "VAL" in resnames


class TestEstimateCrossSection:
    """Tests for estimate_cross_section function."""

    def test_returns_non_negative_float(self, sample_pdb_path):
        pdb_data = parse_pdb(str(sample_pdb_path))
        atoms = pdb_data["atoms"]
        z_center = pdb_data["center"][2]
        area = estimate_cross_section(atoms, z_center=z_center, thickness=10.0)
        assert isinstance(area, float)
        assert area >= 0.0

    def test_returns_positive_area_at_center(self, sample_pdb_path):
        pdb_data = parse_pdb(str(sample_pdb_path))
        atoms = pdb_data["atoms"]
        z_center = pdb_data["center"][2]
        area = estimate_cross_section(atoms, z_center=z_center, thickness=40.0)
        # With a thick enough slice, should capture atoms
        assert area > 0.0

    def test_returns_zero_for_z_far_from_atoms(self, sample_pdb_path):
        pdb_data = parse_pdb(str(sample_pdb_path))
        atoms = pdb_data["atoms"]
        # Place the z_center far away from all atoms
        area = estimate_cross_section(atoms, z_center=9999.0, thickness=1.0)
        assert area == 0.0

    def test_returns_zero_for_empty_atoms(self):
        area = estimate_cross_section([], z_center=0.0, thickness=10.0)
        assert area == 0.0


class TestSuggestMembraneParams:
    """Tests for suggest_membrane_params function."""

    def test_returns_expected_keys(self, sample_pdb_path):
        result = suggest_membrane_params(str(sample_pdb_path))
        expected_keys = {
            "membrane_x_size",
            "membrane_y_size",
            "n_lipids_per_leaflet",
            "protein_center",
            "protein_bbox",
            "cross_section_area",
        }
        assert set(result.keys()) == expected_keys

    def test_membrane_sizes_are_positive(self, sample_pdb_path):
        result = suggest_membrane_params(str(sample_pdb_path))
        assert result["membrane_x_size"] > 0
        assert result["membrane_y_size"] > 0

    def test_n_lipids_per_leaflet_is_non_negative(self, sample_pdb_path):
        result = suggest_membrane_params(str(sample_pdb_path))
        assert result["n_lipids_per_leaflet"] >= 0

    def test_custom_xy_padding(self, sample_pdb_path):
        result_small = suggest_membrane_params(str(sample_pdb_path), xy_padding=5.0)
        result_large = suggest_membrane_params(str(sample_pdb_path), xy_padding=30.0)
        assert result_large["membrane_x_size"] > result_small["membrane_x_size"]
        assert result_large["membrane_y_size"] > result_small["membrane_y_size"]
