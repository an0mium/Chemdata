"""Tests for surface analysis functionality."""

import pytest
import numpy as np
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

from ..base import BaseSurfaceAnalyzer
from ..protein import ProteinSurfaceAnalyzer
from ..binding import BindingSurfaceAnalyzer


def create_test_structure():
    """Create a simple test structure."""
    structure = Structure.Structure("test")
    model = Model.Model(0)
    chain = Chain.Chain("A")

    # Create a few residues
    for i, resname in enumerate(["ALA", "GLY", "SER"]):
        res = Residue.Residue((" ", i, " "), resname, "")
        # Add some atoms
        for name, coord in [
            ("N", [i * 4, 0, 0]),
            ("CA", [i * 4 + 1, 0, 0]),
            ("C", [i * 4 + 2, 0, 0]),
            ("O", [i * 4 + 2, 1, 0]),
        ]:
            atom = Atom.Atom(
                name,
                coord,
                20.0,  # B-factor
                1.0,  # Occupancy
                " ",  # Altloc
                name,  # Fullname
                None,  # Serial number
                element=name[0],
            )
            res.add(atom)
        chain.add(res)

    model.add(chain)
    structure.add(model)
    return structure


class TestBaseSurfaceAnalyzer:
    """Test base surface analysis functionality."""

    def setup_method(self):
        self.analyzer = BaseSurfaceAnalyzer()
        self.structure = create_test_structure()

    def test_get_surface_atoms(self):
        atoms = self.analyzer.get_surface_atoms(self.structure)
        assert len(atoms) > 0
        assert all(isinstance(a, Atom.Atom) for a in atoms)

    def test_calculate_surface_area(self):
        atoms = self.analyzer.get_surface_atoms(self.structure)
        area = self.analyzer.calculate_surface_area(atoms)
        assert area > 0

    def test_calculate_volume(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        volume = self.analyzer.calculate_volume(coords)
        assert volume > 0


class TestProteinSurfaceAnalyzer:
    """Test protein surface analysis functionality."""

    def setup_method(self):
        self.analyzer = ProteinSurfaceAnalyzer()
        self.structure = create_test_structure()

    def test_analyze_surface(self):
        result = self.analyzer.analyze_surface(self.structure)
        assert isinstance(result, dict)
        assert "surface_area" in result
        assert "surface_roughness" in result

    def test_analyze_site(self):
        result = self.analyzer.analyze_site(self.structure, [0, 1])
        assert isinstance(result, dict)
        assert "surface_area" in result
        assert "exposure" in result


class TestBindingSurfaceAnalyzer:
    """Test binding site analysis functionality."""

    def setup_method(self):
        self.analyzer = BindingSurfaceAnalyzer()
        self.structure = create_test_structure()

    def test_analyze_binding_site(self):
        result = self.analyzer.analyze_binding_site(self.structure, [0, 1])
        assert isinstance(result, dict)
        assert "volume" in result
        assert "depth" in result
        assert "sub_pockets" in result

    def test_find_sub_pockets(self):
        atoms = self.analyzer.get_surface_atoms(self.structure)
        pockets = self.analyzer._find_sub_pockets(atoms)
        assert isinstance(pockets, list)
        for pocket in pockets:
            assert "volume" in pocket
            assert "center" in pocket
            assert "residues" in pocket
