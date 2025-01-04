"""Unit tests for structure depiction classes.

Tests cover:
1. Basic 2D and 3D depiction functionality
2. Grid layouts and conformer generation
3. Substructure highlighting
4. Custom coloring and labeling
5. Error handling and edge cases
"""

import os
import unittest
from pathlib import Path
from typing import List, Optional
import numpy as np
from PIL import Image
from rdkit import Chem

from ..depiction_2d import Structure2DDepiction
from ..depiction_3d import Structure3DDepiction


class TestStructureDepiction(unittest.TestCase):
    """Test cases for structure depiction classes."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        # Create test molecules
        cls.smiles = [
            "CC1=CC=CC=C1",  # Toluene
            "CC(=O)O",  # Acetic acid
            "c1ccccc1O",  # Phenol
        ]
        cls.mols = [Chem.MolFromSmiles(s) for s in cls.smiles]

        # Create depiction objects
        cls.depict_2d = Structure2DDepiction()
        cls.depict_3d = Structure3DDepiction()

        # Create test output directory
        cls.test_output = Path("tests/output")
        cls.test_output.mkdir(parents=True, exist_ok=True)

    def test_2d_depiction_basic(self):
        """Test basic 2D depiction functionality."""
        for i, mol in enumerate(self.mols):
            # Test PNG output
            img_data = self.depict_2d.depict_2d(mol)
            self.assertIsNotNone(img_data)
            self.assertIsInstance(img_data, bytes)

            # Test PIL output
            img = self.depict_2d.depict_2d(mol, return_pil=True)
            self.assertIsNotNone(img)
            self.assertIsInstance(img, Image.Image)

            # Save test output
            if img:
                img.save(self.test_output / f"test_2d_{i}.png")

    def test_2d_grid_depiction(self):
        """Test 2D grid depiction."""
        # Test with legends
        legends = ["Toluene", "Acetic acid", "Phenol"]
        grid = self.depict_2d.depict_grid(
            self.mols,
            legends=legends,
            return_pil=True,
        )
        self.assertIsNotNone(grid)
        self.assertIsInstance(grid, Image.Image)

        if grid:
            grid.save(self.test_output / "test_grid.png")

        # Test without legends
        grid = self.depict_2d.depict_grid(self.mols, return_pil=True)
        self.assertIsNotNone(grid)

    def test_2d_highlighting(self):
        """Test substructure highlighting."""
        # Test phenol substructure
        query = Chem.MolFromSmarts("[OH]c1ccccc1")
        self.assertIsNotNone(query)

        if query:
            img = self.depict_2d.depict_with_highlights(
                self.mols[2],  # Phenol
                query,
                return_pil=True,
            )
            self.assertIsNotNone(img)
            self.assertIsInstance(img, Image.Image)

            if img:
                img.save(self.test_output / "test_highlight.png")

    def test_3d_depiction_basic(self):
        """Test basic 3D depiction functionality."""
        for i, mol in enumerate(self.mols):
            # Generate 3D conformer
            mol_3d = self.depict_3d.generate_3d_conformer(mol)
            self.assertIsNotNone(mol_3d)

            if mol_3d:
                # Test PNG output
                img_data = self.depict_3d.depict_3d(mol_3d)
                self.assertIsNotNone(img_data)
                self.assertIsInstance(img_data, bytes)

                # Test PIL output
                img = self.depict_3d.depict_3d(mol_3d, return_pil=True)
                self.assertIsNotNone(img)
                self.assertIsInstance(img, Image.Image)

                if img:
                    img.save(self.test_output / f"test_3d_{i}.png")

    def test_3d_conformer_generation(self):
        """Test 3D conformer generation."""
        mol = self.mols[0]  # Toluene

        # Generate multiple conformers
        n_conf = 5
        mol_3d = self.depict_3d.generate_conformers(mol, n_conf=n_conf)
        self.assertIsNotNone(mol_3d)

        if mol_3d:
            # Check number of conformers
            self.assertEqual(mol_3d.GetNumConformers(), n_conf)

            # Test conformer grid
            grid = self.depict_3d.depict_conformer_grid(
                mol_3d,
                mols_per_row=3,
                return_pil=True,
            )
            self.assertIsNotNone(grid)

            if grid:
                grid.save(self.test_output / "test_conformers.png")

    def test_custom_drawing_options(self):
        """Test custom drawing options."""
        # Test 2D options
        custom_2d = Structure2DDepiction(
            {
                "size": (800, 800),
                "show_atom_numbers": True,
                "annotate_stereo": True,
            }
        )
        img = custom_2d.depict_2d(self.mols[0], return_pil=True)
        self.assertIsNotNone(img)

        if img:
            img.save(self.test_output / "test_custom_2d.png")

        # Test 3D options
        custom_3d = Structure3DDepiction(
            {
                "size": (800, 800),
                "add_hydrogens": True,
                "optimize_geometry": True,
            }
        )
        mol_3d = custom_3d.generate_3d_conformer(self.mols[0])
        if mol_3d:
            img = custom_3d.depict_3d(mol_3d, return_pil=True)
            self.assertIsNotNone(img)

            if img:
                img.save(self.test_output / "test_custom_3d.png")

    def test_error_handling(self):
        """Test error handling."""
        # Test invalid molecule
        invalid_mol = Chem.MolFromSmiles("invalid")
        self.assertIsNone(invalid_mol)

        # Test empty molecule list
        grid = self.depict_2d.depict_grid([], return_pil=True)
        self.assertIsNone(grid)

        # Test invalid drawing options
        with self.assertRaises(Exception):
            Structure2DDepiction({"invalid_option": True})

        # Test invalid substructure query
        invalid_query = Chem.MolFromSmarts("invalid")
        self.assertIsNone(invalid_query)

        if invalid_query:
            img = self.depict_2d.depict_with_highlights(
                self.mols[0],
                invalid_query,
                return_pil=True,
            )
            self.assertIsNone(img)

    def test_pharmacophore_depiction(self):
        """Test pharmacophore feature depiction."""
        mol = self.mols[2]  # Phenol
        mol_3d = self.depict_3d.generate_3d_conformer(mol)
        self.assertIsNotNone(mol_3d)

        if mol_3d:
            # Define pharmacophore features
            features = {
                "donor": [6],  # OH group
                "aromatic": [0, 1, 2, 3, 4, 5],  # Benzene ring
            }

            img = self.depict_3d.depict_pharmacophore(
                mol_3d,
                features,
                return_pil=True,
            )
            self.assertIsNotNone(img)

            if img:
                img.save(self.test_output / "test_pharmacophore.png")

    def test_custom_coloring(self):
        """Test custom atom and bond coloring."""
        mol = self.mols[0]  # Toluene

        # Define custom colors
        atom_colors = {
            0: (1, 0, 0),  # Red
            1: (0, 1, 0),  # Green
        }
        bond_colors = {
            0: (0, 0, 1),  # Blue
        }

        img = self.depict_2d.depict_2d(
            mol,
            atom_colors=atom_colors,
            bond_colors=bond_colors,
            return_pil=True,
        )
        self.assertIsNotNone(img)

        if img:
            img.save(self.test_output / "test_colors.png")

    @classmethod
    def tearDownClass(cls):
        """Clean up test files."""
        # Optionally remove test output directory
        # import shutil
        # shutil.rmtree(cls.test_output)
        pass


if __name__ == "__main__":
    unittest.main()
