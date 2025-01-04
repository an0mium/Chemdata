"""Tests for mock data generator.

Tests:
1. Random molecule generation
2. Mock descriptor generation
3. Mock pharmacophore features
4. Mock similarity calculations
5. Mock activity data
"""

import unittest
import numpy as np
from rdkit import Chem

from .mock_generator import MockGenerator


class TestMockGenerator(unittest.TestCase):
    """Test mock data generator functionality."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.mock_gen = MockGenerator(seed=42)
        cls.test_mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene

    def test_random_mol_generation(self):
        """Test random molecule generation."""
        # Test single molecule
        mol = self.mock_gen.get_random_mol()
        self.assertIsNotNone(mol)
        self.assertIsInstance(mol, Chem.Mol)

        # Test multiple molecules
        mols = self.mock_gen.get_random_mols(n=3)
        self.assertEqual(len(mols), 3)
        self.assertTrue(all(isinstance(m, Chem.Mol) for m in mols))

        # Test with drug molecules
        drug_mol = self.mock_gen.get_random_mol(include_drugs=True)
        self.assertIsNotNone(drug_mol)

    def test_mock_descriptors(self):
        """Test mock descriptor generation."""
        # Test random descriptors
        desc = self.mock_gen.generate_mock_descriptors()
        self.assertIsInstance(desc, dict)
        self.assertTrue(all(isinstance(v, (int, float)) for v in desc.values()))

        # Test with actual molecule
        mol_desc = self.mock_gen.generate_mock_descriptors(self.test_mol)
        self.assertIsInstance(mol_desc, dict)
        self.assertGreater(mol_desc["MolWt"], 0)

    def test_mock_conformer(self):
        """Test 3D conformer generation."""
        # Test single conformer
        mol_3d = self.mock_gen.generate_mock_conformer(self.test_mol)
        self.assertIsNotNone(mol_3d)
        self.assertEqual(mol_3d.GetNumConformers(), 1)

        # Test multiple conformers
        mol_multi = self.mock_gen.generate_mock_conformer(self.test_mol, n_conf=3)
        self.assertEqual(mol_multi.GetNumConformers(), 3)

    def test_mock_pharmacophore(self):
        """Test pharmacophore feature generation."""
        # Test random features
        features = self.mock_gen.generate_mock_pharmacophore()
        self.assertIsInstance(features, dict)
        self.assertTrue(all(isinstance(v, list) for v in features.values()))

        # Test with actual molecule
        mol_features = self.mock_gen.generate_mock_pharmacophore(self.test_mol)
        self.assertIsInstance(mol_features, dict)
        self.assertTrue("aromatic" in mol_features)

    def test_mock_similarity(self):
        """Test similarity matrix generation."""
        # Test matrix shape and properties
        n_mols = 5
        matrix = self.mock_gen.generate_mock_similarity_matrix(n_mols)
        self.assertEqual(matrix.shape, (n_mols, n_mols))
        self.assertTrue(np.allclose(matrix, matrix.T))  # Symmetric
        self.assertTrue(np.all(np.diag(matrix) == 1.0))  # Self-similarity

        # Test similarity range
        sim_range = (0.2, 0.8)
        matrix = self.mock_gen.generate_mock_similarity_matrix(
            n_mols, similarity_range=sim_range
        )
        non_diag = matrix[~np.eye(n_mols, dtype=bool)]
        self.assertTrue(np.all(non_diag >= sim_range[0]))
        self.assertTrue(np.all(non_diag <= sim_range[1]))

    def test_mock_fingerprint(self):
        """Test fingerprint generation."""
        # Test random fingerprint
        fp = self.mock_gen.generate_mock_fingerprint()
        self.assertEqual(len(fp), 1024)
        self.assertTrue(np.all((fp == 0) | (fp == 1)))

        # Test with actual molecule
        mol_fp = self.mock_gen.generate_mock_fingerprint(self.test_mol)
        self.assertEqual(len(mol_fp), 1024)
        self.assertTrue(np.any(mol_fp == 1))  # Should have some bits set

    def test_mock_structure_match(self):
        """Test substructure match generation."""
        # Test random matches
        matches = self.mock_gen.generate_mock_structure_match()
        self.assertIsInstance(matches, dict)
        self.assertTrue(all(isinstance(v, list) for v in matches.values()))

        # Test with actual molecule
        mol_matches = self.mock_gen.generate_mock_structure_match(self.test_mol)
        self.assertIsInstance(mol_matches, dict)
        self.assertTrue(any(len(v) > 0 for v in mol_matches.values()))

    def test_mock_activity_data(self):
        """Test activity data generation."""
        # Test single compound
        activities = self.mock_gen.generate_mock_activity_data(n_compounds=1)
        self.assertEqual(len(activities), 1)
        activity = activities[0]
        self.assertIsInstance(activity["molecule"], Chem.Mol)
        self.assertIsInstance(activity["activity_value"], float)
        self.assertTrue(0.5 <= activity["confidence"] <= 1.0)

        # Test multiple compounds
        activities = self.mock_gen.generate_mock_activity_data(n_compounds=5)
        self.assertEqual(len(activities), 5)

    def test_mock_sar_data(self):
        """Test SAR data generation."""
        sar_data = self.mock_gen.generate_mock_sar_data(n_compounds=10)
        self.assertIn("active", sar_data)
        self.assertIn("inactive", sar_data)
        total_compounds = len(sar_data["active"]) + len(sar_data["inactive"])
        self.assertLessEqual(total_compounds, 10)

        # Check compound features
        for group in sar_data.values():
            for comp in group:
                self.assertIn("features", comp)
                self.assertIsInstance(comp["features"], dict)


if __name__ == "__main__":
    unittest.main()
