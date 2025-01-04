"""Base test class for structure processor tests.

Provides:
1. Common test fixtures
2. Shared test utilities
3. Helper methods for structure testing
4. Mock data access
5. Test cleanup
"""

import unittest
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .mock_generator import MockGenerator
from .test_data import TestMolecules, TestSmarts, TestFiles, TestUtils


class BaseStructureTest(unittest.TestCase):
    """Base class for structure processor tests."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        # Initialize mock generator with fixed seed
        cls.mock_gen = MockGenerator(seed=42)

        # Get test molecules
        cls.test_mols = TestMolecules.get_test_mols()
        cls.drug_mols = TestMolecules.get_drug_mols()
        cls.feature_mols = TestMolecules.get_feature_mols()

        # Get test SMARTS patterns
        cls.patterns = TestSmarts.get_compiled_patterns()

        # Set up test directories
        cls.test_dir = TestFiles.TEST_DIR
        cls.data_dir = TestFiles.TEST_DATA_DIR
        cls.output_dir = TestFiles.TEST_OUTPUT_DIR

        # Create test molecule
        cls.benzene = Chem.MolFromSmiles("c1ccccc1")
        cls.benzene_3d = TestUtils.generate_3d_conformer(cls.benzene)

    @classmethod
    def tearDownClass(cls):
        """Clean up test files."""
        TestFiles.cleanup_test_files()

    def assertMolEqual(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        compare_coords: bool = False,
        msg: Optional[str] = None,
    ):
        """Assert that two molecules are equal."""
        self.assertTrue(
            TestUtils.compare_molecules(mol1, mol2, compare_coords),
            msg or "Molecules are not equal",
        )

    def assertMolValid(self, mol: Chem.Mol, msg: Optional[str] = None):
        """Assert that molecule is valid."""
        self.assertIsNotNone(mol, msg or "Molecule is None")
        self.assertIsInstance(mol, Chem.Mol, msg or "Not a RDKit molecule")
        self.assertGreater(mol.GetNumAtoms(), 0, msg or "Molecule has no atoms")

    def assertDescriptorsValid(
        self, descriptors: Dict[str, float], msg: Optional[str] = None
    ):
        """Assert that molecular descriptors are valid."""
        self.assertIsInstance(descriptors, dict, msg or "Descriptors not a dictionary")
        self.assertTrue(
            all(isinstance(v, (int, float)) for v in descriptors.values()),
            msg or "Invalid descriptor values",
        )
        self.assertGreater(len(descriptors), 0, msg or "No descriptors present")

    def assertFingerprintValid(
        self, fp: np.ndarray, n_bits: int = 1024, msg: Optional[str] = None
    ):
        """Assert that fingerprint is valid."""
        self.assertEqual(len(fp), n_bits, msg or f"Fingerprint not {n_bits} bits")
        self.assertTrue(
            np.all((fp == 0) | (fp == 1)),
            msg or "Invalid fingerprint values",
        )

    def assertSimilarityMatrixValid(
        self,
        matrix: np.ndarray,
        n_mols: int,
        msg: Optional[str] = None,
    ):
        """Assert that similarity matrix is valid."""
        self.assertEqual(
            matrix.shape,
            (n_mols, n_mols),
            msg or "Invalid matrix shape",
        )
        self.assertTrue(
            np.allclose(matrix, matrix.T),
            msg or "Matrix not symmetric",
        )
        self.assertTrue(
            np.all(np.diag(matrix) == 1.0),
            msg or "Invalid self-similarities",
        )
        self.assertTrue(
            np.all((matrix >= 0) & (matrix <= 1)),
            msg or "Similarities not in [0,1]",
        )

    def assertPharmacophoreValid(
        self,
        features: Dict[str, List[int]],
        mol: Optional[Chem.Mol] = None,
        msg: Optional[str] = None,
    ):
        """Assert that pharmacophore features are valid."""
        self.assertIsInstance(features, dict, msg or "Features not a dictionary")
        self.assertTrue(
            all(isinstance(v, list) for v in features.values()),
            msg or "Feature values not lists",
        )
        if mol is not None:
            n_atoms = mol.GetNumAtoms()
            for indices in features.values():
                self.assertTrue(
                    all(0 <= i < n_atoms for i in indices),
                    msg or "Invalid atom indices",
                )

    def assertStructureMatchValid(
        self,
        matches: Dict[str, List[Tuple[int, ...]]],
        mol: Optional[Chem.Mol] = None,
        msg: Optional[str] = None,
    ):
        """Assert that substructure matches are valid."""
        self.assertIsInstance(matches, dict, msg or "Matches not a dictionary")
        self.assertTrue(
            all(isinstance(v, list) for v in matches.values()),
            msg or "Match values not lists",
        )
        if mol is not None:
            n_atoms = mol.GetNumAtoms()
            for match_list in matches.values():
                for match in match_list:
                    self.assertTrue(
                        all(0 <= i < n_atoms for i in match),
                        msg or "Invalid atom indices in match",
                    )

    def get_test_file_path(self, filename: str) -> Path:
        """Get path to test data file."""
        return self.data_dir / filename

    def get_output_path(self, filename: str) -> Path:
        """Get path for test output file."""
        return self.output_dir / filename

    def create_test_mol(
        self,
        smiles: str,
        gen_3d: bool = False,
        n_conf: int = 1,
    ) -> Optional[Chem.Mol]:
        """Create test molecule from SMILES."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None and gen_3d:
            try:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, randomSeed=42)
                for conf_id in range(n_conf):
                    AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
            except Exception:
                return None
        return mol


# Example usage
if __name__ == "__main__":
    # Create test class
    class ExampleTest(BaseStructureTest):
        """Example test using base class."""

        def test_molecule_creation(self):
            """Test molecule creation."""
            mol = self.create_test_mol("c1ccccc1", gen_3d=True)
            self.assertMolValid(mol)
            self.assertEqual(mol.GetNumConformers(), 1)

        def test_mock_data(self):
            """Test mock data generation."""
            # Get mock descriptors
            desc = self.mock_gen.generate_mock_descriptors(self.benzene)
            self.assertDescriptorsValid(desc)

            # Get mock fingerprint
            fp = self.mock_gen.generate_mock_fingerprint(self.benzene)
            self.assertFingerprintValid(fp)

            # Get mock pharmacophore
            features = self.mock_gen.generate_mock_pharmacophore(self.benzene)
            self.assertPharmacophoreValid(features, self.benzene)

    # Run tests
    unittest.main()
