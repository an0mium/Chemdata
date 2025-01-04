"""Mock data generator for structure processor tests.

Provides:
1. Random molecule generation
2. Mock descriptor data
3. Mock pharmacophore features
4. Mock similarity scores
5. Mock structure analysis results
"""

import random
from typing import Any, Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from .test_data import TestMolecules, TestSmarts


class MockGenerator:
    """Generate mock data for testing."""

    def __init__(self, seed: Optional[int] = None):
        """Initialize mock generator with optional seed."""
        self.rng = random.Random(seed)
        np.random.seed(seed)
        self.test_mols = TestMolecules.get_test_mols()
        self.drug_mols = TestMolecules.get_drug_mols()

    def get_random_mol(self, include_drugs: bool = False) -> Optional[Chem.Mol]:
        """Get random test molecule."""
        mols = self.test_mols + (self.drug_mols if include_drugs else [])
        return self.rng.choice(mols) if mols else None

    def get_random_mols(
        self, n: int = 5, include_drugs: bool = False
    ) -> List[Chem.Mol]:
        """Get list of random molecules."""
        mols = self.test_mols + (self.drug_mols if include_drugs else [])
        return self.rng.sample(mols, min(n, len(mols)))

    def generate_mock_descriptors(
        self, mol: Optional[Chem.Mol] = None
    ) -> Dict[str, float]:
        """Generate mock molecular descriptors."""
        base = {
            "MolWt": self.rng.uniform(100, 500),
            "LogP": self.rng.uniform(-2, 5),
            "TPSA": self.rng.uniform(20, 140),
            "HBA": self.rng.randint(0, 10),
            "HBD": self.rng.randint(0, 5),
            "RotBonds": self.rng.randint(0, 10),
            "AromaticRings": self.rng.randint(0, 4),
        }

        if mol is not None:
            # Override with actual values where possible
            try:
                base["MolWt"] = Descriptors.ExactMolWt(mol)
                base["LogP"] = Descriptors.MolLogP(mol)
                base["TPSA"] = Descriptors.TPSA(mol)
                base["HBA"] = Descriptors.NumHAcceptors(mol)
                base["HBD"] = Descriptors.NumHDonors(mol)
                base["RotBonds"] = Descriptors.NumRotatableBonds(mol)
                base["AromaticRings"] = Descriptors.NumAromaticRings(mol)
            except Exception:
                pass

        return base

    def generate_mock_conformer(
        self, mol: Chem.Mol, n_conf: int = 1
    ) -> Optional[Chem.Mol]:
        """Generate mock 3D conformer(s)."""
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMultipleConfs(
                mol,
                numConfs=n_conf,
                randomSeed=self.rng.randint(1, 1000),
            )
            return mol
        except Exception:
            return None

    def generate_mock_pharmacophore(
        self, mol: Optional[Chem.Mol] = None
    ) -> Dict[str, List[int]]:
        """Generate mock pharmacophore features."""
        if mol is None:
            n_atoms = self.rng.randint(5, 20)
            features = {
                "donor": self.rng.sample(range(n_atoms), self.rng.randint(0, 3)),
                "acceptor": self.rng.sample(range(n_atoms), self.rng.randint(0, 3)),
                "aromatic": self.rng.sample(range(n_atoms), self.rng.randint(0, 6)),
                "hydrophobe": self.rng.sample(range(n_atoms), self.rng.randint(0, 4)),
            }
        else:
            # Use actual SMARTS matches where possible
            patterns = TestSmarts.get_compiled_patterns()["pharmacophore"]
            features = {}
            for name, pat in patterns.items():
                if pat is not None:
                    matches = mol.GetSubstructMatches(pat)
                    features[name] = sorted(
                        set(idx for match in matches for idx in match)
                    )

        return features

    def generate_mock_similarity_matrix(
        self, n_mols: int, similarity_range: Tuple[float, float] = (0.1, 0.9)
    ) -> np.ndarray:
        """Generate mock similarity matrix."""
        min_sim, max_sim = similarity_range
        matrix = np.zeros((n_mols, n_mols))

        # Fill upper triangle with random similarities
        for i in range(n_mols):
            for j in range(i + 1, n_mols):
                sim = self.rng.uniform(min_sim, max_sim)
                matrix[i, j] = matrix[j, i] = sim
            matrix[i, i] = 1.0  # Self-similarity

        return matrix

    def generate_mock_fingerprint(
        self, mol: Optional[Chem.Mol] = None, n_bits: int = 1024
    ) -> np.ndarray:
        """Generate mock molecular fingerprint."""
        if mol is not None:
            try:
                # Use actual Morgan fingerprint if possible
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
                return np.array(list(fp.ToBitString()), dtype=int)
            except Exception:
                pass

        # Generate random fingerprint
        return np.random.binomial(1, 0.3, n_bits)

    def generate_mock_structure_match(
        self, mol: Optional[Chem.Mol] = None
    ) -> Dict[str, List[Tuple[int, ...]]]:
        """Generate mock substructure matches."""
        if mol is None:
            n_atoms = self.rng.randint(5, 20)
            matches = {
                "ring": [
                    tuple(sorted(self.rng.sample(range(n_atoms), 6)))
                    for _ in range(self.rng.randint(0, 3))
                ],
                "chain": [
                    tuple(sorted(self.rng.sample(range(n_atoms), 3)))
                    for _ in range(self.rng.randint(0, 4))
                ],
            }
        else:
            # Use actual SMARTS matches where possible
            patterns = TestSmarts.get_compiled_patterns()["ring_systems"]
            matches = {}
            for name, pat in patterns.items():
                if pat is not None:
                    matches[name] = list(mol.GetSubstructMatches(pat))

        return matches

    def generate_mock_descriptor_matrix(
        self, mols: List[Chem.Mol]
    ) -> Tuple[np.ndarray, List[str]]:
        """Generate mock descriptor matrix for multiple molecules."""
        descriptors = []
        for mol in mols:
            desc = self.generate_mock_descriptors(mol)
            descriptors.append(list(desc.values()))

        return np.array(descriptors), list(desc.keys())

    def generate_mock_activity_data(
        self, n_compounds: int = 10
    ) -> List[Dict[str, Any]]:
        """Generate mock activity data."""
        activities = []
        for _ in range(n_compounds):
            mol = self.get_random_mol(include_drugs=True)
            if mol is None:
                continue

            activity = {
                "molecule": mol,
                "smiles": Chem.MolToSmiles(mol),
                "activity_type": self.rng.choice(
                    ["agonist", "antagonist", "inhibitor"]
                ),
                "activity_value": self.rng.lognormvariate(0, 1),
                "confidence": self.rng.uniform(0.5, 1.0),
            }
            activities.append(activity)

        return activities

    def generate_mock_sar_data(
        self, n_compounds: int = 10
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Generate mock structure-activity relationship data."""
        compounds = self.generate_mock_activity_data(n_compounds)

        # Group by activity
        sar_data = {
            "active": [],
            "inactive": [],
        }

        for comp in compounds:
            # Classify as active/inactive
            is_active = comp["activity_value"] < 1.0
            group = "active" if is_active else "inactive"

            # Add structural features
            features = self.generate_mock_pharmacophore(comp["molecule"])
            comp["features"] = features

            sar_data[group].append(comp)

        return sar_data


# Example usage
if __name__ == "__main__":
    # Initialize generator with seed
    mock_gen = MockGenerator(seed=42)

    # Generate mock data
    mol = mock_gen.get_random_mol(include_drugs=True)
    if mol is not None:
        # Generate various mock data
        descriptors = mock_gen.generate_mock_descriptors(mol)
        pharmacophore = mock_gen.generate_mock_pharmacophore(mol)
        fingerprint = mock_gen.generate_mock_fingerprint(mol)
        matches = mock_gen.generate_mock_structure_match(mol)

        # Generate mock data for multiple molecules
        mols = mock_gen.get_random_mols(5, include_drugs=True)
        sim_matrix = mock_gen.generate_mock_similarity_matrix(len(mols))
        desc_matrix, desc_names = mock_gen.generate_mock_descriptor_matrix(mols)

        # Generate activity data
        activities = mock_gen.generate_mock_activity_data()
        sar_data = mock_gen.generate_mock_sar_data()

        print("Mock data generated successfully")
