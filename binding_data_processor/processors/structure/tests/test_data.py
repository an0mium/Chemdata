"""Test data and utilities for structure processor tests.

Provides:
1. Test molecules and SMARTS patterns
2. ML model test data and utilities 
3. Test file paths and utilities
4. Mock data generation
5. Test helper functions
6. Common test fixtures
7. Pharmacological activity prediction
8. Toxicity prediction
9. Abuse potential assessment
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Draw,
    Descriptors,
    rdDecomposition,
    rdMolDescriptors,
    rdMolTransforms,
    rdMolAlign,
    rdFMCS,
)
from rdkit.Chem.Draw import rdDepictor
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
import torch
from torch_geometric.data import Data
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# Configure logging
logger = logging.getLogger(__name__)

# Test molecules as SMILES
TEST_MOLECULES = {
    # Basic molecules
    "benzene": "c1ccccc1",
    "phenol": "Oc1ccccc1",
    "aniline": "Nc1ccccc1",
    "toluene": "Cc1ccccc1",
    "naphthalene": "c1ccc2ccccc2c1",
    # Classical psychedelics
    "lsd": "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3c[nH]c4cccc(C2=C1)c34",
    "psilocybin": "CN1C[C@@H](O[P](=O)(O)O)C=C2c3cccc4[nH]cc(c34)C[C@H]21",
    "dmt": "CN(C)CCc1c[nH]c2ccccc12",
    "mescaline": "COc1cc(CC(N)C)cc(OC)c1OC",
    # Empathogens
    "mdma": "CC(NC)Cc1ccc2c(c1)OCO2",
    "mda": "CC(N)Cc1ccc2c(c1)OCO2",
    "methylone": "CNC(C)C(=O)c1ccc2c(c1)OCO2",
    # Dissociatives
    "ketamine": "CN1C(=O)CCC(C1=O)(c2ccccc2)N(C)C",
    "pcp": "c1ccc(cc1)C2(CCCC2)N3CCCCC3",
    "dxm": "COc1ccc2C[C@H]3[C@H]4C[C@H](C[C@@H](C4)C3)N2c1",
    # Stimulants
    "amphetamine": "CC(N)Cc1ccccc1",
    "cocaine": "COC(=O)[C@H]1[C@@H]2CC[C@H](C[C@@H]1OC(=O)c3ccccc3)N2C",
    "methylphenidate": "COC(=O)C(C1CCCCN1)c2ccccc2",
    # Antidepressants
    "fluoxetine": "CNCCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2",
    "sertraline": "CNC1CCC(c2ccc(Cl)c(Cl)c2)c3ccccc13",
    "paroxetine": "Fc1ccc(cc1)[C@@H]2CCNC[C@H]2OC3=C(O)C=CC(=C3)F",
    # Anxiolytics
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13",
    "alprazolam": "Cc1nnc2CN=C(c3ccccc3)c4cc(Cl)ccc4-n12",
    "buspirone": "O=C1CC2(CCCC2)CC(=O)N1CCCCN3CCN(CC3)c4ncccn4",
}

# Pharmacological classes
PHARM_CLASSES = {
    "psychedelics": ["lsd", "psilocybin", "dmt", "mescaline"],
    "empathogens": ["mdma", "mda", "methylone"],
    "dissociatives": ["ketamine", "pcp", "dxm"],
    "stimulants": ["amphetamine", "cocaine", "methylphenidate"],
    "antidepressants": ["fluoxetine", "sertraline", "paroxetine"],
    "anxiolytics": ["diazepam", "alprazolam", "buspirone"],
}

# Receptor binding profiles
BINDING_PROFILES = {
    "5ht2a_agonists": ["lsd", "psilocybin", "dmt", "mescaline"],
    "sert_inhibitors": ["fluoxetine", "sertraline", "paroxetine", "mdma"],
    "nmda_antagonists": ["ketamine", "pcp", "dxm"],
    "dopamine_active": ["amphetamine", "cocaine", "methylphenidate"],
    "gaba_active": ["diazepam", "alprazolam"],
}

# Known toxicity risks
TOXICITY_RISKS = {
    "cardiotoxicity": ["cocaine", "mdma", "methylone"],
    "hepatotoxicity": ["mdma", "cocaine", "methylphenidate"],
    "neurotoxicity": ["mdma", "methamphetamine", "pcp"],
    "addiction": ["cocaine", "amphetamine", "diazepam", "alprazolam"],
}

# SMARTS patterns
SMARTS_PATTERNS = {
    # Core structures
    "phenethylamine": "c1ccccc1CCN",
    "tryptamine": "c1ccc2c(c1)c(CCN)[nH]2",
    "ergoline": "CN1C[C@@H](C=C2c3cccc4[nH]cc(c34)C[C@H]21)C(=O)N",
    # Functional groups
    "basic_amine": "[NX3;H2,H1,H0][CX4;!$(C(=[O,S,N])[N,O,S])]",
    "phenol": "[OH]c1ccccc1",
    "catechol": "[OH]c1c([OH])cccc1",
    "indole": "[nH]1ccc2ccccc12",
    # Metabolic sites
    "cyp_oxidation": "[CH2,CH3][NX3]",
    "glucuronidation": "[OH,NH2]",
    "acetylation": "[NH2]",
    # Toxicophores
    "quinone_former": "Oc1c([OH,NH2])cccc1",
    "nitro_reduction": "[N+](=O)[O-]",
    "epoxide_former": "C=C[CH2,CH3]",
}


class TestMolecules:
    """Test molecule provider with ML features."""

    @staticmethod
    def get_test_mols() -> Dict[str, Chem.Mol]:
        """Get test molecules."""
        return {
            name: Chem.MolFromSmiles(smiles) for name, smiles in TEST_MOLECULES.items()
        }

    @staticmethod
    def get_class_mols(pharm_class: str) -> Dict[str, Chem.Mol]:
        """Get molecules by pharmacological class."""
        if pharm_class not in PHARM_CLASSES:
            return {}
        return {
            name: Chem.MolFromSmiles(TEST_MOLECULES[name])
            for name in PHARM_CLASSES[pharm_class]
        }

    @staticmethod
    def get_binding_profile_mols(profile: str) -> Dict[str, Chem.Mol]:
        """Get molecules by binding profile."""
        if profile not in BINDING_PROFILES:
            return {}
        return {
            name: Chem.MolFromSmiles(TEST_MOLECULES[name])
            for name in BINDING_PROFILES[profile]
        }

    @staticmethod
    def get_toxicity_risk_mols(risk: str) -> Dict[str, Chem.Mol]:
        """Get molecules by toxicity risk."""
        if risk not in TOXICITY_RISKS:
            return {}
        return {
            name: Chem.MolFromSmiles(TEST_MOLECULES[name])
            for name in TOXICITY_RISKS[risk]
        }

    @staticmethod
    def get_mol_descriptors(mol: Chem.Mol) -> Dict[str, float]:
        """Calculate molecular descriptors."""
        desc_calc = MoleculeDescriptors.MolecularDescriptorCalculator(
            [x[0] for x in Descriptors._descList]
        )
        descriptors = desc_calc.CalcDescriptors(mol)
        return dict(zip(desc_calc.GetDescriptorNames(), descriptors))

    @staticmethod
    def get_graph_features(mol: Chem.Mol) -> Data:
        """Convert molecule to graph features for GNN."""
        # Node features
        atom_features = []
        for atom in mol.GetAtoms():
            features = [
                atom.GetAtomicNum(),
                atom.GetTotalDegree(),
                atom.GetFormalCharge(),
                atom.GetTotalNumHs(),
                atom.GetIsAromatic(),
                atom.GetMass(),
                atom.GetExplicitValence(),
                atom.GetImplicitValence(),
                int(atom.IsInRing()),
                int(atom.IsInRingSize(5)),
                int(atom.IsInRingSize(6)),
            ]
            atom_features.append(features)

        # Edge features
        edge_indices = []
        edge_features = []
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            edge_indices += [[i, j], [j, i]]
            features = [
                bond.GetBondTypeAsDouble(),
                bond.GetIsAromatic(),
                bond.IsInRing(),
                bond.IsInRingSize(5),
                bond.IsInRingSize(6),
                float(bond.GetIsConjugated()),
                float(bond.GetIsRotatable()),
            ]
            edge_features += [features, features]

        return Data(
            x=torch.tensor(atom_features, dtype=torch.float),
            edge_index=torch.tensor(edge_indices, dtype=torch.long).t().contiguous(),
            edge_attr=torch.tensor(edge_features, dtype=torch.float),
        )

    @staticmethod
    def get_fingerprints(
        mol: Chem.Mol,
        radius: int = 2,
        nBits: int = 2048,
        use_features: bool = True,
    ) -> np.ndarray:
        """Generate molecular fingerprints."""
        return np.array(
            AllChem.GetMorganFingerprintAsBitVect(
                mol, radius, nBits=nBits, useFeatures=use_features
            )
        )

    @staticmethod
    def get_scaffolds(mol: Chem.Mol) -> Dict[str, Chem.Mol]:
        """Get molecular scaffolds."""
        return {
            "murcko": MurckoScaffold.GetScaffoldForMol(mol),
            "generic_murcko": MurckoScaffold.MakeScaffoldGeneric(
                MurckoScaffold.GetScaffoldForMol(mol)
            ),
        }

    @staticmethod
    def get_fragments(
        mol: Chem.Mol,
        pattern: str,
    ) -> List[Chem.Mol]:
        """Get molecular fragments matching SMARTS pattern."""
        pattern_mol = Chem.MolFromSmarts(pattern)
        if not pattern_mol:
            return []

        matches = mol.GetSubstructMatches(pattern_mol)
        fragments = []

        for match in matches:
            fragment = Chem.RWMol(mol)
            atoms_to_keep = set(match)
            for atom in reversed(range(fragment.GetNumAtoms())):
                if atom not in atoms_to_keep:
                    fragment.RemoveAtom(atom)
            fragments.append(fragment.GetMol())

        return fragments


class TestSmarts:
    """SMARTS pattern provider."""

    @staticmethod
    def get_compiled_patterns() -> Dict[str, Chem.Mol]:
        """Get compiled SMARTS patterns."""
        return {
            name: Chem.MolFromSmarts(pattern)
            for name, pattern in SMARTS_PATTERNS.items()
        }

    @staticmethod
    def match_pattern(
        mol: Chem.Mol,
        pattern_name: str,
        unique: bool = True,
    ) -> List[Tuple[int, ...]]:
        """Match SMARTS pattern against molecule."""
        pattern = Chem.MolFromSmarts(SMARTS_PATTERNS[pattern_name])
        if pattern is None:
            return []
        if unique:
            return mol.GetSubstructMatches(pattern, uniquify=True)
        return mol.GetSubstructMatches(pattern)

    @staticmethod
    def get_matching_atoms(
        mol: Chem.Mol,
        pattern_name: str,
    ) -> Set[int]:
        """Get atoms matching SMARTS pattern."""
        matches = TestSmarts.match_pattern(mol, pattern_name)
        return set([atom for match in matches for atom in match])

    @staticmethod
    def highlight_matches(
        mol: Chem.Mol,
        pattern_name: str,
    ) -> Tuple[List[int], List[int]]:
        """Get atoms and bonds to highlight for pattern."""
        pattern = Chem.MolFromSmarts(SMARTS_PATTERNS[pattern_name])
        if pattern is None:
            return [], []

        matches = mol.GetSubstructMatches(pattern)
        highlight_atoms = set()
        highlight_bonds = set()

        for match in matches:
            highlight_atoms.update(match)
            for bond in pattern.GetBonds():
                aid1 = match[bond.GetBeginAtomIdx()]
                aid2 = match[bond.GetEndAtomIdx()]
                highlight_bonds.add(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

        return list(highlight_atoms), list(highlight_bonds)


class TestFiles:
    """Test file management."""

    TEST_DIR = Path(__file__).parent
    TEST_DATA_DIR = TEST_DIR / "data"
    TEST_OUTPUT_DIR = TEST_DIR / "output"
    TEST_MODEL_DIR = TEST_DIR / "models"

    @classmethod
    def ensure_dirs(cls):
        """Ensure test directories exist."""
        for directory in [cls.TEST_DATA_DIR, cls.TEST_OUTPUT_DIR, cls.TEST_MODEL_DIR]:
            directory.mkdir(parents=True, exist_ok=True)

    @classmethod
    def cleanup_test_files(cls):
        """Clean up test files."""
        import shutil

        for directory in [cls.TEST_OUTPUT_DIR]:
            if directory.exists():
                shutil.rmtree(directory)
                directory.mkdir(parents=True)

    @classmethod
    def save_test_data(cls, data: Dict, filename: str):
        """Save test data to JSON."""
        path = cls.TEST_DATA_DIR / filename
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load_test_data(cls, filename: str) -> Dict:
        """Load test data from JSON."""
        path = cls.TEST_DATA_DIR / filename
        with open(path) as f:
            return json.load(f)


class TestUtils:
    """Test utility functions."""

    @staticmethod
    def compare_molecules(
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        compare_coords: bool = False,
    ) -> bool:
        """Compare two molecules for equality."""
        if mol1 is None or mol2 is None:
            return False

        # Compare basic properties
        if (
            mol1.GetNumAtoms() != mol2.GetNumAtoms()
            or mol1.GetNumBonds() != mol2.GetNumBonds()
        ):
            return False

        # Compare SMILES
        smiles1 = Chem.MolToSmiles(mol1, canonical=True)
        smiles2 = Chem.MolToSmiles(mol2, canonical=True)
        if smiles1 != smiles2:
            return False

        # Compare 3D coordinates if requested
        if compare_coords:
            if (
                mol1.GetNumConformers() != mol2.GetNumConformers()
                or mol1.GetNumConformers() == 0
            ):
                return False

            conf1 = mol1.GetConformer()
            conf2 = mol2.GetConformer()
            positions1 = conf1.GetPositions()
            positions2 = conf2.GetPositions()
            return np.allclose(positions1, positions2)

        return True

    @staticmethod
    def generate_3d_conformer(
        mol: Chem.Mol,
        random_seed: int = 42,
        num_conformers: int = 1,
        optimize: bool = True,
    ) -> Optional[Chem.Mol]:
        """Generate 3D conformer(s) for molecule."""
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_conformers,
                randomSeed=random_seed,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True,
                enforceChirality=True,
            )

            if optimize:
                for conf_id in range(mol.GetNumConformers()):
                    if AllChem.MMFFHasAllMoleculeParams(mol):
                        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
                    else:
                        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)

            return mol

        except Exception as e:
            logger.error(f"Error generating 3D conformer: {e}")
            return None

    @staticmethod
    def depict_molecule(
        mol: Chem.Mol,
        size: Tuple[int, int] = (300, 300),
        output_path: Optional[Path] = None,
        highlight_atoms: Optional[List[int]] = None,
        highlight_bonds: Optional[List[int]] = None,
        legend: Optional[str] = None,
    ) -> Optional[np.ndarray]:
        """Generate 2D depiction of molecule."""
        try:
            rdDepictor.Compute2DCoords(mol)

            drawer = Draw.rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            drawer.drawOptions().addStereoAnnotation = True
            drawer.drawOptions().addAtomIndices = True

            if highlight_atoms or highlight_bonds:
                colors = [(0.7, 0.7, 0, 0.3)]  # Yellow highlight
                drawer.DrawMolecule(
                    mol,
                    highlightAtoms=highlight_atoms,
                    highlightBonds=highlight_bonds,
                    highlightColor=colors[0],
                    legend=legend,
                )
            else:
                drawer.DrawMolecule(mol, legend=legend)

            drawer.FinishDrawing()
            img = drawer.GetDrawingText()

            if output_path is not None:
                with open(output_path, "wb") as f:
                    f.write(img)

            return np.frombuffer(img, dtype=np.uint8)

        except Exception as e:
            logger.error(f"Error depicting molecule: {e}")
            return None


# Initialize test environment
TestFiles.ensure_dirs()

# Example usage
if __name__ == "__main__":
    # Get test molecules by class
    psychedelics = TestMolecules.get_class_mols("psychedelics")
    sert_inhibitors = TestMolecules.get_binding_profile_mols("sert_inhibitors")
    cardiotoxic = TestMolecules.get_toxicity_risk_mols("cardiotoxicity")

    # Analyze LSD
    lsd = psychedelics["lsd"]

    # Get descriptors and graph features
    descriptors = TestMolecules.get_mol_descriptors(lsd)
    graph = TestMolecules.get_graph_features(lsd)
    fingerprint = TestMolecules.get_fingerprints(lsd, use_features=True)

    # Get scaffolds and fragments
    scaffolds = TestMolecules.get_scaffolds(lsd)
    ergoline_fragments = TestMolecules.get_fragments(lsd, SMARTS_PATTERNS["ergoline"])

    # Generate 3D conformer
    lsd_3d = TestUtils.generate_3d_conformer(lsd, num_conformers=3, optimize=True)

    # Highlight pharmacophores
    basic_amine_atoms, basic_amine_bonds = TestSmarts.highlight_matches(
        lsd, "basic_amine"
    )

    # Generate depiction
    if lsd_3d:
        img = TestUtils.depict_molecule(
            lsd_3d,
            output_path=TestFiles.TEST_OUTPUT_DIR / "lsd.png",
            highlight_atoms=basic_amine_atoms,
            highlight_bonds=basic_amine_bonds,
            legend="LSD with highlighted basic amine",
        )

    print("Test data module loaded successfully")
