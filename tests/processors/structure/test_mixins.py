"""Tests for structure processing mixins."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from binding_data_processor.processors.structure.mixins import (
    StructureValidationMixin,
    ConformerGenerationMixin,
    PharmacophoreFeatureMixin,
    StructureAlignmentMixin,
    StructureVisualizationMixin,
)


class TestStructure(
    StructureValidationMixin,
    ConformerGenerationMixin,
    PharmacophoreFeatureMixin,
    StructureAlignmentMixin,
    StructureVisualizationMixin,
):
    """Test class combining all mixins."""

    pass


@pytest.fixture
def test_structure():
    """Create test structure instance."""
    return TestStructure()


@pytest.fixture
def aspirin():
    """Create aspirin molecule with 3D coordinates."""
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


@pytest.fixture
def nap():
    """Create N-acetyl phenylalanine with 3D coordinates."""
    smiles = "CC(=O)N[C@H](Cc1ccccc1)C(=O)O"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def test_structure_validation(test_structure):
    """Test structure validation."""
    # Valid SMILES
    mol = test_structure.validate_structure("CC(=O)O")
    assert mol is not None
    assert mol.GetNumAtoms() > 0
    assert mol.GetNumConformers() == 0  # No 3D coords yet

    # Invalid SMILES
    mol = test_structure.validate_structure("invalid")
    assert mol is None

    # Empty input
    mol = test_structure.validate_structure("")
    assert mol is None


def test_conformer_generation(test_structure):
    """Test conformer generation."""
    mol = Chem.MolFromSmiles("CC(=O)O")
    mol = Chem.AddHs(mol)

    # Generate single conformer
    mol_conf = test_structure.generate_conformers(mol, n_confs=1)
    assert mol_conf is not None
    assert mol_conf.GetNumConformers() == 1

    # Generate multiple conformers
    mol_conf = test_structure.generate_conformers(mol, n_confs=3)
    assert mol_conf is not None
    assert mol_conf.GetNumConformers() == 3

    # Test optimization flag
    mol_conf = test_structure.generate_conformers(mol, n_confs=1, optimize=False)
    assert mol_conf is not None
    assert mol_conf.GetNumConformers() == 1


def test_detailed_feature_detection(test_structure, nap):
    """Test detailed feature detection."""
    features = test_structure.detect_features(nap)

    # Group features by type
    feature_dict = {}
    for feature in features:
        if feature["type"] not in feature_dict:
            feature_dict[feature["type"]] = []
        feature_dict[feature["type"]].append(feature)

    # Check hydrogen bond donors (NH)
    assert "hbd" in feature_dict
    hbd_features = feature_dict["hbd"]
    assert len(hbd_features) >= 1

    # Check hydrogen bond acceptors (C=O)
    assert "hba" in feature_dict
    hba_features = feature_dict["hba"]
    assert len(hba_features) >= 2  # Amide and carboxyl C=O

    # Check aromatic features
    assert "aromatic" in feature_dict
    ar_features = feature_dict["aromatic"]
    assert len(ar_features) >= 1  # Phenyl ring

    # Check feature properties
    for feature in features:
        assert "type" in feature
        assert "atoms" in feature
        assert "position" in feature
        if feature["type"] in ["hbd", "hba"]:
            assert "vector" in feature


def test_feature_based_alignment(test_structure, aspirin, nap):
    """Test feature-based alignment."""
    # Get features
    ref_features = test_structure.detect_features(aspirin)
    probe_features = test_structure.detect_features(nap)

    # Align using features
    rmsd, aligned_mol = test_structure.align_structures(
        aspirin,
        nap,
        ref_features,
        probe_features,
    )
    assert isinstance(rmsd, float)
    assert rmsd >= 0
    assert aligned_mol is not None


def test_visualization_features(test_structure, aspirin):
    """Test visualization capabilities."""
    features = test_structure.detect_features(aspirin)

    # Basic structure visualization
    svg = test_structure.visualize_structure(aspirin)
    assert svg is not None
    assert isinstance(svg, str)
    assert svg.startswith("<?xml")

    # Highlight specific atoms
    highlight_atoms = [0, 1, 2]
    svg = test_structure.visualize_structure(
        aspirin,
        highlight_atoms=highlight_atoms,
    )
    assert svg is not None

    # Feature visualization
    svg = test_structure.visualize_features(aspirin, features)
    assert svg is not None
    assert isinstance(svg, str)
    assert svg.startswith("<?xml")


def test_feature_vectors(test_structure, aspirin):
    """Test feature vector generation."""
    features = test_structure.detect_features(aspirin, include_vectors=True)

    for feature in features:
        if feature["type"] in ["hbd", "hba", "donor", "acceptor"]:
            if "vector" in feature:
                vector = np.array(feature["vector"])
                # Check vector normalization
                norm = np.linalg.norm(vector)
                assert abs(norm - 1.0) < 1e-6


def test_comprehensive_error_handling(test_structure):
    """Test comprehensive error handling."""
    # Structure validation
    assert test_structure.validate_structure(None) is None
    assert test_structure.validate_structure("invalid") is None

    # Conformer generation
    mol = Chem.MolFromSmiles("[Xe]")  # Molecule that can't form conformers
    assert test_structure.generate_conformers(mol) is None

    # Feature detection
    assert test_structure.detect_features(None) == []

    # Structure alignment
    rmsd, mol = test_structure.align_structures(None, None, [], [])
    assert rmsd == float("inf")
    assert mol is None

    # Visualization
    assert test_structure.visualize_structure(None) is None
    assert test_structure.visualize_features(None, []) is None


def test_feature_definitions(test_structure):
    """Test pharmacophore feature definitions."""
    assert hasattr(test_structure, "FEATURE_DEFINITIONS")
    assert isinstance(test_structure.FEATURE_DEFINITIONS, dict)
    assert len(test_structure.FEATURE_DEFINITIONS) > 0

    # Test each SMARTS pattern
    for name, smarts in test_structure.FEATURE_DEFINITIONS.items():
        pattern = Chem.MolFromSmarts(smarts)
        assert pattern is not None, f"Invalid SMARTS pattern for {name}"


def test_transformation_calculation(test_structure):
    """Test coordinate transformation calculations."""
    # Create test coordinates
    coords1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    coords2 = np.array([[1, 1, 1], [2, 1, 1], [1, 2, 1]])

    # Calculate transformation
    rotation, translation = test_structure._get_transformation(coords1, coords2)

    # Check rotation matrix properties
    assert rotation.shape == (3, 3)
    assert np.allclose(rotation.dot(rotation.T), np.eye(3))  # Orthogonal
    assert np.abs(np.linalg.det(rotation) - 1.0) < 1e-6  # Proper rotation

    # Check translation vector
    assert translation.shape == (3,)

    # Test transformation
    transformed = coords1.dot(rotation.T) + translation
    assert np.allclose(transformed, coords2, atol=1e-6)


def test_color_handling(test_structure):
    """Test color handling in visualization."""
    # Test color generation for different feature types
    feature_types = ["hbd", "hba", "aromatic", "hydrophobe"]
    for feature_type in feature_types:
        color = test_structure._get_feature_color(feature_type)
        assert len(color) == 3  # RGB tuple
        assert all(0 <= c <= 1 for c in color)  # Valid color values

    # Test default color for unknown feature type
    color = test_structure._get_feature_color("unknown")
    assert color == (0.5, 0.5, 0.5)  # Default gray
