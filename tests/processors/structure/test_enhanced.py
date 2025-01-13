"""Tests for enhanced structure processing functionality."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from binding_data_processor.processors.structure.pharmacophore.enhanced import EnhancedPharmacophoreGenerator
from binding_data_processor.processors.structure.pharmacophore.features import PharmacophoreFeature


@pytest.fixture
def generator():
    """Create enhanced pharmacophore generator instance."""
    return EnhancedPharmacophoreGenerator()


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
def ibuprofen():
    """Create ibuprofen molecule with 3D coordinates."""
    smiles = "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"
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


def test_generate_from_smiles(generator):
    """Test generating features from SMILES."""
    # Valid SMILES
    features = generator.generate_from_smiles("CC(=O)N[C@H](Cc1ccccc1)C(=O)O")
    assert features
    assert len(features) > 0
    assert all(isinstance(f, PharmacophoreFeature) for f in features)

    # Check for expected feature types
    feature_types = {f.feature_type for f in features}
    assert "hbd" in feature_types  # NH group
    assert "hba" in feature_types  # C=O groups
    assert "aromatic" in feature_types  # Phenyl ring

    # Invalid SMILES
    features = generator.generate_from_smiles("invalid")
    assert features == []


def test_detailed_feature_detection(generator, nap):
    """Test detailed feature detection."""
    features = generator.detect_features(nap)

    # Group features by type
    feature_dict = {}
    for feature in features:
        if feature.feature_type not in feature_dict:
            feature_dict[feature.feature_type] = []
        feature_dict[feature.feature_type].append(feature)

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
        assert hasattr(feature, "atoms")
        assert hasattr(feature, "position")
        if feature.feature_type in ["hbd", "hba"]:
            assert hasattr(feature, "vector")


def test_pharmacophore_alignment(generator, aspirin, ibuprofen):
    """Test pharmacophore-based alignment."""
    rmsd, ref_features, aligned_features = generator.align_pharmacophores(
        aspirin,
        ibuprofen,
    )

    # Basic checks
    assert isinstance(rmsd, float)
    assert rmsd >= 0
    assert ref_features
    assert aligned_features
    assert len(ref_features) > 0
    assert len(aligned_features) > 0

    # Compare pharmacophores before and after alignment
    metrics_before = generator.compare_pharmacophores(
        generator.generate(aspirin),
        generator.generate(ibuprofen),
    )
    metrics_after = generator.compare_pharmacophores(ref_features, aligned_features)

    # Alignment should improve similarity
    assert metrics_after["similarity"] >= metrics_before["similarity"]


def test_feature_comparison(generator, aspirin, ibuprofen, nap):
    """Test pharmacophore comparison functionality."""
    # Generate features
    features1 = generator.generate(aspirin)
    features2 = generator.generate(ibuprofen)
    features3 = generator.generate(nap)

    # Compare similar molecules
    metrics = generator.compare_pharmacophores(features1, features2)
    assert isinstance(metrics, dict)
    assert "overlap" in metrics
    assert "coverage" in metrics
    assert "similarity" in metrics
    assert all(0 <= v <= 1 for v in metrics.values())

    # Compare with different distance cutoffs
    metrics_strict = generator.compare_pharmacophores(features1, features2, distance_cutoff=1.0)
    metrics_loose = generator.compare_pharmacophores(features1, features2, distance_cutoff=3.0)
    assert metrics_loose["similarity"] >= metrics_strict["similarity"]

    # Compare different molecule types
    metrics_diff = generator.compare_pharmacophores(features1, features3)
    assert metrics_diff["similarity"] < metrics["similarity"]  # Should be less similar


def test_visualization_features(generator, aspirin):
    """Test pharmacophore visualization capabilities."""
    features = generator.generate(aspirin)

    # Basic visualization
    svg1 = generator.visualize_pharmacophore(aspirin, features, show_vectors=False)
    assert svg1 is not None
    assert isinstance(svg1, str)
    assert svg1.startswith("<?xml")

    # Visualization with vectors
    svg2 = generator.visualize_pharmacophore(aspirin, features, show_vectors=True)
    assert svg2 is not None
    assert isinstance(svg2, str)
    assert svg2.startswith("<?xml")

    # Vectors should add more content to SVG
    assert len(svg2) > len(svg1)


def test_color_schemes(generator, nap):
    """Test color scheme handling."""
    features = generator.generate(nap)

    # Test default color scheme
    vis = generator.visualize_pharmacophore(nap, features)
    assert vis is not None

    # Ensure each feature type has a color
    feature_types = {f.feature_type for f in features}
    for feature_type in feature_types:
        color = generator._get_feature_color(feature_type)
        assert len(color) == 3  # RGB tuple
        assert all(0 <= c <= 1 for c in color)  # Valid color values


def test_comprehensive_error_handling(generator):
    """Test comprehensive error handling."""
    # Structure validation
    features = generator.generate_from_smiles("invalid")
    assert features == []

    # Conformer generation
    mol = Chem.MolFromSmiles("[Xe]")  # Molecule that can't form conformers
    features = generator.generate(mol)
    assert features == []

    # Feature detection with invalid molecule
    features = generator.detect_features(None)
    assert features == []

    # Alignment with invalid molecules
    rmsd, ref_features, aligned_features = generator.align_pharmacophores(None, None)
    assert rmsd == float("inf")
    assert ref_features == []
    assert aligned_features == []

    # Comparison with empty features
    metrics = generator.compare_pharmacophores([], [])
    assert metrics["overlap"] == 0
    assert metrics["coverage"] == 0
    assert metrics["similarity"] == 0

    # Visualization with invalid input
    svg = generator.visualize_pharmacophore(None, [])
    assert svg is None


def test_feature_vectors(generator, nap):
    """Test generation and properties of feature vectors."""
    features = generator.detect_features(nap, include_vectors=True)

    for feature in features:
        if feature.feature_type in ["hbd", "hba", "donor", "acceptor"]:
            if hasattr(feature, "vector") and feature.vector is not None:
                vector = np.array(feature.vector)
                # Check vector normalization
                norm = np.linalg.norm(vector)
                assert abs(norm - 1.0) < 1e-6


def test_conformer_handling(generator):
    """Test handling of molecules with/without conformers."""
    # Molecule without conformers
    mol = Chem.MolFromSmiles("CC(=O)O")
    features = generator.generate(mol)
    assert features == []

    # Generate conformers
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    features = generator.generate(mol)
    assert features
    assert len(features) > 0

    # Multiple conformers
    mol_confs = generator.generate_conformers(mol, n_confs=3)
    assert mol_confs is not None
    assert mol_confs.GetNumConformers() == 3
