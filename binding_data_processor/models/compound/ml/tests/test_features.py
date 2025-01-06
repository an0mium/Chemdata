"""Tests for chemical feature extraction."""

import pytest
from pathlib import Path
from typing import List, Dict, Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from ....core import CompoundData
from ..features import (
    FeatureExtractor,
    MorganFeatures,
    MACCSFeatures,
    RDKitFeatures,
    PharmacophoreFeatures,
)


@pytest.fixture
def test_compounds() -> List[CompoundData]:
    """Create test compounds."""
    compounds = []
    for i in range(5):
        compound = CompoundData(
            name=f"Test Compound {i}",
            smiles=f"CC{i}CC",  # Simple alkanes
            cas_number=f"123-45-{i}",
        )
        compounds.append(compound)
    return compounds


@pytest.fixture
def feature_extractor(tmp_path: Path) -> FeatureExtractor:
    """Create feature extractor."""
    return FeatureExtractor(
        cache_dir=str(tmp_path / "cache"),
        feature_types=["morgan", "maccs", "rdkit", "pharmacophore"],
    )


def test_morgan_features():
    """Test Morgan fingerprint extraction."""
    extractor = MorganFeatures()
    mol = Chem.MolFromSmiles("CCO")  # Ethanol
    features = extractor.extract(mol)

    assert isinstance(features, np.ndarray)
    assert features.shape == (2048,)  # Default Morgan length
    assert features.dtype == np.int8
    assert features.sum() > 0  # Should have some bits set


def test_maccs_features():
    """Test MACCS key extraction."""
    extractor = MACCSFeatures()
    mol = Chem.MolFromSmiles("CCO")  # Ethanol
    features = extractor.extract(mol)

    assert isinstance(features, np.ndarray)
    assert features.shape == (167,)  # MACCS keys length
    assert features.dtype == np.int8
    assert features.sum() > 0  # Should have some bits set


def test_rdkit_features():
    """Test RDKit descriptor extraction."""
    extractor = RDKitFeatures()
    mol = Chem.MolFromSmiles("CCO")  # Ethanol
    features = extractor.extract(mol)

    assert isinstance(features, np.ndarray)
    assert len(features) > 0
    assert features.dtype == np.float32
    assert not np.any(np.isnan(features))  # No NaN values


def test_pharmacophore_features():
    """Test pharmacophore feature extraction."""
    extractor = PharmacophoreFeatures()
    mol = Chem.MolFromSmiles("CCO")  # Ethanol
    features = extractor.extract(mol)

    assert isinstance(features, np.ndarray)
    assert len(features) > 0
    assert features.dtype == np.float32
    assert not np.any(np.isnan(features))  # No NaN values


def test_feature_extraction(
    feature_extractor: FeatureExtractor,
    test_compounds: List[CompoundData],
):
    """Test full feature extraction pipeline."""
    # Extract features
    features = feature_extractor.extract_features(test_compounds[0])

    # Check feature structure
    assert isinstance(features, Dict)
    assert all(k in features for k in ["morgan", "maccs", "rdkit", "pharmacophore"])
    assert all(isinstance(v, np.ndarray) for v in features.values())
    assert all(len(v) > 0 for v in features.values())
    assert all(not np.any(np.isnan(v)) for v in features.values())


def test_feature_caching(
    feature_extractor: FeatureExtractor,
    test_compounds: List[CompoundData],
):
    """Test feature caching."""
    compound = test_compounds[0]

    # First extraction
    features1 = feature_extractor.extract_features(compound)

    # Second extraction (should use cache)
    features2 = feature_extractor.extract_features(compound)

    # Check features match
    assert all(np.array_equal(features1[k], features2[k]) for k in features1.keys())


def test_batch_extraction(
    feature_extractor: FeatureExtractor,
    test_compounds: List[CompoundData],
):
    """Test batch feature extraction."""
    # Extract features
    feature_matrix = feature_extractor.extract_batch(test_compounds)

    # Check matrix structure
    assert isinstance(feature_matrix, np.ndarray)
    assert feature_matrix.shape[0] == len(test_compounds)
    assert feature_matrix.shape[1] > 0
    assert not np.any(np.isnan(feature_matrix))


def test_feature_persistence(
    feature_extractor: FeatureExtractor,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test feature persistence."""
    compound = test_compounds[0]

    # Extract and save features
    features1 = feature_extractor.extract_features(compound)
    feature_extractor.save_features(str(tmp_path / "features"))

    # Create new extractor and load features
    new_extractor = FeatureExtractor(
        cache_dir=str(tmp_path / "new_cache"),
        feature_types=["morgan", "maccs", "rdkit", "pharmacophore"],
    )
    new_extractor.load_features(str(tmp_path / "features"))

    # Extract features with new extractor
    features2 = new_extractor.extract_features(compound)

    # Check features match
    assert all(np.array_equal(features1[k], features2[k]) for k in features1.keys())


def test_invalid_smiles(feature_extractor: FeatureExtractor):
    """Test handling of invalid SMILES."""
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    with pytest.raises(ValueError):
        feature_extractor.extract_features(invalid_compound)


def test_custom_parameters():
    """Test custom feature parameters."""
    # Morgan with custom radius and length
    morgan = MorganFeatures(radius=3, n_bits=1024)
    mol = Chem.MolFromSmiles("CCO")
    features = morgan.extract(mol)
    assert features.shape == (1024,)

    # RDKit with custom descriptor set
    rdkit = RDKitFeatures(descriptors=["MolWt", "NumRotatableBonds"])
    features = rdkit.extract(mol)
    assert len(features) == 2

    # Pharmacophore with custom features
    pharm = PharmacophoreFeatures(features=["Donor", "Acceptor", "Aromatic"])
    features = pharm.extract(mol)
    assert len(features) == 3


def test_feature_statistics(
    feature_extractor: FeatureExtractor,
    test_compounds: List[CompoundData],
):
    """Test feature statistics calculation."""
    # Extract batch features
    feature_matrix = feature_extractor.extract_batch(test_compounds)

    # Calculate statistics
    stats = feature_extractor.calculate_statistics(feature_matrix)

    assert isinstance(stats, Dict)
    assert "mean" in stats
    assert "std" in stats
    assert "min" in stats
    assert "max" in stats
    assert all(isinstance(v, np.ndarray) for v in stats.values())
    assert all(v.shape == (feature_matrix.shape[1],) for v in stats.values())


def test_feature_importance(
    feature_extractor: FeatureExtractor,
    test_compounds: List[CompoundData],
):
    """Test feature importance calculation."""
    # Create random labels
    labels = np.random.choice([0, 1], size=len(test_compounds))

    # Calculate importance
    importance = feature_extractor.calculate_importance(
        test_compounds,
        labels,
        method="random_forest",
    )

    assert isinstance(importance, Dict)
    assert all(k in importance for k in feature_extractor.feature_types)
    assert all(isinstance(v, np.ndarray) for v in importance.values())
    assert all(len(v) > 0 for v in importance.values())
    assert all(not np.any(np.isnan(v)) for v in importance.values())
