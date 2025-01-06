"""Tests for ML ensemble models."""

import logging
import pytest
from pathlib import Path
from typing import List

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import accuracy_score, mean_squared_error

from ....core import CompoundData
from ..ensemble import ClassifierEnsemble, RegressorEnsemble


@pytest.fixture
def test_compounds() -> List[CompoundData]:
    """Create test compounds."""
    compounds = []
    for i in range(10):
        compound = CompoundData(
            name=f"Test Compound {i}",
            smiles=f"CC{i}CC",
            cas_number=f"123-45-{i}",
        )
        compounds.append(compound)
    return compounds


@pytest.fixture
def classifier_ensemble(tmp_path: Path) -> ClassifierEnsemble:
    """Create test classifier ensemble."""
    return ClassifierEnsemble(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        feature_types=["morgan", "maccs", "rdkit"],
        log_level=logging.DEBUG,
    )


@pytest.fixture
def regressor_ensemble(tmp_path: Path) -> RegressorEnsemble:
    """Create test regressor ensemble."""
    return RegressorEnsemble(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        feature_types=["morgan", "maccs", "rdkit"],
        log_level=logging.DEBUG,
    )


def test_classifier_ensemble_creation(classifier_ensemble: ClassifierEnsemble):
    """Test classifier ensemble creation."""
    ensemble = classifier_ensemble.create_ensemble(
        "test",
        ensemble_type="voting",
        base_estimator=RandomForestClassifier(n_estimators=10),
        n_models=3,
    )

    assert len(ensemble.estimators_) == 3
    for model in ensemble.estimators_:
        assert isinstance(model, RandomForestClassifier)


def test_regressor_ensemble_creation(regressor_ensemble: RegressorEnsemble):
    """Test regressor ensemble creation."""
    ensemble = regressor_ensemble.create_ensemble(
        "test",
        ensemble_type="stacking",
        base_estimator=RandomForestRegressor(n_estimators=10),
        n_models=3,
    )

    assert len(ensemble.estimators_) == 3
    for model in ensemble.estimators_:
        assert isinstance(model, RandomForestRegressor)


def test_classifier_ensemble_prediction(
    classifier_ensemble: ClassifierEnsemble,
    test_compounds: List[CompoundData],
):
    """Test classifier ensemble prediction."""
    # Create and fit ensemble
    labels = np.random.choice(["A", "B"], size=len(test_compounds))
    ensemble = classifier_ensemble.fit("test", test_compounds, labels)

    # Make predictions
    for compound in test_compounds:
        pred, conf = classifier_ensemble.predict("test", compound)
        assert isinstance(pred, str)
        assert isinstance(conf, float)
        assert 0 <= conf <= 1


def test_regressor_ensemble_prediction(
    regressor_ensemble: RegressorEnsemble,
    test_compounds: List[CompoundData],
):
    """Test regressor ensemble prediction."""
    # Create and fit ensemble
    values = np.random.random(len(test_compounds))
    ensemble = regressor_ensemble.fit("test", test_compounds, values)

    # Make predictions
    for compound in test_compounds:
        pred, conf = regressor_ensemble.predict("test", compound)
        assert isinstance(pred, float)
        assert isinstance(conf, float)
        assert 0 <= conf <= 1


def test_classifier_ensemble_persistence(
    classifier_ensemble: ClassifierEnsemble,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test classifier ensemble persistence."""
    # Create and fit ensemble
    labels = np.random.choice(["A", "B"], size=len(test_compounds))
    ensemble = classifier_ensemble.fit("test", test_compounds, labels)

    # Save models
    save_dir = tmp_path / "saved_models"
    classifier_ensemble.save_models(save_dir)

    # Create new ensemble and load models
    new_ensemble = ClassifierEnsemble(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
        feature_types=["morgan", "maccs", "rdkit"],
    )

    # Compare predictions
    for compound in test_compounds:
        pred1, conf1 = classifier_ensemble.predict("test", compound)
        pred2, conf2 = new_ensemble.predict("test", compound)
        assert pred1 == pred2
        assert abs(conf1 - conf2) < 1e-6


def test_regressor_ensemble_persistence(
    regressor_ensemble: RegressorEnsemble,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test regressor ensemble persistence."""
    # Create and fit ensemble
    values = np.random.random(len(test_compounds))
    ensemble = regressor_ensemble.fit("test", test_compounds, values)

    # Save models
    save_dir = tmp_path / "saved_models"
    regressor_ensemble.save_models(save_dir)

    # Create new ensemble and load models
    new_ensemble = RegressorEnsemble(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
        feature_types=["morgan", "maccs", "rdkit"],
    )

    # Compare predictions
    for compound in test_compounds:
        pred1, conf1 = regressor_ensemble.predict("test", compound)
        pred2, conf2 = new_ensemble.predict("test", compound)
        assert abs(pred1 - pred2) < 1e-6
        assert abs(conf1 - conf2) < 1e-6


def test_classifier_ensemble_metrics(
    classifier_ensemble: ClassifierEnsemble,
    test_compounds: List[CompoundData],
):
    """Test classifier ensemble metrics."""
    # Create and fit ensemble
    labels = np.random.choice(["A", "B"], size=len(test_compounds))
    metrics = classifier_ensemble.fit("test", test_compounds, labels)

    assert isinstance(metrics, dict)
    assert "accuracy" in metrics
    assert "f1" in metrics
    assert "precision" in metrics
    assert "recall" in metrics
    assert all(isinstance(v, float) for v in metrics.values())
    assert all(0 <= v <= 1 for v in metrics.values())


def test_regressor_ensemble_metrics(
    regressor_ensemble: RegressorEnsemble,
    test_compounds: List[CompoundData],
):
    """Test regressor ensemble metrics."""
    # Create and fit ensemble
    values = np.random.random(len(test_compounds))
    metrics = regressor_ensemble.fit("test", test_compounds, values)

    assert isinstance(metrics, dict)
    assert "mse" in metrics
    assert "rmse" in metrics
    assert "mae" in metrics
    assert "r2" in metrics
    assert all(isinstance(v, float) for v in metrics.values())
    assert all(v >= 0 for v in metrics.values())


def test_classifier_ensemble_error_handling(classifier_ensemble: ClassifierEnsemble):
    """Test classifier ensemble error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    with pytest.raises(ValueError):
        classifier_ensemble.predict("test", invalid_compound)

    # Test with missing model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    with pytest.raises(KeyError):
        classifier_ensemble.predict("missing_model", valid_compound)


def test_regressor_ensemble_error_handling(regressor_ensemble: RegressorEnsemble):
    """Test regressor ensemble error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    with pytest.raises(ValueError):
        regressor_ensemble.predict("test", invalid_compound)

    # Test with missing model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    with pytest.raises(KeyError):
        regressor_ensemble.predict("missing_model", valid_compound)
