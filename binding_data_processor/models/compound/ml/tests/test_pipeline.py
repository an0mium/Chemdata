"""Tests for ML pipeline integration."""

import logging
import pytest
from pathlib import Path
from typing import List, Dict

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from ....core import CompoundData
from ..pipeline import MLPipeline
from ..features import FeatureExtractor
from ..ensemble import ClassifierEnsemble, RegressorEnsemble


@pytest.fixture
def test_compounds() -> List[CompoundData]:
    """Create test compounds."""
    compounds = []
    for i in range(10):
        compound = CompoundData(
            name=f"Test Compound {i}",
            smiles=f"CC{i}CC",  # Simple alkanes
            cas_number=f"123-45-{i}",
        )
        compounds.append(compound)
    return compounds


@pytest.fixture
def ml_pipeline(tmp_path: Path) -> MLPipeline:
    """Create ML pipeline."""
    return MLPipeline(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        feature_types=["morgan", "maccs", "rdkit"],
        log_level=logging.DEBUG,
    )


def test_pipeline_initialization(ml_pipeline: MLPipeline):
    """Test pipeline initialization."""
    assert isinstance(ml_pipeline.feature_extractor, FeatureExtractor)
    assert isinstance(ml_pipeline.classifiers, Dict)
    assert isinstance(ml_pipeline.regressors, Dict)


def test_classification_training(
    ml_pipeline: MLPipeline,
    test_compounds: List[CompoundData],
):
    """Test classification model training."""
    # Create random labels
    labels = np.random.choice(["A", "B"], size=len(test_compounds))

    # Train model
    metrics = ml_pipeline.train_classifier(
        "test_clf",
        test_compounds,
        labels,
        model_type="random_forest",
        n_models=3,
    )

    # Check metrics
    assert isinstance(metrics, Dict)
    assert "accuracy" in metrics
    assert "f1" in metrics
    assert "precision" in metrics
    assert "recall" in metrics
    assert all(isinstance(v, float) for v in metrics.values())
    assert all(0 <= v <= 1 for v in metrics.values())


def test_regression_training(
    ml_pipeline: MLPipeline,
    test_compounds: List[CompoundData],
):
    """Test regression model training."""
    # Create random values
    values = np.random.random(len(test_compounds))

    # Train model
    metrics = ml_pipeline.train_regressor(
        "test_reg",
        test_compounds,
        values,
        model_type="random_forest",
        n_models=3,
    )

    # Check metrics
    assert isinstance(metrics, Dict)
    assert "mse" in metrics
    assert "rmse" in metrics
    assert "mae" in metrics
    assert "r2" in metrics
    assert all(isinstance(v, float) for v in metrics.values())
    assert all(v >= 0 for v in metrics.values())


def test_classification_prediction(
    ml_pipeline: MLPipeline,
    test_compounds: List[CompoundData],
):
    """Test classification prediction."""
    # Train model
    labels = np.random.choice(["A", "B"], size=len(test_compounds))
    ml_pipeline.train_classifier(
        "test_clf",
        test_compounds,
        labels,
        model_type="random_forest",
        n_models=3,
    )

    # Make predictions
    for compound in test_compounds:
        pred, conf = ml_pipeline.predict_classifier("test_clf", compound)
        assert isinstance(pred, str)
        assert isinstance(conf, float)
        assert 0 <= conf <= 1


def test_regression_prediction(
    ml_pipeline: MLPipeline,
    test_compounds: List[CompoundData],
):
    """Test regression prediction."""
    # Train model
    values = np.random.random(len(test_compounds))
    ml_pipeline.train_regressor(
        "test_reg",
        test_compounds,
        values,
        model_type="random_forest",
        n_models=3,
    )

    # Make predictions
    for compound in test_compounds:
        pred, conf = ml_pipeline.predict_regressor("test_reg", compound)
        assert isinstance(pred, float)
        assert isinstance(conf, float)
        assert 0 <= conf <= 1


def test_model_persistence(
    ml_pipeline: MLPipeline,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test model persistence."""
    # Train models
    labels = np.random.choice(["A", "B"], size=len(test_compounds))
    values = np.random.random(len(test_compounds))

    ml_pipeline.train_classifier(
        "test_clf",
        test_compounds,
        labels,
        model_type="random_forest",
        n_models=3,
    )
    ml_pipeline.train_regressor(
        "test_reg",
        test_compounds,
        values,
        model_type="random_forest",
        n_models=3,
    )

    # Save models
    save_dir = tmp_path / "saved_models"
    ml_pipeline.save_models(save_dir)

    # Create new pipeline and load models
    new_pipeline = MLPipeline(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
        feature_types=["morgan", "maccs", "rdkit"],
    )

    # Compare predictions
    for compound in test_compounds:
        # Classification
        pred1, conf1 = ml_pipeline.predict_classifier("test_clf", compound)
        pred2, conf2 = new_pipeline.predict_classifier("test_clf", compound)
        assert pred1 == pred2
        assert abs(conf1 - conf2) < 1e-6

        # Regression
        pred1, conf1 = ml_pipeline.predict_regressor("test_reg", compound)
        pred2, conf2 = new_pipeline.predict_regressor("test_reg", compound)
        assert abs(pred1 - pred2) < 1e-6
        assert abs(conf1 - conf2) < 1e-6


def test_feature_caching(
    ml_pipeline: MLPipeline,
    test_compounds: List[CompoundData],
):
    """Test feature caching in pipeline."""
    compound = test_compounds[0]

    # First extraction
    features1 = ml_pipeline.feature_extractor.extract_features(compound)

    # Second extraction (should use cache)
    features2 = ml_pipeline.feature_extractor.extract_features(compound)

    # Check features match
    assert all(np.array_equal(features1[k], features2[k]) for k in features1.keys())


def test_error_handling(ml_pipeline: MLPipeline):
    """Test error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    with pytest.raises(ValueError):
        ml_pipeline.predict_classifier("test_clf", invalid_compound)

    with pytest.raises(ValueError):
        ml_pipeline.predict_regressor("test_reg", invalid_compound)

    # Test with missing model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    with pytest.raises(KeyError):
        ml_pipeline.predict_classifier("missing_clf", valid_compound)

    with pytest.raises(KeyError):
        ml_pipeline.predict_regressor("missing_reg", valid_compound)


def test_model_metrics(
    ml_pipeline: MLPipeline,
    test_compounds: List[CompoundData],
):
    """Test model metrics calculation."""
    # Train models
    labels = np.random.choice(["A", "B"], size=len(test_compounds))
    values = np.random.random(len(test_compounds))

    # Classification metrics
    clf_metrics = ml_pipeline.train_classifier(
        "test_clf",
        test_compounds,
        labels,
        model_type="random_forest",
        n_models=3,
    )

    assert isinstance(clf_metrics, Dict)
    assert all(k in clf_metrics for k in ["accuracy", "f1", "precision", "recall"])
    assert all(isinstance(v, float) for v in clf_metrics.values())
    assert all(0 <= v <= 1 for v in clf_metrics.values())

    # Regression metrics
    reg_metrics = ml_pipeline.train_regressor(
        "test_reg",
        test_compounds,
        values,
        model_type="random_forest",
        n_models=3,
    )

    assert isinstance(reg_metrics, Dict)
    assert all(k in reg_metrics for k in ["mse", "rmse", "mae", "r2"])
    assert all(isinstance(v, float) for v in reg_metrics.values())
    assert all(v >= 0 for v in reg_metrics.values())
