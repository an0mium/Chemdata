"""Tests for ML predictors module."""

import pytest
from binding_data_processor.models.compound.ml import predictors


@pytest.fixture
def mock_predictor():
    """Create a mock predictor for testing."""
    return predictors.BasePredictor()


def test_predictor_initialization(mock_predictor):
    """Test predictor initialization."""
    assert isinstance(mock_predictor, predictors.BasePredictor)
    assert hasattr(mock_predictor, "predict")


def test_prediction(mock_predictor):
    """Test prediction functionality."""
    with pytest.raises(NotImplementedError):
        mock_predictor.predict(None)


def test_model_training(mock_predictor):
    """Test model training."""
    with pytest.raises(NotImplementedError):
        mock_predictor.train(None, None)


def test_feature_extraction(mock_predictor):
    """Test feature extraction."""
    with pytest.raises(NotImplementedError):
        mock_predictor.extract_features(None)
