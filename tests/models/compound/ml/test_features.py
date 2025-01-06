"""Tests for ML features module."""

import pytest
import numpy as np
from binding_data_processor.models.compound.ml import features


@pytest.fixture
def mock_feature_extractor():
    """Create a mock feature extractor for testing."""
    return features.FeatureExtractor()


def test_feature_extractor_initialization(mock_feature_extractor):
    """Test feature extractor initialization."""
    assert isinstance(mock_feature_extractor, features.FeatureExtractor)
    assert hasattr(mock_feature_extractor, "extract")


def test_feature_extraction(mock_feature_extractor):
    """Test feature extraction."""
    with pytest.raises(NotImplementedError):
        mock_feature_extractor.extract(None)


def test_feature_validation(mock_feature_extractor):
    """Test feature validation."""
    with pytest.raises(NotImplementedError):
        mock_feature_extractor.validate_features(np.array([]))


def test_feature_scaling(mock_feature_extractor):
    """Test feature scaling."""
    with pytest.raises(NotImplementedError):
        mock_feature_extractor.scale_features(np.array([]))


def test_feature_selection(mock_feature_extractor):
    """Test feature selection."""
    with pytest.raises(NotImplementedError):
        mock_feature_extractor.select_features(np.array([]), np.array([]))


def test_feature_importance(mock_feature_extractor):
    """Test feature importance calculation."""
    with pytest.raises(NotImplementedError):
        mock_feature_extractor.calculate_importance(np.array([]), np.array([]))


def test_feature_engineering(mock_feature_extractor):
    """Test feature engineering."""
    with pytest.raises(NotImplementedError):
        mock_feature_extractor.engineer_features(np.array([]))
