"""Tests for analysis base module."""

import pytest
from binding_data_processor.models.compound.analysis import base


@pytest.fixture
def mock_analyzer():
    """Create a mock analyzer for testing."""
    return base.BaseAnalyzer()


def test_analyzer_initialization(mock_analyzer):
    """Test analyzer initialization."""
    assert isinstance(mock_analyzer, base.BaseAnalyzer)
    assert hasattr(mock_analyzer, "analyze")


def test_analysis(mock_analyzer):
    """Test analysis functionality."""
    with pytest.raises(NotImplementedError):
        mock_analyzer.analyze(None)


def test_result_validation(mock_analyzer):
    """Test result validation."""
    with pytest.raises(NotImplementedError):
        mock_analyzer.validate_results(None)


def test_result_aggregation(mock_analyzer):
    """Test result aggregation."""
    with pytest.raises(NotImplementedError):
        mock_analyzer.aggregate_results([])
