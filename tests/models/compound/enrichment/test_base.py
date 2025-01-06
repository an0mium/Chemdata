"""Tests for enrichment base module."""

import pytest
from binding_data_processor.models.compound.enrichment import base


@pytest.fixture
def mock_enricher():
    """Create a mock enricher for testing."""
    return base.BaseEnricher()


def test_enricher_initialization(mock_enricher):
    """Test enricher initialization."""
    assert isinstance(mock_enricher, base.BaseEnricher)
    assert hasattr(mock_enricher, "enrich")


def test_enrichment(mock_enricher):
    """Test enrichment functionality."""
    with pytest.raises(NotImplementedError):
        mock_enricher.enrich(None)


def test_data_validation(mock_enricher):
    """Test data validation."""
    with pytest.raises(NotImplementedError):
        mock_enricher.validate_data(None)


def test_data_integration(mock_enricher):
    """Test data integration."""
    with pytest.raises(NotImplementedError):
        mock_enricher.integrate_data({})


def test_source_configuration(mock_enricher):
    """Test source configuration."""
    with pytest.raises(NotImplementedError):
        mock_enricher.configure_sources([])
