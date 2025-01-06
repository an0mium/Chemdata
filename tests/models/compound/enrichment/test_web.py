"""Tests for web enrichment module."""

import pytest
from binding_data_processor.models.compound.enrichment import web


@pytest.fixture
def mock_enricher():
    """Create a mock web enricher for testing."""
    return web.WebEnricher()


def test_enricher_initialization(mock_enricher):
    """Test enricher initialization."""
    assert isinstance(mock_enricher, web.WebEnricher)
    assert hasattr(mock_enricher, "enrich")


def test_web_enrichment(mock_enricher):
    """Test web enrichment functionality."""
    with pytest.raises(NotImplementedError):
        mock_enricher.enrich(None)


def test_data_validation(mock_enricher):
    """Test data validation."""
    with pytest.raises(NotImplementedError):
        mock_enricher.validate_data(None)


def test_data_integration(mock_enricher):
    """Test data integration."""
    with pytest.raises(NotImplementedError):
        mock_enricher.integrate_data(None, None)
