"""Tests for enhanced web enrichment manager."""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch

from ..manager_enhanced import (
    WebEnrichmentManagerEnhanced,
    EnrichmentConfig,
)
from ...models.compound import Compound
from ...pipeline.infrastructure.circuit_breaker import CircuitConfig


@pytest.fixture
def config(tmp_path):
    """Create test configuration."""
    return EnrichmentConfig(
        reddit_client_id="test_id",
        reddit_client_secret="test_secret",
        twitter_bearer_token="test_token",
        model_dir=tmp_path / "models",
        cache_dir=tmp_path / "cache",
        circuit_config=CircuitConfig(
            failure_threshold=2,
            failure_timeout=1,
            reset_timeout=1,
        ),
    )


@pytest.fixture
def manager(config):
    """Create test manager."""
    return WebEnrichmentManagerEnhanced(config)


@pytest.fixture
def compound():
    """Create test compound."""
    return Compound(
        name="Test Compound",
        smiles="CC1=CC=CC=C1",
        cas_number="100-00-0",
    )


def test_initialization(manager):
    """Test manager initialization."""
    assert manager.http is not None
    assert manager.swiss_client is not None
    assert manager.community_client is not None
    assert manager.social_client is not None
    assert manager.executor is not None


def test_enrich_compounds(manager, compound):
    """Test compound enrichment."""
    # Mock client process methods
    manager.swiss_client.process_compounds = Mock()
    manager.community_client.process_compounds = Mock()
    manager.social_client.process_compounds = Mock()

    # Enrich compound
    manager.enrich_compounds([compound])

    # Check client calls
    manager.swiss_client.process_compounds.assert_called_once()
    manager.community_client.process_compounds.assert_called_once()
    manager.social_client.process_compounds.assert_called_once()

    # Check metadata
    assert compound.enrichment_metadata is not None
    assert "timestamp" in compound.enrichment_metadata


def test_skip_predictions(manager, compound):
    """Test skipping predictions."""
    # Mock client process methods
    manager.swiss_client.process_compounds = Mock()
    manager.community_client.process_compounds = Mock()
    manager.social_client.process_compounds = Mock()

    # Enrich compound with skip_predictions
    manager.enrich_compounds([compound], skip_predictions=True)

    # Check client calls
    manager.swiss_client.process_compounds.assert_not_called()
    manager.community_client.process_compounds.assert_called_once()
    manager.social_client.process_compounds.assert_called_once()


def test_skip_web_data(manager, compound):
    """Test skipping web data."""
    # Mock client process methods
    manager.swiss_client.process_compounds = Mock()
    manager.community_client.process_compounds = Mock()
    manager.social_client.process_compounds = Mock()

    # Enrich compound with skip_web_data
    manager.enrich_compounds([compound], skip_web_data=True)

    # Check client calls
    manager.swiss_client.process_compounds.assert_called_once()
    manager.community_client.process_compounds.assert_not_called()
    manager.social_client.process_compounds.assert_not_called()


def test_batch_processing(manager):
    """Test batch processing."""
    # Create test compounds
    compounds = [
        Compound(name=f"Compound {i}", smiles="CC", cas_number=str(i))
        for i in range(10)
    ]

    # Mock client process methods
    manager.swiss_client.process_compounds = Mock()
    manager.community_client.process_compounds = Mock()
    manager.social_client.process_compounds = Mock()

    # Set small batch size
    manager.config.batch_size = 3

    # Enrich compounds
    manager.enrich_compounds(compounds)

    # Should have processed in 4 batches
    assert manager.swiss_client.process_compounds.call_count == 4
    assert manager.community_client.process_compounds.call_count == 4
    assert manager.social_client.process_compounds.call_count == 4


def test_error_handling(manager, compound):
    """Test error handling."""
    # Mock client to raise error
    manager.swiss_client.process_compounds = Mock(
        side_effect=Exception("Test error")
    )

    # Should log error but not raise
    manager.enrich_compounds([compound])

    # Metadata should still be present
    assert compound.enrichment_metadata is not None
    assert "timestamp" in compound.enrichment_metadata


def test_metrics(manager):
    """Test metrics collection."""
    metrics = manager.get_metrics()
    assert "http_client" in metrics
    assert "swiss_client" in metrics
    assert "community_client" in metrics
    assert "social_client" in metrics


def test_cleanup(manager):
    """Test cleanup on close."""
    # Mock client close methods
    manager.swiss_client.close = Mock()
    manager.community_client.close = Mock()
    manager.social_client.close = Mock()
    manager.http.close = Mock()

    # Close manager
    manager.close()

    # Check client cleanup
    manager.swiss_client.close.assert_called_once()
    manager.community_client.close.assert_called_once()
    manager.social_client.close.assert_called_once()
    manager.http.close.assert_called_once()
