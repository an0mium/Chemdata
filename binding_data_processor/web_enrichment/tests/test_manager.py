"""Tests for web enrichment manager."""

import logging
from datetime import datetime, timedelta
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from ..manager import (
    WebEnrichmentManager,
    WebEnrichmentError,
    EnrichmentConfig,
)
from ..clients.swiss import SwissClient
from ..clients.community import CommunityClient
from ..clients.social import SocialClient
from ...models.compound import Compound
from ...pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig
from ...pipeline.infrastructure.monitoring import MetricsCollector


@pytest.fixture
def test_compounds():
    """Create test compounds with realistic data."""
    compounds = []

    # Caffeine
    caffeine = Compound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    compounds.append(caffeine)

    # Amphetamine
    amphetamine = Compound(
        name="Amphetamine",
        smiles="CC(N)CC1=CC=CC=C1",
        cas_number="300-62-9",
    )
    compounds.append(amphetamine)

    # Aspirin
    aspirin = Compound(
        name="Aspirin",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        cas_number="50-78-2",
    )
    compounds.append(aspirin)

    return compounds


@pytest.fixture
def config():
    """Create test config."""
    return EnrichmentConfig(
        reddit_client_id="test_id",
        reddit_client_secret="test_secret",
        twitter_bearer_token="test_token",
        skip_predictions=False,
        skip_web_data=False,
        use_cache=True,
        n_workers=2,
        batch_size=10,
        model_dir=Path("/tmp/models"),
        cache_dir=Path("/tmp/cache"),
        circuit_config=CircuitBreakerConfig(
            failure_threshold=3,
            recovery_timeout=60,
        ),
        client_configs={
            "swiss": {"base_url": "http://test.swiss.com"},
            "community": {"base_url": "http://test.community.com"},
            "social": {"base_url": "http://test.social.com"},
        },
    )


@pytest.fixture
def manager(config):
    """Create test manager."""
    return WebEnrichmentManager(
        config=config,
        logger=logging.getLogger("test"),
        metrics_collector=Mock(spec=MetricsCollector),
    )


def test_manager_initialization(manager, config):
    """Test manager initialization."""
    # Check client initialization
    assert len(manager.clients) == 3
    assert isinstance(manager.clients["swiss"], SwissClient)
    assert isinstance(manager.clients["community"], CommunityClient)
    assert isinstance(manager.clients["social"], SocialClient)

    # Check config
    assert manager.config == config
    assert manager.executor._max_workers == config.n_workers

    # Check initial state
    assert len(manager.processed_compounds) == 0
    assert len(manager.failed_compounds) == 0


def test_register_client(manager):
    """Test client registration."""
    # Create mock client class
    mock_client = Mock()
    mock_client.name = "test"

    # Register client
    with patch("..base.WebClient", return_value=mock_client):
        manager.register_client(
            mock_client.__class__,
            base_url="http://test.com",
        )

    assert "test" in manager.clients
    assert manager.clients["test"] == mock_client


def test_get_client(manager):
    """Test getting client."""
    # Get existing client
    client = manager.get_client("swiss")
    assert isinstance(client, SwissClient)

    # Get non-existent client
    with pytest.raises(WebEnrichmentError) as exc:
        manager.get_client("invalid")
    assert "Client not found" in str(exc.value)


def test_enrich_compounds_success(manager, test_compounds):
    """Test successful compound enrichment."""
    # Mock client responses
    swiss_data = {
        "targets": [
            {
                "target": "Adenosine A2a receptor",
                "probability": 0.95,
            },
        ],
    }
    community_data = {
        "reports": [
            {
                "source": "PsychonautWiki",
                "text": "Test report",
            },
        ],
    }
    social_data = {
        "posts": [
            {
                "platform": "Reddit",
                "text": "Test post",
            },
        ],
    }

    # Record start time
    start_time = datetime.now()

    # Mock client methods
    with patch.object(
        manager.clients["swiss"],
        "predict_targets",
        return_value=swiss_data,
    ), patch.object(
        manager.clients["community"],
        "get_data",
        return_value=community_data,
    ), patch.object(
        manager.clients["social"],
        "get_data",
        return_value=social_data,
    ):
        manager.enrich_compounds(test_compounds)

    # Check enrichment results
    for compound in test_compounds:
        # Check data
        assert compound.swiss_data == swiss_data
        assert compound.community_data == community_data
        assert compound.social_data == social_data

        # Check metadata
        assert "timestamp" in compound.enrichment_metadata
        timestamp = datetime.fromisoformat(compound.enrichment_metadata["timestamp"])
        assert start_time <= timestamp <= datetime.now()
        assert timestamp - start_time < timedelta(seconds=10)

        assert compound.enrichment_metadata["sources"] == [
            "swiss",
            "community",
            "social",
        ]

        # Check tracking
        assert compound.smiles in manager.processed_compounds


def test_enrich_compounds_skip_options(manager, test_compounds):
    """Test compound enrichment with skip options."""
    # Test skip predictions
    manager.enrich_compounds(
        compounds=test_compounds,
        skip_predictions=True,
        skip_web_data=False,
    )
    for compound in test_compounds:
        assert not compound.swiss_data
        assert compound.community_data
        assert compound.social_data

    # Test skip web data
    test_compounds_2 = test_compounds.copy()
    manager.enrich_compounds(
        compounds=test_compounds_2,
        skip_predictions=False,
        skip_web_data=True,
    )
    for compound in test_compounds_2:
        assert compound.swiss_data
        assert not compound.community_data
        assert not compound.social_data


def test_enrich_compounds_batch_processing(manager):
    """Test batch processing of compounds."""
    # Create test compounds
    compounds = [
        Compound(
            name=f"Compound {i}",
            smiles=f"C{'C' * i}",
            cas_number=f"123-{i:02d}-{i:1d}",
        )
        for i in range(25)  # More than batch size
    ]

    # Mock client responses
    swiss_data = {"targets": [{"target": "test", "probability": 0.9}]}
    community_data = {"reports": [{"source": "test", "text": "test"}]}
    social_data = {"posts": [{"platform": "test", "text": "test"}]}

    with patch.object(
        manager.clients["swiss"],
        "predict_targets",
        return_value=swiss_data,
    ), patch.object(
        manager.clients["community"],
        "get_data",
        return_value=community_data,
    ), patch.object(
        manager.clients["social"],
        "get_data",
        return_value=social_data,
    ):
        manager.enrich_compounds(compounds)

    # Check all compounds were processed
    assert len(manager.processed_compounds) == len(compounds)
    for compound in compounds:
        assert compound.smiles in manager.processed_compounds


def test_enrich_compounds_error_handling(manager, test_compounds):
    """Test error handling during enrichment."""
    # Mock client error
    error = Exception("Test error")
    with patch.object(
        manager.clients["swiss"],
        "predict_targets",
        side_effect=error,
    ):
        with pytest.raises(Exception) as exc:
            manager.enrich_compounds([test_compounds[0]])

    assert str(exc.value) == "Test error"
    assert test_compounds[0].smiles in manager.failed_compounds


def test_metrics(manager, test_compounds):
    """Test metrics collection."""
    # Mock successful enrichment
    with patch.object(
        manager.clients["swiss"],
        "predict_targets",
        return_value={},
    ), patch.object(
        manager.clients["community"],
        "get_data",
        return_value={},
    ), patch.object(
        manager.clients["social"],
        "get_data",
        return_value={},
    ):
        manager.enrich_compounds(test_compounds[:2])

    # Mock failed enrichment
    with patch.object(
        manager.clients["swiss"],
        "predict_targets",
        side_effect=Exception("Test error"),
    ):
        with pytest.raises(Exception):
            manager.enrich_compounds([test_compounds[2]])

    # Check metrics
    metrics = manager.get_metrics()
    assert metrics["processed_compounds"] == 2
    assert metrics["failed_compounds"] == 1
    assert "swiss" in metrics["clients"]
    assert "community" in metrics["clients"]
    assert "social" in metrics["clients"]


def test_cleanup(manager):
    """Test manager cleanup."""
    with patch.object(manager.executor, "shutdown") as mock_executor_close:
        with patch.object(manager.clients["swiss"], "close") as mock_swiss_close:
            with patch.object(manager.clients["community"], "close") as mock_community_close:
                with patch.object(manager.clients["social"], "close") as mock_social_close:
                    manager.close()

    mock_executor_close.assert_called_once()
    mock_swiss_close.assert_called_once()
    mock_community_close.assert_called_once()
    mock_social_close.assert_called_once()
