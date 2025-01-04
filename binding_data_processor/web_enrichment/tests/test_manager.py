"""Tests for web enrichment manager."""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

from ..manager import WebEnrichmentManager, EnrichmentConfig
from ..http_client import HTTPClient
from ..swiss_client import SwissClient
from ..community_client import CommunityClient
from ..social_client import SocialClient
from ...models.compound import Compound


@pytest.fixture
def test_config():
    """Create test enrichment configuration."""
    return EnrichmentConfig(
        reddit_client_id="test_reddit_id",
        reddit_client_secret="test_reddit_secret",
        twitter_bearer_token="test_twitter_token",
        model_dir=Path("/tmp/models"),
        cache_dir=Path("/tmp/cache"),
        n_workers=2,
        batch_size=10,
    )


@pytest.fixture
def test_compounds():
    """Create test compounds."""
    compounds = []
    
    # Create test compounds
    caffeine = Compound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    compounds.append(caffeine)
    
    amphetamine = Compound(
        name="Amphetamine",
        smiles="CC(N)CC1=CC=CC=C1",
        cas_number="300-62-9",
    )
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def mock_clients():
    """Create mock web clients."""
    with patch("binding_data_processor.web_enrichment.manager.HTTPClient") as mock_http, \
         patch("binding_data_processor.web_enrichment.manager.SwissClient") as mock_swiss, \
         patch("binding_data_processor.web_enrichment.manager.CommunityClient") as mock_community, \
         patch("binding_data_processor.web_enrichment.manager.SocialClient") as mock_social:
        
        yield {
            "http": mock_http,
            "swiss": mock_swiss,
            "community": mock_community,
            "social": mock_social,
        }


class TestWebEnrichmentManager:
    """Tests for WebEnrichmentManager class."""

    def test_initialization(self, test_config, mock_clients):
        """Test manager initialization."""
        manager = WebEnrichmentManager(test_config)

        # Check client initialization
        assert isinstance(manager.http, HTTPClient)
        assert isinstance(manager.swiss_client, SwissClient)
        assert isinstance(manager.community_client, CommunityClient)
        assert isinstance(manager.social_client, SocialClient)

        # Check config
        assert manager.config == test_config
        assert manager.executor._max_workers == test_config.n_workers

    def test_enrich_compounds(self, test_config, test_compounds, mock_clients):
        """Test compound enrichment."""
        manager = WebEnrichmentManager(test_config)

        # Mock client process_compounds methods
        manager.swiss_client.process_compounds = Mock()
        manager.community_client.process_compounds = Mock()
        manager.social_client.process_compounds = Mock()

        # Record start time
        start_time = datetime.now()

        # Enrich compounds
        manager.enrich_compounds(test_compounds)

        # Check client calls
        assert manager.swiss_client.process_compounds.call_count == len(test_compounds)
        assert manager.community_client.process_compounds.call_count == len(test_compounds)
        assert manager.social_client.process_compounds.call_count == len(test_compounds)

        # Check compound data
        for compound in test_compounds:
            assert isinstance(compound.swiss_data, dict)
            assert isinstance(compound.community_data, dict)
            assert isinstance(compound.social_data, dict)
            assert isinstance(compound.enrichment_metadata, dict)
            
            # Check timestamp
            timestamp = datetime.fromisoformat(compound.enrichment_metadata["timestamp"])
            assert start_time <= timestamp <= datetime.now()
            assert timestamp - start_time < timedelta(seconds=10)
            
            assert "sources" in compound.enrichment_metadata

    def test_skip_predictions(self, test_config, test_compounds, mock_clients):
        """Test skipping predictions."""
        manager = WebEnrichmentManager(test_config)

        # Mock client process_compounds methods
        manager.swiss_client.process_compounds = Mock()
        manager.community_client.process_compounds = Mock()
        manager.social_client.process_compounds = Mock()

        # Enrich compounds with skip_predictions=True
        manager.enrich_compounds(test_compounds, skip_predictions=True)

        # Check client calls
        assert manager.swiss_client.process_compounds.call_count == 0
        assert manager.community_client.process_compounds.call_count == len(test_compounds)
        assert manager.social_client.process_compounds.call_count == len(test_compounds)

    def test_skip_web_data(self, test_config, test_compounds, mock_clients):
        """Test skipping web data."""
        manager = WebEnrichmentManager(test_config)

        # Mock client process_compounds methods
        manager.swiss_client.process_compounds = Mock()
        manager.community_client.process_compounds = Mock()
        manager.social_client.process_compounds = Mock()

        # Enrich compounds with skip_web_data=True
        manager.enrich_compounds(test_compounds, skip_web_data=True)

        # Check client calls
        assert manager.swiss_client.process_compounds.call_count == len(test_compounds)
        assert manager.community_client.process_compounds.call_count == 0
        assert manager.social_client.process_compounds.call_count == 0

    def test_batch_processing(self, test_config, test_compounds, mock_clients):
        """Test batch processing."""
        # Set small batch size
        test_config.batch_size = 1
        manager = WebEnrichmentManager(test_config)

        # Mock client process_compounds methods
        manager.swiss_client.process_compounds = Mock()
        manager.community_client.process_compounds = Mock()
        manager.social_client.process_compounds = Mock()

        # Enrich compounds
        manager.enrich_compounds(test_compounds)

        # Check batch processing
        assert manager.swiss_client.process_compounds.call_count == len(test_compounds)
        assert manager.community_client.process_compounds.call_count == len(test_compounds)
        assert manager.social_client.process_compounds.call_count == len(test_compounds)

    def test_error_handling(self, test_config, test_compounds, mock_clients):
        """Test error handling."""
        manager = WebEnrichmentManager(test_config)

        # Mock client process_compounds methods to raise errors
        manager.swiss_client.process_compounds = Mock(side_effect=Exception("Swiss error"))
        manager.community_client.process_compounds = Mock(side_effect=Exception("Community error"))
        manager.social_client.process_compounds = Mock(side_effect=Exception("Social error"))

        # Enrich compounds (should not raise errors)
        manager.enrich_compounds(test_compounds)

        # Check compound data still initialized
        for compound in test_compounds:
            assert isinstance(compound.swiss_data, dict)
            assert isinstance(compound.community_data, dict)
            assert isinstance(compound.social_data, dict)
            assert isinstance(compound.enrichment_metadata, dict)

    def test_cleanup(self, test_config, mock_clients):
        """Test cleanup on close."""
        manager = WebEnrichmentManager(test_config)

        # Mock close methods
        manager.http.close = Mock()
        manager.swiss_client.close = Mock()
        manager.community_client.close = Mock()
        manager.social_client.close = Mock()

        # Close manager
        manager.close()

        # Check close calls
        assert manager.http.close.called
        assert manager.swiss_client.close.called
        assert manager.community_client.close.called
        assert manager.social_client.close.called
