"""Tests for enhanced web enrichment manager.

This module contains tests for the WebEnrichmentManagerEnhanced class, which coordinates
web enrichment operations across multiple clients:

- HTTP client for making web requests
- Swiss client for predictions (SwissTargetPrediction, SwissADME)
- Community client for web scraping (PsychonautWiki, Erowid, TripSit)
- Social client for social media data (Reddit, Twitter)
- Crawl4AI client for research papers and patents

The tests verify:
1. Basic initialization and configuration
2. Single compound enrichment with all data sources
3. Batch processing with configurable batch size
4. Error handling and resilience
5. Circuit breaker integration for fault tolerance
6. Metrics collection for monitoring
7. Proper cleanup of resources
"""

import pytest
from unittest.mock import Mock, AsyncMock, MagicMock, patch

from ..manager_enhanced import (
    WebEnrichmentManagerEnhanced,
    EnrichmentConfig,
    EnrichmentResult,
)
from ...models.compound import Compound
from ...pipeline.infrastructure.circuit_breaker import CircuitConfig


# Import paths
BASE_PATH = "binding_data_processor.web_enrichment"
HTTP_CLIENT_PATH = f"{BASE_PATH}.http_client_enhanced.HTTPClientEnhanced"
SWISS_CLIENT_PATH = f"{BASE_PATH}.swiss_client.SwissClient"
COMMUNITY_CLIENT_PATH = f"{BASE_PATH}.community_client.CommunityClient"
SOCIAL_CLIENT_PATH = f"{BASE_PATH}.social_client.SocialClient"
CRAWL4AI_CLIENT_PATH = f"{BASE_PATH}.crawl4ai_client.Crawl4AIClient"


@pytest.fixture
def mock_http():
    """Mock HTTP client."""
    with patch(HTTP_CLIENT_PATH) as mock:
        client = MagicMock()
        client.get_metrics = MagicMock(return_value={"requests": 10})
        mock.return_value = client
        yield client


@pytest.fixture
def mock_swiss():
    """Mock Swiss client."""
    with patch(SWISS_CLIENT_PATH) as mock:
        client = MagicMock()
        client.get_compound_data = AsyncMock()
        client.process_compounds = Mock()
        client.processed_compounds = ["compound1"]
        client.failed_compounds = []
        mock.return_value = client
        yield client


@pytest.fixture
def mock_community():
    """Mock community client."""
    with patch(COMMUNITY_CLIENT_PATH) as mock:
        client = MagicMock()
        client.process_compounds = Mock()
        client.processed_compounds = ["compound1"]
        client.failed_compounds = []
        mock.return_value = client
        yield client


@pytest.fixture
def mock_social():
    """Mock social client."""
    with patch(SOCIAL_CLIENT_PATH) as mock:
        client = MagicMock()
        client.get_compound_mentions = AsyncMock()
        client.process_compounds = Mock()
        client.processed_compounds = ["compound1"]
        client.failed_compounds = []
        mock.return_value = client
        yield client


@pytest.fixture
def mock_crawl4ai():
    """Mock Crawl4AI client."""
    with patch(CRAWL4AI_CLIENT_PATH) as mock:
        client = MagicMock()
        client.scrape_research = AsyncMock()
        client.scrape_patents = AsyncMock()
        client.scrape_community = AsyncMock()
        mock.return_value = client
        yield client


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
def manager(config, mock_http, mock_swiss, mock_community, mock_social, mock_crawl4ai):
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
    assert manager.crawl4ai_client is not None
    assert manager.executor is not None


@pytest.mark.asyncio
async def test_enrich_compound(manager, mock_swiss, mock_social, mock_crawl4ai):
    """Test enriching a single compound."""
    # Mock client responses
    mock_swiss.get_compound_data.return_value = {"predictions": [0.8]}
    mock_social.get_compound_mentions.return_value = {"mentions": 10}
    mock_crawl4ai.scrape_research.return_value = [{"title": "Paper 1"}]
    mock_crawl4ai.scrape_patents.return_value = [{"title": "Patent 1"}]
    mock_crawl4ai.scrape_community.return_value = [{"title": "Post 1"}]

    # Test enrichment
    result = await manager.enrich_compound("test_compound", smiles="C1=CC=CC=C1")

    # Verify result
    assert isinstance(result, EnrichmentResult)
    assert result.research_papers == [{"title": "Paper 1"}]
    assert result.patents == [{"title": "Patent 1"}]
    assert result.community_posts == [{"title": "Post 1"}]
    assert result.swiss_data == {"predictions": [0.8]}
    assert result.social_data == {"mentions": 10}
    assert "test_compound C1=CC=CC=C1" in result.metadata["query"]


@pytest.mark.asyncio
async def test_enrich_compound_skip_web(manager, mock_swiss):
    """Test enriching a compound with web data skipped."""
    # Mock Swiss response
    mock_swiss.get_compound_data.return_value = {"predictions": [0.8]}

    # Test enrichment with web data skipped
    result = await manager.enrich_compound(
        "test_compound",
        skip_web_data=True,
    )

    # Verify result
    assert isinstance(result, EnrichmentResult)
    assert not result.research_papers
    assert not result.patents
    assert not result.community_posts
    assert result.swiss_data == {"predictions": [0.8]}
    assert not result.social_data


@pytest.mark.asyncio
async def test_enrich_compound_skip_predictions(manager, mock_crawl4ai, mock_social):
    """Test enriching a compound with predictions skipped."""
    # Mock client responses
    mock_crawl4ai.scrape_research.return_value = [{"title": "Paper 1"}]
    mock_crawl4ai.scrape_patents.return_value = [{"title": "Patent 1"}]
    mock_crawl4ai.scrape_community.return_value = [{"title": "Post 1"}]
    mock_social.get_compound_mentions.return_value = {"mentions": 10}

    # Test enrichment with predictions skipped
    result = await manager.enrich_compound(
        "test_compound",
        skip_predictions=True,
    )

    # Verify result
    assert isinstance(result, EnrichmentResult)
    assert result.research_papers == [{"title": "Paper 1"}]
    assert result.patents == [{"title": "Patent 1"}]
    assert result.community_posts == [{"title": "Post 1"}]
    assert not result.swiss_data
    assert result.social_data == {"mentions": 10}


def test_batch_processing(manager, compound):
    """Test batch processing."""
    # Create test compounds
    compounds = [Compound(name=f"Compound {i}", smiles="CC", cas_number=str(i)) for i in range(10)]

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
    manager.swiss_client.process_compounds = Mock(side_effect=Exception("Test error"))

    # Should log error but not raise
    manager.enrich_compounds([compound])

    # Metadata should still be present
    assert compound.enrichment_metadata is not None
    assert "timestamp" in compound.enrichment_metadata


@pytest.mark.asyncio
async def test_circuit_breaker(manager, mock_swiss):
    """Test circuit breaker integration."""
    # Mock Swiss client to fail
    mock_swiss.get_compound_data = AsyncMock(side_effect=Exception("Test error"))

    # Make multiple requests to trigger circuit breaker
    for _ in range(3):
        await manager.enrich_compound("test")

    # Circuit should be open
    assert manager.http.circuit_breaker.is_open()


def test_metrics(manager):
    """Test metrics collection."""
    metrics = manager.get_metrics()
    assert "http_client" in metrics
    assert "swiss_client" in metrics
    assert "community_client" in metrics
    assert "social_client" in metrics


def test_cleanup(manager, mock_http, mock_swiss, mock_community, mock_social):
    """Test cleanup on close."""
    manager.close()
    mock_http.close.assert_called_once()
    mock_swiss.close.assert_called_once()
    mock_community.close.assert_called_once()
    mock_social.close.assert_called_once()
