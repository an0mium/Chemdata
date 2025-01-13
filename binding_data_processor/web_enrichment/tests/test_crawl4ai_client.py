"""Tests for Crawl4AI client integration."""

import pytest
from unittest.mock import AsyncMock, MagicMock, Mock, patch
from pathlib import Path
from ..crawl4ai_client import (
    Crawl4AIClient,
    ResearchData,
    PatentData,
    CommunityData,
    ResearchSchema,
    PatentSchema,
    CommunitySchema,
    CircuitBreaker,
    CircuitBreakerConfig,
)


@pytest.fixture
def mock_crawler():
    """Mock Crawl4AI Crawler."""
    with patch("crawl4ai.AsyncWebCrawler") as _:
        crawler = MagicMock()
        crawler.arun = AsyncMock()
        crawler.return_value.__aenter__.return_value = crawler
        crawler.return_value.__aexit__.return_value = None
        crawler.arun.return_value = Mock(
            success=True,
            extracted_content=[
                {
                    "title": "Test Title",
                    "abstract": "Test Abstract",
                    "authors": ["Author 1"],
                    "compounds": ["Compound A"],
                    "effects": ["Effect 1"],
                    "mechanisms": ["Mechanism 1"],
                    "safety_notes": ["Safety note 1"],
                    "confidence": 0.95,
                }
            ],
        )
        yield crawler


@pytest.fixture
def mock_sync_crawler():
    """Mock synchronous Crawl4AI Crawler for legacy support."""
    with patch("crawl4ai.Crawler") as mock:
        crawler = MagicMock()
        crawler.crawl = AsyncMock()
        mock.return_value = crawler
        yield crawler


@pytest.fixture
def mock_llm():
    """Mock LLM extraction strategy."""
    with patch("crawl4ai.extraction_strategy.LLMExtractionStrategy") as mock:
        llm = MagicMock()
        mock.return_value = llm
        yield llm


@pytest.fixture
def mock_circuit_breaker():
    """Mock circuit breaker."""
    with patch("..pipeline.infrastructure.circuit_breaker.CircuitBreaker") as mock:
        breaker = MagicMock()
        mock.return_value = breaker
        yield breaker


@pytest.fixture
def mock_validator():
    """Create a mock validator."""
    with patch("..validation.schema.DataValidator") as mock:
        mock.validate.return_value = Mock(is_valid=True)
        yield mock


@pytest.fixture
def client():
    """Create test client."""
    return Crawl4AIClient(
        api_key="test_key",
        llm_provider="ollama/llama2",  # Use local model for testing
        api_token="no_token",
    )


@pytest.fixture
def client_with_circuit_config():
    """Create client with custom circuit breaker config."""
    config = CircuitBreakerConfig(
        error_threshold=3,
        recovery_timeout=60,
        half_open_timeout=30,
    )
    return Crawl4AIClient(circuit_config=config)


@pytest.mark.asyncio
async def test_scrape_research_with_css(client, mock_crawler):
    """Test scraping research papers with CSS extraction."""
    # Mock crawl results
    mock_crawler.arun.return_value = MagicMock(
        extracted_content=[
            {
                "title": "Test Paper",
                "abstract": "Test abstract",
                "authors": ["Author 1", "Author 2"],
                "journal": "Test Journal",
                "year": 2023,
                "doi": "10.1234/test",
                "full_text": "Full paper text",
                "references": ["ref1", "ref2"],
                "citations": ["cite1", "cite2"],
                "metadata": {"source": "scholar"},
            }
        ]
    )

    # Execute scrape
    results = await client.scrape_research("test query")

    # Verify results
    assert len(results) == 1
    paper = results[0]
    assert isinstance(paper, ResearchData)
    assert paper.title == "Test Paper"
    assert paper.abstract == "Test abstract"
    assert paper.authors == ["Author 1", "Author 2"]
    assert paper.journal == "Test Journal"
    assert paper.year == 2023
    assert paper.doi == "10.1234/test"
    assert paper.full_text == "Full paper text"
    assert paper.references == ["ref1", "ref2"]
    assert paper.citations == ["cite1", "cite2"]
    assert paper.metadata == {"source": "scholar"}


@pytest.mark.asyncio
async def test_scrape_research_with_llm_fallback(client, mock_crawler, mock_llm):
    """Test LLM fallback when CSS extraction fails."""
    # Mock CSS extraction failure then LLM success
    mock_crawler.arun.side_effect = [
        MagicMock(extracted_content=None),  # CSS fails
        MagicMock(
            extracted_content=[
                {  # LLM succeeds
                    "title": "LLM Paper",
                    "abstract": "LLM abstract",
                    "authors": ["LLM Author"],
                    "metadata": {"source": "llm"},
                }
            ]
        ),
    ]

    # Execute scrape
    results = await client.scrape_research("test query")

    # Verify LLM results
    assert len(results) == 1
    paper = results[0]
    assert paper.title == "LLM Paper"
    assert paper.abstract == "LLM abstract"
    assert paper.authors == ["LLM Author"]
    assert paper.metadata == {"source": "llm"}


@pytest.mark.asyncio
async def test_scrape_patents(client, mock_crawler):
    """Test scraping patents."""
    # Mock crawl results
    mock_crawler.arun.return_value = MagicMock(
        extracted_content=[
            {
                "title": "Test Patent",
                "abstract": "Test abstract",
                "inventors": ["Inventor 1", "Inventor 2"],
                "assignee": "Test Company",
                "filing_date": "2023-01-01",
                "publication_date": "2023-06-01",
                "patent_number": "US123456",
                "claims": ["claim1", "claim2"],
                "description": "Full description",
                "metadata": {"source": "google"},
            }
        ]
    )

    # Execute scrape
    results = await client.scrape_patents("test query")

    # Verify results
    assert len(results) == 1
    patent = results[0]
    assert isinstance(patent, PatentData)
    assert patent.title == "Test Patent"
    assert patent.abstract == "Test abstract"
    assert patent.inventors == ["Inventor 1", "Inventor 2"]
    assert patent.assignee == "Test Company"
    assert patent.filing_date == "2023-01-01"
    assert patent.publication_date == "2023-06-01"
    assert patent.patent_number == "US123456"
    assert patent.claims == ["claim1", "claim2"]
    assert patent.description == "Full description"
    assert patent.metadata == {"source": "google"}


@pytest.mark.asyncio
async def test_scrape_community(client, mock_crawler):
    """Test scraping community forums."""
    # Mock crawl results
    mock_crawler.arun.return_value = MagicMock(
        extracted_content=[
            {
                "source": "reddit",
                "title": "Test Post",
                "content": "Test content",
                "author": "user123",
                "date": "2023-06-01",
                "url": "https://reddit.com/r/test/123",
                "replies": [{"text": "reply1"}, {"text": "reply2"}],
                "metadata": {"subreddit": "test"},
            }
        ]
    )

    # Execute scrape
    results = await client.scrape_community("test query")

    # Verify results
    assert len(results) == 1
    post = results[0]
    assert isinstance(post, CommunityData)
    assert post.source == "reddit"
    assert post.title == "Test Post"
    assert post.content == "Test content"
    assert post.author == "user123"
    assert post.date == "2023-06-01"
    assert post.url == "https://reddit.com/r/test/123"
    assert post.replies == [{"text": "reply1"}, {"text": "reply2"}]
    assert post.metadata == {"subreddit": "test"}


@pytest.mark.asyncio
async def test_scrape_all(client, mock_crawler):
    """Test scraping all sources."""
    results = await client.scrape_all("psychoactive compounds")
    assert isinstance(results, dict)
    assert "research" in results
    assert "patents" in results
    assert "community" in results
    assert isinstance(results["research"], list)
    assert isinstance(results["patents"], list)
    assert isinstance(results["community"], list)


def test_research_schema():
    """Test research data schema."""
    schema = ResearchSchema(
        title="Test Paper",
        abstract="Test abstract",
        authors=["Author 1", "Author 2"],
        journal="Test Journal",
        year=2024,
        doi="10.1234/test",
        compounds=["Compound A", "Compound B"],
        effects=["Effect 1", "Effect 2"],
        mechanisms=["Mechanism 1"],
        safety_notes=["Safety note 1"],
        confidence=0.95,
    )
    assert schema.title == "Test Paper"
    assert len(schema.authors) == 2
    assert len(schema.compounds) == 2
    assert schema.confidence > 0.9


def test_patent_schema():
    """Test patent data schema."""
    schema = PatentSchema(
        patent_number="US123456",
        title="Test Patent",
        abstract="Test abstract",
        inventors=["Inventor 1"],
        assignee="Test Company",
        filing_date="2024-01-01",
        compounds=["Compound A"],
        synthesis_methods=["Method 1"],
        uses=["Use 1"],
        claims=["Claim 1"],
        confidence=0.9,
    )
    assert schema.patent_number == "US123456"
    assert schema.assignee == "Test Company"
    assert len(schema.compounds) == 1
    assert schema.confidence > 0.8


def test_community_schema():
    """Test community data schema."""
    schema = CommunitySchema(
        source="Test Forum",
        thread_title="Test Thread",
        post_content="Test content",
        timestamp="2024-01-01T12:00:00",
        compounds=["Compound A"],
        effects=["Effect 1"],
        dosage_info="10mg",
        safety_notes=["Safety note 1"],
        confidence=0.85,
    )
    assert schema.source == "Test Forum"
    assert schema.dosage_info == "10mg"
    assert len(schema.compounds) == 1
    assert schema.confidence > 0.8


@pytest.mark.asyncio
async def test_error_handling(client):
    """Test error handling in scraping functions."""
    with pytest.raises(Exception):
        await client.scrape_research("")  # Empty query should raise error

    with pytest.raises(Exception):
        await client.scrape_patents("")  # Empty query should raise error

    with pytest.raises(Exception):
        await client.scrape_community("")  # Empty query should raise error


@pytest.mark.asyncio
async def test_circuit_breaker_custom_config(client_with_circuit_config, mock_crawler):
    """Test circuit breaker with custom configuration."""
    # Mock failures up to threshold
    mock_crawler.arun.side_effect = [Exception("Failed")] * 3

    # Execute scrapes until circuit opens
    with pytest.raises(CircuitBreaker.OpenError):
        for _ in range(3):  # Custom threshold of 3
            try:
                await client_with_circuit_config.scrape_research("test query")
            except Exception:
                continue

    # Verify circuit opened after 3 errors
    assert mock_crawler.arun.call_count == 3


@pytest.mark.asyncio
async def test_circuit_breaker_open(client, mock_circuit_breaker):
    """Test circuit breaker prevents requests when open."""
    # Mock open circuit
    mock_circuit_breaker.__enter__.side_effect = CircuitBreaker.OpenError()

    # Verify request is blocked
    with pytest.raises(CircuitBreaker.OpenError):
        await client.scrape_research("test query")


@pytest.mark.asyncio
async def test_circuit_breaker_half_open(client, mock_circuit_breaker, mock_crawler):
    """Test circuit breaker allows test requests when half-open."""
    # Mock successful test request
    mock_crawler.arun.return_value = MagicMock(
        extracted_content=[{"title": "Test"}],
    )

    # Execute scrape
    results = await client.scrape_research("test query")

    # Verify circuit closed after success
    assert len(results) == 1
    mock_circuit_breaker.__exit__.assert_called_once()


@pytest.mark.asyncio
async def test_screenshot_capture(client, mock_crawler):
    """Test screenshot capture functionality."""
    # Mock successful screenshot
    screenshot_path = Path("screenshots/test.png")
    mock_crawler.arun.return_value = MagicMock(
        screenshot=screenshot_path,
        extracted_content=[{"title": "Test"}],
    )

    # Execute scrape
    results = await client.scrape_research("test query")

    # Verify screenshot was captured
    assert mock_crawler.arun.call_args[1]["config"].screenshot.enabled
    assert mock_crawler.arun.call_args[1]["config"].screenshot.full_page
    assert results[0].metadata.get("screenshot") == str(screenshot_path)


@pytest.mark.asyncio
async def test_javascript_execution(client, mock_crawler):
    """Test JavaScript execution functionality."""
    # Mock successful JS execution
    mock_crawler.arun.return_value = MagicMock(
        javascript_logs=["Executed script"],
        extracted_content=[{"title": "Test"}],
    )

    # Execute scrape
    results = await client.scrape_research("test query")

    # Verify JS config
    config = mock_crawler.arun.call_args[1]["config"]
    assert config.javascript.enabled
    assert config.javascript.wait_for_network
    assert config.javascript.wait_for_selectors
    assert results[0].metadata.get("javascript_logs") == ["Executed script"]


@pytest.mark.asyncio
async def test_media_extraction(client):
    """Test media extraction functionality."""
    with patch.object(client, "request") as mock_request:
        mock_request.return_value = {
            "images": ["image1.jpg", "image2.jpg"],
            "videos": ["video1.mp4"],
        }

        result = client.extract_media("https://example.com")
        assert "images" in result
        assert "videos" in result
        assert len(result["images"]) == 2
        assert len(result["videos"]) == 1


@pytest.mark.asyncio
async def test_link_extraction(client):
    """Test link extraction functionality."""
    with patch.object(client, "request") as mock_request:
        mock_request.return_value = [
            {"url": "https://example.com/page1", "text": "Link 1"},
            {"url": "https://example.com/page2", "text": "Link 2"},
        ]

        result = client.extract_links("https://example.com", depth=2)
        assert len(result) == 2
        assert result[0]["url"] == "https://example.com/page1"


@pytest.mark.asyncio
async def test_metadata_extraction(client):
    """Test metadata extraction functionality."""
    with patch.object(client, "request") as mock_request:
        mock_request.return_value = {
            "title": "Test Page",
            "description": "Test Description",
            "schema_org": {"@type": "Article"},
        }

        result = client.extract_metadata("https://example.com")
        assert result["title"] == "Test Page"
        assert "schema_org" in result


@pytest.mark.asyncio
async def test_legacy_crawler_fallback(client, mock_sync_crawler):
    """Test fallback to synchronous crawler if async not available."""
    # Mock async crawler import error
    with patch("crawl4ai.AsyncWebCrawler", side_effect=ImportError):
        # Mock successful sync crawl
        mock_sync_crawler.crawl.return_value = [{"title": "Test"}]

        # Execute crawl
        results = await client._execute_crawl("test", {"urls": [], "selectors": {}})

        # Verify sync crawler was used
        assert results == [{"title": "Test"}]
        mock_sync_crawler.crawl.assert_called_once()


def test_get_headers(client):
    """Test header generation."""
    headers = client._get_headers()
    assert "Authorization" in headers
    assert headers["Authorization"] == "Bearer test_key"
    assert headers["Content-Type"] == "application/json"


@pytest.mark.asyncio
async def test_hooks(client, caplog):
    """Test request lifecycle hooks."""
    url = "https://example.com"

    await client._before_request_hook(url)
    assert "Starting request" in caplog.text

    await client._after_request_hook(url, {"status": "success"})
    assert "Completed request" in caplog.text

    await client._on_error_hook(url, Exception("Test error"))
    assert "Error scraping" in caplog.text


def test_timestamp_format(client):
    """Test timestamp generation."""
    timestamp = client._get_timestamp()
    assert isinstance(timestamp, str)
    assert "T" in timestamp  # ISO format contains T between date and time
