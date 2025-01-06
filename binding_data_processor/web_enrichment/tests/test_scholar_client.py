"""Tests for enhanced Google Scholar client."""

import pytest
from unittest.mock import AsyncMock, MagicMock, patch
from datetime import datetime, timedelta
from pathlib import Path

from ..clients.scholar import EnhancedScholarClient, ScholarResearchData, SessionManager
from ...models.compound import CompoundData


@pytest.fixture
def mock_crawler():
    """Mock AsyncWebCrawler."""
    with patch("crawl4ai.AsyncWebCrawler") as mock:
        crawler = AsyncMock()
        crawler.arun.return_value = MagicMock(
            success=True,
            extracted_content=[
                {
                    "title": [
                        "Effects of caffeine on cognitive performance",
                        "Adenosine receptor antagonism by caffeine",
                    ],
                    "abstract": [
                        "Study examining caffeine's effects on memory and attention...",
                        "Investigation of caffeine's mechanism as an adenosine antagonist...",
                    ],
                    "authors": [
                        "J Smith, A Jones - Journal of Psychopharmacology, 2024",
                        "B Wilson, C Brown - Neuropharmacology, 2023",
                    ],
                    "citations": [
                        "Cited by 42",
                        "Cited by 28",
                    ],
                    "doi": [
                        "https://doi.org/10.1234/test.1",
                        "https://doi.org/10.1234/test.2",
                    ],
                    "compounds": ["caffeine", "adenosine"],
                    "effects": ["improved memory", "increased focus"],
                    "mechanisms": ["adenosine receptor antagonism"],
                    "safety": ["well tolerated", "mild side effects"],
                    "url": ["https://example.com/paper1", "https://example.com/paper2"],
                }
            ],
            screenshot=Path("test.png"),
            javascript_logs=["Page loaded successfully"],
        )
        mock.return_value.__aenter__.return_value = crawler
        mock.return_value.__aexit__.return_value = None
        yield crawler


@pytest.fixture
def mock_validator():
    """Mock validator."""
    with patch("..validation.enhanced.EnhancedValidator") as mock:
        validator = AsyncMock()
        validator.validate.return_value = MagicMock(is_valid=True)
        mock.return_value = validator
        yield validator


@pytest.fixture
def client(mock_crawler, mock_validator):
    """Create test client."""
    return EnhancedScholarClient(
        validator=mock_validator,
        llm_provider="test_llm",
        api_token="test_token",
        cache_dir=None,
        proxy_enabled=True,
        proxy_rotation=True,
        proxy_retry_count=3,
        rate_limit=10,
        rate_limit_delay=60,
    )


@pytest.fixture
def test_compound():
    """Create test compound."""
    return CompoundData(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
        iupac_name="1,3,7-trimethylpurine-2,6-dione",
        common_name_1="caffeine",
        common_name_2="trimethylxanthine",
        common_name_3="N/A",
    )


@pytest.mark.asyncio
async def test_session_manager():
    """Test session management."""
    manager = SessionManager(
        max_retries=2,
        backoff_factor=1.0,
        max_concurrent=2,
        rate_limit=2,
        rate_window=1,
    )

    # Test successful execution
    mock_func = AsyncMock(return_value="success")
    success, result = await manager.execute(mock_func)
    assert success
    assert result == "success"

    # Test rate limiting
    success1, _ = await manager.execute(mock_func)
    success2, _ = await manager.execute(mock_func)
    assert success1 and success2
    assert len(manager.request_times) == 2

    # Test retry on failure
    mock_func = AsyncMock(side_effect=[Exception("error"), "success"])
    success, result = await manager.execute(mock_func)
    assert success
    assert result == "success"

    # Test max retries
    mock_func = AsyncMock(side_effect=Exception("error"))
    success, result = await manager.execute(mock_func)
    assert not success
    assert result is None


@pytest.mark.asyncio
async def test_search_articles_success(client, mock_crawler):
    """Test successful article search."""
    # Search articles
    articles = await client.search_articles(
        query="caffeine",
        max_results=10,
        from_year=2020,
        to_year=2024,
        sort="relevance",
    )

    # Verify results
    assert len(articles) == 2
    article = articles[0]
    assert isinstance(article, ScholarResearchData)
    assert article.title == "Effects of caffeine on cognitive performance"
    assert "caffeine" in article.abstract
    assert len(article.authors) == 2
    assert article.journal == "Journal of Psychopharmacology"
    assert article.year == 2024
    assert article.doi == "10.1234/test.1"
    assert article.citations == 42
    assert len(article.compounds) == 2
    assert len(article.effects) == 2
    assert len(article.mechanisms) == 1
    assert len(article.safety_notes) == 2
    assert article.confidence > 0.5
    assert article.metadata["source"] == "google_scholar"
    assert "extraction_date" in article.metadata
    assert "screenshot" in article.metadata
    assert "javascript_logs" in article.metadata
    assert "url" in article.metadata


@pytest.mark.asyncio
async def test_search_articles_pagination(client, mock_crawler):
    """Test article search pagination."""
    # Mock multiple pages
    mock_crawler.arun.side_effect = [
        # First page
        MagicMock(
            success=True,
            extracted_content=[
                {
                    "title": ["Article 1"],
                    "abstract": ["Abstract 1"],
                    "authors": ["Author 1 - Journal 1, 2024"],
                    "citations": ["Cited by 10"],
                    "compounds": ["compound1"],
                    "effects": ["effect1"],
                    "mechanisms": ["mechanism1"],
                    "url": ["https://example.com/1"],
                }
            ],
        ),
        # Second page
        MagicMock(
            success=True,
            extracted_content=[
                {
                    "title": ["Article 2"],
                    "abstract": ["Abstract 2"],
                    "authors": ["Author 2 - Journal 2, 2023"],
                    "citations": ["Cited by 20"],
                    "compounds": ["compound2"],
                    "effects": ["effect2"],
                    "mechanisms": ["mechanism2"],
                    "url": ["https://example.com/2"],
                }
            ],
        ),
    ]

    # Search articles
    articles = await client.search_articles("test", max_results=30)

    # Verify results
    assert len(articles) == 2
    assert articles[0].title == "Article 1"
    assert articles[1].title == "Article 2"
    assert articles[0].metadata["url"] == "https://example.com/1"
    assert articles[1].metadata["url"] == "https://example.com/2"


@pytest.mark.asyncio
async def test_search_articles_partial_data(client, mock_crawler):
    """Test article search with partial data."""
    # Mock crawler response with missing fields
    mock_crawler.arun.return_value = MagicMock(
        success=True,
        extracted_content=[
            {
                "title": ["Test Title"],
                "abstract": [],
                "authors": [],
                "citations": [],
                "compounds": [],
                "effects": [],
                "mechanisms": [],
                "safety": [],
                "url": [],
            }
        ],
        screenshot=Path("test.png"),
        javascript_logs=[],
    )

    # Search articles
    articles = await client.search_articles("test")

    # Verify results
    assert len(articles) == 1
    article = articles[0]
    assert article.title == "Test Title"
    assert article.abstract is None
    assert article.authors == []
    assert article.journal is None
    assert article.year is None
    assert article.citations is None
    assert article.compounds == []
    assert article.effects == []
    assert article.mechanisms == []
    assert article.safety_notes == []
    assert article.confidence < 0.5  # Lower confidence with missing fields


@pytest.mark.asyncio
async def test_search_articles_caching(client, mock_crawler):
    """Test article search caching."""
    # First search
    articles1 = await client.search_articles("caffeine")
    assert len(articles1) > 0

    # Mock crawler failure
    mock_crawler.arun.side_effect = Exception("Crawl error")

    # Second search should use cache
    articles2 = await client.search_articles("caffeine")
    assert len(articles2) > 0
    assert articles2[0].title == articles1[0].title


@pytest.mark.asyncio
async def test_get_compound_articles_success(client, mock_crawler, test_compound):
    """Test successful compound article search."""
    # Get articles
    articles = await client.get_compound_articles(
        compound=test_compound,
        max_results=10,
    )

    # Verify results
    assert len(articles) > 0
    article = articles[0]
    assert isinstance(article, ScholarResearchData)
    assert article.title == "Effects of caffeine on cognitive performance"
    assert "caffeine" in article.compounds
    assert len(article.effects) > 0
    assert len(article.mechanisms) > 0
    assert article.confidence > 0.5


@pytest.mark.asyncio
async def test_get_compound_articles_caching(client, mock_crawler, test_compound):
    """Test compound article search caching."""
    # First search
    articles1 = await client.get_compound_articles(test_compound)
    assert len(articles1) > 0

    # Mock crawler failure
    mock_crawler.arun.side_effect = Exception("Crawl error")

    # Second search should use cache
    articles2 = await client.get_compound_articles(test_compound)
    assert len(articles2) > 0
    assert articles2[0].title == articles1[0].title


@pytest.mark.asyncio
async def test_track_citations(client, mock_crawler):
    """Test citation tracking."""
    # Create test article
    article = ScholarResearchData(
        title="Test Article",
        abstract="Test abstract",
        authors=["Test Author"],
        journal="Test Journal",
        year=2024,
        citations=10,
        metadata={
            "url": "https://example.com/article",
            "extraction_date": (datetime.now() - timedelta(days=10)).isoformat(),
        },
    )

    # Mock current citations
    mock_crawler.arun.return_value = MagicMock(
        success=True,
        extracted_content=[{"citations": ["Cited by 15"]}],
    )

    # Track citations
    data = await client.track_citations(article)

    # Verify results
    assert data["initial_citations"] == 10
    assert data["current_citations"] == 15
    assert data["citation_change"] == 5
    assert data["citation_velocity"] == 0.5  # 5 citations / 10 days
    assert data["tracking_period_days"] == 10
    assert "last_updated" in data


@pytest.mark.asyncio
async def test_error_handling(client, mock_crawler):
    """Test error handling."""
    # Mock crawler failure
    mock_crawler.arun.side_effect = Exception("Crawl error")

    # Test search error handling
    articles = await client.search_articles("caffeine")
    assert articles == []

    # Test compound search error handling
    articles = await client.get_compound_articles(test_compound())
    assert articles == []


@pytest.mark.asyncio
async def test_validation(client, mock_crawler, mock_validator):
    """Test data validation."""
    # Mock validation failure
    mock_validator.validate.return_value = MagicMock(is_valid=False)

    # Test search validation
    articles = await client.search_articles("caffeine")
    assert len(articles) == 0

    # Test compound search validation
    articles = await client.get_compound_articles(test_compound())
    assert len(articles) == 0


def test_build_search_url(client):
    """Test search URL building."""
    # Basic URL
    url = client._build_search_url("caffeine")
    assert "scholar.google.com/scholar" in url
    assert "q=caffeine" in url

    # URL with filters
    url = client._build_search_url(
        query="caffeine",
        from_year=2020,
        to_year=2024,
        sort="date",
        start=10,
    )
    assert "q=caffeine" in url
    assert "as_ylo=2020" in url
    assert "as_yhi=2024" in url
    assert "scisbd=1" in url
    assert "start=10" in url


def test_build_compound_query(client, test_compound):
    """Test compound query building."""
    query = client._build_compound_query(test_compound)
    assert '"Caffeine"' in query
    assert '"58-08-2"' in query
    assert '"1,3,7-trimethylpurine-2,6-dione"' in query
    assert '"caffeine"' in query
    assert '"trimethylxanthine"' in query
    assert "N/A" not in query


@pytest.mark.asyncio
async def test_relevance_check(client, test_compound):
    """Test article relevance checking."""
    # Create test article
    article = ScholarResearchData(
        title="Test Article",
        abstract="Study of caffeine effects",
        authors=["Author 1"],
        journal="Test Journal",
        year=2024,
        compounds=["caffeine"],
        effects=["improved memory"],
        mechanisms=["adenosine antagonism"],
    )

    # Check relevance
    assert await client._is_relevant(article, test_compound)

    # Test irrelevant article
    irrelevant = ScholarResearchData(
        title="Unrelated Article",
        abstract="Study of glucose metabolism",
        authors=["Author 2"],
        journal="Test Journal",
        year=2024,
        compounds=[],
        effects=[],
        mechanisms=[],
    )
    assert not await client._is_relevant(irrelevant, test_compound)


def test_parse_year(client):
    """Test year parsing."""
    assert client._parse_year("2024") == 2024
    assert client._parse_year("Published in 2024") == 2024
    assert client._parse_year("2024 - Journal") == 2024
    assert client._parse_year("Invalid") is None
    assert client._parse_year(None) is None


def test_parse_citations(client):
    """Test citation parsing."""
    assert client._parse_citations("Cited by 42") == 42
    assert client._parse_citations("Cited by 1000") == 1000
    assert client._parse_citations("42 citations") is None
    assert client._parse_citations("No citations") is None
    assert client._parse_citations(None) is None


def test_extract_url(client):
    """Test URL extraction."""
    html = '<a href="https://example.com">Title</a>'
    assert client._extract_url(html) == "https://example.com"
    assert client._extract_url("<a>No URL</a>") is None
    assert client._extract_url("Invalid HTML") is None
    assert client._extract_url(None) is None


def test_extract_doi(client):
    """Test DOI extraction."""
    assert client._extract_doi("https://doi.org/10.1234/test") == "10.1234/test"
    assert client._extract_doi("doi:10.1234/test") == "10.1234/test"
    assert client._extract_doi("https://example.com") is None
    assert client._extract_doi(None) is None


def test_calculate_confidence(client):
    """Test confidence scoring."""
    # Full data
    full_data = {
        "title": "Title",
        "abstract": "Abstract",
        "authors": ["Author"],
        "journal": "Journal",
        "year": 2024,
        "doi": "10.1234/test",
        "compounds": ["compound"],
        "effects": ["effect"],
        "mechanisms": ["mechanism"],
        "safety": ["safe"],
    }
    assert client._calculate_confidence(full_data) == 1.0

    # Partial data
    partial_data = {
        "title": "Title",
        "abstract": "Abstract",
        "authors": ["Author"],
    }
    assert 0.2 < client._calculate_confidence(partial_data) < 0.5

    # No data
    empty_data = {}
    assert client._calculate_confidence(empty_data) == 0.0


def test_research_data_schema():
    """Test ScholarResearchData schema."""
    data = ScholarResearchData(
        title="Test Title",
        abstract="Test Abstract",
        authors=["Test Author"],
        journal="Test Journal",
        year=2024,
        doi="10.1234/test",
        citations=42,
        compounds=["Test Compound"],
        effects=["Test Effect"],
        mechanisms=["Test Mechanism"],
        safety_notes=["Test Safety"],
        binding_data=[{"type": "Ki", "value": 1.0}],
        confidence=0.9,
        metadata={
            "source": "test",
            "url": "https://example.com",
            "extraction_date": datetime.now().isoformat(),
            "screenshot": "test.png",
            "javascript_logs": ["test log"],
        },
    )

    assert data.title == "Test Title"
    assert data.abstract == "Test Abstract"
    assert data.authors == ["Test Author"]
    assert data.journal == "Test Journal"
    assert data.year == 2024
    assert data.doi == "10.1234/test"
    assert data.citations == 42
    assert data.compounds == ["Test Compound"]
    assert data.effects == ["Test Effect"]
    assert data.mechanisms == ["Test Mechanism"]
    assert data.safety_notes == ["Test Safety"]
    assert data.binding_data == [{"type": "Ki", "value": 1.0}]
    assert data.confidence == 0.9
    assert data.metadata["source"] == "test"
    assert data.metadata["url"] == "https://example.com"
    assert "extraction_date" in data.metadata
    assert data.metadata["screenshot"] == "test.png"
    assert data.metadata["javascript_logs"] == ["test log"]
