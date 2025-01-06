"""Tests for enhanced ScienceDirect client."""

import pytest
from unittest.mock import AsyncMock, MagicMock, Mock, patch
from datetime import datetime
from pathlib import Path
from bs4 import BeautifulSoup

from ..clients.sciencedirect import EnhancedScienceDirectClient, ScienceDirectResearchData
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
                        ["J Smith", "A Jones"],
                        ["B Wilson", "C Brown"],
                    ],
                    "journal": [
                        "Journal of Psychopharmacology",
                        "Neuropharmacology",
                    ],
                    "year": [
                        "2024",
                        "2023",
                    ],
                    "citations": [
                        "42",
                        "28",
                    ],
                    "doi": [
                        "10.1234/test.1",
                        "10.1234/test.2",
                    ],
                    "pii": [
                        "S123456",
                        "S789012",
                    ],
                    "compounds": ["caffeine", "adenosine"],
                    "effects": ["improved memory", "increased focus"],
                    "mechanisms": ["adenosine receptor antagonism"],
                    "safety": ["well tolerated", "mild side effects"],
                    "url": [
                        "/science/article/pii/S123456",
                        "/science/article/pii/S789012",
                    ],
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
    return EnhancedScienceDirectClient(
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
    assert isinstance(article, ScienceDirectResearchData)
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
    assert article.metadata["source"] == "sciencedirect"
    assert article.metadata["url"] == "https://www.sciencedirect.com/science/article/pii/S123456"
    assert "extraction_date" in article.metadata
    assert "screenshot" in article.metadata
    assert "javascript_logs" in article.metadata


@pytest.mark.asyncio
async def test_search_articles_pagination(client, mock_crawler):
    """Test article search pagination."""
    # Mock multiple pages
    mock_crawler.arun.side_effect = [
        MagicMock(
            success=True,
            extracted_content=[
                {
                    "title": [f"Article {i}" for i in range(25)],
                    "abstract": [""] * 25,
                    "authors": [[]] * 25,
                    "pii": [f"S{i}" for i in range(25)],
                    "url": [f"/article/{i}" for i in range(25)],
                }
            ],
            screenshot=Path("test1.png"),
            javascript_logs=["Page 1 loaded"],
        )
        for _ in range(3)
    ]

    # Search articles with pagination
    articles = await client.search_articles(query="test", max_results=60)

    # Verify pagination
    assert len(articles) == 60
    assert articles[0].title == "Article 0"
    assert articles[25].title == "Article 0"
    assert articles[50].title == "Article 0"


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
                "journal": [],
                "year": [],
                "citations": [],
                "doi": [],
                "pii": [],
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
    mock_crawler.arun.side_effect = Exception("API error")

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
    assert isinstance(article, ScienceDirectResearchData)
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
    mock_crawler.arun.side_effect = Exception("API error")

    # Second search should use cache
    articles2 = await client.get_compound_articles(test_compound)
    assert len(articles2) > 0
    assert articles2[0].title == articles1[0].title


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
    assert "sciencedirect.com/search" in url
    assert "qs=caffeine" in url

    # URL with filters
    url = client._build_search_url(
        query="caffeine",
        from_year=2020,
        to_year=2024,
        sort="date",
        start=25,
    )
    assert "qs=caffeine" in url
    assert "date=2020" in url
    assert "end=2024" in url
    assert "sortBy=date" in url
    assert "offset=25" in url


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
    article = ScienceDirectResearchData(
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
    irrelevant = ScienceDirectResearchData(
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
    assert client._parse_citations("42") == 42
    assert client._parse_citations("Cited by 42") == 42
    assert client._parse_citations("42 citations") == 42
    assert client._parse_citations("No citations") is None
    assert client._parse_citations(None) is None


def test_extract_url(client):
    """Test URL extraction."""
    # Test relative URL
    assert (
        client._extract_url("/science/article/pii/S123456")
        == "https://www.sciencedirect.com/science/article/pii/S123456"
    )

    # Test absolute URL
    assert client._extract_url("https://www.sciencedirect.com/article") == "https://www.sciencedirect.com/article"

    # Test invalid URLs
    assert client._extract_url("<a>No URL</a>") is None
    assert client._extract_url("Invalid HTML") is None
    assert client._extract_url(None) is None


def test_extract_doi(client):
    """Test DOI extraction."""
    # Test DOI formats
    assert client._extract_doi("10.1234/test") == "10.1234/test"
    assert client._extract_doi("doi.org/10.1234/test") == "10.1234/test"

    # Test PII extraction
    assert client._extract_doi("/pii/S123456") == "S123456"
    assert client._extract_doi("pii/S123456") == "S123456"

    # Test invalid inputs
    assert client._extract_doi("Invalid") is None
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
        "pii": "S123456",
        "url": "/article/123",
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
    """Test ScienceDirectResearchData schema."""
    data = ScienceDirectResearchData(
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
            "source": "sciencedirect",
            "url": "https://www.sciencedirect.com/science/article/pii/S123456",
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
    assert data.metadata["source"] == "sciencedirect"
    assert data.metadata["url"] == "https://www.sciencedirect.com/science/article/pii/S123456"
    assert "extraction_date" in data.metadata
    assert data.metadata["screenshot"] == "test.png"
    assert data.metadata["javascript_logs"] == ["test log"]
