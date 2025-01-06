"""Tests for enhanced patent client with Google Patents, USPTO and Espacenet integration."""

import pytest
from unittest.mock import AsyncMock, MagicMock, Mock, patch
from pathlib import Path

from ..clients.patents import (
    EnhancedPatentClient,
    PatentData,
    PatentFamily,
    PatentCitation,
    PatentLegalStatus,
)


@pytest.fixture
def mock_crawler():
    """Mock Crawl4AI crawler."""
    with patch("crawl4ai.AsyncWebCrawler") as mock:
        crawler = MagicMock()
        crawler.arun = AsyncMock()
        crawler.return_value.__aenter__.return_value = crawler
        crawler.return_value.__aexit__.return_value = None
        mock.return_value = crawler
        yield crawler


@pytest.fixture
def mock_session():
    """Mock aiohttp session."""
    with patch("aiohttp.ClientSession") as mock:
        session = MagicMock()
        session.get = AsyncMock()
        mock.return_value.__aenter__.return_value = session
        mock.return_value.__aexit__.return_value = None
        yield session


@pytest.fixture
def client():
    """Create test client with all API keys."""
    return EnhancedPatentClient(
        llm_provider="test_llm",
        api_token="test_token",
        uspto_api_key="test_uspto",
        espacenet_api_key="test_espacenet",
    )


@pytest.mark.asyncio
async def test_search_patents_success(client, mock_crawler):
    """Test successful patent search."""
    # Mock crawler response
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=["Executed script"],
        extracted_content=[
            {
                "title": ["Test Patent"],
                "abstract": ["Test Abstract"],
                "inventors": [["John Smith"]],
                "assignee": ["Test Company"],
                "filing_date": ["2024-01-01"],
                "publication_date": ["2024-01-05"],
                "patent_number": ["US12345678"],
                "compounds": ["Compound A"],
                "effects": ["Effect 1"],
                "mechanisms": ["Mechanism 1"],
                "safety": ["Safety note 1"],
            }
        ],
    )

    # Execute search
    results = await client.search_patents("test query")

    # Verify results
    assert len(results) == 1
    patent = results[0]
    assert isinstance(patent, PatentData)
    assert patent.title == "Test Patent"
    assert patent.abstract == "Test Abstract"
    assert patent.inventors == ["John Smith"]
    assert patent.assignee == "Test Company"
    assert patent.filing_date == "2024-01-01"
    assert patent.publication_date == "2024-01-05"
    assert patent.patent_number == "US12345678"
    assert patent.compounds == ["Compound A"]
    assert patent.effects == ["Effect 1"]
    assert patent.mechanisms == ["Mechanism 1"]
    assert patent.safety_notes == ["Safety note 1"]
    assert patent.confidence > 0.8  # High confidence with all fields present
    assert "google_patents" in patent.metadata["sources"]
    assert "espacenet" in patent.metadata["sources"]
    assert "uspto" in patent.metadata["sources"]
    assert patent.metadata["screenshot"] == str(Path("test.png"))
    assert patent.metadata["javascript_logs"] == ["Executed script"]


@pytest.mark.asyncio
async def test_search_patents_with_structure(client, mock_crawler):
    """Test patent search with chemical structure."""
    # Mock crawler response
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=[],
        extracted_content=[{"title": ["Test Patent"]}],
    )

    # Execute search with structure
    results = await client.search_patents(
        query="test",
        structure="CC(=O)Oc1ccccc1C(=O)O",  # Aspirin SMILES
    )

    # Verify results and search URL
    assert len(results) == 1
    url = mock_crawler.arun.call_args[1]["urls"][0]
    assert "structure:(CC(=O)Oc1ccccc1C(=O)O)" in url


@pytest.mark.asyncio
async def test_search_patents_with_date_range(client, mock_crawler):
    """Test patent search with date range."""
    # Mock crawler response
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=[],
        extracted_content=[{"title": ["Test Patent"]}],
    )

    # Execute search with date range
    results = await client.search_patents(
        query="test",
        date_range=("2020-01-01", "2024-01-01"),
    )

    # Verify results and search URL
    assert len(results) == 1
    url = mock_crawler.arun.call_args[1]["urls"][0]
    assert "after:2020-01-01" in url
    assert "before:2024-01-01" in url


@pytest.mark.asyncio
async def test_search_patents_partial_data(client, mock_crawler):
    """Test patent search with partial data."""
    # Mock crawler response with missing fields
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=[],
        extracted_content=[
            {
                "title": ["Test Patent"],
                "abstract": [],
                "inventors": [],
                "assignee": [],
                "filing_date": [],
                "publication_date": [],
                "patent_number": [],
                "compounds": [],
                "effects": [],
                "mechanisms": [],
                "safety": [],
            }
        ],
    )

    # Execute search
    results = await client.search_patents("test query")

    # Verify results
    assert len(results) == 1
    patent = results[0]
    assert patent.title == "Test Patent"
    assert patent.abstract is None
    assert patent.inventors == []
    assert patent.assignee is None
    assert patent.filing_date is None
    assert patent.publication_date is None
    assert patent.patent_number is None
    assert patent.compounds == []
    assert patent.effects == []
    assert patent.mechanisms == []
    assert patent.safety_notes == []
    assert patent.confidence < 0.5  # Lower confidence with missing fields


@pytest.mark.asyncio
async def test_search_patents_no_results(client, mock_crawler):
    """Test patent search with no results."""
    # Mock crawler response with no content
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=[],
        extracted_content=[{}],
    )

    # Execute search
    results = await client.search_patents("test query")

    # Verify results
    assert len(results) == 0


@pytest.mark.asyncio
async def test_search_patents_error(client, mock_crawler):
    """Test patent search with error."""
    # Mock crawler failure
    mock_crawler.arun.return_value = Mock(success=False)

    # Execute search
    results = await client.search_patents("test query")

    # Verify results
    assert len(results) == 0


@pytest.mark.asyncio
async def test_espacenet_family(client, mock_session):
    """Test getting patent family from Espacenet."""
    # Mock Espacenet API response
    mock_session.get.return_value.__aenter__.return_value = Mock(
        status=200,
        json=AsyncMock(
            return_value={
                "ops:world-patent-data": {
                    "ops:patent-family": {
                        "@family-id": "12345",
                        "family-member": [
                            {
                                "publication-reference": {
                                    "document-id": {
                                        "doc-number": "US12345678",
                                        "country": "US",
                                        "date": "2024-01-01",
                                    }
                                }
                            },
                            {
                                "publication-reference": {
                                    "document-id": {
                                        "doc-number": "EP12345678",
                                        "country": "EP",
                                        "date": "2024-01-02",
                                    }
                                }
                            },
                        ],
                    }
                }
            }
        ),
    )

    # Get family data
    family = await client._get_espacenet_family("US12345678")

    # Verify family data
    assert isinstance(family, PatentFamily)
    assert family.family_id == "12345"
    assert "US12345678" in family.members
    assert "EP12345678" in family.members
    assert family.priority_date == "2024-01-01"
    assert family.countries == {"US", "EP"}
    assert family.metadata["source"] == "espacenet"


@pytest.mark.asyncio
async def test_espacenet_citations(client, mock_session):
    """Test getting citations from Espacenet."""
    # Mock Espacenet API response
    mock_session.get.return_value.__aenter__.return_value = Mock(
        status=200,
        json=AsyncMock(
            return_value={
                "ops:world-patent-data": {
                    "ops:citation-list": {
                        "citation": [
                            {
                                "patcit": {
                                    "document-id": {
                                        "doc-number": "US12345678",
                                        "date": "2024-01-01",
                                    }
                                },
                                "passage": {"text": "Test Citation"},
                                "@cited-phase": "examination",
                                "@cited-category": "X",
                            }
                        ]
                    }
                }
            }
        ),
    )

    # Get citations
    citations = await client._get_espacenet_citations("US12345678")

    # Verify citations
    assert len(citations) == 1
    citation = citations[0]
    assert isinstance(citation, PatentCitation)
    assert citation.patent_number == "US12345678"
    assert citation.title == "Test Citation"
    assert citation.filing_date == "2024-01-01"
    assert citation.relevance == 1.0  # X category = highest relevance
    assert citation.citation_type == "examination"
    assert citation.metadata["source"] == "espacenet"
    assert citation.metadata["category"] == "X"


@pytest.mark.asyncio
async def test_espacenet_legal_status(client, mock_session):
    """Test getting legal status from Espacenet."""
    # Mock Espacenet API response
    mock_session.get.return_value.__aenter__.return_value = Mock(
        status=200,
        json=AsyncMock(
            return_value={
                "ops:world-patent-data": {
                    "ops:legal-status": [
                        {
                            "@status-code": "granted",
                            "date": "2024-01-01",
                            "country": "US",
                            "text": "Patent granted",
                        }
                    ]
                }
            }
        ),
    )

    # Get legal status
    status = await client._get_espacenet_legal_status("US12345678")

    # Verify legal status
    assert isinstance(status, PatentLegalStatus)
    assert status.status == "granted"
    assert status.date == "2024-01-01"
    assert status.country == "US"
    assert status.description == "Patent granted"
    assert status.metadata["source"] == "espacenet"


def test_parse_citation_relevance():
    """Test citation relevance scoring."""
    client = EnhancedPatentClient()

    # Test all citation categories
    categories = {
        "X": 1.0,  # Particularly relevant if taken alone
        "Y": 0.8,  # Particularly relevant if combined
        "A": 0.5,  # Background art
        "O": 0.3,  # Non-written disclosure
        "P": 0.4,  # Intermediate document
        "T": 0.6,  # Theory/principle
        "E": 0.7,  # Earlier patent document
        "D": 0.4,  # Document cited in application
    }

    for category, expected_score in categories.items():
        citation = {"@cited-category": category}
        assert client._parse_citation_relevance(citation) == expected_score

    # Test unknown category
    assert client._parse_citation_relevance({"@cited-category": "Z"}) is None
    assert client._parse_citation_relevance({}) is None


def test_build_patent_url():
    """Test patent URL building."""
    client = EnhancedPatentClient()
    assert client._build_patent_url("US12345678") == "https://patents.google.com/patent/US12345678"
    assert client._build_patent_url(None) is None


def test_calculate_confidence():
    """Test confidence score calculation."""
    client = EnhancedPatentClient()

    # All fields present
    full_data = {
        "title": "Title",
        "abstract": "Abstract",
        "inventors": ["Inventor"],
        "assignee": "Company",
        "filing_date": "2024-01-01",
        "publication_date": "2024-01-05",
        "patent_number": "US12345678",
        "url": "https://example.com",
        "compounds": ["Compound"],
        "effects": ["Effect"],
        "mechanisms": ["Mechanism"],
        "safety": ["Safety"],
        "family": {"id": "12345"},
        "citations": [{"patent": "US87654321"}],
        "classifications": [{"code": "A61K"}],
    }
    assert client._calculate_confidence(full_data) == 1.0

    # Some fields missing
    partial_data = {
        "title": "Title",
        "abstract": "Abstract",
        "inventors": ["Inventor"],
    }
    assert 0.2 < client._calculate_confidence(partial_data) < 0.5

    # No fields
    empty_data = {}
    assert client._calculate_confidence(empty_data) == 0.0
