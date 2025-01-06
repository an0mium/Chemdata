"""Tests for Reddit Crawl4AI client."""

import pytest
from unittest.mock import AsyncMock, MagicMock, Mock, patch
from pathlib import Path

from ..clients.reddit import RedditCrawl4AIClient, RedditPostData


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
def client():
    """Create test client."""
    return RedditCrawl4AIClient(
        llm_provider="test_llm",
        api_token="test_token",
    )


@pytest.mark.asyncio
async def test_search_posts_success(client, mock_crawler):
    """Test successful post search."""
    # Mock crawler response
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=["Executed script"],
        extracted_content=[
            {
                "title": ["Test Post"],
                "content": ["Test Content"],
                "author": ["TestUser"],
                "subreddit": ["nootropics"],
                "score": ["42"],
                "num_comments": ["10"],
                "compounds": ["Compound A"],
                "effects": ["Effect 1"],
                "mechanisms": ["Mechanism 1"],
                "safety": ["Safety note 1"],
            }
        ],
    )

    # Execute search
    results = await client.search_posts("test query")

    # Verify results
    assert len(results) == 1
    post = results[0]
    assert isinstance(post, RedditPostData)
    assert post.title == "Test Post"
    assert post.content == "Test Content"
    assert post.author == "TestUser"
    assert post.subreddit == "nootropics"
    assert post.score == 42
    assert post.num_comments == 10
    assert post.compounds == ["Compound A"]
    assert post.effects == ["Effect 1"]
    assert post.mechanisms == ["Mechanism 1"]
    assert post.safety_notes == ["Safety note 1"]
    assert post.confidence > 0.8  # High confidence with all fields present
    assert post.metadata["source"] == "reddit"
    assert post.metadata["screenshot"] == str(Path("test.png"))
    assert post.metadata["javascript_logs"] == ["Executed script"]


@pytest.mark.asyncio
async def test_search_posts_custom_subreddits(client, mock_crawler):
    """Test post search with custom subreddits."""
    # Mock crawler response
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=[],
        extracted_content=[{"title": ["Test Post"]}],
    )

    # Execute search with custom subreddits
    results = await client.search_posts(
        query="test",
        subreddits=["DrugNerds", "Psychonaut"],
    )

    # Verify results and search URLs
    assert len(results) == 1
    urls = [call[1]["urls"][0] for call in mock_crawler.arun.call_args_list]
    assert any("r/DrugNerds" in url for url in urls)
    assert any("r/Psychonaut" in url for url in urls)


@pytest.mark.asyncio
async def test_search_posts_partial_data(client, mock_crawler):
    """Test post search with partial data."""
    # Mock crawler response with missing fields
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=[],
        extracted_content=[
            {
                "title": ["Test Post"],
                "content": [],
                "author": [],
                "subreddit": [],
                "score": [],
                "num_comments": [],
                "compounds": [],
                "effects": [],
                "mechanisms": [],
                "safety": [],
            }
        ],
    )

    # Execute search
    results = await client.search_posts("test query")

    # Verify results
    assert len(results) == 1
    post = results[0]
    assert post.title == "Test Post"
    assert post.content is None
    assert post.author is None
    assert post.subreddit == "researchchemicals"  # Default first subreddit
    assert post.score is None
    assert post.num_comments is None
    assert post.compounds == []
    assert post.effects == []
    assert post.mechanisms == []
    assert post.safety_notes == []
    assert post.confidence < 0.5  # Lower confidence with missing fields


@pytest.mark.asyncio
async def test_search_posts_no_results(client, mock_crawler):
    """Test post search with no results."""
    # Mock crawler response with no content
    mock_crawler.arun.return_value = Mock(
        success=True,
        screenshot=Path("test.png"),
        javascript_logs=[],
        extracted_content=[{}],
    )

    # Execute search
    results = await client.search_posts("test query")

    # Verify results
    assert len(results) == 0


@pytest.mark.asyncio
async def test_search_posts_error(client, mock_crawler):
    """Test post search with error."""
    # Mock crawler failure
    mock_crawler.arun.return_value = Mock(success=False)

    # Execute search
    results = await client.search_posts("test query")

    # Verify results
    assert len(results) == 0


def test_parse_score():
    """Test score parsing."""
    client = RedditCrawl4AIClient()
    assert client._parse_score("42") == 42
    assert client._parse_score("Score: 1000") == 1000
    assert client._parse_score("Invalid") is None
    assert client._parse_score(None) is None


def test_parse_num_comments():
    """Test comment count parsing."""
    client = RedditCrawl4AIClient()
    assert client._parse_num_comments("10") == 10
    assert client._parse_num_comments("100 comments") == 100
    assert client._parse_num_comments("Invalid") is None
    assert client._parse_num_comments(None) is None


def test_extract_url():
    """Test URL extraction from title element."""
    client = RedditCrawl4AIClient()
    html = '<a href="/r/nootropics/comments/123">Title</a>'
    assert client._extract_url(html) == "https://old.reddit.com/r/nootropics/comments/123"
    assert client._extract_url("<a>No URL</a>") is None
    assert client._extract_url("Invalid HTML") is None
    assert client._extract_url(None) is None


def test_calculate_confidence():
    """Test confidence score calculation."""
    client = RedditCrawl4AIClient()

    # All fields present
    full_data = {
        "title": "Title",
        "content": "Content",
        "author": "Author",
        "subreddit": "Subreddit",
        "score": 42,
        "num_comments": 10,
        "url": "https://example.com",
        "compounds": ["Compound"],
        "effects": ["Effect"],
        "mechanisms": ["Mechanism"],
        "safety": ["Safety"],
    }
    assert client._calculate_confidence(full_data) == 1.0

    # Some fields missing
    partial_data = {
        "title": "Title",
        "content": "Content",
        "author": "Author",
    }
    assert 0.2 < client._calculate_confidence(partial_data) < 0.5

    # No fields
    empty_data = {}
    assert client._calculate_confidence(empty_data) == 0.0
