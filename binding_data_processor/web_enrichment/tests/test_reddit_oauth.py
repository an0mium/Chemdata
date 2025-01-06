"""Tests for Reddit OAuth client."""

import pytest
from datetime import datetime, timedelta
from unittest.mock import AsyncMock, patch

from ..clients.reddit_oauth import RedditOAuthClient, RedditAPIConfig, RedditPostData


@pytest.fixture
def mock_config():
    """Create mock Reddit API config."""
    return RedditAPIConfig(
        client_id="test_client_id",
        client_secret="test_client_secret",
        user_agent="test_user_agent",
        username="test_user",
        password="test_pass",
    )


@pytest.fixture
def mock_post_data():
    """Create mock Reddit post data."""
    return {
        "data": {
            "title": "Test Compound Analysis",
            "selftext": "Analysis of test compound shows interesting effects",
            "author": "test_user",
            "subreddit": "DrugNerds",
            "score": 100,
            "num_comments": 10,
            "permalink": "/r/DrugNerds/test_post",
            "created_utc": datetime.now().timestamp(),
            "upvote_ratio": 0.95,
            "is_self": True,
            "domain": "self.DrugNerds",
        }
    }


@pytest.fixture
async def mock_client(mock_config):
    """Create mock Reddit OAuth client."""
    with patch("praw.Reddit"), patch("aiohttp.ClientSession"):
        client = RedditOAuthClient(config=mock_config)
        client.access_token = "test_token"
        client.token_expiry = datetime.now() + timedelta(hours=1)
        yield client


@pytest.mark.asyncio
async def test_search_posts_api(mock_client, mock_post_data):
    """Test searching posts via API."""
    # Mock API response
    mock_response = AsyncMock()
    mock_response.status = 200
    mock_response.json.return_value = {
        "data": {
            "children": [mock_post_data],
        }
    }

    # Mock session get
    mock_client.session = AsyncMock()
    mock_client.session.get.return_value.__aenter__.return_value = mock_response

    # Mock content analysis
    with patch("..llm_utils.analyze_text") as mock_analyze:
        mock_analyze.return_value = {
            "compounds": ["Test Compound"],
            "effects": ["Test Effect"],
            "mechanisms": ["Test Mechanism"],
            "safety": ["Test Safety Note"],
        }

        results = await mock_client._search_posts_api(
            query="test",
            subreddits=["DrugNerds"],
            max_results=1,
        )

    assert len(results) == 1
    result = results[0]
    assert isinstance(result, RedditPostData)
    assert result.title == mock_post_data["data"]["title"]
    assert result.compounds == ["Test Compound"]
    assert result.effects == ["Test Effect"]
    assert result.confidence > 0


@pytest.mark.asyncio
async def test_search_posts_fallback(mock_client):
    """Test fallback to scraping when API fails."""
    # Mock API failure
    mock_client._search_posts_api = AsyncMock(side_effect=Exception("API Error"))

    # Mock scraper success
    mock_client.scraper.search_posts = AsyncMock(
        return_value=[
            RedditPostData(
                title="Test Post",
                content="Test Content",
                author="test_user",
                subreddit="DrugNerds",
                score=10,
                num_comments=5,
                compounds=["Test Compound"],
                effects=["Test Effect"],
                mechanisms=["Test Mechanism"],
                safety_notes=["Test Safety"],
                confidence=0.8,
                metadata={},
            )
        ]
    )

    results = await mock_client.search_posts(
        query="test",
        subreddits=["DrugNerds"],
        max_results=1,
    )

    assert len(results) == 1
    assert results[0].title == "Test Post"
    assert results[0].compounds == ["Test Compound"]


@pytest.mark.asyncio
async def test_monitor_subreddits(mock_client, mock_post_data):
    """Test subreddit monitoring."""
    # Mock getting new posts
    mock_client._get_new_posts = AsyncMock(return_value=[mock_post_data["data"]])

    # Mock content analysis
    with patch("..llm_utils.analyze_text") as mock_analyze:
        mock_analyze.return_value = {
            "compounds": ["Test Compound"],
            "safety": ["Test Safety Note"],
        }

        # Start monitoring (will run indefinitely, so we'll break after first iteration)
        with patch("asyncio.sleep", AsyncMock()) as mock_sleep:
            mock_sleep.side_effect = [None, Exception("Stop monitoring")]
            with pytest.raises(Exception, match="Stop monitoring"):
                await mock_client.monitor_subreddits(
                    subreddits=["DrugNerds"],
                    interval=1,
                )

    # Verify monitoring behavior
    assert mock_client._get_new_posts.called
    assert mock_analyze.called


@pytest.mark.asyncio
async def test_token_refresh(mock_client):
    """Test access token refresh."""
    # Mock expired token
    mock_client.access_token = "old_token"
    mock_client.token_expiry = datetime.now() - timedelta(hours=1)

    # Mock token refresh response
    mock_response = AsyncMock()
    mock_response.status = 200
    mock_response.json.return_value = {
        "access_token": "new_token",
        "expires_in": 3600,
        "refresh_token": "new_refresh_token",
    }

    # Mock session post
    mock_client.session = AsyncMock()
    mock_client.session.post.return_value.__aenter__.return_value = mock_response

    await mock_client.refresh_token()

    assert mock_client.access_token == "new_token"
    assert mock_client.config.refresh_token == "new_refresh_token"
    assert mock_client.token_expiry > datetime.now()


@pytest.mark.asyncio
async def test_error_handling(mock_client):
    """Test error handling in various scenarios."""
    # Test API error handling
    mock_response = AsyncMock()
    mock_response.status = 403
    mock_response.text = AsyncMock(return_value="Forbidden")

    mock_client.session = AsyncMock()
    mock_client.session.get.return_value.__aenter__.return_value = mock_response

    results = await mock_client._search_posts_api(
        query="test",
        subreddits=["DrugNerds"],
        max_results=1,
    )
    assert len(results) == 0

    # Test token refresh error handling
    mock_client.session.post.return_value.__aenter__.return_value = mock_response
    with pytest.raises(Exception, match="Token refresh failed"):
        await mock_client.refresh_token()

    # Test content analysis error handling
    with patch("..llm_utils.analyze_text", side_effect=Exception("LLM Error")):
        result = await mock_client._analyze_content("test content")
        assert result == {
            "compounds": [],
            "effects": [],
            "mechanisms": [],
            "safety": [],
        }


def test_confidence_calculation(mock_client, mock_post_data):
    """Test confidence score calculation."""
    # Test high quality post
    score = mock_client._calculate_confidence(mock_post_data["data"])
    assert 0.8 <= score <= 1.0

    # Test low quality post
    low_quality = {
        "title": "Test",
        "score": 1,
        "num_comments": 0,
        "selftext": "",
        "upvote_ratio": 0.5,
    }
    score = mock_client._calculate_confidence(low_quality)
    assert 0 <= score <= 0.3
