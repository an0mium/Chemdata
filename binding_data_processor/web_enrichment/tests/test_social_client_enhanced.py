"""Tests for enhanced social media client."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

import responses

from ..social_client_enhanced import SocialClientEnhanced
from ...models.compound import Compound
from ...pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig


@pytest.fixture
def client(tmp_path):
    """Create test client."""
    return SocialClientEnhanced(
        reddit_client_id="test_id",
        reddit_client_secret="test_secret",
        twitter_bearer_token="test_token",
        cache_dir=tmp_path / "cache",
        circuit_config=CircuitBreakerConfig(
            failure_threshold=2,
            failure_timeout=1,
            reset_timeout=1,
        ),
    )


@pytest.fixture
def compound():
    """Create test compound."""
    return Compound(
        name="Test Compound",
        smiles="CC1=CC=CC=C1",
        cas_number="100-00-0",
    )


@pytest.fixture
def mock_reddit_post():
    """Create mock Reddit post."""
    post = Mock()
    post.id = "test_post"
    post.title = "Test Post"
    post.selftext = "Test compound experience report"
    post.permalink = "/r/researchchemicals/comments/test"
    post.author = "test_user"
    post.subreddit = Mock(display_name="researchchemicals")
    post.score = 10
    post.created_utc = datetime.now().timestamp()
    post.comments = Mock()
    post.comments.list = lambda: []
    return post


@pytest.fixture
def mock_tweet():
    """Create mock tweet."""
    tweet = Mock()
    tweet.id = "123456"
    tweet.text = "Test compound research"
    tweet.author = Mock(username="test_user")
    tweet.created_at = datetime.now()
    tweet.public_metrics = {"retweet_count": 5, "like_count": 10}
    tweet.entities = {"hashtags": [{"tag": "research"}]}
    return tweet


@patch("praw.Reddit")
def test_reddit_search(mock_reddit, client, compound, mock_reddit_post):
    """Test Reddit search functionality."""
    # Mock subreddit search
    mock_subreddit = Mock()
    mock_subreddit.search.return_value = [mock_reddit_post]
    mock_reddit.return_value.subreddit.return_value = mock_subreddit

    # Get Reddit data
    data = client._get_reddit_data(compound.name)

    # Check data
    assert data is not None
    assert len(data["posts"]) == 1
    assert data["posts"][0]["id"] == "test_post"
    assert data["posts"][0]["title"] == "Test Post"
    assert data["subreddits"]["researchchemicals"] == 1


@patch("tweepy.Paginator")
def test_twitter_search(mock_paginator, client, compound, mock_tweet):
    """Test Twitter search functionality."""
    # Mock tweet search
    mock_paginator.return_value.flatten.return_value = [mock_tweet]

    # Get Twitter data
    data = client._get_twitter_data(compound.name)

    # Check data
    assert data is not None
    assert len(data["tweets"]) == 1
    assert data["tweets"][0]["id"] == "123456"
    assert data["tweets"][0]["text"] == "Test compound research"
    assert data["users"]["test_user"] == 1
    assert data["hashtags"]["research"] == 1


def test_get_cached_reddit_data(client, compound, tmp_path):
    """Test getting cached Reddit data."""
    # Create cache file
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_file = cache_dir / f"reddit_{compound.name}.json"

    data = {
        "posts": [
            {
                "id": "test_post",
                "title": "Test Post",
                "text": "Test text",
            }
        ],
        "comments": [],
        "subreddits": {"researchchemicals": 1},
    }
    with cache_file.open("w") as f:
        json.dump(data, f)

    # Get cached data
    cached = client._get_cached_reddit_data(compound.name)

    # Check cached data
    assert cached is not None
    assert len(cached["posts"]) == 1
    assert cached["posts"][0]["id"] == "test_post"
    assert cached["subreddits"]["researchchemicals"] == 1


def test_get_cached_twitter_data(client, compound, tmp_path):
    """Test getting cached Twitter data."""
    # Create cache file
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_file = cache_dir / f"twitter_{compound.name}.json"

    data = {
        "tweets": [
            {
                "id": "123456",
                "text": "Test tweet",
                "author": "test_user",
            }
        ],
        "users": {"test_user": 1},
        "hashtags": {"research": 1},
    }
    with cache_file.open("w") as f:
        json.dump(data, f)

    # Get cached data
    cached = client._get_cached_twitter_data(compound.name)

    # Check cached data
    assert cached is not None
    assert len(cached["tweets"]) == 1
    assert cached["tweets"][0]["id"] == "123456"
    assert cached["users"]["test_user"] == 1


@patch("praw.Reddit")
@patch("tweepy.Paginator")
def test_process_compounds(
    mock_paginator,
    mock_reddit,
    client,
    compound,
    mock_reddit_post,
    mock_tweet,
):
    """Test processing compounds."""
    # Mock Reddit search
    mock_subreddit = Mock()
    mock_subreddit.search.return_value = [mock_reddit_post]
    mock_reddit.return_value.subreddit.return_value = mock_subreddit

    # Mock Twitter search
    mock_paginator.return_value.flatten.return_value = [mock_tweet]

    # Process compound
    client.process_compounds([compound])

    # Check results
    assert compound.social_data is not None
    assert "reddit" in compound.social_data
    assert "twitter" in compound.social_data
    assert len(client.processed_compounds) == 1
    assert len(client.failed_compounds) == 0


def test_metrics(client, compound):
    """Test metrics collection."""
    # Process successful compound
    client.processed_compounds.append(compound.name)
    client.source_stats["reddit"]["success"] += 1
    client.source_stats["reddit"]["posts"] += 2
    client.source_stats["twitter"]["success"] += 1
    client.source_stats["twitter"]["tweets"] += 3

    # Process failed compound
    client.failed_compounds.append("Failed Compound")

    # Get metrics
    metrics = client.get_metrics()

    # Check metrics
    assert metrics["processed_compounds"] == 1
    assert metrics["failed_compounds"] == 1
    assert metrics["success_rate"] == 0.5
    assert metrics["source_stats"]["reddit"]["success"] == 1
    assert metrics["source_stats"]["reddit"]["posts"] == 2
    assert metrics["source_stats"]["twitter"]["success"] == 1
    assert metrics["source_stats"]["twitter"]["tweets"] == 3


@patch("praw.Reddit")
@patch("tweepy.Paginator")
def test_error_handling(mock_paginator, mock_reddit, client, compound):
    """Test error handling."""
    # Mock failing endpoints
    mock_reddit.side_effect = Exception("Reddit API error")
    mock_paginator.side_effect = Exception("Twitter API error")

    # Process compound
    client.process_compounds([compound])

    # Check error handling
    assert len(client.processed_compounds) == 0
    assert len(client.failed_compounds) == 1
    assert client.failed_compounds[0] == compound.name
    assert client.source_stats["reddit"]["failure"] == 1
    assert client.source_stats["twitter"]["failure"] == 1
