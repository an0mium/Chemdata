"""Tests for social media data enrichment client.

This module tests the SocialClient which integrates with social media APIs
and ML models for sentiment analysis and entity recognition.
"""

from unittest.mock import Mock, patch

import pytest

from ..social import SocialClient
from ..base import WebClientError


@pytest.fixture
def mock_http_client():
    """Mock HTTP client fixture."""
    return Mock()


@pytest.fixture
def mock_logger():
    """Mock logger fixture."""
    return Mock()


@pytest.fixture
def mock_text_classifier():
    """Mock text classifier fixture."""
    mock = Mock()
    mock.return_value = [{
        "label": "positive",
        "score": 0.95
    }]
    return mock


@pytest.fixture
def mock_ner_model():
    """Mock NER model fixture."""
    mock = Mock()
    mock.return_value = [{
        "word": "caffeine",
        "entity": "COMPOUND",
        "score": 0.9
    }]
    return mock


@pytest.fixture
def mock_reddit_response():
    """Mock Reddit API response fixture."""
    return {
        "data": {
            "subreddits": [
                "Nootropics",
                "DrugNerds",
                "Psychonaut",
            ],
            "posts": [
                {
                    "id": "abc123",
                    "title": "Experience Report: Caffeine + L-Theanine",
                    "body": "Detailed experience report...",
                    "url": "/r/nootropics/abc123",
                    "author": "user123",
                    "subreddit": "nootropics",
                    "score": 150,
                    "created_utc": "2023-01-01T00:00:00Z",
                    "num_comments": 45,
                    "classification": {
                        "label": "positive",
                        "score": 0.95,
                    },
                    "entities": [
                        {
                            "word": "caffeine",
                            "entity": "COMPOUND",
                            "score": 0.9,
                        },
                    ],
                    "novel_compounds": ["caffeine"],
                },
            ],
            "comments": [
                {
                    "id": "def456",
                    "body": "Interesting experience, thanks for sharing!",
                    "author": "user456",
                    "score": 25,
                    "created_utc": "2023-01-01T01:00:00Z",
                    "post_id": "abc123",
                    "classification": {
                        "label": "positive",
                        "score": 0.9,
                    },
                },
            ],
            "stats": {
                "total_posts": 1,
                "total_comments": 1,
                "unique_users": 2,
                "novel_mentions": 1,
                "sentiment_breakdown": {
                    "positive": 2,
                    "neutral": 0,
                    "negative": 0,
                },
            },
        },
        "metadata": {
            "query_time": 0.15,
            "api_version": "2023.1",
        },
    }


@pytest.fixture
def mock_twitter_response():
    """Mock Twitter API response fixture."""
    return {
        "data": {
            "tweets": [
                {
                    "id": "123456789",
                    "text": "Just tried caffeine + l-theanine stack. Great focus!",
                    "author": "nootropics_user",
                    "created_at": "2023-01-01T00:00:00Z",
                    "retweet_count": 10,
                    "like_count": 50,
                    "reply_count": 5,
                    "quote_count": 2,
                    "url": "https://twitter.com/nootropics_user/status/123456789",
                    "classification": {
                        "label": "positive",
                        "score": 0.95,
                    },
                    "entities": [
                        {
                            "word": "caffeine",
                            "entity": "COMPOUND",
                            "score": 0.9,
                        },
                    ],
                    "novel_compounds": ["caffeine"],
                },
            ],
            "hashtags": [
                "#Nootropics",
                "#BrainHealth",
                "#Productivity",
            ],
            "mentions": [
                "@nootropics_expert",
                "@science_account",
            ],
            "stats": {
                "total_tweets": 1,
                "unique_users": 1,
                "total_likes": 50,
                "total_retweets": 10,
                "novel_mentions": 1,
                "sentiment_breakdown": {
                    "positive": 1,
                    "neutral": 0,
                    "negative": 0,
                },
            },
        },
        "metadata": {
            "query_time": 0.12,
            "api_version": "2023.1",
        },
    }


@pytest.fixture
def mock_bluesky_response():
    """Mock Bluesky API response fixture."""
    return {
        "data": {
            "posts": [
                {
                    "uri": "at://user.bsky.social/posts/123",
                    "text": "Exploring the cognitive effects of caffeine",
                    "author": {
                        "did": "did:plc:123",
                        "handle": "researcher.bsky.social",
                    },
                    "createdAt": "2023-01-01T00:00:00Z",
                    "likeCount": 30,
                    "repostCount": 5,
                    "replyCount": 8,
                    "classification": {
                        "label": "neutral",
                        "score": 0.8,
                    },
                    "entities": [
                        {
                            "word": "caffeine",
                            "entity": "COMPOUND",
                            "score": 0.9,
                        },
                    ],
                    "novel_compounds": ["caffeine"],
                },
            ],
            "threads": [
                {
                    "uri": "at://user.bsky.social/posts/123",
                    "replies": [
                        {
                            "uri": "at://user.bsky.social/posts/456",
                            "text": "Interesting findings!",
                            "author": {
                                "did": "did:plc:456",
                                "handle": "scientist.bsky.social",
                            },
                            "createdAt": "2023-01-01T01:00:00Z",
                            "classification": {
                                "label": "positive",
                                "score": 0.85,
                            },
                        },
                    ],
                },
            ],
            "stats": {
                "total_posts": 1,
                "total_replies": 1,
                "unique_users": 2,
                "total_likes": 30,
                "total_reposts": 5,
                "novel_mentions": 1,
                "sentiment_breakdown": {
                    "positive": 1,
                    "neutral": 1,
                    "negative": 0,
                },
            },
        },
        "metadata": {
            "query_time": 0.18,
            "api_version": "2023.1",
        },
    }


@pytest.fixture
def social_client(
    mock_http_client,
    mock_logger,
    mock_text_classifier,
    mock_ner_model,
    tmp_path,
):
    """Social client fixture with mocked dependencies."""
    with patch("transformers.pipeline") as mock_pipeline:
        mock_pipeline.return_value = mock_text_classifier
        
        client = SocialClient(
            http_client=mock_http_client,
            logger=mock_logger,
            model_dir=tmp_path,
            cache_dir=tmp_path,
        )
        
        # Override NER model
        client.ner_model = mock_ner_model
        
        return client


def test_init(mock_http_client, mock_logger):
    """Test client initialization."""
    client = SocialClient(
        http_client=mock_http_client,
        logger=mock_logger,
    )
    
    assert client.name == "social"
    assert client.http == mock_http_client
    assert client.logger == mock_logger
    assert client.reddit_url == "https://api.reddit.com/v1"
    assert client.twitter_url == "https://api.twitter.com/v2"
    assert client.bluesky_url == "https://api.bsky.app/v1"


def test_init_defaults():
    """Test initialization with default values."""
    client = SocialClient()
    
    assert client.name == "social"
    assert client.http is not None
    assert client.logger is not None
    assert client.text_classifier is None
    assert client.ner_model is None


def test_get_reddit_data_success(
    social_client,
    mock_http_client,
    mock_reddit_response,
):
    """Test successful Reddit data retrieval."""
    mock_http_client.get.return_value = mock_reddit_response
    
    response = social_client.get_reddit_data(
        compound_name="caffeine",
        subreddits=["Nootropics", "DrugNerds"],
        time_range="month",
    )
    
    assert response == mock_reddit_response
    mock_http_client.get.assert_called_with(
        url=f"{social_client.reddit_url}/search",
        params={
            "q": "caffeine",
            "subreddit": "Nootropics+DrugNerds",
            "t": "month",
            "sort": "relevance",
            "limit": 100,
        },
        headers={"Accept": "application/json"},
    )


def test_get_reddit_data_error(social_client, mock_http_client):
    """Test Reddit data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "API error",
        status_code=500,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        social_client.get_reddit_data(
            compound_name="invalid",
            subreddits=["invalid"],
            time_range="invalid",
        )
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_get_twitter_data_success(
    social_client,
    mock_http_client,
    mock_twitter_response,
):
    """Test successful Twitter data retrieval."""
    mock_http_client.get.return_value = mock_twitter_response
    
    response = social_client.get_twitter_data(
        compound_name="caffeine",
        hashtags=["#Nootropics"],
        time_range="month",
    )
    
    assert response == mock_twitter_response
    mock_http_client.get.assert_called_with(
        url=f"{social_client.twitter_url}/search",
        params={
            "q": "caffeine #Nootropics",
            "time": "month",
            "sort": "relevance",
            "limit": 100,
        },
        headers={"Accept": "application/json"},
    )


def test_get_twitter_data_error(social_client, mock_http_client):
    """Test Twitter data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "API error",
        status_code=500,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        social_client.get_twitter_data(
            compound_name="invalid",
            hashtags=["invalid"],
            time_range="invalid",
        )
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_get_bluesky_data_success(
    social_client,
    mock_http_client,
    mock_bluesky_response,
):
    """Test successful Bluesky data retrieval."""
    mock_http_client.get.return_value = mock_bluesky_response
    
    response = social_client.get_bluesky_data(
        compound_name="caffeine",
        time_range="month",
    )
    
    assert response == mock_bluesky_response
    mock_http_client.get.assert_called_with(
        url=f"{social_client.bluesky_url}/search",
        params={
            "q": "caffeine",
            "time": "month",
            "sort": "relevance",
            "limit": 100,
        },
        headers={"Accept": "application/json"},
    )


def test_get_bluesky_data_error(social_client, mock_http_client):
    """Test Bluesky data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "API error",
        status_code=500,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        social_client.get_bluesky_data(
            compound_name="invalid",
            time_range="invalid",
        )
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_get_all_data_success(
    social_client,
    mock_http_client,
    mock_reddit_response,
    mock_twitter_response,
    mock_bluesky_response,
):
    """Test successful retrieval of all social data."""
    mock_http_client.get.side_effect = [
        mock_reddit_response,
        mock_twitter_response,
        mock_bluesky_response,
    ]
    
    response = social_client.get_all_data(
        compound_name="caffeine",
        subreddits=["Nootropics", "DrugNerds"],
        hashtags=["#Nootropics"],
        time_range="month",
    )
    
    assert response == {
        "reddit": mock_reddit_response,
        "twitter": mock_twitter_response,
        "bluesky": mock_bluesky_response,
    }


def test_get_all_data_partial_success(
    social_client,
    mock_http_client,
    mock_reddit_response,
    mock_twitter_response,
):
    """Test partial success in retrieving social data."""
    mock_http_client.get.side_effect = [
        mock_reddit_response,
        mock_twitter_response,
        WebClientError("API error", status_code=500),
    ]
    
    response = social_client.get_all_data(
        compound_name="caffeine",
        subreddits=["Nootropics", "DrugNerds"],
        hashtags=["#Nootropics"],
        time_range="month",
    )
    
    assert response == {
        "reddit": mock_reddit_response,
        "twitter": mock_twitter_response,
        "bluesky": None,
    }


def test_analyze_text_success(social_client):
    """Test successful text analysis."""
    text = "Caffeine provides great focus and energy"
    
    classification = social_client.analyze_text(text)
    
    assert classification["label"] == "positive"
    assert classification["score"] == 0.95


def test_analyze_text_no_classifier(social_client):
    """Test text analysis with no classifier."""
    social_client.text_classifier = None
    text = "Test text"
    
    classification = social_client.analyze_text(text)
    
    assert classification is None


def test_extract_entities_success(social_client):
    """Test successful entity extraction."""
    text = "Caffeine is a stimulant"
    
    entities = social_client.extract_entities(text)
    
    assert len(entities) == 1
    assert entities[0]["word"] == "caffeine"
    assert entities[0]["entity"] == "COMPOUND"
    assert entities[0]["score"] == 0.9


def test_extract_entities_no_model(social_client):
    """Test entity extraction with no model."""
    social_client.ner_model = None
    text = "Test text"
    
    entities = social_client.extract_entities(text)
    
    assert entities == []


def test_validate_compound_name_success(social_client):
    """Test successful compound name validation."""
    valid_name = "caffeine"
    
    # Should not raise any errors
    social_client._validate_compound_name(valid_name)


def test_validate_compound_name_empty(social_client):
    """Test compound name validation with empty string."""
    with pytest.raises(ValueError) as exc_info:
        social_client._validate_compound_name("")
    
    assert "Compound name cannot be empty" in str(exc_info.value)


def test_validate_compound_name_invalid_type(social_client):
    """Test compound name validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        social_client._validate_compound_name(123)
    
    assert "Compound name must be a string" in str(exc_info.value)


def test_validate_subreddits_success(social_client):
    """Test successful subreddits validation."""
    valid_subreddits = ["Nootropics", "DrugNerds"]
    
    # Should not raise any errors
    social_client._validate_subreddits(valid_subreddits)


def test_validate_subreddits_empty(social_client):
    """Test subreddits validation with empty list."""
    with pytest.raises(ValueError) as exc_info:
        social_client._validate_subreddits([])
    
    assert "Subreddits list cannot be empty" in str(exc_info.value)


def test_validate_subreddits_invalid_type(social_client):
    """Test subreddits validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        social_client._validate_subreddits("not a list")
    
    assert "Subreddits must be a list" in str(exc_info.value)


def test_validate_hashtags_success(social_client):
    """Test successful hashtags validation."""
    valid_hashtags = ["#Nootropics", "#BrainHealth"]
    
    # Should not raise any errors
    social_client._validate_hashtags(valid_hashtags)


def test_validate_hashtags_empty(social_client):
    """Test hashtags validation with empty list."""
    with pytest.raises(ValueError) as exc_info:
        social_client._validate_hashtags([])
    
    assert "Hashtags list cannot be empty" in str(exc_info.value)


def test_validate_hashtags_invalid_type(social_client):
    """Test hashtags validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        social_client._validate_hashtags("not a list")
    
    assert "Hashtags must be a list" in str(exc_info.value)


def test_validate_hashtags_invalid_format(social_client):
    """Test hashtags validation with invalid format."""
    with pytest.raises(ValueError) as exc_info:
        social_client._validate_hashtags(["invalid"])
    
    assert "Hashtags must start with #" in str(exc_info.value)


def test_validate_time_range_success(social_client):
    """Test successful time range validation."""
    valid_ranges = ["hour", "day", "week", "month", "year", "all"]
    
    for time_range in valid_ranges:
        # Should not raise any errors
        social_client._validate_time_range(time_range)


def test_validate_time_range_invalid(social_client):
    """Test time range validation with invalid value."""
    with pytest.raises(ValueError) as exc_info:
        social_client._validate_time_range("invalid")
    
    assert "Invalid time range" in str(exc_info.value)


def test_validate_time_range_invalid_type(social_client):
    """Test time range validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        social_client._validate_time_range(123)
    
    assert "Time range must be a string" in str(exc_info.value)


def test_validate_reddit_data(social_client, mock_reddit_response):
    """Test Reddit data validation."""
    # Should not raise any errors
    social_client._validate_reddit_data(mock_reddit_response)
    
    # Missing required fields
    invalid_data = {"data": {"subreddits": []}}
    with pytest.raises(WebClientError) as exc_info:
        social_client._validate_reddit_data(invalid_data)
    assert "Missing required fields" in str(exc_info.value)
    
    # Invalid types
    invalid_data = {
        "data": {
            "subreddits": "not a list",  # Should be list
            "posts": "not a list",  # Should be list
            "comments": "not a list",  # Should be list
        }
    }
    with pytest.raises(WebClientError) as exc_info:
        social_client._validate_reddit_data(invalid_data)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_twitter_data(social_client, mock_twitter_response):
    """Test Twitter data validation."""
    # Should not raise any errors
    social_client._validate_twitter_data(mock_twitter_response)
    
    # Missing required fields
    invalid_data = {"data": {"tweets": []}}
    with pytest.raises(WebClientError) as exc_info:
        social_client._validate_twitter_data(invalid_data)
    assert "Missing required fields" in str(exc_info.value)
    
    # Invalid types
    invalid_data = {
        "data": {
            "tweets": "not a list",  # Should be list
            "hashtags": "not a list",  # Should be list
            "mentions": "not a list",  # Should be list
        }
    }
    with pytest.raises(WebClientError) as exc_info:
        social_client._validate_twitter_data(invalid_data)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_bluesky_data(social_client, mock_bluesky_response):
    """Test Bluesky data validation."""
    # Should not raise any errors
    social_client._validate_bluesky_data(mock_bluesky_response)
    
    # Missing required fields
    invalid_data = {"data": {"posts": []}}
    with pytest.raises(WebClientError) as exc_info:
        social_client._validate_bluesky_data(invalid_data)
    assert "Missing required fields" in str(exc_info.value)
    
    # Invalid types
    invalid_data = {
        "data": {
            "posts": "not a list",  # Should be list
            "threads": "not a list",  # Should be list
        }
    }
    with pytest.raises(WebClientError) as exc_info:
        social_client._validate_bluesky_data(invalid_data)
    assert "Invalid field types" in str(exc_info.value)
