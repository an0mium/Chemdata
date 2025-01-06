"""Tests for social media data validation."""

import pytest

from ..data import validate_social_data
from ..validation import (
    validate_reddit_data,
    validate_twitter_data,
    validate_bluesky_data,
    ValidationError,
)


@pytest.fixture
def valid_reddit_post():
    """Valid Reddit post fixture."""
    return {
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
    }


@pytest.fixture
def valid_reddit_data(valid_reddit_post):
    """Valid Reddit data fixture."""
    return {
        "subreddits": [
            "Nootropics",
            "DrugNerds",
            "Psychonaut",
        ],
        "posts": [valid_reddit_post],
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
            }
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
    }


@pytest.fixture
def valid_twitter_tweet():
    """Valid Twitter tweet fixture."""
    return {
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
    }


@pytest.fixture
def valid_twitter_data(valid_twitter_tweet):
    """Valid Twitter data fixture."""
    return {
        "tweets": [valid_twitter_tweet],
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
    }


@pytest.fixture
def valid_bluesky_post():
    """Valid Bluesky post fixture."""
    return {
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
    }


@pytest.fixture
def valid_bluesky_data(valid_bluesky_post):
    """Valid Bluesky data fixture."""
    return {
        "posts": [valid_bluesky_post],
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
                    }
                ],
            }
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
    }


@pytest.fixture
def valid_social_data(valid_reddit_data, valid_twitter_data, valid_bluesky_data):
    """Valid social data fixture."""
    return {
        "reddit": valid_reddit_data,
        "twitter": valid_twitter_data,
        "bluesky": valid_bluesky_data,
    }


def test_validate_social_data_success(valid_social_data):
    """Test successful validation of social data."""
    # Should not raise any errors
    validate_social_data(valid_social_data)


def test_validate_social_data_missing_sources():
    """Test validation with missing data sources."""
    invalid_data = {
        "reddit": {},  # Missing required fields
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_social_data(invalid_data)
    assert "Missing required fields in reddit data" in str(exc_info.value)


def test_validate_reddit_data_success(valid_reddit_data):
    """Test successful Reddit data validation."""
    # Should not raise any errors
    validate_reddit_data(valid_reddit_data)


def test_validate_reddit_data_missing_fields():
    """Test Reddit data validation with missing fields."""
    invalid_data = {
        "subreddits": [],
        # Missing required fields
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_reddit_data(invalid_data)
    assert "Missing required field: posts" in str(exc_info.value)


def test_validate_reddit_data_invalid_types():
    """Test Reddit data validation with invalid types."""
    invalid_data = {
        "subreddits": "not a list",  # Should be list
        "posts": "not a list",  # Should be list
        "comments": "not a list",  # Should be list
        "stats": "not an object",  # Should be dict
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_reddit_data(invalid_data)
    assert "Invalid type for field: subreddits" in str(exc_info.value)


def test_validate_reddit_data_invalid_post():
    """Test Reddit data validation with invalid post."""
    invalid_data = {
        "subreddits": ["Nootropics"],
        "posts": [
            {
                "id": 123,  # Should be string
                "score": "not a number",  # Should be number
                # Missing required fields
            }
        ],
        "comments": [],
        "stats": {"total_posts": 1},
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_reddit_data(invalid_data)
    assert "Invalid post data" in str(exc_info.value)


def test_validate_twitter_data_success(valid_twitter_data):
    """Test successful Twitter data validation."""
    # Should not raise any errors
    validate_twitter_data(valid_twitter_data)


def test_validate_twitter_data_missing_fields():
    """Test Twitter data validation with missing fields."""
    invalid_data = {
        "tweets": [],
        # Missing required fields
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_twitter_data(invalid_data)
    assert "Missing required field: hashtags" in str(exc_info.value)


def test_validate_twitter_data_invalid_types():
    """Test Twitter data validation with invalid types."""
    invalid_data = {
        "tweets": "not a list",  # Should be list
        "hashtags": "not a list",  # Should be list
        "mentions": "not a list",  # Should be list
        "stats": "not an object",  # Should be dict
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_twitter_data(invalid_data)
    assert "Invalid type for field: tweets" in str(exc_info.value)


def test_validate_twitter_data_invalid_tweet():
    """Test Twitter data validation with invalid tweet."""
    invalid_data = {
        "tweets": [
            {
                "id": 123,  # Should be string
                "like_count": "not a number",  # Should be number
                # Missing required fields
            }
        ],
        "hashtags": [],
        "mentions": [],
        "stats": {"total_tweets": 1},
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_twitter_data(invalid_data)
    assert "Invalid tweet data" in str(exc_info.value)


def test_validate_bluesky_data_success(valid_bluesky_data):
    """Test successful Bluesky data validation."""
    # Should not raise any errors
    validate_bluesky_data(valid_bluesky_data)


def test_validate_bluesky_data_missing_fields():
    """Test Bluesky data validation with missing fields."""
    invalid_data = {
        "posts": [],
        # Missing required fields
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_bluesky_data(invalid_data)
    assert "Missing required field: threads" in str(exc_info.value)


def test_validate_bluesky_data_invalid_types():
    """Test Bluesky data validation with invalid types."""
    invalid_data = {
        "posts": "not a list",  # Should be list
        "threads": "not a list",  # Should be list
        "stats": "not an object",  # Should be dict
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_bluesky_data(invalid_data)
    assert "Invalid type for field: posts" in str(exc_info.value)


def test_validate_bluesky_data_invalid_post():
    """Test Bluesky data validation with invalid post."""
    invalid_data = {
        "posts": [
            {
                "uri": 123,  # Should be string
                "likeCount": "not a number",  # Should be number
                # Missing required fields
            }
        ],
        "threads": [],
        "stats": {"total_posts": 1},
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_bluesky_data(invalid_data)
    assert "Invalid post data" in str(exc_info.value)


def test_validate_bluesky_data_invalid_thread():
    """Test Bluesky data validation with invalid thread."""
    invalid_data = {
        "posts": [],
        "threads": [
            {
                "uri": 123,  # Should be string
                "replies": "not a list",  # Should be list
            }
        ],
        "stats": {"total_posts": 1},
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_bluesky_data(invalid_data)
    assert "Invalid thread data" in str(exc_info.value)
