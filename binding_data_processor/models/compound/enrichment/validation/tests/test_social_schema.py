"""Tests for social media data schema validation.

This module tests both JSON schema validation and object schema validation
for social media data structures from Reddit, Twitter, and Bluesky.
"""

import pytest
from jsonschema import validate as validate_json_schema
from jsonschema.exceptions import ValidationError as JsonSchemaError

from ..schema import (
    SOCIAL_DATA_SCHEMA,
    RedditPostSchema,
    RedditDataSchema,
    TwitterTweetSchema,
    TwitterDataSchema,
    BlueskyPostSchema,
    BlueskyThreadSchema,
    BlueskyDataSchema,
    SocialDataSchema,
    SchemaValidationError,
)


# Common Fixtures
@pytest.fixture
def valid_classification():
    """Valid classification data fixture."""
    return {
        "label": "positive",
        "score": 0.95,
    }


@pytest.fixture
def valid_entity():
    """Valid entity data fixture."""
    return {
        "word": "caffeine",
        "entity": "COMPOUND",
        "score": 0.9,
    }


@pytest.fixture
def valid_stats():
    """Valid stats data fixture."""
    return {
        "total_posts": 1,
        "total_comments": 1,
        "unique_users": 2,
        "novel_mentions": 1,
        "sentiment_breakdown": {
            "positive": 2,
            "neutral": 0,
            "negative": 0,
        },
    }


# Reddit Fixtures
@pytest.fixture
def valid_reddit_post(valid_classification, valid_entity):
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
        "classification": valid_classification,
        "entities": [valid_entity],
        "novel_compounds": ["caffeine"],
    }


@pytest.fixture
def valid_reddit_comment(valid_classification):
    """Valid Reddit comment fixture."""
    return {
        "id": "def456",
        "body": "Interesting experience, thanks for sharing!",
        "author": "user456",
        "score": 25,
        "created_utc": "2023-01-01T01:00:00Z",
        "post_id": "abc123",
        "classification": valid_classification,
    }


@pytest.fixture
def valid_reddit_data(valid_reddit_post, valid_reddit_comment, valid_stats):
    """Valid Reddit data fixture."""
    return {
        "subreddits": [
            "Nootropics",
            "DrugNerds",
            "Psychonaut",
        ],
        "posts": [valid_reddit_post],
        "comments": [valid_reddit_comment],
        "stats": valid_stats,
    }


# Twitter Fixtures
@pytest.fixture
def valid_twitter_tweet(valid_classification, valid_entity):
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
        "classification": valid_classification,
        "entities": [valid_entity],
        "novel_compounds": ["caffeine"],
    }


@pytest.fixture
def valid_twitter_data(valid_twitter_tweet, valid_stats):
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
        "stats": valid_stats,
    }


# Bluesky Fixtures
@pytest.fixture
def valid_bluesky_post(valid_classification, valid_entity):
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
        "classification": valid_classification,
        "entities": [valid_entity],
        "novel_compounds": ["caffeine"],
    }


@pytest.fixture
def valid_bluesky_thread(valid_classification):
    """Valid Bluesky thread fixture."""
    return {
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
                "classification": valid_classification,
            }
        ],
    }


@pytest.fixture
def valid_bluesky_data(valid_bluesky_post, valid_bluesky_thread, valid_stats):
    """Valid Bluesky data fixture."""
    return {
        "posts": [valid_bluesky_post],
        "threads": [valid_bluesky_thread],
        "stats": valid_stats,
    }


@pytest.fixture
def valid_social_data(valid_reddit_data, valid_twitter_data, valid_bluesky_data):
    """Valid social data fixture."""
    return {
        "reddit": valid_reddit_data,
        "twitter": valid_twitter_data,
        "bluesky": valid_bluesky_data,
    }


# JSON Schema Tests
def test_json_schema_valid(valid_social_data):
    """Test JSON schema validation with valid data."""
    # Should not raise any errors
    validate_json_schema(instance=valid_social_data, schema=SOCIAL_DATA_SCHEMA)


def test_json_schema_missing_required():
    """Test JSON schema validation with missing required fields."""
    # Missing reddit posts
    invalid_data = {
        "reddit": {
            "stats": {
                "total_posts": 0,
                "total_comments": 0,
                "unique_users": 0,
            },
        },
        "twitter": {
            "tweets": [],
            "stats": {
                "total_tweets": 0,
                "unique_users": 0,
                "total_likes": 0,
                "total_retweets": 0,
            },
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate_json_schema(instance=invalid_data, schema=SOCIAL_DATA_SCHEMA)


def test_json_schema_invalid_types():
    """Test JSON schema validation with invalid types."""
    invalid_data = {
        "reddit": {
            "posts": [
                {
                    "id": 123,  # Should be string
                    "title": "Test Post",
                    "body": "Test body",
                    "url": "/r/test/123",
                    "author": "testuser",
                    "subreddit": "test",
                    "score": "10",  # Should be number
                    "created_utc": "2023-01-01T00:00:00Z",
                    "comments": [],
                }
            ],
            "stats": {
                "total_posts": "1",  # Should be number
                "total_comments": 0,
                "unique_users": 1,
            },
        },
        "twitter": {
            "tweets": [],
            "stats": {
                "total_tweets": 0,
                "unique_users": 0,
                "total_likes": 0,
                "total_retweets": 0,
            },
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate_json_schema(instance=invalid_data, schema=SOCIAL_DATA_SCHEMA)


# Object Schema Tests
def test_reddit_post_schema_success(valid_reddit_post):
    """Test successful Reddit post schema validation."""
    schema = RedditPostSchema()
    # Should not raise any errors
    schema.validate(valid_reddit_post)


def test_reddit_post_schema_missing_fields():
    """Test Reddit post schema with missing fields."""
    schema = RedditPostSchema()
    invalid_post = {
        "id": "123",
        # Missing required fields
    }
    
    with pytest.raises(SchemaValidationError) as exc_info:
        schema.validate(invalid_post)
    assert "Missing required field: title" in str(exc_info.value)


def test_reddit_data_schema_success(valid_reddit_data):
    """Test successful Reddit data schema validation."""
    schema = RedditDataSchema()
    # Should not raise any errors
    schema.validate(valid_reddit_data)


def test_reddit_data_schema_missing_fields():
    """Test Reddit data schema with missing fields."""
    schema = RedditDataSchema()
    invalid_data = {
        "subreddits": [],
        # Missing required fields
    }
    
    with pytest.raises(SchemaValidationError) as exc_info:
        schema.validate(invalid_data)
    assert "Missing required field: posts" in str(exc_info.value)


def test_twitter_tweet_schema_success(valid_twitter_tweet):
    """Test successful Twitter tweet schema validation."""
    schema = TwitterTweetSchema()
    # Should not raise any errors
    schema.validate(valid_twitter_tweet)


def test_twitter_tweet_schema_missing_fields():
    """Test Twitter tweet schema with missing fields."""
    schema = TwitterTweetSchema()
    invalid_tweet = {
        "id": "123",
        # Missing required fields
    }
    
    with pytest.raises(SchemaValidationError) as exc_info:
        schema.validate(invalid_tweet)
    assert "Missing required field: text" in str(exc_info.value)


def test_twitter_data_schema_success(valid_twitter_data):
    """Test successful Twitter data schema validation."""
    schema = TwitterDataSchema()
    # Should not raise any errors
    schema.validate(valid_twitter_data)


def test_twitter_data_schema_missing_fields():
    """Test Twitter data schema with missing fields."""
    schema = TwitterDataSchema()
    invalid_data = {
        "tweets": [],
        # Missing required fields
    }
    
    with pytest.raises(SchemaValidationError) as exc_info:
        schema.validate(invalid_data)
    assert "Missing required field: hashtags" in str(exc_info.value)


def test_bluesky_post_schema_success(valid_bluesky_post):
    """Test successful Bluesky post schema validation."""
    schema = BlueskyPostSchema()
    # Should not raise any errors
    schema.validate(valid_bluesky_post)


def test_bluesky_post_schema_missing_fields():
    """Test Bluesky post schema with missing fields."""
    schema = BlueskyPostSchema()
    invalid_post = {
        "uri": "123",
        # Missing required fields
    }
    
    with pytest.raises(SchemaValidationError) as exc_info:
        schema.validate(invalid_post)
    assert "Missing required field: text" in str(exc_info.value)


def test_bluesky_thread_schema_success(valid_bluesky_thread):
    """Test successful Bluesky thread schema validation."""
    schema = BlueskyThreadSchema()
    # Should not raise any errors
    schema.validate(valid_bluesky_thread)


def test_bluesky_thread_schema_missing_fields():
    """Test Bluesky thread schema with missing fields."""
    schema = BlueskyThreadSchema()
    invalid_thread = {
        "uri": "123",
        # Missing required fields
    }
    
    with pytest.raises(SchemaValidationError) as exc_info:
        schema.validate(invalid_thread)
    assert "Missing required field: replies" in str(exc_info.value)


def test_bluesky_data_schema_success(valid_bluesky_data):
    """Test successful Bluesky data schema validation."""
    schema = BlueskyDataSchema()
    # Should not raise any errors
    schema.validate(valid_bluesky_data)


def test_bluesky_data_schema_missing_fields():
    """Test Bluesky data schema with missing fields."""
    schema = BlueskyDataSchema()
    invalid_data = {
        "posts": [],
        # Missing required fields
    }
    
    with pytest.raises(SchemaValidationError) as exc_info:
        schema.validate(invalid_data)
    assert "Missing required field: threads" in str(exc_info.value)


def test_social_data_schema_success(valid_social_data):
    """Test successful social data schema validation."""
    schema = SocialDataSchema()
    # Should not raise any errors
    schema.validate(valid_social_data)


def test_social_data_schema_missing_sources():
    """Test social data schema with missing sources."""
    schema = SocialDataSchema()
    invalid_data = {
        "reddit": {},  # Missing required fields
    }
    
    with pytest.raises(SchemaValidationError) as exc_info:
        schema.validate(invalid_data)
    assert "Missing required fields in reddit data" in str(exc_info.value)
