"""Shared test fixtures for web enrichment data validation."""

import pytest


@pytest.fixture
def valid_swiss_data():
    """Valid Swiss data fixture."""
    return {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 0.95,
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": {
            "mw": 194.19,
            "logp": -0.07,
            "hbd": 0,
            "hba": 6,
            "tpsa": 58.44,
        },
    }


@pytest.fixture
def valid_community_data():
    """Valid community data fixture."""
    return {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": ["Stimulation"],
            "dosage": {"common": "100-200 mg"},
            "duration": {"onset": "15-45 minutes"},
        },
        "erowid": {
            "name": "Caffeine",
            "url": "https://erowid.org/chemicals/caffeine/",
            "experiences": [],
            "health": [],
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": {"dose_common": "100-200mg"},
            "interactions": [],
        },
    }


@pytest.fixture
def valid_social_data():
    """Valid social data fixture."""
    return {
        "reddit": {
            "posts": [
                {
                    "id": "123",
                    "title": "Test Post",
                    "body": "Test caffeine effects",
                    "url": "/r/nootropics/123",
                    "author": "testuser",
                    "subreddit": "nootropics",
                    "score": 10,
                    "created_utc": "2023-01-01T00:00:00Z",
                    "comments": [],
                }
            ],
            "stats": {
                "total_posts": 1,
                "total_comments": 0,
                "unique_users": 1,
            },
        },
        "twitter": {
            "tweets": [
                {
                    "id": "456",
                    "text": "Testing caffeine",
                    "author": "testuser",
                    "created_at": "2023-01-01T00:00:00Z",
                    "metrics": {"likes": 5},
                    "hashtags": ["nootropics"],
                }
            ],
            "stats": {
                "total_tweets": 1,
                "unique_users": 1,
                "total_likes": 5,
            },
        },
    }


@pytest.fixture
def invalid_swiss_data():
    """Invalid Swiss data fixture."""
    return {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 1.5,  # Invalid probability
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": {
            "mw": "not a number",  # Invalid numeric value
            "logp": -0.07,
            "hbd": 0,
            "hba": 6,
            "tpsa": 58.44,
        },
    }


@pytest.fixture
def invalid_community_data():
    """Invalid community data fixture."""
    return {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "not a url",  # Invalid URL format
            "effects": "not an array",  # Invalid array type
            "dosage": {"common": "100-200 mg"},
            "duration": {"onset": "15-45 minutes"},
        },
        "erowid": {
            "name": "Caffeine",
            # Missing required fields
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": {"dose_common": "100-200mg"},
            "interactions": "not an array",  # Invalid array type
        },
    }


@pytest.fixture
def invalid_social_data():
    """Invalid social data fixture."""
    return {
        "reddit": {
            "posts": [
                {
                    "id": "123",
                    "title": "Test Post",
                    # Missing required fields
                }
            ],
            "stats": {
                "total_posts": "not a number",  # Invalid numeric value
                "total_comments": 0,
                "unique_users": 1,
            },
        },
        "twitter": {
            "tweets": [
                {
                    "id": "456",
                    "text": "Testing caffeine",
                    "author": "testuser",
                    "created_at": "not a date",  # Invalid date format
                    "metrics": {"likes": "not a number"},  # Invalid numeric value
                    "hashtags": "not an array",  # Invalid array type
                }
            ],
            "stats": {
                "total_tweets": 1,
                "unique_users": 1,
                "total_likes": "not a number",  # Invalid numeric value
            },
        },
    }
