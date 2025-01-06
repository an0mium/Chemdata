"""Tests for web enrichment data schemas."""

import pytest
from jsonschema import validate, ValidationError as JsonSchemaError

from ..schema import (
    SWISS_DATA_SCHEMA,
    COMMUNITY_DATA_SCHEMA,
    SOCIAL_DATA_SCHEMA,
)


def test_swiss_data_schema():
    """Test Swiss data schema validation."""
    # Valid data
    data = {
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
    validate(instance=data, schema=SWISS_DATA_SCHEMA)
    
    # Missing required field
    invalid_data = {
        "adme_properties": {
            "mw": 194.19,
        }
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=SWISS_DATA_SCHEMA)
    
    # Invalid probability type
    invalid_data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": "0.95",  # Should be number
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": {
            "mw": 194.19,
        },
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=SWISS_DATA_SCHEMA)
    
    # Invalid probability range
    invalid_data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 1.5,  # Should be between 0 and 1
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": {
            "mw": 194.19,
        },
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=SWISS_DATA_SCHEMA)


def test_community_data_schema():
    """Test community data schema validation."""
    # Valid data
    data = {
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
    validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)
    
    # Missing required field in PsychonautWiki data
    invalid_data = {
        "psychonautwiki": {
            "name": "Caffeine",
        },
        "erowid": data["erowid"],
        "tripsit": data["tripsit"],
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=COMMUNITY_DATA_SCHEMA)
    
    # Invalid URL format
    invalid_data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "not a url",  # Should be URL format
            "effects": ["Stimulation"],
            "dosage": {"common": "100-200 mg"},
            "duration": {"onset": "15-45 minutes"},
        },
        "erowid": data["erowid"],
        "tripsit": data["tripsit"],
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=COMMUNITY_DATA_SCHEMA)
    
    # Invalid array type
    invalid_data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": "Stimulation",  # Should be array
            "dosage": {"common": "100-200 mg"},
            "duration": {"onset": "15-45 minutes"},
        },
        "erowid": data["erowid"],
        "tripsit": data["tripsit"],
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=COMMUNITY_DATA_SCHEMA)


def test_social_data_schema():
    """Test social data schema validation."""
    # Valid data
    data = {
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
    validate(instance=data, schema=SOCIAL_DATA_SCHEMA)
    
    # Missing required field in Reddit post
    invalid_data = {
        "reddit": {
            "posts": [
                {
                    "id": "123",
                    "title": "Test Post",
                }
            ],
            "stats": data["reddit"]["stats"],
        },
        "twitter": data["twitter"],
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=SOCIAL_DATA_SCHEMA)
    
    # Invalid numeric type
    invalid_data = {
        "reddit": {
            "posts": data["reddit"]["posts"],
            "stats": {
                "total_posts": "1",  # Should be number
                "total_comments": 0,
                "unique_users": 1,
            },
        },
        "twitter": data["twitter"],
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=SOCIAL_DATA_SCHEMA)
    
    # Invalid date format
    invalid_data = {
        "reddit": data["reddit"],
        "twitter": {
            "tweets": [
                {
                    "id": "456",
                    "text": "Testing caffeine",
                    "author": "testuser",
                    "created_at": "not a date",  # Should be ISO date
                    "metrics": {"likes": 5},
                    "hashtags": ["nootropics"],
                }
            ],
            "stats": data["twitter"]["stats"],
        },
    }
    with pytest.raises(JsonSchemaError):
        validate(instance=invalid_data, schema=SOCIAL_DATA_SCHEMA)
