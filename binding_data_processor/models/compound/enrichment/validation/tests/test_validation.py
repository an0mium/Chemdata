"""Tests for web enrichment data validation."""

import pytest

from ..data import (
    validate_swiss_data,
    validate_community_data,
    validate_social_data,
    ValidationError,
)


def test_validate_swiss_data():
    """Test Swiss data validation."""
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
    validate_swiss_data(data)
    
    # Missing target predictions
    invalid_data = {
        "adme_properties": {
            "mw": 194.19,
        }
    }
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(invalid_data)
    assert "Missing target predictions" in str(exc_info.value)
    
    # Invalid probability value
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
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(invalid_data)
    assert "Invalid probability value" in str(exc_info.value)
    
    # Invalid numeric property
    invalid_data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 0.95,
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": {
            "mw": "not a number",
        },
    }
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(invalid_data)
    assert "Invalid numeric value" in str(exc_info.value)


def test_validate_community_data():
    """Test community data validation."""
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
    validate_community_data(data)
    
    # Missing required field in PsychonautWiki data
    invalid_data = {
        "psychonautwiki": {
            "name": "Caffeine",
        },
        "erowid": data["erowid"],
        "tripsit": data["tripsit"],
    }
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(invalid_data)
    assert "Missing required field" in str(exc_info.value)
    
    # Missing required field in Erowid data
    invalid_data = {
        "psychonautwiki": data["psychonautwiki"],
        "erowid": {
            "name": "Caffeine",
        },
        "tripsit": data["tripsit"],
    }
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(invalid_data)
    assert "Missing required field" in str(exc_info.value)
    
    # Missing required field in TripSit data
    invalid_data = {
        "psychonautwiki": data["psychonautwiki"],
        "erowid": data["erowid"],
        "tripsit": {
            "name": "Caffeine",
        },
    }
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(invalid_data)
    assert "Missing required field" in str(exc_info.value)


def test_validate_social_data():
    """Test social data validation."""
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
    validate_social_data(data)
    
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
    with pytest.raises(ValidationError) as exc_info:
        validate_social_data(invalid_data)
    assert "Missing required field" in str(exc_info.value)
    
    # Missing required field in Twitter tweet
    invalid_data = {
        "reddit": data["reddit"],
        "twitter": {
            "tweets": [
                {
                    "id": "456",
                    "text": "Testing caffeine",
                }
            ],
            "stats": data["twitter"]["stats"],
        },
    }
    with pytest.raises(ValidationError) as exc_info:
        validate_social_data(invalid_data)
    assert "Missing required field" in str(exc_info.value)
    
    # Invalid numeric value in Reddit stats
    invalid_data = {
        "reddit": {
            "posts": data["reddit"]["posts"],
            "stats": {
                "total_posts": "not a number",
                "total_comments": 0,
                "unique_users": 1,
            },
        },
        "twitter": data["twitter"],
    }
    with pytest.raises(ValidationError) as exc_info:
        validate_social_data(invalid_data)
    assert "Invalid numeric value" in str(exc_info.value)
    
    # Invalid numeric value in Twitter stats
    invalid_data = {
        "reddit": data["reddit"],
        "twitter": {
            "tweets": data["twitter"]["tweets"],
            "stats": {
                "total_tweets": "not a number",
                "unique_users": 1,
                "total_likes": 5,
            },
        },
    }
    with pytest.raises(ValidationError) as exc_info:
        validate_social_data(invalid_data)
    assert "Invalid numeric value" in str(exc_info.value)
