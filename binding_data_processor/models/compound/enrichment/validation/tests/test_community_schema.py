"""Tests for community data JSON schema."""

import pytest
from jsonschema import validate, ValidationError as JsonSchemaError

from ..schema import COMMUNITY_DATA_SCHEMA


@pytest.fixture
def valid_psychonautwiki_data():
    """Valid PsychonautWiki data fixture."""
    return {
        "name": "Caffeine",
        "url": "https://psychonautwiki.org/wiki/Caffeine",
        "effects": [
            "Stimulation",
            "Anxiety",
            "Increased heart rate",
        ],
        "dosage": {
            "threshold": "20 mg",
            "light": "50-100 mg",
            "common": "100-200 mg",
            "strong": "200-500 mg",
            "heavy": "500+ mg",
        },
        "duration": {
            "onset": "15-45 minutes",
            "duration": "2-5 hours",
            "after_effects": "2-4 hours",
        },
        "risks": [
            "Anxiety",
            "Insomnia",
            "Dependence",
        ],
        "interactions": [
            {
                "substance": "Alcohol",
                "note": "May reduce perception of alcohol intoxication",
            },
        ],
    }


@pytest.fixture
def valid_erowid_data():
    """Valid Erowid data fixture."""
    return {
        "name": "Caffeine",
        "url": "https://erowid.org/chemicals/caffeine/",
        "experiences": [
            {
                "title": "First Time Experience",
                "date": "2023-01-01",
                "rating": 4,
                "body": "Sample experience report",
                "author": "Anonymous",
                "tags": ["stimulant", "positive"],
            }
        ],
        "health": [
            "May cause anxiety in sensitive individuals",
            "Can lead to physical dependence",
        ],
        "effects": {
            "positive": [
                "Increased alertness",
                "Enhanced focus",
            ],
            "neutral": [
                "Mild stimulation",
            ],
            "negative": [
                "Anxiety",
                "Insomnia",
            ],
        },
        "dosage": {
            "notes": "Dosage information from user reports",
            "ranges": {
                "light": "50-100mg",
                "common": "100-200mg",
                "strong": "200-500mg",
            },
        },
    }


@pytest.fixture
def valid_tripsit_data():
    """Valid TripSit data fixture."""
    return {
        "name": "Caffeine",
        "url": "https://drugs.tripsit.me/caffeine",
        "properties": {
            "dose_light": "50-100mg",
            "dose_common": "100-200mg",
            "dose_strong": "200-500mg",
            "duration_total": "2-5 hours",
            "roas": ["oral", "sublingual"],
        },
        "interactions": [
            {
                "substance": "Alcohol",
                "status": "Caution",
                "note": "May reduce perception of alcohol intoxication",
            },
        ],
        "categories": [
            "stimulant",
            "common",
            "legal",
        ],
        "summary": "Common stimulant found in coffee and tea",
        "effects": {
            "physical": [
                "Increased heart rate",
                "Increased blood pressure",
            ],
            "cognitive": [
                "Alertness",
                "Focus",
            ],
        },
    }


def test_community_data_schema_valid(
    valid_psychonautwiki_data,
    valid_erowid_data,
    valid_tripsit_data,
):
    """Test schema validation with valid data."""
    data = {
        "psychonautwiki": valid_psychonautwiki_data,
        "erowid": valid_erowid_data,
        "tripsit": valid_tripsit_data,
    }
    
    # Should not raise any errors
    validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)


def test_community_data_schema_missing_required():
    """Test schema validation with missing required fields."""
    # Missing psychonautwiki data
    data = {
        "erowid": {
            "name": "Caffeine",
            "url": "https://erowid.org/chemicals/caffeine/",
            "experiences": [],
            "health": [],
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": {},
            "interactions": [],
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)


def test_community_data_schema_invalid_types():
    """Test schema validation with invalid types."""
    data = {
        "psychonautwiki": {
            "name": 123,  # Should be string
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": ["Stimulation"],
            "dosage": {"common": 100},  # Should be string
            "duration": {"onset": 15},  # Should be string
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
            "properties": {},
            "interactions": [],
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)


def test_community_data_schema_invalid_arrays():
    """Test schema validation with invalid arrays."""
    data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": "not an array",  # Should be array
            "dosage": {},
            "duration": {},
        },
        "erowid": {
            "name": "Caffeine",
            "url": "https://erowid.org/chemicals/caffeine/",
            "experiences": "not an array",  # Should be array
            "health": "not an array",  # Should be array
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": {},
            "interactions": "not an array",  # Should be array
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)


def test_community_data_schema_invalid_objects():
    """Test schema validation with invalid objects."""
    data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": [],
            "dosage": "not an object",  # Should be object
            "duration": "not an object",  # Should be object
        },
        "erowid": {
            "name": "Caffeine",
            "url": "https://erowid.org/chemicals/caffeine/",
            "experiences": [],
            "health": [],
            "effects": "not an object",  # Should be object
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": "not an object",  # Should be object
            "interactions": [],
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)


def test_community_data_schema_invalid_experience():
    """Test schema validation with invalid experience data."""
    data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": [],
            "dosage": {},
            "duration": {},
        },
        "erowid": {
            "name": "Caffeine",
            "url": "https://erowid.org/chemicals/caffeine/",
            "experiences": [
                {
                    "title": 123,  # Should be string
                    "date": "not a date",  # Should be date format
                    "rating": "not a number",  # Should be number
                    "body": 123,  # Should be string
                }
            ],
            "health": [],
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": {},
            "interactions": [],
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)


def test_community_data_schema_invalid_interactions():
    """Test schema validation with invalid interaction data."""
    data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": [],
            "dosage": {},
            "duration": {},
            "interactions": [
                {
                    "substance": 123,  # Should be string
                    "note": 123,  # Should be string
                }
            ],
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
            "properties": {},
            "interactions": [
                {
                    "substance": 123,  # Should be string
                    "status": 123,  # Should be string
                    "note": 123,  # Should be string
                }
            ],
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)


def test_community_data_schema_invalid_effects():
    """Test schema validation with invalid effects data."""
    data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": [123],  # Should be strings
            "dosage": {},
            "duration": {},
        },
        "erowid": {
            "name": "Caffeine",
            "url": "https://erowid.org/chemicals/caffeine/",
            "experiences": [],
            "health": [],
            "effects": {
                "positive": [123],  # Should be strings
                "neutral": [123],  # Should be strings
                "negative": [123],  # Should be strings
            },
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": {},
            "interactions": [],
            "effects": {
                "physical": [123],  # Should be strings
                "cognitive": [123],  # Should be strings
            },
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=COMMUNITY_DATA_SCHEMA)
