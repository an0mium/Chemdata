"""Tests for community data validation functions."""

import pytest

from ..data import validate_community_data, ValidationError


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


def test_validate_community_data_success(
    valid_psychonautwiki_data,
    valid_erowid_data,
    valid_tripsit_data,
):
    """Test successful validation of community data."""
    data = {
        "psychonautwiki": valid_psychonautwiki_data,
        "erowid": valid_erowid_data,
        "tripsit": valid_tripsit_data,
    }
    
    # Should not raise any errors
    validate_community_data(data)


def test_validate_community_data_missing_required_fields():
    """Test validation with missing required fields."""
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
    
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(data)
    assert "Missing required source 'psychonautwiki'" in str(exc_info.value)


def test_validate_psychonautwiki_data_invalid(valid_erowid_data, valid_tripsit_data):
    """Test validation with invalid PsychonautWiki data."""
    # Invalid effects type
    data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": "not a list",  # Should be list
            "dosage": {"common": "100-200 mg"},
            "duration": {"onset": "15-45 minutes"},
        },
        "erowid": valid_erowid_data,
        "tripsit": valid_tripsit_data,
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(data)
    assert "Invalid effects format in psychonautwiki data" in str(exc_info.value)


def test_validate_erowid_data_invalid(valid_psychonautwiki_data, valid_tripsit_data):
    """Test validation with invalid Erowid data."""
    # Invalid experience rating
    data = {
        "psychonautwiki": valid_psychonautwiki_data,
        "erowid": {
            "name": "Caffeine",
            "url": "https://erowid.org/chemicals/caffeine/",
            "experiences": [
                {
                    "title": "Test Experience",
                    "date": "2023-01-01",
                    "rating": "not a number",  # Should be number
                    "body": "Test report",
                }
            ],
            "health": [],
        },
        "tripsit": valid_tripsit_data,
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(data)
    assert "Invalid experience rating in erowid data" in str(exc_info.value)


def test_validate_tripsit_data_invalid(valid_psychonautwiki_data, valid_erowid_data):
    """Test validation with invalid TripSit data."""
    # Invalid properties format
    data = {
        "psychonautwiki": valid_psychonautwiki_data,
        "erowid": valid_erowid_data,
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": "not an object",  # Should be object
            "interactions": [],
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(data)
    assert "Invalid properties format in tripsit data" in str(exc_info.value)


def test_validate_community_data_invalid_urls():
    """Test validation with invalid URLs."""
    data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "not a url",  # Invalid URL
            "effects": [],
            "dosage": {},
            "duration": {},
        },
        "erowid": {
            "name": "Caffeine",
            "url": "not a url",  # Invalid URL
            "experiences": [],
            "health": [],
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "not a url",  # Invalid URL
            "properties": {},
            "interactions": [],
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(data)
    assert "Invalid URL format" in str(exc_info.value)


def test_validate_community_data_invalid_dates():
    """Test validation with invalid dates."""
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
                    "title": "Test Experience",
                    "date": "not a date",  # Invalid date
                    "rating": 4,
                    "body": "Test report",
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
    
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(data)
    assert "Invalid date format" in str(exc_info.value)


def test_validate_community_data_invalid_dosages():
    """Test validation with invalid dosage formats."""
    data = {
        "psychonautwiki": {
            "name": "Caffeine",
            "url": "https://psychonautwiki.org/wiki/Caffeine",
            "effects": [],
            "dosage": {
                "common": 100,  # Should be string
            },
            "duration": {},
        },
        "erowid": {
            "name": "Caffeine",
            "url": "https://erowid.org/chemicals/caffeine/",
            "experiences": [],
            "health": [],
            "dosage": {
                "ranges": {
                    "common": 100,  # Should be string
                },
            },
        },
        "tripsit": {
            "name": "Caffeine",
            "url": "https://drugs.tripsit.me/caffeine",
            "properties": {
                "dose_common": 100,  # Should be string
            },
            "interactions": [],
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_community_data(data)
    assert "Invalid dosage format" in str(exc_info.value)
