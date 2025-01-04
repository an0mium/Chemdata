"""Tests for activity analysis functionality."""

import pytest

from ..activity_analysis import ActivityAnalysisMixin


class TestCompound(ActivityAnalysisMixin):
    """Test class implementing ActivityAnalysisMixin."""

    def __init__(self, primary_activity=None, mechanism=None, effect_profile=None):
        self.primary_activity = primary_activity
        self.mechanism_of_action = mechanism
        self.effect_profile = effect_profile or {}


@pytest.fixture
def test_activity_data():
    """Create test activity data."""
    return {
        "primary_activity": 1.2,  # nM
        "mechanism_of_action": "serotonin 5-HT2A receptor agonist",
        "effect_profile": {
            "psychedelic": (0.9, 0.95),
            "mood": (0.8, 0.9),
            "cognition": (0.6, 0.7),
            "perception": (0.85, 0.9),
            "anxiety": (-0.3, 0.6),  # Negative effect with low confidence
        },
    }


def test_primary_activity(test_activity_data):
    """Test analyzing primary activity."""
    compound = TestCompound(
        primary_activity=test_activity_data["primary_activity"],
        mechanism=test_activity_data["mechanism_of_action"]
    )
    result = compound._analyze_primary_activity()
    
    assert result["value"] == 1.2
    assert result["unit"] == "nM"
    assert result["mechanism"] == "serotonin 5-HT2A receptor agonist"
    assert result["potency"] == "high"  # nM range is high potency


def test_effect_profile(test_activity_data):
    """Test analyzing effect profile."""
    compound = TestCompound(effect_profile=test_activity_data["effect_profile"])
    profile = compound._analyze_effect_profile()
    
    # Check positive effects
    assert "psychedelic" in profile["positive_effects"]
    assert profile["positive_effects"]["psychedelic"]["score"] >= 0.9
    assert profile["positive_effects"]["psychedelic"]["confidence"] >= 0.95

    # Check negative effects
    assert "anxiety" in profile["negative_effects"]
    assert profile["negative_effects"]["anxiety"]["score"] <= -0.3
    assert profile["negative_effects"]["anxiety"]["confidence"] <= 0.6

    # Check overall profile
    assert profile["primary_domain"] == "psychedelic"
    assert profile["confidence_range"] == (0.6, 0.95)


def test_activity_patterns(test_activity_data):
    """Test analyzing activity patterns."""
    compound = TestCompound(
        primary_activity=test_activity_data["primary_activity"],
        mechanism=test_activity_data["mechanism_of_action"],
        effect_profile=test_activity_data["effect_profile"]
    )
    patterns = compound._analyze_activity_patterns()
    
    # Should identify psychedelic pattern
    psychedelic_pattern = next(
        p for p in patterns if p["domain"] == "psychedelic"
    )
    assert psychedelic_pattern["effects"] == ["psychedelic", "perception"]
    assert psychedelic_pattern["confidence"] >= 0.9

    # Should identify mood pattern
    mood_pattern = next(
        p for p in patterns if p["domain"] == "mood"
    )
    assert "mood" in mood_pattern["effects"]
    assert mood_pattern["confidence"] >= 0.8


def test_activity_confidence(test_activity_data):
    """Test analyzing activity confidence."""
    compound = TestCompound(
        primary_activity=test_activity_data["primary_activity"],
        mechanism=test_activity_data["mechanism_of_action"],
        effect_profile=test_activity_data["effect_profile"]
    )
    confidence = compound._analyze_activity_confidence()
    
    assert 0.8 <= confidence["overall"] <= 0.95
    assert confidence["mechanism"] > 0
    assert len(confidence["by_effect"]) == len(test_activity_data["effect_profile"])
    assert confidence["by_effect"]["psychedelic"] >= 0.95


def test_no_activity_data():
    """Test analysis with no activity data."""
    compound = TestCompound()
    
    # Test primary activity analysis
    primary = compound._analyze_primary_activity()
    assert primary["value"] == 0
    assert primary["mechanism"] == "unknown"
    assert primary["potency"] == "unknown"

    # Test effect profile analysis
    profile = compound._analyze_effect_profile()
    assert not profile["positive_effects"]
    assert not profile["negative_effects"]
    assert profile["primary_domain"] == "unknown"
    assert profile["confidence_range"] == (0, 0)

    # Test activity patterns analysis
    patterns = compound._analyze_activity_patterns()
    assert not patterns

    # Test confidence analysis
    confidence = compound._analyze_activity_confidence()
    assert confidence["overall"] == 0
    assert confidence["mechanism"] == 0
    assert not confidence["by_effect"]
