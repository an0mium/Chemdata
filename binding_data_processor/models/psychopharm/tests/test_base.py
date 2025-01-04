"""Tests for base psychopharm functionality."""

import pytest
from datetime import datetime

from ..base import (
    PsychoactiveClass,
    NootropicMechanism,
    RiskLevel,
    DoseRange,
    TimeRange,
    EffectScore,
    ReceptorBinding,
    BaseMixin,
)


def test_psychoactive_class_enum():
    """Test PsychoactiveClass enumeration."""
    assert PsychoactiveClass.UNKNOWN.value == "unknown"
    assert PsychoactiveClass.PSYCHEDELIC.value == "psychedelic"
    assert PsychoactiveClass.STIMULANT.value == "stimulant"
    assert PsychoactiveClass.DEPRESSANT.value == "depressant"


def test_nootropic_mechanism_enum():
    """Test NootropicMechanism enumeration."""
    assert NootropicMechanism.CHOLINERGIC.value == "cholinergic"
    assert NootropicMechanism.GLUTAMATERGIC.value == "glutamatergic"
    assert NootropicMechanism.DOPAMINERGIC.value == "dopaminergic"


def test_risk_level_enum():
    """Test RiskLevel enumeration."""
    assert RiskLevel.UNKNOWN.value == 0
    assert RiskLevel.LOW.value == 1
    assert RiskLevel.MODERATE.value == 2
    assert RiskLevel.HIGH.value == 3
    assert RiskLevel.SEVERE.value == 4


def test_dose_range():
    """Test DoseRange type."""
    dose = (10.0, 50.0, 25.0)  # min, max, recommended
    assert isinstance(dose, DoseRange)
    assert len(dose) == 3
    assert all(isinstance(d, float) for d in dose)


def test_time_range():
    """Test TimeRange type."""
    time = (30.0, 240.0)  # onset, duration
    assert isinstance(time, TimeRange)
    assert len(time) == 2
    assert all(isinstance(t, float) for t in time)


def test_effect_score():
    """Test EffectScore type."""
    score = (0.8, 0.9)  # magnitude, confidence
    assert isinstance(score, EffectScore)
    assert len(score) == 2
    assert all(isinstance(s, float) for s in score)
    assert 0.0 <= score[0] <= 1.0  # magnitude
    assert 0.0 <= score[1] <= 1.0  # confidence


def test_receptor_binding():
    """Test ReceptorBinding type."""
    binding = (1.0, 0.8, "agonist")  # affinity, confidence, activity
    assert isinstance(binding, ReceptorBinding)
    assert len(binding) == 3
    assert isinstance(binding[0], float)  # affinity
    assert isinstance(binding[1], float)  # confidence
    assert isinstance(binding[2], str)    # activity


class TestBaseMixin:
    """Tests for BaseMixin class."""

    def test_initialization(self):
        """Test initialization of BaseMixin."""
        base = BaseMixin()
        assert base.confidence_scores == {}
        assert base.last_updated is None

    def test_get_base_dict(self):
        """Test get_base_dict method."""
        base = BaseMixin()
        base.confidence_scores = {"test": 0.8}
        base.last_updated = datetime(2023, 1, 1)
        
        data = base.get_base_dict()
        assert isinstance(data, dict)
        assert "confidence_scores" in data
        assert "last_updated" in data
        assert data["confidence_scores"]["test"] == 0.8
        assert data["last_updated"] == "2023-01-01T00:00:00"

    def test_merge_base_data(self):
        """Test merge_base_data method."""
        base1 = BaseMixin()
        base2 = BaseMixin()
        
        # Set up test data
        base1.confidence_scores = {"score1": 0.7}
        base1.last_updated = datetime(2023, 1, 1)
        
        base2.confidence_scores = {"score1": 0.8, "score2": 0.9}
        base2.last_updated = datetime(2023, 1, 2)
        
        # Merge data
        base1.merge_base_data(base2)
        
        # Check results
        assert base1.confidence_scores["score1"] == 0.8  # Higher confidence wins
        assert base1.confidence_scores["score2"] == 0.9  # New score added
        assert base1.last_updated == datetime(2023, 1, 2)  # Later date wins

    def test_update_confidence(self):
        """Test update_confidence method."""
        base = BaseMixin()
        
        # Initial update
        base.update_confidence("test", 0.7)
        assert base.confidence_scores["test"] == 0.7
        
        # Update with higher confidence
        base.update_confidence("test", 0.8)
        assert base.confidence_scores["test"] == 0.8
        
        # Update with lower confidence (should not change)
        base.update_confidence("test", 0.6)
        assert base.confidence_scores["test"] == 0.8

    def test_validate_confidence_scores(self):
        """Test validate_confidence_scores method."""
        base = BaseMixin()
        
        # Valid scores
        base.confidence_scores = {
            "score1": 0.7,
            "score2": 0.9,
        }
        assert base.validate_confidence_scores()
        
        # Invalid score (> 1.0)
        base.confidence_scores["invalid"] = 1.1
        assert not base.validate_confidence_scores()
        
        # Invalid score (< 0.0)
        base.confidence_scores["invalid"] = -0.1
        assert not base.validate_confidence_scores()
        
        # Invalid score (wrong type)
        base.confidence_scores["invalid"] = "not a number"
        assert not base.validate_confidence_scores()


if __name__ == "__main__":
    pytest.main([__file__])
