"""Tests for binding analysis functionality."""

import pytest
from rdkit import Chem

from ..binding_analysis import BindingAnalysisMixin


class TestCompound(BindingAnalysisMixin):
    """Test class implementing BindingAnalysisMixin."""

    def __init__(self, targets=None):
        self.targets = targets or []


@pytest.fixture
def test_targets():
    """Create test binding data."""
    return [
        {
            "common_name": "5-HT2A",
            "affinity_value": 1.2,
            "affinity_type": "Ki",
            "affinity_unit": "nM",
            "activity_type": "agonist",
            "confidence": 0.9,
            "assay_details": {
                "family": "serotonin",
                "effects": ["psychedelic", "mood"],
                "mechanisms": ["receptor activation"],
            },
            "reference_dois": ["10.1234/example.1"],
            "experimental_conditions": {"pH": 7.4, "temp": 37},
        },
        {
            "common_name": "5-HT2B",
            "affinity_value": 15.0,
            "affinity_type": "Ki",
            "affinity_unit": "nM",
            "activity_type": "agonist",
            "confidence": 0.85,
            "assay_details": {
                "family": "serotonin",
                "effects": ["cardiovascular"],
                "mechanisms": ["receptor activation"],
            },
            "reference_dois": ["10.1234/example.2"],
            "experimental_conditions": {"pH": 7.4, "temp": 37},
        },
        {
            "common_name": "D2",
            "affinity_value": 120.0,
            "affinity_type": "Ki",
            "affinity_unit": "nM",
            "activity_type": "antagonist",
            "confidence": 0.8,
            "assay_details": {
                "family": "dopamine",
                "effects": ["antipsychotic"],
                "mechanisms": ["receptor blockade"],
            },
            "reference_dois": ["10.1234/example.3"],
            "experimental_conditions": {"pH": 7.4, "temp": 37},
        },
    ]


def test_strongest_binding(test_targets):
    """Test finding strongest binding interaction."""
    compound = TestCompound(test_targets)
    result = compound._find_strongest_binding()
    
    assert result is not None
    assert result["target"] == "5-HT2A"
    assert result["affinity"] == 1.2
    assert result["type"] == "Ki"
    assert result["unit"] == "nM"
    assert result["activity"] == "agonist"
    assert result["confidence"] == 0.9


def test_binding_patterns(test_targets):
    """Test analyzing binding patterns."""
    compound = TestCompound(test_targets)
    patterns = compound._analyze_binding_patterns()
    
    # Should find serotonin receptor pattern
    serotonin_pattern = next(
        p for p in patterns if p["family"] == "serotonin"
    )
    assert serotonin_pattern["target_count"] == 2
    assert serotonin_pattern["common_activity"] == "agonist"
    assert 0.85 <= serotonin_pattern["confidence"] <= 0.9

    # Should find dopamine receptor pattern
    dopamine_pattern = next(
        p for p in patterns if p["family"] == "dopamine"
    )
    assert dopamine_pattern["target_count"] == 1
    assert dopamine_pattern["common_activity"] == "antagonist"
    assert dopamine_pattern["confidence"] == 0.8


def test_target_selectivity(test_targets):
    """Test analyzing target selectivity."""
    compound = TestCompound(test_targets)
    selectivity = compound._analyze_target_selectivity()
    
    assert selectivity["selectivity"] in ["highly selective", "moderately selective", "non-selective"]
    assert isinstance(selectivity["score"], float)
    assert len(selectivity["ratios"]) > 0

    # Should be selective between serotonin and dopamine families
    assert selectivity["score"] > 10  # At least moderately selective


def test_receptor_families(test_targets):
    """Test analyzing receptor families."""
    compound = TestCompound(test_targets)
    families = compound._analyze_receptor_families()
    
    assert "serotonin" in families
    assert "dopamine" in families

    # Check serotonin family
    serotonin = families["serotonin"]
    assert serotonin["count"] == 2
    assert serotonin["activities"]["agonist"] == 2
    assert 0.85 <= serotonin["confidence"] <= 0.9

    # Check dopamine family
    dopamine = families["dopamine"]
    assert dopamine["count"] == 1
    assert dopamine["activities"]["antagonist"] == 1
    assert dopamine["confidence"] == 0.8


def test_binding_confidence(test_targets):
    """Test analyzing binding confidence."""
    compound = TestCompound(test_targets)
    confidence = compound._analyze_binding_confidence()
    
    assert 0.8 <= confidence["overall"] <= 0.9
    assert len(confidence["by_target"]) == 3

    # Check individual target confidence
    assert confidence["by_target"]["5-HT2A"]["confidence"] == 0.9
    assert confidence["by_target"]["5-HT2B"]["confidence"] == 0.85
    assert confidence["by_target"]["D2"]["confidence"] == 0.8

    # Check experimental support
    assert confidence["by_target"]["5-HT2A"]["experimental"]
    assert confidence["by_target"]["5-HT2B"]["experimental"]
    assert confidence["by_target"]["D2"]["experimental"]


def test_no_targets():
    """Test analysis with no targets."""
    compound = TestCompound([])
    
    assert compound._find_strongest_binding() is None
    assert compound._analyze_binding_patterns() == []
    assert compound._analyze_target_selectivity() == {"selectivity": "unknown", "score": 0.0}
    assert compound._analyze_receptor_families() == {}
    assert compound._analyze_binding_confidence() == {"overall": 0.0, "by_target": {}}


def test_single_target(test_targets):
    """Test analysis with single target."""
    compound = TestCompound([test_targets[0]])  # Just 5-HT2A
    
    # Should still find strongest binding
    strongest = compound._find_strongest_binding()
    assert strongest is not None
    assert strongest["target"] == "5-HT2A"

    # Should have limited patterns
    patterns = compound._analyze_binding_patterns()
    assert len(patterns) == 1
    assert patterns[0]["family"] == "serotonin"

    # Can't calculate selectivity with single target
    selectivity = compound._analyze_target_selectivity()
    assert selectivity["selectivity"] == "unknown"

    # Should still analyze family
    families = compound._analyze_receptor_families()
    assert len(families) == 1
    assert "serotonin" in families


def test_invalid_targets():
    """Test analysis with invalid target data."""
    # Missing required fields
    invalid_targets = [
        {
            "common_name": "TEST",
            # Missing affinity_value
            "affinity_type": "Ki",
        }
    ]
    compound = TestCompound(invalid_targets)
    
    # Should handle missing data gracefully
    assert compound._find_strongest_binding() is None
    assert compound._analyze_binding_patterns() == []
    assert compound._analyze_target_selectivity() == {"selectivity": "unknown", "score": 0.0}
    assert compound._analyze_receptor_families() == {}
    assert compound._analyze_binding_confidence() == {"overall": 0.0, "by_target": {}}
