"""Tests for safety analysis functionality."""

import pytest

from ..safety_analysis import SafetyAnalysisMixin


class TestCompound(SafetyAnalysisMixin):
    """Test class implementing SafetyAnalysisMixin."""

    def __init__(self, safety_profile=None, property_data=None):
        self.safety_profile = safety_profile or {}
        self.molecular_weight = property_data.get("molecular_weight") if property_data else None
        self.logp = property_data.get("logp") if property_data else None
        self.hbd = property_data.get("hbd") if property_data else None
        self.hba = property_data.get("hba") if property_data else None
        self.tpsa = property_data.get("tpsa") if property_data else None
        self.rotatable_bonds = property_data.get("rotatable_bonds") if property_data else None


@pytest.fixture
def test_safety_data():
    """Create test safety data."""
    return {
        "safety_profile": {
            "known_risks": [
                {
                    "description": "Cardiovascular stimulation",
                    "severity": "moderate",
                    "confidence": 0.9,
                    "mechanism": "beta-adrenergic activation",
                    "conditions": ["high dose", "pre-existing conditions"],
                },
                {
                    "description": "Anxiety",
                    "severity": "mild",
                    "confidence": 0.85,
                    "mechanism": "CNS stimulation",
                    "conditions": ["sensitive individuals"],
                },
            ],
            "known_interactions": [
                {
                    "description": "MAO inhibitors",
                    "severity": "severe",
                    "confidence": 0.95,
                    "mechanism": "synergistic effects",
                    "recommendations": ["contraindicated"],
                },
                {
                    "description": "Other stimulants",
                    "severity": "moderate",
                    "confidence": 0.9,
                    "mechanism": "additive effects",
                    "recommendations": ["dose reduction"],
                },
            ],
        },
        "property_data": {
            "molecular_weight": 450.0,
            "logp": 4.2,
            "hbd": 3,
            "hba": 6,
            "tpsa": 90.0,
            "rotatable_bonds": 8,
        },
    }


def test_risk_assessment(test_safety_data):
    """Test risk assessment analysis."""
    compound = TestCompound(
        safety_profile=test_safety_data["safety_profile"],
        property_data=test_safety_data["property_data"]
    )
    assessment = compound._analyze_risk_assessment()
    
    # Check overall risk
    assert assessment["overall_risk"]["level"] in ["low", "moderate", "high", "severe"]
    assert 0 <= assessment["overall_risk"]["score"] <= 1.0
    assert 0.8 <= assessment["overall_risk"]["confidence"] <= 1.0

    # Check risk categories
    assert "cardiovascular" in assessment["risk_categories"]
    cv_risk = assessment["risk_categories"]["cardiovascular"]
    assert cv_risk["severity"] == "moderate"
    assert cv_risk["confidence"] == 0.9

    # Check risk factors
    assert len(assessment["risk_factors"]) >= 2
    assert any(f["type"] == "pharmacological" for f in assessment["risk_factors"])
    assert any(f["type"] == "chemical" for f in assessment["risk_factors"])


def test_safety_alerts(test_safety_data):
    """Test safety alerts analysis."""
    compound = TestCompound(
        safety_profile=test_safety_data["safety_profile"],
        property_data=test_safety_data["property_data"]
    )
    alerts = compound._analyze_safety_alerts()
    
    # Check property-based alerts
    assert any(a["type"] == "molecular_weight" for a in alerts)
    assert any(a["type"] == "logp" for a in alerts)

    # Check risk-based alerts
    assert any(
        a["description"] == "Cardiovascular stimulation"
        for a in alerts if a["type"] == "known_risk"
    )

    # Check interaction alerts
    assert any(
        a["description"] == "MAO inhibitors"
        for a in alerts if a["type"] == "interaction"
    )


def test_interaction_analysis(test_safety_data):
    """Test interaction analysis."""
    compound = TestCompound(safety_profile=test_safety_data["safety_profile"])
    interactions = compound._analyze_interactions()
    
    # Check severe interactions
    severe = interactions["severe"]
    assert len(severe) == 1
    assert severe[0]["description"] == "MAO inhibitors"
    assert severe[0]["confidence"] == 0.95

    # Check moderate interactions
    moderate = interactions["moderate"]
    assert len(moderate) == 1
    assert moderate[0]["description"] == "Other stimulants"
    assert moderate[0]["confidence"] == 0.9

    # Check overall assessment
    assert interactions["overall_risk"] == "high"
    assert interactions["recommendations"]


def test_safety_confidence(test_safety_data):
    """Test safety confidence analysis."""
    compound = TestCompound(
        safety_profile=test_safety_data["safety_profile"],
        property_data=test_safety_data["property_data"]
    )
    confidence = compound._analyze_safety_confidence()
    
    # Check overall confidence
    assert 0.8 <= confidence["overall"] <= 1.0

    # Check risk confidence
    assert len(confidence["by_risk"]) == 2
    assert confidence["by_risk"]["Cardiovascular stimulation"] == 0.9
    assert confidence["by_risk"]["Anxiety"] == 0.85

    # Check interaction confidence
    assert len(confidence["by_interaction"]) == 2
    assert confidence["by_interaction"]["MAO inhibitors"] == 0.95
    assert confidence["by_interaction"]["Other stimulants"] == 0.9


def test_no_safety_data():
    """Test analysis with no safety data."""
    compound = TestCompound()
    
    # Test risk assessment
    assessment = compound._analyze_risk_assessment()
    assert assessment["overall_risk"]["level"] == "unknown"
    assert assessment["overall_risk"]["score"] == 0
    assert not assessment["risk_categories"]
    assert not assessment["risk_factors"]

    # Test safety alerts
    alerts = compound._analyze_safety_alerts()
    assert not alerts

    # Test interactions
    interactions = compound._analyze_interactions()
    assert not interactions["severe"]
    assert not interactions["moderate"]
    assert interactions["overall_risk"] == "unknown"
    assert not interactions["recommendations"]

    # Test confidence
    confidence = compound._analyze_safety_confidence()
    assert confidence["overall"] == 0
    assert not confidence["by_risk"]
    assert not confidence["by_interaction"]


def test_partial_safety_data(test_safety_data):
    """Test analysis with partial safety data."""
    # Only known risks, no interactions
    partial_profile = {
        "known_risks": test_safety_data["safety_profile"]["known_risks"]
    }
    compound = TestCompound(safety_profile=partial_profile)
    
    # Should still analyze risks
    assessment = compound._analyze_risk_assessment()
    assert assessment["overall_risk"]["level"] != "unknown"
    assert assessment["risk_categories"]

    # Should have risk alerts but no interaction alerts
    alerts = compound._analyze_safety_alerts()
    assert any(a["type"] == "known_risk" for a in alerts)
    assert not any(a["type"] == "interaction" for a in alerts)

    # Should handle missing interactions gracefully
    interactions = compound._analyze_interactions()
    assert not interactions["severe"]
    assert not interactions["moderate"]
    assert interactions["overall_risk"] == "unknown"

    # Should have risk confidence but no interaction confidence
    confidence = compound._analyze_safety_confidence()
    assert confidence["by_risk"]
    assert not confidence["by_interaction"]


def test_property_based_alerts(test_safety_data):
    """Test property-based safety alerts."""
    compound = TestCompound(property_data=test_safety_data["property_data"])
    alerts = compound._analyze_safety_alerts()
    
    # Check molecular weight alert
    mw_alert = next(a for a in alerts if a["type"] == "molecular_weight")
    assert mw_alert["value"] == 450.0
    assert mw_alert["threshold"] == 500.0
    assert mw_alert["severity"] == "info"

    # Check LogP alert
    logp_alert = next(a for a in alerts if a["type"] == "logp")
    assert logp_alert["value"] == 4.2
    assert logp_alert["threshold"] == 5.0
    assert logp_alert["severity"] == "warning"

    # Check TPSA alert
    tpsa_alert = next(a for a in alerts if a["type"] == "tpsa")
    assert tpsa_alert["value"] == 90.0
    assert tpsa_alert["threshold"] == 140.0
    assert tpsa_alert["severity"] == "info"
