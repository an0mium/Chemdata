"""Tests for safety assessment and risk analysis functionality."""

import pytest

from ..base import RiskLevel, PsychoactiveClass
from ..safety import SafetyProfileMixin


class TestSafetyProfileMixin:
    """Tests for SafetyProfileMixin class."""

    def setup_method(self):
        """Set up test instance."""
        self.safety = SafetyProfileMixin()
        
        # Set up test data
        self.safety.safety_alerts = {
            "cardiotoxicity": RiskLevel.HIGH,
            "psychosis": RiskLevel.MODERATE,
        }
        self.safety.contraindications = {
            "cardiovascular_disease",
            "psychotic_disorders",
        }
        self.safety.interaction_risks = {
            "SSRI": ("serotonin_syndrome", RiskLevel.SEVERE),
            "MAOI": ("hypertensive_crisis", RiskLevel.SEVERE),
        }

    def test_initialization(self):
        """Test initialization of SafetyProfileMixin."""
        safety = SafetyProfileMixin()
        assert safety.safety_alerts == {}
        assert safety.contraindications == set()
        assert safety.interaction_risks == {}

    def test_get_safety_dict(self):
        """Test get_safety_dict method."""
        data = self.safety.get_safety_dict()
        
        assert isinstance(data, dict)
        assert "alerts" in data
        assert "contraindications" in data
        assert "interactions" in data
        assert "risk_assessment" in data
        assert "recommendations" in data
        
        # Check alerts
        alerts = data["alerts"]
        assert alerts["cardiotoxicity"] == RiskLevel.HIGH.value
        assert alerts["psychosis"] == RiskLevel.MODERATE.value
        
        # Check contraindications
        contraindications = data["contraindications"]
        assert "cardiovascular_disease" in contraindications
        assert "psychotic_disorders" in contraindications
        
        # Check interactions
        interactions = data["interactions"]
        assert interactions["SSRI"]["mechanism"] == "serotonin_syndrome"
        assert interactions["SSRI"]["risk_level"] == RiskLevel.SEVERE.value

    def test_analyze_binding_risks(self):
        """Test _analyze_binding_risks method."""
        # Set up receptor data
        self.safety.receptor_profiles = {
            "5-HT2B": (1.2, 0.9, "agonist"),    # Cardiotoxicity risk
            "D2": (5.6, 0.8, "agonist"),        # Addiction risk
        }
        
        risks = self.safety._analyze_binding_risks()
        
        assert "cardiotoxicity" in risks
        cardiotoxicity = risks["cardiotoxicity"]
        assert cardiotoxicity["risk_level"] == RiskLevel.SEVERE.value
        assert cardiotoxicity["mechanism"] == "5-HT2B agonism"
        assert cardiotoxicity["affinity"] == 1.2
        assert cardiotoxicity["confidence"] == 0.9
        
        assert "addiction" in risks
        addiction = risks["addiction"]
        assert addiction["risk_level"] == RiskLevel.HIGH.value
        assert addiction["mechanism"] == "D2 agonism"

    def test_analyze_activity_risks(self):
        """Test _analyze_activity_risks method."""
        # Set up activity data
        self.safety.psychoactive_class = PsychoactiveClass.PSYCHEDELIC
        self.safety.nootropic_mechanisms = {
            "GLUTAMATERGIC",
        }
        
        risks = self.safety._analyze_activity_risks()
        
        assert "class_risks" in risks
        class_risks = risks["class_risks"]
        assert "psychosis" in class_risks
        assert class_risks["psychosis"]["risk_level"] == RiskLevel.HIGH.value
        assert "hppd" in class_risks
        
        assert "excitotoxicity" in risks
        assert risks["excitotoxicity"]["risk_level"] == RiskLevel.MODERATE.value

    def test_analyze_interaction_risks(self):
        """Test _analyze_interaction_risks method."""
        # Set up receptor data for drug class detection
        self.safety.receptor_profiles = {
            "SERT": (1.2, 0.9, "inhibitor"),    # SSRI
            "MAO": (5.6, 0.8, "inhibitor"),     # MAOI
        }
        
        risks = self.safety._analyze_interaction_risks()
        
        assert "SSRI" in risks
        ssri_risks = risks["SSRI"]["interactions"]
        assert "MAOI" in ssri_risks
        assert ssri_risks["MAOI"]["risk_level"] == RiskLevel.SEVERE.value
        
        assert "MAOI" in risks
        maoi_risks = risks["MAOI"]["interactions"]
        assert "SSRI" in maoi_risks
        assert maoi_risks["SSRI"]["risk_level"] == RiskLevel.SEVERE.value

    def test_generate_safety_recommendations(self):
        """Test _generate_safety_recommendations method."""
        # Set up risk data
        self.safety.receptor_profiles = {
            "5-HT2B": (1.2, 0.9, "agonist"),    # Cardiotoxicity
            "NMDA": (3.4, 0.8, "antagonist"),   # Neurotoxicity
        }
        self.safety.psychoactive_class = PsychoactiveClass.PSYCHEDELIC
        
        recommendations = self.safety._generate_safety_recommendations()
        
        # Check contraindications
        assert (
            "Contraindicated in patients with cardiovascular conditions"
            in recommendations["contraindications"]
        )
        
        # Check monitoring recommendations
        assert (
            "Regular cardiovascular monitoring required"
            in recommendations["monitoring"]
        )
        assert (
            "Monitor for cognitive and neurological effects"
            in recommendations["monitoring"]
        )
        
        # Check precautions
        assert (
            "Use with caution due to potential neurotoxicity"
            in recommendations["precautions"]
        )
        
        # Check risk patterns
        assert "risk_patterns" in recommendations
        patterns = recommendations["risk_patterns"]
        assert any("severe risks" in pattern.lower() for pattern in patterns)

    def test_merge_safety_data(self):
        """Test merge_safety_data method."""
        other = SafetyProfileMixin()
        
        # Set up other instance data
        other.safety_alerts = {
            "cardiotoxicity": RiskLevel.SEVERE,    # Higher risk
            "neurotoxicity": RiskLevel.HIGH,       # New alert
        }
        other.contraindications = {
            "cardiovascular_disease",              # Duplicate
            "seizure_disorders",                   # New
        }
        other.interaction_risks = {
            "SSRI": ("serotonin_toxicity", RiskLevel.HIGH),     # Lower risk (shouldn't update)
            "TCA": ("serotonin_syndrome", RiskLevel.HIGH),      # New interaction
        }
        
        # Merge data
        self.safety.merge_safety_data(other)
        
        # Check alerts merging (highest risk level wins)
        assert self.safety.safety_alerts["cardiotoxicity"] == RiskLevel.SEVERE
        assert self.safety.safety_alerts["neurotoxicity"] == RiskLevel.HIGH
        assert self.safety.safety_alerts["psychosis"] == RiskLevel.MODERATE
        
        # Check contraindications merging (union)
        assert "cardiovascular_disease" in self.safety.contraindications
        assert "seizure_disorders" in self.safety.contraindications
        assert "psychotic_disorders" in self.safety.contraindications
        
        # Check interaction risks merging (highest risk level wins)
        assert self.safety.interaction_risks["SSRI"] == ("serotonin_syndrome", RiskLevel.SEVERE)
        assert self.safety.interaction_risks["TCA"] == ("serotonin_syndrome", RiskLevel.HIGH)
        assert self.safety.interaction_risks["MAOI"] == ("hypertensive_crisis", RiskLevel.SEVERE)


if __name__ == "__main__":
    pytest.main([__file__])
