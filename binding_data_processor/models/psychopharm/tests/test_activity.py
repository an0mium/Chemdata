"""Tests for psychoactive and nootropic activity functionality."""

import pytest

from ..base import PsychoactiveClass, NootropicMechanism
from ..activity import ActivityProfileMixin


class TestActivityProfileMixin:
    """Tests for ActivityProfileMixin class."""

    def setup_method(self):
        """Set up test instance."""
        self.activity = ActivityProfileMixin()
        
        # Set up test data
        self.activity.psychoactive_class = PsychoactiveClass.PSYCHEDELIC
        self.activity.secondary_classes = {PsychoactiveClass.STIMULANT}
        self.activity.effect_profile = {
            "visual": (0.8, 0.9),  # magnitude, confidence
            "euphoria": (0.6, 0.8),
        }
        self.activity.nootropic_mechanisms = {
            NootropicMechanism.GLUTAMATERGIC,
            NootropicMechanism.DOPAMINERGIC,
        }
        self.activity.cognitive_effects = {
            "focus": (0.7, 0.8),
            "memory": (0.5, 0.7),
        }
        self.activity.side_effects = {
            "anxiety": (0.3, 0.6),
            "insomnia": (0.4, 0.7),
        }

    def test_initialization(self):
        """Test initialization of ActivityProfileMixin."""
        activity = ActivityProfileMixin()
        assert activity.psychoactive_class == PsychoactiveClass.UNKNOWN
        assert activity.secondary_classes == set()
        assert activity.effect_profile == {}
        assert activity.nootropic_mechanisms == set()
        assert activity.cognitive_effects == {}
        assert activity.side_effects == {}

    def test_get_activity_dict(self):
        """Test get_activity_dict method."""
        data = self.activity.get_activity_dict()
        
        assert isinstance(data, dict)
        assert "psychoactive" in data
        assert "nootropic" in data
        assert "effects" in data
        assert "mechanisms" in data
        assert "predictions" in data
        
        # Check psychoactive data
        psychoactive = data["psychoactive"]
        assert psychoactive["primary_class"] == PsychoactiveClass.PSYCHEDELIC.value
        assert PsychoactiveClass.STIMULANT.value in psychoactive["secondary_classes"]
        
        # Check nootropic data
        nootropic = data["nootropic"]
        assert NootropicMechanism.GLUTAMATERGIC.value in nootropic["mechanisms"]
        assert "focus" in nootropic["cognitive_effects"]
        assert "anxiety" in nootropic["side_effects"]

    def test_analyze_effects(self):
        """Test _analyze_effects method."""
        # Add receptor data to test effect analysis
        self.activity.receptor_profiles = {
            "5-HT2A": (1.2, 0.9, "agonist"),
            "D2": (5.6, 0.8, "agonist"),
        }
        
        effects = self.activity._analyze_effects()
        
        assert "psychedelic" in effects
        psychedelic = effects["psychedelic"]
        assert psychedelic["probability"] > 0
        assert psychedelic["confidence"] > 0
        assert "5-HT2A_agonist" in psychedelic["mechanisms"]

    def test_analyze_mechanisms(self):
        """Test _analyze_mechanisms method."""
        # Add receptor data to test mechanism analysis
        self.activity.receptor_profiles = {
            "AMPA": (2.3, 0.8, "modulator"),
            "D1": (3.4, 0.7, "agonist"),
        }
        
        mechanisms = self.activity._analyze_mechanisms()
        
        assert NootropicMechanism.GLUTAMATERGIC.value in mechanisms
        glutamatergic = mechanisms[NootropicMechanism.GLUTAMATERGIC.value]
        assert glutamatergic["strength"] > 0
        assert glutamatergic["confidence"] > 0
        assert "AMPA_modulator" in glutamatergic["active_pathways"]

    def test_predict_psychoactive_class(self):
        """Test _predict_psychoactive_class method."""
        # Add receptor data to test class prediction
        self.activity.receptor_profiles = {
            "5-HT2A": (1.2, 0.9, "agonist"),
            "5-HT2C": (3.4, 0.8, "agonist"),
        }
        
        prediction = self.activity._predict_psychoactive_class()
        
        assert "primary" in prediction
        assert prediction["primary"] == "psychedelic"
        assert "scores" in prediction
        assert "psychedelic" in prediction["scores"]

    def test_predict_nootropic_effects(self):
        """Test _predict_nootropic_effects method."""
        # Add receptor data to test effect prediction
        self.activity.receptor_profiles = {
            "AMPA": (2.3, 0.8, "modulator"),
            "NMDA": (4.5, 0.7, "modulator"),
        }
        
        predictions = self.activity._predict_nootropic_effects()
        
        assert NootropicMechanism.GLUTAMATERGIC.value in predictions
        prediction = predictions[NootropicMechanism.GLUTAMATERGIC.value]
        assert prediction["score"] > 0
        assert prediction["confidence"] > 0
        assert len(prediction["active_pathways"]) > 0

    def test_predict_activity_risks(self):
        """Test _predict_activity_risks method."""
        risks = self.activity._predict_activity_risks()
        
        assert "class_risks" in risks
        class_risks = risks["class_risks"]
        assert "psychosis" in class_risks
        assert "hppd" in class_risks

    def test_merge_activity_data(self):
        """Test merge_activity_data method."""
        other = ActivityProfileMixin()
        
        # Set up other instance data
        other.psychoactive_class = PsychoactiveClass.STIMULANT
        other.effect_profile = {
            "visual": (0.9, 0.95),    # Higher confidence
            "energy": (0.8, 0.9),     # New effect
        }
        other.nootropic_mechanisms = {
            NootropicMechanism.DOPAMINERGIC,
            NootropicMechanism.CHOLINERGIC,    # New mechanism
        }
        other.cognitive_effects = {
            "focus": (0.8, 0.9),      # Higher confidence
            "alertness": (0.7, 0.8),  # New effect
        }
        
        # Merge data
        self.activity.merge_activity_data(other)
        
        # Check psychoactive class merging
        assert self.activity.psychoactive_class == PsychoactiveClass.STIMULANT
        assert PsychoactiveClass.PSYCHEDELIC in self.activity.secondary_classes
        
        # Check effect profile merging
        assert self.activity.effect_profile["visual"] == (0.9, 0.95)    # Updated
        assert self.activity.effect_profile["energy"] == (0.8, 0.9)     # Added
        assert self.activity.effect_profile["euphoria"] == (0.6, 0.8)   # Unchanged
        
        # Check nootropic mechanism merging
        assert NootropicMechanism.CHOLINERGIC in self.activity.nootropic_mechanisms
        assert NootropicMechanism.GLUTAMATERGIC in self.activity.nootropic_mechanisms
        assert NootropicMechanism.DOPAMINERGIC in self.activity.nootropic_mechanisms
        
        # Check cognitive effects merging
        assert self.activity.cognitive_effects["focus"] == (0.8, 0.9)      # Updated
        assert self.activity.cognitive_effects["alertness"] == (0.7, 0.8)  # Added
        assert self.activity.cognitive_effects["memory"] == (0.5, 0.7)     # Unchanged


if __name__ == "__main__":
    pytest.main([__file__])
