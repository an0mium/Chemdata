"""Tests for receptor binding functionality."""

import pytest
from ..binding import ReceptorProfileMixin


class TestReceptorProfileMixin:
    """Tests for ReceptorProfileMixin class."""

    def setup_method(self):
        """Set up test instance."""
        self.binding = ReceptorProfileMixin()
        
        # Set up test data
        self.test_profiles = {
            "5-HT2A": (1.2, 0.9, "agonist"),      # affinity, confidence, activity
            "D2": (5.6, 0.8, "antagonist"),
            "NMDA": (10.3, 0.7, "antagonist"),
        }
        self.binding.receptor_profiles = self.test_profiles.copy()

    def test_initialization(self):
        """Test initialization of ReceptorProfileMixin."""
        binding = ReceptorProfileMixin()
        assert binding.receptor_profiles == {}
        assert binding.primary_targets == set()
        assert binding.secondary_targets == set()

    def test_get_binding_dict(self):
        """Test get_binding_dict method."""
        data = self.binding.get_binding_dict()
        
        assert isinstance(data, dict)
        assert "profiles" in data
        assert "primary_targets" in data
        assert "secondary_targets" in data
        
        # Check profile data
        for receptor, (affinity, confidence, activity) in self.test_profiles.items():
            assert receptor in data["profiles"]
            profile = data["profiles"][receptor]
            assert profile["affinity"] == affinity
            assert profile["confidence"] == confidence
            assert profile["activity"] == activity

    def test_add_receptor_binding(self):
        """Test add_receptor_binding method."""
        # Add new receptor
        self.binding.add_receptor_binding(
            "5-HT1A",
            affinity=2.3,
            confidence=0.85,
            activity="partial_agonist"
        )
        
        assert "5-HT1A" in self.binding.receptor_profiles
        binding = self.binding.receptor_profiles["5-HT1A"]
        assert binding[0] == 2.3  # affinity
        assert binding[1] == 0.85  # confidence
        assert binding[2] == "partial_agonist"  # activity
        
        # Update existing with higher confidence
        self.binding.add_receptor_binding(
            "5-HT1A",
            affinity=2.0,
            confidence=0.9,
            activity="agonist"
        )
        
        binding = self.binding.receptor_profiles["5-HT1A"]
        assert binding[0] == 2.0  # New affinity
        assert binding[1] == 0.9  # Higher confidence
        assert binding[2] == "agonist"  # New activity
        
        # Try update with lower confidence (should not change)
        self.binding.add_receptor_binding(
            "5-HT1A",
            affinity=1.8,
            confidence=0.8,
            activity="antagonist"
        )
        
        binding = self.binding.receptor_profiles["5-HT1A"]
        assert binding[0] == 2.0  # Unchanged
        assert binding[1] == 0.9  # Unchanged
        assert binding[2] == "agonist"  # Unchanged

    def test_get_strongest_targets(self):
        """Test get_strongest_targets method."""
        # Add more test data
        self.binding.add_receptor_binding("5-HT1A", 15.0, 0.7, "agonist")
        self.binding.add_receptor_binding("SERT", 0.8, 0.95, "inhibitor")
        
        # Get strongest targets
        primary, secondary = self.binding.get_strongest_targets()
        
        # SERT should be primary (lowest affinity)
        assert "SERT" in primary
        
        # 5-HT2A should be secondary (second lowest affinity)
        assert "5-HT2A" in secondary
        
        # 5-HT1A should not be in either (high affinity)
        assert "5-HT1A" not in primary
        assert "5-HT1A" not in secondary

    def test_get_binding_profile(self):
        """Test get_binding_profile method."""
        profile = self.binding.get_binding_profile("5-HT2A")
        assert profile is not None
        assert profile[0] == 1.2  # affinity
        assert profile[1] == 0.9  # confidence
        assert profile[2] == "agonist"  # activity
        
        # Test non-existent receptor
        assert self.binding.get_binding_profile("not_a_receptor") is None

    def test_has_strong_binding(self):
        """Test has_strong_binding method."""
        # Test strong binding (low affinity)
        assert self.binding.has_strong_binding("5-HT2A")
        
        # Test weak binding (high affinity)
        self.binding.add_receptor_binding("weak", 500.0, 0.9, "agonist")
        assert not self.binding.has_strong_binding("weak")
        
        # Test non-existent receptor
        assert not self.binding.has_strong_binding("not_a_receptor")

    def test_merge_binding_data(self):
        """Test merge_binding_data method."""
        other = ReceptorProfileMixin()
        
        # Add some overlapping and new data
        other.add_receptor_binding("5-HT2A", 1.0, 0.95, "partial_agonist")  # Higher confidence
        other.add_receptor_binding("D2", 6.0, 0.7, "antagonist")  # Lower confidence
        other.add_receptor_binding("NEW", 3.0, 0.8, "agonist")  # New receptor
        
        # Merge data
        self.binding.merge_binding_data(other)
        
        # Check results
        assert self.binding.receptor_profiles["5-HT2A"][0] == 1.0  # Updated affinity
        assert self.binding.receptor_profiles["5-HT2A"][1] == 0.95  # Higher confidence
        assert self.binding.receptor_profiles["5-HT2A"][2] == "partial_agonist"  # Updated activity
        
        assert self.binding.receptor_profiles["D2"][0] == 5.6  # Unchanged
        assert self.binding.receptor_profiles["D2"][1] == 0.8  # Unchanged
        assert self.binding.receptor_profiles["D2"][2] == "antagonist"  # Unchanged
        
        assert "NEW" in self.binding.receptor_profiles  # New receptor added

    def test_validate_binding_data(self):
        """Test validate_binding_data method."""
        # Valid data
        assert self.binding.validate_binding_data()
        
        # Invalid affinity (negative)
        self.binding.receptor_profiles["invalid"] = (-1.0, 0.8, "agonist")
        assert not self.binding.validate_binding_data()
        
        # Invalid confidence (> 1.0)
        self.binding.receptor_profiles["invalid"] = (1.0, 1.1, "agonist")
        assert not self.binding.validate_binding_data()
        
        # Invalid activity (not a string)
        self.binding.receptor_profiles["invalid"] = (1.0, 0.8, 123)
        assert not self.binding.validate_binding_data()
        
        # Invalid tuple length
        self.binding.receptor_profiles["invalid"] = (1.0, 0.8)
        assert not self.binding.validate_binding_data()


if __name__ == "__main__":
    pytest.main([__file__])
