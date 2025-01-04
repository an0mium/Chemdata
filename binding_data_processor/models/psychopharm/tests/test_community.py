"""Tests for community data integration functionality."""

import pytest
from datetime import datetime

from ..base import RiskLevel
from ..community import CommunityDataMixin


class TestCommunityDataMixin:
    """Tests for CommunityDataMixin class."""

    def setup_method(self):
        """Set up test instance."""
        self.community = CommunityDataMixin()
        
        # Set up test data
        self.community.experience_reports = {
            "psychonaut_0": {
                "dosage": {
                    "route": "oral",
                    "amount": 25.0,
                },
                "timeline": {
                    "onset": {"onset": 30.0, "duration": 240.0},
                    "peak": {"onset": 90.0, "duration": 180.0},
                },
                "effects": {
                    "visual": {"magnitude": 0.8, "confidence": 0.9},
                    "euphoria": {"magnitude": 0.7, "confidence": 0.8},
                },
                "risks": {
                    "anxiety": "MODERATE",
                    "insomnia": "LOW",
                },
            },
        }
        self.community.report_sources = {"psychonaut"}
        self.community.total_reports = 1
        self.community.last_report_date = datetime(2023, 1, 1)
        
        # Set up dosage data
        self.community.dosage_data = {
            "oral": (20.0, 30.0, 25.0),  # min, max, recommended
        }
        self.community.duration_data = {
            "onset": (30.0, 240.0),  # onset, duration
            "peak": (90.0, 180.0),
        }
        
        # Set up effect data
        self.community.reported_effects = {
            "visual": [(0.8, 0.9)],
            "euphoria": [(0.7, 0.8)],
        }
        self.community.effect_frequencies = {
            "visual": 1.0,
            "euphoria": 1.0,
        }
        
        # Set up safety data
        self.community.reported_risks = {
            "anxiety": RiskLevel.MODERATE,
            "insomnia": RiskLevel.LOW,
        }
        self.community.contraindications = {
            "anxiety_disorders",
        }
        self.community.interactions = {
            "alcohol": {
                "risk_level": RiskLevel.HIGH.value,
                "mechanism": "CNS_depression",
            },
        }

    def test_initialization(self):
        """Test initialization of CommunityDataMixin."""
        community = CommunityDataMixin()
        assert community.experience_reports == {}
        assert community.report_sources == set()
        assert community.total_reports == 0
        assert community.last_report_date is None
        assert community.dosage_data == {}
        assert community.duration_data == {}
        assert community.reported_effects == {}
        assert community.effect_frequencies == {}
        assert community.reported_risks == {}
        assert community.contraindications == set()
        assert community.interactions == {}

    def test_get_community_dict(self):
        """Test get_community_dict method."""
        data = self.community.get_community_dict()
        
        assert isinstance(data, dict)
        assert "reports" in data
        assert "dosage" in data
        assert "duration" in data
        assert "effects" in data
        assert "safety" in data
        
        # Check report data
        reports = data["reports"]
        assert reports["total"] == 1
        assert "psychonaut" in reports["sources"]
        assert reports["last_date"] == "2023-01-01T00:00:00"
        
        # Check dosage data
        dosage = data["dosage"]["oral"]
        assert dosage["min"] == 20.0
        assert dosage["max"] == 30.0
        assert dosage["recommended"] == 25.0
        
        # Check duration data
        duration = data["duration"]["onset"]
        assert duration["onset"] == 30.0
        assert duration["duration"] == 240.0
        
        # Check effect data
        effects = data["effects"]["visual"]
        assert len(effects["scores"]) == 1
        assert effects["frequency"] == 1.0
        
        # Check safety data
        safety = data["safety"]
        assert safety["risks"]["anxiety"] == RiskLevel.MODERATE.value
        assert "anxiety_disorders" in safety["contraindications"]
        assert "alcohol" in safety["interactions"]

    def test_add_experience_report(self):
        """Test add_experience_report method."""
        report_data = {
            "dosage": {
                "route": "oral",
                "amount": 22.0,
            },
            "timeline": {
                "onset": {"onset": 35.0, "duration": 220.0},
            },
            "effects": {
                "visual": {"magnitude": 0.7, "confidence": 0.8},
                "focus": {"magnitude": 0.6, "confidence": 0.7},
            },
            "risks": {
                "anxiety": "HIGH",
                "headache": "LOW",
            },
            "interactions": {
                "caffeine": {
                    "risk_level": "MODERATE",
                    "mechanism": "stimulation",
                },
            },
        }
        report_date = datetime(2023, 2, 1)
        
        # Add report
        self.community.add_experience_report("erowid", report_data, report_date)
        
        # Check report tracking
        assert "erowid_0" in self.community.experience_reports
        assert "erowid" in self.community.report_sources
        assert self.community.total_reports == 2
        assert self.community.last_report_date == report_date
        
        # Check dosage updates
        oral_dose = self.community.dosage_data["oral"]
        assert oral_dose[0] == 20.0  # min unchanged
        assert oral_dose[1] == 30.0  # max unchanged
        assert oral_dose[2] == 23.5  # new average
        
        # Check effect updates
        assert len(self.community.reported_effects["visual"]) == 2
        assert "focus" in self.community.reported_effects
        assert self.community.effect_frequencies["visual"] == 1.0
        assert self.community.effect_frequencies["focus"] == 0.5
        
        # Check risk updates
        assert self.community.reported_risks["anxiety"] == RiskLevel.HIGH
        assert self.community.reported_risks["headache"] == RiskLevel.LOW
        assert "caffeine" in self.community.interactions

    def test_merge_community_data(self):
        """Test merge_community_data method."""
        other = CommunityDataMixin()
        
        # Set up other instance data
        other.experience_reports = {
            "tripsit_0": {
                "dosage": {"route": "oral", "amount": 28.0},
                "effects": {"visual": {"magnitude": 0.9, "confidence": 0.95}},
            },
        }
        other.report_sources = {"tripsit"}
        other.total_reports = 1
        other.last_report_date = datetime(2023, 3, 1)
        
        other.dosage_data = {
            "oral": (25.0, 35.0, 30.0),  # Higher doses
            "sublingual": (15.0, 25.0, 20.0),  # New ROA
        }
        
        other.reported_effects = {
            "visual": [(0.9, 0.95)],  # Higher confidence
            "creativity": [(0.8, 0.9)],  # New effect
        }
        
        other.reported_risks = {
            "anxiety": RiskLevel.HIGH,  # Higher risk
            "nausea": RiskLevel.LOW,    # New risk
        }
        
        # Merge data
        self.community.merge_community_data(other)
        
        # Check report merging
        assert "tripsit_0" in self.community.experience_reports
        assert "tripsit" in self.community.report_sources
        assert self.community.total_reports == 2
        assert self.community.last_report_date == datetime(2023, 3, 1)
        
        # Check dosage merging
        oral = self.community.dosage_data["oral"]
        assert oral[0] == 20.0  # Keep lower min
        assert oral[1] == 35.0  # Take higher max
        assert oral[2] == 27.5  # Average recommended
        assert "sublingual" in self.community.dosage_data
        
        # Check effect merging
        assert len(self.community.reported_effects["visual"]) == 2
        assert "creativity" in self.community.reported_effects
        assert self.community.effect_frequencies["visual"] == 1.0
        assert self.community.effect_frequencies["creativity"] == 0.5
        
        # Check risk merging
        assert self.community.reported_risks["anxiety"] == RiskLevel.HIGH
        assert self.community.reported_risks["nausea"] == RiskLevel.LOW
        assert self.community.reported_risks["insomnia"] == RiskLevel.LOW


if __name__ == "__main__":
    pytest.main([__file__])
