"""Tests for compound functionality."""

import pytest
from datetime import datetime

from ..base import (
    PsychoactiveClass,
    NootropicMechanism,
    RiskLevel,
)
from ..compound import PsychoactiveCompound


class TestPsychoactiveCompound:
    """Tests for PsychoactiveCompound class."""

    def setup_method(self):
        """Set up test instance."""
        self.compound = PsychoactiveCompound(
            name="Test Compound",
            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            cas_number="58-08-2",
        )
        
        # Set up receptor data
        self.compound.add_receptor_binding(
            "5-HT2A",
            affinity=1.2,
            confidence=0.9,
            activity="agonist"
        )
        self.compound.add_receptor_binding(
            "D2",
            affinity=5.6,
            confidence=0.8,
            activity="antagonist"
        )
        
        # Set up activity data
        self.compound.psychoactive_class = PsychoactiveClass.STIMULANT
        self.compound.effect_profile = {
            "stimulation": (0.8, 0.9),
            "focus": (0.7, 0.8),
        }
        self.compound.nootropic_mechanisms = {
            NootropicMechanism.DOPAMINERGIC,
        }
        
        # Set up safety data
        self.compound.safety_alerts = {
            "anxiety": RiskLevel.MODERATE,
            "insomnia": RiskLevel.HIGH,
        }
        self.compound.contraindications = {
            "anxiety_disorders",
            "cardiovascular_disease",
        }
        
        # Set up community data
        self.compound.experience_reports = {
            "erowid_0": {
                "dosage": {
                    "route": "oral",
                    "amount": 200.0,  # mg
                },
                "effects": {
                    "stimulation": {"magnitude": 0.8, "confidence": 0.9},
                    "anxiety": {"magnitude": 0.4, "confidence": 0.7},
                },
            },
        }
        self.compound.dosage_data = {
            "oral": (100.0, 400.0, 200.0),  # min, max, recommended
        }
        
        # Set up enrichment data
        self.compound.external_ids = {
            "pubchem_cid": "2519",
            "chembl_id": "CHEMBL113",
        }
        self.compound.literature_data = {
            "pubmed_ids": ["12345678"],
            "citations": {
                "12345678": {
                    "title": "Caffeine Study",
                    "findings": {
                        "binding_data": {
                            "A2A": {"ki": 0.8, "confidence": 0.95},
                        },
                    },
                },
            },
        }

    def test_initialization(self):
        """Test initialization of PsychoactiveCompound."""
        compound = PsychoactiveCompound(
            name="New Compound",
            smiles="CCO",  # Ethanol
            cas_number="64-17-5",
        )
        
        # Check basic attributes
        assert compound.name == "New Compound"
        assert compound.smiles == "CCO"
        assert compound.cas_number == "64-17-5"
        
        # Check mixin initializations
        assert compound.receptor_profiles == {}
        assert compound.psychoactive_class == PsychoactiveClass.UNKNOWN
        assert compound.safety_alerts == {}
        assert compound.experience_reports == {}
        assert compound.external_ids == {}

    @pytest.mark.parametrize("attribute,value,expected", [
        ("name", "Updated Name", "Updated Name"),
        ("smiles", "C1=CC=CC=C1", "C1=CC=CC=C1"),  # Benzene
        ("cas_number", "71-43-2", "71-43-2"),
    ])
    def test_basic_attributes(self, attribute, value, expected):
        """Test setting basic compound attributes."""
        setattr(self.compound, attribute, value)
        assert getattr(self.compound, attribute) == expected

    def test_get_compound_dict(self):
        """Test get_compound_dict method."""
        data = self.compound.get_compound_dict()
        
        assert isinstance(data, dict)
        assert "basic" in data
        assert "binding" in data
        assert "activity" in data
        assert "safety" in data
        assert "community" in data
        assert "enrichment" in data
        
        # Check basic data
        basic = data["basic"]
        assert basic["name"] == "Test Compound"
        assert basic["smiles"] == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        assert basic["cas_number"] == "58-08-2"
        
        # Check binding data
        binding = data["binding"]
        assert "5-HT2A" in binding["profiles"]
        assert binding["profiles"]["5-HT2A"]["affinity"] == 1.2
        
        # Check activity data
        activity = data["activity"]
        assert activity["primary_class"] == PsychoactiveClass.STIMULANT.value
        assert "stimulation" in activity["effects"]
        
        # Check safety data
        safety = data["safety"]
        assert safety["alerts"]["anxiety"] == RiskLevel.MODERATE.value
        assert "anxiety_disorders" in safety["contraindications"]
        
        # Check community data
        community = data["community"]
        assert "erowid_0" in community["reports"]
        assert community["dosage"]["oral"]["recommended"] == 200.0
        
        # Check enrichment data
        enrichment = data["enrichment"]
        assert enrichment["identifiers"]["pubchem_cid"] == "2519"
        assert "12345678" in enrichment["literature"]["pubmed_ids"]

    def test_update_timestamp(self):
        """Test timestamp updates."""
        # Set initial timestamp
        initial_time = datetime(2023, 1, 1)
        self.compound.last_updated = initial_time
        
        # Update compound data
        self.compound.add_receptor_binding(
            "NEW_RECEPTOR",
            affinity=1.0,
            confidence=0.9,
            activity="agonist"
        )
        
        # Check timestamp was updated
        assert self.compound.last_updated > initial_time
        assert isinstance(self.compound.last_updated, datetime)

    def test_merge_compound_data(self):
        """Test merge_compound_data method."""
        other = PsychoactiveCompound(
            name="Test Compound",
            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            cas_number="58-08-2",
        )
        
        # Add some overlapping and new data
        other.add_receptor_binding(
            "5-HT2A",
            affinity=1.0,  # Better affinity
            confidence=0.95,  # Higher confidence
            activity="partial_agonist"  # Different activity
        )
        other.add_receptor_binding(
            "NMDA",  # New receptor
            affinity=10.3,
            confidence=0.7,
            activity="antagonist"
        )
        
        other.psychoactive_class = PsychoactiveClass.STIMULANT
        other.effect_profile = {
            "stimulation": (0.9, 0.95),  # Higher confidence
            "alertness": (0.8, 0.9),  # New effect
        }
        
        other.safety_alerts = {
            "anxiety": RiskLevel.HIGH,  # Higher risk
            "tachycardia": RiskLevel.MODERATE,  # New alert
        }
        
        # Merge data
        self.compound.merge_compound_data(other)
        
        # Check binding merging
        assert self.compound.receptor_profiles["5-HT2A"][0] == 1.0  # Better affinity
        assert self.compound.receptor_profiles["5-HT2A"][1] == 0.95  # Higher confidence
        assert "NMDA" in self.compound.receptor_profiles  # New receptor
        
        # Check activity merging
        assert self.compound.effect_profile["stimulation"] == (0.9, 0.95)  # Higher confidence
        assert "alertness" in self.compound.effect_profile  # New effect
        
        # Check safety merging
        assert self.compound.safety_alerts["anxiety"] == RiskLevel.HIGH  # Higher risk
        assert self.compound.safety_alerts["tachycardia"] == RiskLevel.MODERATE  # New alert


if __name__ == "__main__":
    pytest.main([__file__])
