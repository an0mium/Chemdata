"""Tests for property analysis functionality."""

import pytest

from ..property_analysis import PropertyAnalysisMixin


class TestCompound(PropertyAnalysisMixin):
    """Test class implementing PropertyAnalysisMixin."""

    def __init__(self, property_data=None, experimental_data=None):
        self._property_analysis = {}
        self.molecular_weight = property_data.get("molecular_weight") if property_data else None
        self.logp = property_data.get("logp") if property_data else None
        self.hbd = property_data.get("hbd") if property_data else None
        self.hba = property_data.get("hba") if property_data else None
        self.tpsa = property_data.get("tpsa") if property_data else None
        self.rotatable_bonds = property_data.get("rotatable_bonds") if property_data else None
        self.charge = property_data.get("charge") if property_data else None
        self.stereocenter_count = property_data.get("stereocenter_count") if property_data else None
        self.ring_count = property_data.get("ring_count") if property_data else None
        self.atom_count = property_data.get("atom_count") if property_data else None
        self.experimental_data = experimental_data


@pytest.fixture
def test_property_data():
    """Create test property data."""
    return {
        # Drug-like compound properties
        "drug_like": {
            "molecular_weight": 320.0,
            "logp": 2.8,
            "hbd": 2,
            "hba": 5,
            "tpsa": 90.0,
            "rotatable_bonds": 6,
            "charge": 0,
            "stereocenter_count": 2,
            "ring_count": 3,
            "atom_count": 45,
        },
        # Non-drug-like compound properties
        "non_drug_like": {
            "molecular_weight": 650.0,
            "logp": 6.5,
            "hbd": 8,
            "hba": 12,
            "tpsa": 180.0,
            "rotatable_bonds": 15,
            "charge": 2,
            "stereocenter_count": 0,
            "ring_count": 8,
            "atom_count": 85,
        },
    }


def test_physicochemical_analysis(test_property_data):
    """Test physicochemical property analysis."""
    compound = TestCompound(property_data=test_property_data["drug_like"])
    properties = compound._analyze_physicochemical()
    
    # Check all properties are present
    assert properties["molecular_weight"] == 320.0
    assert properties["logp"] == 2.8
    assert properties["hbd"] == 2
    assert properties["hba"] == 5
    assert properties["tpsa"] == 90.0
    assert properties["rotatable_bonds"] == 6
    assert properties["charge"] == 0
    assert properties["stereocenter_count"] == 2
    assert properties["ring_count"] == 3


def test_drug_likeness_analysis(test_property_data):
    """Test drug-likeness analysis."""
    # Test drug-like compound
    drug_like = TestCompound(property_data=test_property_data["drug_like"])
    dl_analysis = drug_like._analyze_drug_likeness()
    
    # Check Lipinski's Rule of 5
    assert dl_analysis["lipinski"]["mw_ok"]
    assert dl_analysis["lipinski"]["logp_ok"]
    assert dl_analysis["lipinski"]["hbd_ok"]
    assert dl_analysis["lipinski"]["hba_ok"]
    assert dl_analysis["lipinski"]["violations"] == 0

    # Check Veber rules
    assert dl_analysis["veber"]["rotatable_ok"]
    assert dl_analysis["veber"]["tpsa_ok"]

    # Check overall classification
    assert dl_analysis["overall"]["classification"] == "highly drug-like"
    assert dl_analysis["overall"]["score"] >= 0.8

    # Test non-drug-like compound
    non_drug_like = TestCompound(property_data=test_property_data["non_drug_like"])
    ndl_analysis = non_drug_like._analyze_drug_likeness()
    
    # Should fail Lipinski rules
    assert not ndl_analysis["lipinski"]["mw_ok"]
    assert not ndl_analysis["lipinski"]["logp_ok"]
    assert not ndl_analysis["lipinski"]["hbd_ok"]
    assert not ndl_analysis["lipinski"]["hba_ok"]
    assert ndl_analysis["lipinski"]["violations"] >= 3

    # Should fail Veber rules
    assert not ndl_analysis["veber"]["rotatable_ok"]
    assert not ndl_analysis["veber"]["tpsa_ok"]

    # Check overall classification
    assert ndl_analysis["overall"]["classification"] == "non-drug-like"
    assert ndl_analysis["overall"]["score"] <= 0.4


def test_property_alerts(test_property_data):
    """Test property-based alerts."""
    # Test non-drug-like compound (should have many alerts)
    compound = TestCompound(property_data=test_property_data["non_drug_like"])
    alerts = compound._analyze_property_alerts()
    
    # Check molecular weight alert
    mw_alert = next(a for a in alerts if a["type"] == "molecular_weight")
    assert mw_alert["description"] == "High molecular weight may reduce bioavailability"
    assert mw_alert["value"] == 650.0
    assert mw_alert["threshold"] == 500.0
    assert mw_alert["severity"] == "moderate"

    # Check LogP alert
    logp_alert = next(a for a in alerts if a["type"] == "logp")
    assert logp_alert["description"] == "High LogP may cause poor solubility"
    assert logp_alert["value"] == 6.5
    assert logp_alert["threshold"] == 5.0
    assert logp_alert["severity"] == "moderate"

    # Check rotatable bonds alert
    rot_alert = next(a for a in alerts if a["type"] == "rotatable_bonds")
    assert "flexibility" in rot_alert["description"].lower()
    assert rot_alert["value"] == 15
    assert rot_alert["threshold"] == 10
    assert rot_alert["severity"] == "moderate"


def test_bioavailability_analysis(test_property_data):
    """Test bioavailability analysis."""
    # Test drug-like compound
    drug_like = TestCompound(property_data=test_property_data["drug_like"])
    drug_like_bio = drug_like._analyze_bioavailability()
    
    # Should have good bioavailability
    assert drug_like_bio["classification"] in ["high", "moderate"]
    assert drug_like_bio["score"] >= 0.7
    assert not drug_like_bio["limiting_factors"]
    assert drug_like_bio["predictions"]["oral"]
    assert drug_like_bio["predictions"]["intestinal"]

    # Test non-drug-like compound
    non_drug_like = TestCompound(property_data=test_property_data["non_drug_like"])
    non_drug_like_bio = non_drug_like._analyze_bioavailability()
    
    # Should have poor bioavailability
    assert non_drug_like_bio["classification"] in ["low", "very low"]
    assert non_drug_like_bio["score"] <= 0.5
    assert len(non_drug_like_bio["limiting_factors"]) >= 3
    assert not non_drug_like_bio["predictions"]["oral"]
    assert not non_drug_like_bio["predictions"]["intestinal"]


def test_bbb_permeability(test_property_data):
    """Test BBB permeability analysis."""
    # Test drug-like compound
    drug_like = TestCompound(property_data=test_property_data["drug_like"])
    drug_like_bbb = drug_like._analyze_bioavailability()["bbb_permeability"]
    
    # Should have good BBB permeability
    assert drug_like_bbb["classification"] in ["high", "moderate"]
    assert drug_like_bbb["score"] >= 0.6

    # Test non-drug-like compound
    non_drug_like = TestCompound(property_data=test_property_data["non_drug_like"])
    non_drug_like_bbb = non_drug_like._analyze_bioavailability()["bbb_permeability"]
    
    # Should have poor BBB permeability
    assert non_drug_like_bbb["classification"] == "low"
    assert non_drug_like_bbb["score"] <= 0.4


def test_no_property_data():
    """Test analysis with no property data."""
    compound = TestCompound()
    
    # Should handle missing data gracefully
    properties = compound._analyze_physicochemical()
    assert all(v is None for v in properties.values())

    drug_likeness = compound._analyze_drug_likeness()
    assert drug_likeness["overall"]["classification"] == "unknown"
    assert drug_likeness["overall"]["score"] == 0

    alerts = compound._analyze_property_alerts()
    assert not alerts

    bioavailability = compound._analyze_bioavailability()
    assert bioavailability["classification"] == "unknown"
    assert bioavailability["score"] == 0
    assert not bioavailability["limiting_factors"]


def test_experimental_data_confidence(test_property_data):
    """Test confidence adjustment with experimental data."""
    # Test with experimental data
    compound_with_exp = TestCompound(
        property_data=test_property_data["drug_like"],
        experimental_data={"bioavailability": 0.8}
    )
    bio_with_exp = compound_with_exp._analyze_bioavailability()
    
    # Test without experimental data
    compound_no_exp = TestCompound(
        property_data=test_property_data["drug_like"]
    )
    bio_no_exp = compound_no_exp._analyze_bioavailability()
    
    # Confidence should be higher with experimental data
    assert bio_with_exp["confidence"] > bio_no_exp["confidence"]


def test_partial_property_data(test_property_data):
    """Test analysis with partial property data."""
    # Only molecular weight and LogP
    partial_data = {
        "molecular_weight": test_property_data["drug_like"]["molecular_weight"],
        "logp": test_property_data["drug_like"]["logp"],
    }
    compound = TestCompound(property_data=partial_data)
    
    # Should still analyze available properties
    properties = compound._analyze_physicochemical()
    assert properties["molecular_weight"] == partial_data["molecular_weight"]
    assert properties["logp"] == partial_data["logp"]
    assert properties["hbd"] is None

    # Should make limited drug-likeness assessment
    drug_likeness = compound._analyze_drug_likeness()
    assert drug_likeness["lipinski"]["mw_ok"]
    assert drug_likeness["lipinski"]["logp_ok"]
    assert drug_likeness["overall"]["classification"] != "unknown"

    # Should generate relevant alerts
    alerts = compound._analyze_property_alerts()
    alert_types = {a["type"] for a in alerts}
    assert "molecular_weight" in alert_types
    assert "logp" in alert_types
    assert "tpsa" not in alert_types
