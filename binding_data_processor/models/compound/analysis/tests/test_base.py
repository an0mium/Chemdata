"""Tests for base compound analysis functionality."""

import pytest
from rdkit import Chem

from ..base import CompoundAnalysis, AnalysisError


@pytest.fixture
def caffeine_compound():
    """Create test compound with Caffeine data."""
    compound = CompoundAnalysis()
    compound.smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
    compound.name = "Caffeine"
    compound.cas_number = "58-08-2"
    
    # Add test binding data
    compound.targets = [
        {
            "common_name": "A2A",
            "affinity_value": 0.8,
            "affinity_type": "Ki",
            "affinity_unit": "nM",
            "activity_type": "antagonist",
            "confidence": 0.95,
            "assay_details": {
                "family": "adenosine",
                "effects": ["stimulation", "wakefulness"],
                "mechanisms": ["competitive antagonism"],
            },
            "reference_dois": ["10.1234/example.1", "10.1234/example.2"],
            "experimental_conditions": {"pH": 7.4, "temp": 37},
        },
        {
            "common_name": "A1",
            "affinity_value": 12.0,
            "affinity_type": "Ki",
            "affinity_unit": "nM",
            "activity_type": "antagonist",
            "confidence": 0.9,
            "assay_details": {
                "family": "adenosine",
                "effects": ["stimulation"],
                "mechanisms": ["competitive antagonism"],
            },
            "reference_dois": ["10.1234/example.3"],
            "experimental_conditions": {"pH": 7.4, "temp": 37},
        },
    ]

    # Add test activity data
    compound.primary_activity = 0.8  # nM
    compound.mechanism_of_action = "adenosine receptor antagonist"
    compound.effect_profile = {
        "stimulation": (0.8, 0.9),
        "wakefulness": (0.7, 0.8),
        "focus": (0.6, 0.7),
    }

    # Add test safety data
    compound.safety_profile = {
        "known_risks": [
            {
                "description": "Anxiety at high doses",
                "severity": "moderate",
                "confidence": 0.8,
            }
        ],
        "known_interactions": [
            {
                "description": "May interact with other stimulants",
                "severity": "moderate",
                "confidence": 0.7,
            }
        ],
    }

    # Add test property data using RDKit
    mol = Chem.MolFromSmiles(compound.smiles)
    compound.molecular_weight = Chem.Descriptors.ExactMolWt(mol)
    compound.logp = Chem.Descriptors.MolLogP(mol)
    compound.hbd = Chem.Descriptors.NumHDonors(mol)
    compound.hba = Chem.Descriptors.NumHAcceptors(mol)
    compound.tpsa = Chem.Descriptors.TPSA(mol)
    compound.rotatable_bonds = Chem.Descriptors.NumRotatableBonds(mol)
    compound.charge = Chem.GetFormalCharge(mol)
    compound.atom_count = mol.GetNumAtoms()
    compound.ring_count = Chem.Descriptors.RingCount(mol)
    compound.stereocenter_count = len(Chem.FindMolChiralCenters(mol))

    # Add test reference compounds
    compound.reference_compounds = [
        {
            "name": "Theophylline",
            "smiles": "CN1C=NC2=C1C(=O)NC(=O)N2C",
            "primary_activity": 1.2,
        },
        {
            "name": "Theobromine",
            "smiles": "CN1C=NC2=C1C(=O)NC(=O)N2",
            "primary_activity": 2.5,
        },
    ]

    return compound


def test_initialization(lsd_data):
    """Test compound analysis initialization."""
    # Test with complete data
    compound = CompoundAnalysis()
    for key, value in lsd_data.items():
        setattr(compound, key, value)
    
    assert compound.smiles == lsd_data["smiles"]
    assert compound.name == lsd_data["name"]
    assert compound.targets == lsd_data["targets"]
    assert compound._binding_analysis == {}
    assert compound._activity_analysis == {}
    assert compound._safety_analysis == {}
    assert compound._property_analysis == {}
    assert compound._sar_analysis == {}


def test_data_validation(lsd_data):
    """Test input data validation."""
    compound = CompoundAnalysis()
    
    # Test invalid SMILES
    with pytest.raises(AnalysisError) as exc:
        compound.smiles = "invalid_smiles"
        compound.validate_structure()
    assert "Invalid SMILES" in str(exc.value)

    # Test invalid target data
    with pytest.raises(AnalysisError) as exc:
        compound.targets = [{"invalid": "data"}]
        compound.validate_targets()
    assert "Invalid target data" in str(exc.value)

    # Test valid data
    compound = CompoundAnalysis()
    for key, value in lsd_data.items():
        setattr(compound, key, value)
    compound.validate_all()  # Should not raise


def test_analyze_all(caffeine_compound):
    """Test running all analyses."""
    results = caffeine_compound.analyze_all()
    
    # Check all analysis types are present
    assert "binding" in results
    assert "activity" in results
    assert "safety" in results
    assert "properties" in results
    assert "sar" in results

    # Check binding analysis
    binding = results["binding"]
    assert binding["strongest_binding"]["target"] == "A2A"
    assert binding["target_selectivity"]["selectivity"] == "moderately selective"

    # Check activity analysis
    activity = results["activity"]
    assert activity["primary_activity"]["mechanism"] == "adenosine receptor antagonist"
    assert "stimulation" in activity["effect_profile"]

    # Check safety analysis
    safety = results["safety"]
    assert len(safety["property_alerts"]) > 0
    assert safety["risk_assessment"]["overall_risk"]["level"] in ["low", "moderate"]

    # Check property analysis
    properties = results["properties"]
    assert "drug_likeness" in properties
    assert "bioavailability" in properties

    # Check SAR analysis
    sar = results["sar"]
    assert len(sar["pharmacophores"]) > 0
    assert len(sar["activity_cliffs"]) > 0


def test_analysis_caching(lsd_data):
    """Test analysis result caching."""
    compound = CompoundAnalysis()
    for key, value in lsd_data.items():
        setattr(compound, key, value)
    
    # Initial analysis
    result1 = compound.analyze_binding()
    assert compound._binding_analysis == result1

    # Cached result
    result2 = compound.analyze_binding()
    assert result2 is result1  # Same object in memory

    # Cache invalidation
    compound.targets = compound.targets + [{
        "common_name": "D2",
        "affinity_value": 50.0,
        "affinity_type": "Ki",
        "affinity_unit": "nM",
        "activity_type": "antagonist",
        "confidence": 0.8,
    }]
    result3 = compound.analyze_binding()
    assert result3 is not result1  # New analysis


def test_analysis_dependencies(lsd_data):
    """Test handling of analysis dependencies."""
    compound = CompoundAnalysis()
    
    # Test with minimal data
    compound.smiles = lsd_data["smiles"]
    
    # Property analysis should work (only needs SMILES)
    properties = compound.analyze_properties()
    assert properties["drug_likeness"]["overall"]["classification"] != "unknown"

    # Activity analysis should be limited
    activity = compound.analyze_activity()
    assert activity["primary_activity"]["mechanism"] == "unknown"

    # Add activity data
    compound.primary_activity = lsd_data["primary_activity"]
    compound.mechanism_of_action = lsd_data["mechanism_of_action"]
    
    # Activity analysis should now work
    activity = compound.analyze_activity()
    assert activity["primary_activity"]["mechanism"] == lsd_data["mechanism_of_action"]


def test_error_handling():
    """Test error handling during analysis."""
    compound = CompoundAnalysis()
    
    # Missing required data
    with pytest.raises(AnalysisError) as exc:
        compound.analyze_all()
    assert "Missing required data" in str(exc.value)

    # Invalid property data
    compound.smiles = "CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1"
    compound.molecular_weight = "invalid"  # Should be float
    with pytest.raises(AnalysisError) as exc:
        compound.analyze_properties()
    assert "Invalid property data" in str(exc.value)


def test_result_formatting(lsd_data):
    """Test analysis result formatting."""
    compound = CompoundAnalysis()
    for key, value in lsd_data.items():
        setattr(compound, key, value)
    
    # Get full analysis
    results = compound.analyze_all()
    
    # Check structure
    assert isinstance(results, dict)
    assert all(k in results for k in [
        "binding",
        "activity",
        "safety",
        "properties",
        "sar",
    ])

    # Check content types
    assert isinstance(results["binding"], dict)
    assert isinstance(results["activity"], dict)
    assert isinstance(results["safety"], dict)
    assert isinstance(results["properties"], dict)
    assert isinstance(results["sar"], dict)

    # Check nested structure
    binding = results["binding"]
    assert "strongest_binding" in binding
    assert "target_selectivity" in binding
    assert isinstance(binding["strongest_binding"], dict)
    assert isinstance(binding["target_selectivity"], dict)


def test_analysis_summary(caffeine_compound, lsd_data):
    """Test getting analysis summary."""
    # Test with Caffeine
    summary1 = caffeine_compound.get_analysis_summary()
    assert summary1["binding"]["strongest"]["target"] == "A2A"
    assert summary1["binding"]["selectivity"] == "moderately selective"
    assert summary1["activity"]["primary"]["mechanism"] == "adenosine receptor antagonist"
    assert summary1["activity"]["confidence"] > 0

    # Test with LSD
    compound = CompoundAnalysis()
    for key, value in lsd_data.items():
        setattr(compound, key, value)
    summary2 = compound.get_analysis_summary()
    assert summary2["compound_info"]["name"] == "LSD"
    assert summary2["binding_summary"]["strongest_target"] == "5-HT2A"
    assert summary2["activity_summary"]["primary_mechanism"] == "serotonin 5-HT2A receptor agonist"


def test_get_key_findings(caffeine_compound):
    """Test getting key findings."""
    findings = caffeine_compound.get_key_findings()
    
    # Check findings structure
    assert "highlights" in findings
    assert "concerns" in findings
    assert "opportunities" in findings

    # Check findings content
    assert any("A2A" in h for h in findings["highlights"])
    assert any("safety" in c.lower() for c in findings["concerns"])
    assert any("activity cliff" in o.lower() for o in findings["opportunities"])


def test_get_optimization_suggestions(caffeine_compound):
    """Test getting optimization suggestions."""
    suggestions = caffeine_compound.get_optimization_suggestions()
    
    # Check suggestions structure
    assert "binding" in suggestions
    assert "properties" in suggestions
    assert "safety" in suggestions

    # Check suggestions content
    assert any("selectivity" in s.lower() for s in suggestions["binding"])
    assert any(isinstance(s, str) for s in suggestions["properties"])
    assert any(isinstance(s, str) for s in suggestions["safety"])


def test_partial_analysis(lsd_data):
    """Test running partial analysis."""
    compound = CompoundAnalysis()
    for key, value in lsd_data.items():
        setattr(compound, key, value)
    
    # Run specific analyses
    results = compound.analyze_selected([
        "binding",
        "activity",
    ])
    
    # Should only include requested analyses
    assert "binding" in results
    assert "activity" in results
    assert "safety" not in results
    assert "properties" not in results
    assert "sar" not in results

    # Results should be complete
    binding = results["binding"]
    assert binding["strongest_binding"]["target"] == "5-HT2A"
    activity = results["activity"]
    assert activity["primary_activity"]["mechanism"] == "serotonin 5-HT2A receptor agonist"


def test_analysis_flags(lsd_data):
    """Test analysis control flags."""
    compound = CompoundAnalysis()
    for key, value in lsd_data.items():
        setattr(compound, key, value)
    
    # Skip validation
    results = compound.analyze_all(validate=False)
    assert results  # Should complete without validation

    # Skip cache
    result1 = compound.analyze_binding()
    result2 = compound.analyze_binding(use_cache=False)
    assert result2 is not result1  # Should be new analysis

    # Force recompute
    result3 = compound.analyze_binding()
    result4 = compound.analyze_binding(force_recompute=True)
    assert result4 is not result3  # Should be new analysis


def test_analysis_metadata(lsd_data):
    """Test analysis metadata handling."""
    compound = CompoundAnalysis()
    for key, value in lsd_data.items():
        setattr(compound, key, value)
    
    # Run analysis
    compound.analyze_all()
    
    # Check metadata
    metadata = compound.get_analysis_metadata()
    assert "timestamp" in metadata
    assert "version" in metadata
    assert "analysis_types" in metadata
    assert set(metadata["analysis_types"]) == {
        "binding",
        "activity",
        "safety",
        "properties",
        "sar",
    }
