"""Integration tests for compound analysis functionality."""

import pytest

from ..base import CompoundAnalysis


@pytest.fixture
def test_compound_data():
    """Create test compound data."""
    return {
        # Main compound data
        "smiles": "CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1",  # LSD
        "name": "LSD",
        "cas_number": "50-37-3",
        
        # Binding data
        "targets": [
            {
                "common_name": "5-HT2A",
                "affinity_value": 1.2,
                "affinity_type": "Ki",
                "affinity_unit": "nM",
                "activity_type": "agonist",
                "confidence": 0.95,
                "assay_details": {
                    "family": "serotonin",
                    "effects": ["psychedelic", "mood"],
                    "mechanisms": ["receptor activation"],
                },
            },
            {
                "common_name": "5-HT1A",
                "affinity_value": 8.5,
                "affinity_type": "Ki",
                "affinity_unit": "nM",
                "activity_type": "partial agonist",
                "confidence": 0.9,
                "assay_details": {
                    "family": "serotonin",
                    "effects": ["mood", "anxiety"],
                    "mechanisms": ["receptor modulation"],
                },
            },
        ],
        
        # Activity data
        "primary_activity": 1.2,  # nM at 5-HT2A
        "mechanism_of_action": "serotonin 5-HT2A receptor agonist",
        "effect_profile": {
            "psychedelic": (0.9, 0.95),
            "mood": (0.8, 0.9),
            "cognition": (0.6, 0.7),
            "perception": (0.85, 0.9),
            "anxiety": (-0.3, 0.6),
        },
        
        # Safety data
        "safety_profile": {
            "known_risks": [
                {
                    "description": "Psychological distress",
                    "severity": "moderate",
                    "confidence": 0.9,
                    "mechanism": "5-HT2A activation",
                    "conditions": ["predisposed individuals"],
                },
            ],
            "known_interactions": [
                {
                    "description": "Lithium",
                    "severity": "severe",
                    "confidence": 0.95,
                    "mechanism": "serotonergic effects",
                    "recommendations": ["contraindicated"],
                },
            ],
        },
        
        # Reference compounds
        "reference_compounds": [
            {
                "name": "Psilocin",
                "smiles": "CN(C)CCc1c[nH]c2cccc(O)c12",
                "primary_activity": 10.0,
            },
            {
                "name": "DMT",
                "smiles": "CN(C)CCc1c[nH]c2ccccc12",
                "primary_activity": 100.0,
            },
        ],
        
        # Experimental data
        "experimental_data": {
            "bioavailability": 0.7,
            "half_life": 3.6,  # hours
            "volume_distribution": 0.8,  # L/kg
        },
    }


def test_full_analysis(test_compound_data):
    """Test running all analyses together."""
    # Create compound with all data
    compound = CompoundAnalysis()
    for key, value in test_compound_data.items():
        setattr(compound, key, value)
    
    # Run all analyses
    results = compound.analyze_all()
    
    # Check all analysis types are present
    assert "binding" in results
    assert "activity" in results
    assert "safety" in results
    assert "properties" in results
    assert "sar" in results

    # Check binding analysis
    binding = results["binding"]
    assert binding["strongest_binding"]["target"] == "5-HT2A"
    assert binding["target_selectivity"]["selectivity"] == "moderately selective"

    # Check activity analysis
    activity = results["activity"]
    assert activity["primary_activity"]["mechanism"] == "serotonin 5-HT2A receptor agonist"
    assert "psychedelic" in activity["effect_profile"]["positive_effects"]

    # Check safety analysis
    safety = results["safety"]
    assert any(
        "Psychological distress" in r["description"]
        for r in safety["risk_assessment"]["risk_factors"]
    )
    assert safety["interactions"]["severe"]

    # Check property analysis
    properties = results["properties"]
    assert properties["drug_likeness"]["overall"]["classification"] in [
        "highly drug-like",
        "moderately drug-like",
    ]
    assert properties["bioavailability"]["classification"] in ["high", "moderate"]

    # Check SAR analysis
    sar = results["sar"]
    assert len(sar["pharmacophores"]) > 0
    assert len(sar["activity_cliffs"]) > 0


def test_analysis_dependencies():
    """Test that analyses handle missing dependencies correctly."""
    # Create compound with minimal data
    compound = CompoundAnalysis()
    compound.smiles = "CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1"
    
    # Run all analyses
    results = compound.analyze_all()
    
    # Binding analysis should be limited
    assert results["binding"]["target_selectivity"]["selectivity"] == "unknown"
    
    # Activity analysis should be limited
    assert results["activity"]["primary_activity"]["mechanism"] == "unknown"
    
    # Safety analysis should be limited
    assert results["safety"]["risk_assessment"]["overall_risk"]["level"] == "unknown"
    
    # Property analysis should still work (uses SMILES)
    assert results["properties"]["drug_likeness"]["overall"]["classification"] != "unknown"
    
    # SAR analysis should be limited
    assert not results["sar"]["activity_cliffs"]


def test_analysis_caching():
    """Test that analysis results are cached and invalidated properly."""
    compound = CompoundAnalysis()
    for key, value in test_compound_data.items():
        setattr(compound, key, value)
    
    # First analysis should compute and cache
    compound.analyze_all()
    
    # Get cached results
    binding1 = compound._binding_analysis
    activity1 = compound._activity_analysis
    safety1 = compound._safety_analysis
    property1 = compound._property_analysis
    sar1 = compound._sar_analysis
    
    # Second analysis should use cache
    compound.analyze_all()
    
    # Should be same objects in memory
    assert compound._binding_analysis is binding1
    assert compound._activity_analysis is activity1
    assert compound._safety_analysis is safety1
    assert compound._property_analysis is property1
    assert compound._sar_analysis is sar1

    # Modify data
    compound.targets.append({
        "common_name": "D2",
        "affinity_value": 50.0,
        "affinity_type": "Ki",
        "affinity_unit": "nM",
        "activity_type": "antagonist",
        "confidence": 0.8,
    })
    
    # Analysis should recompute
    compound.analyze_all()
    
    # Should be new objects
    assert compound._binding_analysis is not binding1
    assert compound._activity_analysis is not activity1
    assert compound._safety_analysis is not safety1
    assert compound._property_analysis is not property1
    assert compound._sar_analysis is not sar1


def test_analysis_summary(test_compound_data):
    """Test getting analysis summary."""
    compound = CompoundAnalysis()
    for key, value in test_compound_data.items():
        setattr(compound, key, value)
    
    # Get summary
    summary = compound.get_analysis_summary()
    
    # Check summary structure
    assert "binding" in summary
    assert "activity" in summary
    assert "safety" in summary
    assert "properties" in summary
    assert "sar" in summary

    # Check summary content
    assert summary["binding"]["strongest"]["target"] == "5-HT2A"
    assert summary["activity"]["primary"]["mechanism"] == "serotonin 5-HT2A receptor agonist"
    assert summary["safety"]["risk_level"] in ["low", "moderate", "high", "severe"]
    assert summary["properties"]["drug_likeness"] in [
        "highly drug-like",
        "moderately drug-like",
        "weakly drug-like",
        "non-drug-like",
    ]
    assert isinstance(summary["sar"]["pharmacophores"], int)


def test_key_findings(test_compound_data):
    """Test getting key findings."""
    compound = CompoundAnalysis()
    for key, value in test_compound_data.items():
        setattr(compound, key, value)
    
    # Get findings
    findings = compound.get_key_findings()
    
    # Check findings structure
    assert "highlights" in findings
    assert "concerns" in findings
    assert "opportunities" in findings

    # Check findings content
    assert any("5-HT2A" in h for h in findings["highlights"])
    assert any("safety" in c.lower() for c in findings["concerns"])
    assert any("activity cliff" in o.lower() for o in findings["opportunities"])


def test_optimization_suggestions(test_compound_data):
    """Test getting optimization suggestions."""
    compound = CompoundAnalysis()
    for key, value in test_compound_data.items():
        setattr(compound, key, value)
    
    # Get suggestions
    suggestions = compound.get_optimization_suggestions()
    
    # Check suggestions structure
    assert "binding" in suggestions
    assert "properties" in suggestions
    assert "safety" in suggestions

    # Check suggestions content
    assert any("selectivity" in s.lower() for s in suggestions["binding"])
    assert any(isinstance(s, str) for s in suggestions["properties"])
    assert any(isinstance(s, str) for s in suggestions["safety"])
