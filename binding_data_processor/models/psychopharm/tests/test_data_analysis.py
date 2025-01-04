"""Tests for data analysis functionality."""

import pytest
import numpy as np
from collections import Counter

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..analysis import (
    DataAnalyzer,
    AnalysisResult,
    BindingAnalysis,
    ActivityAnalysis,
    SafetyAnalysis,
)


@pytest.fixture
def test_compounds():
    """Create test compounds fixture."""
    compounds = []
    
    # Caffeine
    caffeine = PsychoactiveCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    caffeine.psychoactive_class = PsychoactiveClass.STIMULANT
    caffeine.add_receptor_binding(
        "A2A",
        affinity=0.8,
        confidence=0.95,
        activity="antagonist"
    )
    caffeine.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
    }
    compounds.append(caffeine)
    
    # Amphetamine
    amphetamine = PsychoactiveCompound(
        name="Amphetamine",
        smiles="CC(N)CC1=CC=CC=C1",
        cas_number="300-62-9",
    )
    amphetamine.psychoactive_class = PsychoactiveClass.STIMULANT
    amphetamine.add_receptor_binding(
        "DAT",
        affinity=0.05,
        confidence=0.95,
        activity="inhibitor"
    )
    amphetamine.add_receptor_binding(
        "NET",
        affinity=0.07,
        confidence=0.90,
        activity="inhibitor"
    )
    amphetamine.safety_alerts = {
        "addiction": RiskLevel.HIGH,
        "cardiovascular": RiskLevel.HIGH,
    }
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def analyzer():
    """Create a test data analyzer fixture."""
    return DataAnalyzer()


class TestDataAnalyzer:
    """Tests for DataAnalyzer class."""

    def test_initialization(self, analyzer):
        """Test initialization of DataAnalyzer."""
        assert isinstance(analyzer.binding_analyzer, BindingAnalysis)
        assert isinstance(analyzer.activity_analyzer, ActivityAnalysis)
        assert isinstance(analyzer.safety_analyzer, SafetyAnalysis)
        assert analyzer.stats == {}

    def test_binding_profile_analysis(self, test_compounds, analyzer):
        """Test binding profile analysis."""
        # Analyze binding profiles
        results = analyzer.analyze_binding_profiles(test_compounds)
        
        # Check receptor coverage
        assert "receptor_coverage" in results
        assert "A2A" in results["receptor_coverage"]
        assert "DAT" in results["receptor_coverage"]
        assert "NET" in results["receptor_coverage"]
        
        # Check binding patterns
        assert "binding_patterns" in results
        assert results["binding_patterns"]["stimulants"]["primary_targets"] == ["DAT", "NET"]
        
        # Check selectivity analysis
        assert "selectivity_analysis" in results
        assert results["selectivity_analysis"]["Amphetamine"]["DAT_vs_NET"] < 2.0

    def test_activity_pattern_analysis(self, test_compounds, analyzer):
        """Test activity pattern analysis."""
        # Analyze activity patterns
        results = analyzer.analyze_activity_patterns(test_compounds)
        
        # Check class distribution
        assert "class_distribution" in results
        assert results["class_distribution"][PsychoactiveClass.STIMULANT] == 2
        
        # Check mechanism patterns
        assert "mechanism_patterns" in results
        assert "antagonist" in results["mechanism_patterns"]
        assert "inhibitor" in results["mechanism_patterns"]
        
        # Check activity correlations
        assert "activity_correlations" in results
        assert "DAT_NET_correlation" in results["activity_correlations"]

    def test_safety_trend_analysis(self, test_compounds, analyzer):
        """Test safety trend analysis."""
        # Analyze safety trends
        results = analyzer.analyze_safety_trends(test_compounds)
        
        # Check risk distribution
        assert "risk_distribution" in results
        assert results["risk_distribution"][RiskLevel.HIGH] == 3  # Total high risks
        
        # Check common alerts
        assert "common_alerts" in results
        assert results["common_alerts"]["addiction"]["count"] == 1
        assert results["common_alerts"]["anxiety"]["count"] == 1
        
        # Check risk patterns
        assert "risk_patterns" in results
        assert "stimulant_risks" in results["risk_patterns"]
        assert "cardiovascular" in results["risk_patterns"]["stimulant_risks"]

    def test_structure_activity_analysis(self, test_compounds, analyzer):
        """Test structure-activity relationship analysis."""
        # Analyze SAR
        results = analyzer.analyze_structure_activity(test_compounds)
        
        # Check structural features
        assert "structural_features" in results
        assert "aromatic_rings" in results["structural_features"]
        assert "amine_groups" in results["structural_features"]
        
        # Check activity correlations
        assert "feature_activity_correlations" in results
        assert "aromatic_vs_potency" in results["feature_activity_correlations"]
        
        # Check scaffold analysis
        assert "common_scaffolds" in results
        assert len(results["common_scaffolds"]) > 0

    def test_trend_detection(self, test_compounds, analyzer):
        """Test trend detection in compound data."""
        # Add temporal data
        for compound in test_compounds:
            compound.discovery_date = "2020-01-01"
            compound.last_updated = "2023-01-01"
        
        # Analyze trends
        results = analyzer.detect_trends(test_compounds)
        
        # Check temporal trends
        assert "temporal_trends" in results
        assert "discovery_rate" in results["temporal_trends"]
        assert "update_frequency" in results["temporal_trends"]
        
        # Check property trends
        assert "property_trends" in results
        assert "binding_trends" in results["property_trends"]
        assert "safety_trends" in results["property_trends"]

    def test_similarity_analysis(self, test_compounds, analyzer):
        """Test compound similarity analysis."""
        # Analyze similarities
        results = analyzer.analyze_similarities(test_compounds)
        
        # Check structural similarities
        assert "structural_similarities" in results
        assert isinstance(results["structural_similarities"], np.ndarray)
        assert results["structural_similarities"].shape == (2, 2)
        
        # Check pharmacological similarities
        assert "pharmacological_similarities" in results
        assert "binding_similarity" in results["pharmacological_similarities"]
        assert "effect_similarity" in results["pharmacological_similarities"]

    def test_statistical_analysis(self, test_compounds, analyzer):
        """Test statistical analysis of compound properties."""
        # Analyze statistics
        results = analyzer.analyze_statistics(test_compounds)
        
        # Check basic stats
        assert "property_stats" in results
        assert "affinity_stats" in results["property_stats"]
        assert "risk_stats" in results["property_stats"]
        
        # Check distributions
        assert "distributions" in results
        assert "affinity_distribution" in results["distributions"]
        assert "risk_distribution" in results["distributions"]
        
        # Check correlations
        assert "correlations" in results
        assert "affinity_risk_correlation" in results["correlations"]

    def test_batch_analysis(self, test_compounds, analyzer):
        """Test batch analysis functionality."""
        # Run batch analysis
        results = analyzer.analyze_batch(
            compounds=test_compounds,
            analyses=["binding", "activity", "safety"]
        )
        
        # Check batch results
        assert "binding_analysis" in results
        assert "activity_analysis" in results
        assert "safety_analysis" in results
        assert analyzer.stats["total_analyzed"] == 2

    def test_error_handling(self, test_compounds, analyzer):
        """Test error handling during analysis."""
        # Add invalid data
        test_compounds[0].receptor_profiles["invalid"] = {
            "affinity": "invalid",
            "confidence": None,
        }
        
        # Run analysis with invalid data
        results = analyzer.analyze_binding_profiles(test_compounds)
        
        # Check error handling
        assert "analysis_errors" in results
        assert "invalid_data_points" in analyzer.stats
        assert analyzer.stats["error_count"] > 0

    def test_analysis_performance(self, test_compounds, analyzer):
        """Test analysis performance monitoring."""
        # Run analysis with performance monitoring
        analyzer.analyze_batch(
            compounds=test_compounds,
            analyses=["binding", "activity", "safety"]
        )
        
        # Check performance stats
        assert "analysis_time" in analyzer.stats
        assert isinstance(analyzer.stats["analysis_time"], float)
        assert analyzer.stats["analysis_time"] >= 0


if __name__ == "__main__":
    pytest.main([__file__])
