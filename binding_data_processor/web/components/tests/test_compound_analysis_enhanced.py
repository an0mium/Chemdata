"""Tests for enhanced compound analysis component."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Fragments

from ...components.compound_analysis_enhanced import CompoundAnalysisEnhanced
from ....models.compound import Compound


@pytest.fixture
def component(tmp_path):
    """Create test component."""
    return CompoundAnalysisEnhanced(
        template_dir=tmp_path / "templates",
        static_dir=tmp_path / "static",
    )


@pytest.fixture
def compound():
    """Create test compound."""
    compound = Compound(
        name="Test Compound",
        smiles="CC1=CC=CC=C1",
        cas_number="100-00-0",
    )

    # Add binding data
    compound.binding_data = [
        {
            "target": "5-HT2A",
            "affinity": "1.2",
            "confidence": 0.9,
        },
        {
            "target": "D2",
            "affinity": "2.5",
            "confidence": 0.8,
        },
    ]

    # Add predictions
    compound.predictions = {
        "activity": {
            "stimulant": 0.7,
            "psychedelic": 0.3,
            "sedative": 0.2,
        },
        "safety": {
            "toxicity": 0.2,
            "addiction": 0.4,
            "overdose": 0.3,
        },
    }

    # Add social data
    now = datetime.now()
    compound.social_data = {
        "reddit": {
            "posts": [
                {
                    "id": "1",
                    "title": "Test Post 1",
                    "created_utc": (now - timedelta(days=2)).isoformat(),
                },
                {
                    "id": "2",
                    "title": "Test Post 2",
                    "created_utc": (now - timedelta(days=1)).isoformat(),
                },
            ],
        },
        "twitter": {
            "tweets": [
                {
                    "id": "1",
                    "text": "Test Tweet 1",
                    "created_at": (now - timedelta(days=2)).isoformat(),
                },
                {
                    "id": "2",
                    "text": "Test Tweet 2",
                    "created_at": (now - timedelta(days=1)).isoformat(),
                },
            ],
        },
    }

    return compound


@patch("flask.render_template")
def test_render_analysis(mock_render, component, compound):
    """Test rendering analysis interface."""
    # Mock template rendering
    mock_render.return_value = "<html>Test</html>"

    # Render analysis
    result = component.render_analysis(compound)

    # Check result
    assert result.success
    assert result.data["html"] == "<html>Test</html>"
    assert "analyses" in result.data
    assert "analysis_stats" in result.data

    # Check analyses
    analyses = result.data["analyses"]
    assert "binding_analysis" in analyses
    assert "activity_analysis" in analyses
    assert "safety_analysis" in analyses
    assert "sar_analysis" in analyses
    assert "property_analysis" in analyses
    assert "community_analysis" in analyses

    # Check stats
    assert component.analysis_stats["total_analyses"] == 1
    assert component.analysis_stats["binding_analyses"] == 1
    assert component.analysis_stats["activity_analyses"] == 1
    assert component.analysis_stats["safety_analyses"] == 1
    assert component.analysis_stats["sar_analyses"] == 1
    assert component.analysis_stats["property_analyses"] == 1
    assert component.analysis_stats["community_analyses"] == 1


def test_binding_analysis(component, compound):
    """Test binding analysis."""
    # Perform analysis
    result = component.render_analysis(
        compound,
        analysis_types=["binding"],
    )

    # Check result
    assert result.success
    analysis = result.data["analyses"]["binding_analysis"]

    # Check statistics
    stats = analysis["statistics"]
    assert stats["strongest_target"] == "5-HT2A"
    assert stats["strongest_affinity"] == 1.2
    assert stats["mean_affinity"] == pytest.approx(1.85)
    assert stats["mean_confidence"] == pytest.approx(0.85)
    assert stats["target_count"] == 2

    # Check selectivity
    selectivity = analysis["selectivity"]
    assert "5-HT2A vs D2" in selectivity
    assert "D2 vs 5-HT2A" in selectivity
    assert selectivity["5-HT2A vs D2"] == pytest.approx(2.5 / 1.2)


def test_activity_analysis(component, compound):
    """Test activity analysis."""
    # Perform analysis
    result = component.render_analysis(
        compound,
        analysis_types=["activity"],
    )

    # Check result
    assert result.success
    analysis = result.data["analyses"]["activity_analysis"]

    # Check statistics
    stats = analysis["statistics"]
    assert stats["primary_activity"] == "stimulant"
    assert stats["primary_score"] == 0.7
    assert stats["mean_score"] == pytest.approx(0.4)
    assert stats["activity_count"] == 3

    # Check patterns
    patterns = analysis["patterns"]
    assert "stimulant" in patterns["high"]
    assert "psychedelic" in patterns["moderate"]
    assert "sedative" in patterns["low"]


def test_safety_analysis(component, compound):
    """Test safety analysis."""
    # Perform analysis
    result = component.render_analysis(
        compound,
        analysis_types=["safety"],
    )

    # Check result
    assert result.success
    analysis = result.data["analyses"]["safety_analysis"]

    # Check statistics
    stats = analysis["statistics"]
    assert stats["highest_risk"] == "addiction"
    assert stats["highest_score"] == 0.4
    assert stats["mean_score"] == pytest.approx(0.3)
    assert stats["risk_count"] == 3

    # Check patterns
    patterns = analysis["patterns"]
    assert "addiction" in patterns["moderate"]
    assert "overdose" in patterns["moderate"]
    assert "toxicity" in patterns["low"]


def test_sar_analysis(component, compound):
    """Test SAR analysis."""
    # Perform analysis
    result = component.render_analysis(
        compound,
        analysis_types=["sar"],
    )

    # Check result
    assert result.success
    analysis = result.data["analyses"]["sar_analysis"]

    # Check fingerprint
    assert isinstance(analysis["fingerprint"], list)
    assert len(analysis["fingerprint"]) > 0

    # Check descriptors
    descriptors = analysis["descriptors"]
    assert descriptors["MW"] == pytest.approx(92.14)
    assert descriptors["LogP"] == pytest.approx(2.19)
    assert descriptors["TPSA"] == pytest.approx(0.0)
    assert descriptors["HBA"] == 0
    assert descriptors["HBD"] == 0
    assert descriptors["Rings"] == 1
    assert descriptors["AromaticRings"] == 1

    # Check fragments
    fragments = analysis["fragments"]
    assert fragments["Aromatic"] == 0
    assert fragments["Alkyl"] == 0
    assert fragments["Carbonyl"] == 0
    assert fragments["Carboxyl"] == 0
    assert fragments["Amine"] == 0
    assert fragments["Amide"] == 0


def test_property_analysis(component, compound):
    """Test property analysis."""
    # Perform analysis
    result = component.render_analysis(
        compound,
        analysis_types=["properties"],
    )

    # Check result
    assert result.success
    analysis = result.data["analyses"]["property_analysis"]

    # Check properties
    properties = analysis["properties"]
    assert properties["MW"] == pytest.approx(92.14)
    assert properties["LogP"] == pytest.approx(2.19)
    assert properties["TPSA"] == pytest.approx(0.0)
    assert properties["HBA"] == 0
    assert properties["HBD"] == 0
    assert properties["Rings"] == 1
    assert properties["AromaticRings"] == 1
    assert properties["HeavyAtoms"] == 7
    assert properties["Charge"] == 0

    # Check Lipinski
    lipinski = analysis["lipinski"]
    assert lipinski["MW_ok"] is True
    assert lipinski["LogP_ok"] is True
    assert lipinski["HBA_ok"] is True
    assert lipinski["HBD_ok"] is True
    assert lipinski["RotBonds_ok"] is True
    assert analysis["violations"] == 0


def test_community_analysis(component, compound):
    """Test community analysis."""
    # Perform analysis
    result = component.render_analysis(
        compound,
        analysis_types=["community"],
    )

    # Check result
    assert result.success
    analysis = result.data["analyses"]["community_analysis"]

    # Check statistics
    stats = analysis["statistics"]
    assert stats["reddit_posts"] == 2
    assert stats["twitter_mentions"] == 2
    assert stats["total_mentions"] == 4

    # Check date stats
    date_stats = analysis["date_stats"]
    assert "first_mention" in date_stats
    assert "last_mention" in date_stats
    assert date_stats["mention_days"] == 1


def test_analysis_history(component, compound):
    """Test analysis history tracking."""
    # Perform analyses
    component.render_analysis(compound, analysis_types=["binding"])
    component.render_analysis(compound, analysis_types=["activity"])
    component.render_analysis(compound, analysis_types=["safety"])

    # Check history
    history = component.analysis_stats["analysis_history"]
    assert len(history) == 3
    assert history[0]["analysis_types"] == ["binding"]
    assert history[1]["analysis_types"] == ["activity"]
    assert history[2]["analysis_types"] == ["safety"]


def test_error_handling(component, compound):
    """Test error handling."""
    # Test invalid SMILES
    compound.smiles = "invalid"
    result = component.render_analysis(
        compound,
        analysis_types=["sar"],
    )
    assert not result.success
    assert "Invalid SMILES" in result.error

    # Test template error
    with patch("flask.render_template", side_effect=Exception("Template error")):
        result = component.render_analysis(compound)
        assert not result.success
        assert "Template error" in result.error
