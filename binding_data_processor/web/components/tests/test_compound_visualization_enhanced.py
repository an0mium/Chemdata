"""Tests for enhanced compound visualization component."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem

from ...components.compound_visualization_enhanced import CompoundVisualizationEnhanced
from ....models.compound import Compound


@pytest.fixture
def component(tmp_path):
    """Create test component."""
    return CompoundVisualizationEnhanced(
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
def test_render_visualization(mock_render, component, compound):
    """Test rendering visualization interface."""
    # Mock template rendering
    mock_render.return_value = "<html>Test</html>"

    # Render visualization
    result = component.render_visualization(compound)

    # Check result
    assert result.success
    assert result.data["html"] == "<html>Test</html>"
    assert "visualizations" in result.data
    assert "visualization_stats" in result.data

    # Check visualizations
    visualizations = result.data["visualizations"]
    assert "structure" in visualizations
    assert "binding_plot" in visualizations
    assert "activity_plot" in visualizations
    assert "safety_plot" in visualizations
    assert "social_plot" in visualizations

    # Check stats
    assert component.visualization_stats["total_visualizations"] == 1
    assert component.visualization_stats["structure_views"] == 1
    assert component.visualization_stats["binding_plots"] == 1
    assert component.visualization_stats["activity_plots"] == 1
    assert component.visualization_stats["safety_plots"] == 1
    assert component.visualization_stats["social_plots"] == 1


def test_structure_rendering(component, compound):
    """Test structure visualization."""
    # Render structure
    result = component.render_visualization(
        compound,
        plot_types=["structure"],
    )

    # Check result
    assert result.success
    assert "structure" in result.data["visualizations"]
    svg = result.data["visualizations"]["structure"]
    assert isinstance(svg, str)
    assert svg.startswith("<?xml")
    assert "svg" in svg

    # Test invalid SMILES
    compound.smiles = "invalid"
    result = component.render_visualization(compound)
    assert not result.success
    assert "Invalid SMILES" in result.error


def test_binding_plot(component, compound):
    """Test binding profile plot."""
    # Create plot
    result = component.render_visualization(
        compound,
        plot_types=["binding"],
    )

    # Check result
    assert result.success
    assert "binding_plot" in result.data["visualizations"]
    plot = result.data["visualizations"]["binding_plot"]

    # Check interactive plot
    assert isinstance(plot, dict)
    assert "data" in plot
    assert "layout" in plot
    assert plot["data"][0]["type"] == "bar"
    assert len(plot["data"][0]["x"]) == 2  # Two targets
    assert plot["layout"]["title"]["text"] == "Receptor Binding Profile"

    # Check static plot
    result = component.render_visualization(
        compound,
        plot_types=["binding"],
        interactive=False,
    )
    plot = result.data["visualizations"]["binding_plot"]
    assert isinstance(plot, str)
    assert plot.startswith("<?xml")


def test_activity_plot(component, compound):
    """Test activity profile plot."""
    # Create plot
    result = component.render_visualization(
        compound,
        plot_types=["activity"],
    )

    # Check result
    assert result.success
    assert "activity_plot" in result.data["visualizations"]
    plot = result.data["visualizations"]["activity_plot"]

    # Check interactive plot
    assert isinstance(plot, dict)
    assert "data" in plot
    assert "layout" in plot
    assert plot["data"][0]["type"] == "scatterpolar"
    assert len(plot["data"][0]["theta"]) == 3  # Three activities
    assert plot["layout"]["title"]["text"] == "Activity Profile"

    # Check static plot
    result = component.render_visualization(
        compound,
        plot_types=["activity"],
        interactive=False,
    )
    plot = result.data["visualizations"]["activity_plot"]
    assert isinstance(plot, str)
    assert plot.startswith("<?xml")


def test_safety_plot(component, compound):
    """Test safety profile plot."""
    # Create plot
    result = component.render_visualization(
        compound,
        plot_types=["safety"],
    )

    # Check result
    assert result.success
    assert "safety_plot" in result.data["visualizations"]
    plot = result.data["visualizations"]["safety_plot"]

    # Check interactive plot
    assert isinstance(plot, dict)
    assert "data" in plot
    assert "layout" in plot
    assert plot["data"][0]["type"] == "heatmap"
    assert len(plot["data"][0]["x"]) == 3  # Three risk categories
    assert plot["layout"]["title"]["text"] == "Safety Profile"

    # Check static plot
    result = component.render_visualization(
        compound,
        plot_types=["safety"],
        interactive=False,
    )
    plot = result.data["visualizations"]["safety_plot"]
    assert isinstance(plot, str)
    assert plot.startswith("<?xml")


def test_social_plot(component, compound):
    """Test social data visualization."""
    # Create plot
    result = component.render_visualization(
        compound,
        plot_types=["social"],
    )

    # Check result
    assert result.success
    assert "social_plot" in result.data["visualizations"]
    plot = result.data["visualizations"]["social_plot"]

    # Check interactive plot
    assert isinstance(plot, dict)
    assert "data" in plot
    assert "layout" in plot
    assert plot["data"][0]["type"] == "scatter"
    assert len(plot["data"][0]["x"]) == 4  # Four social mentions
    assert plot["layout"]["title"]["text"] == "Social Media Timeline"

    # Check static plot
    result = component.render_visualization(
        compound,
        plot_types=["social"],
        interactive=False,
    )
    plot = result.data["visualizations"]["social_plot"]
    assert isinstance(plot, str)
    assert plot.startswith("<?xml")


def test_visualization_history(component, compound):
    """Test visualization history tracking."""
    # Perform visualizations
    component.render_visualization(compound, plot_types=["structure"])
    component.render_visualization(compound, plot_types=["binding"])
    component.render_visualization(compound, plot_types=["activity"])

    # Check history
    history = component.visualization_stats["visualization_history"]
    assert len(history) == 3
    assert history[0]["plot_types"] == ["structure"]
    assert history[1]["plot_types"] == ["binding"]
    assert history[2]["plot_types"] == ["activity"]


def test_error_handling(component, compound):
    """Test error handling."""
    # Test invalid plot type
    result = component.render_visualization(
        compound,
        plot_types=["invalid"],
    )
    assert result.success  # Should still succeed but skip invalid plot
    assert "invalid_plot" not in result.data["visualizations"]

    # Test template error
    with patch("flask.render_template", side_effect=Exception("Template error")):
        result = component.render_visualization(compound)
        assert not result.success
        assert "Template error" in result.error
