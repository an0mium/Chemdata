"""Tests for enhanced compound detail view component."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

import pandas as pd
import plotly.graph_objects as go

from ...components.compound_detail_enhanced import CompoundDetailEnhanced
from ....models.compound import Compound


@pytest.fixture
def component(tmp_path):
    """Create test component."""
    return CompoundDetailEnhanced(
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

    # Add social data
    now = datetime.now()
    compound.social_data = {
        "reddit": {
            "posts": [
                {
                    "id": "1",
                    "title": "Test Post 1",
                    "created_utc": (now - timedelta(days=2)).timestamp(),
                },
                {
                    "id": "2",
                    "title": "Test Post 2",
                    "created_utc": (now - timedelta(days=1)).timestamp(),
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
def test_render_detail(mock_render, component, compound):
    """Test rendering compound detail."""
    # Mock template rendering
    mock_render.return_value = "<html>Test</html>"

    # Render detail
    result = component.render_detail(compound)

    # Check result
    assert result.success
    assert result.data["html"] == "<html>Test</html>"
    assert result.data["compound"] == compound
    assert "visualizations" in result.data

    # Check visualizations
    visualizations = result.data["visualizations"]
    assert "structure" in visualizations
    assert "binding" in visualizations
    assert "social" in visualizations

    # Check stats
    assert component.view_stats["total_views"] == 1
    assert component.view_stats["structure_views"] == 1
    assert component.view_stats["binding_views"] == 1
    assert component.view_stats["social_views"] == 1


@patch("rdkit.Chem.MolFromSmiles")
@patch("rdkit.Chem.Draw.MolToImage")
def test_generate_structure(mock_draw, mock_mol, component, compound):
    """Test structure visualization."""
    # Mock RDKit
    mock_mol.return_value = "mock_mol"
    mock_draw.return_value.tostring.return_value = "mock_svg"

    # Generate structure
    svg = component._generate_structure(compound)

    # Check result
    assert svg == "mock_svg"
    mock_mol.assert_called_once_with(compound.smiles)
    mock_draw.assert_called_once_with("mock_mol")


def test_generate_binding_plot(component, compound):
    """Test binding plot generation."""
    # Generate plot
    plot_json = component._generate_binding_plot(compound)

    # Parse plot
    plot = json.loads(plot_json)

    # Check plot
    assert plot["data"][0]["type"] == "bar"
    assert len(plot["data"][0]["x"]) == 2  # Two targets
    assert plot["data"][0]["name"] == "Binding Affinity"
    assert "error_y" in plot["data"][0]
    assert "Binding Affinity by Target" in plot["layout"]["title"]["text"]


def test_generate_social_plots(component, compound):
    """Test social plot generation."""
    # Generate plots
    plots = component._generate_social_plots(compound)

    # Check Reddit plot
    reddit_plot = json.loads(plots["reddit_activity"])
    assert reddit_plot["data"][0]["type"] == "scatter"
    assert len(reddit_plot["data"][0]["x"]) == 2  # Two posts
    assert reddit_plot["data"][0]["name"] == "Reddit Posts"
    assert "Reddit Activity" in reddit_plot["layout"]["title"]["text"]

    # Check Twitter plot
    twitter_plot = json.loads(plots["twitter_activity"])
    assert twitter_plot["data"][0]["type"] == "scatter"
    assert len(twitter_plot["data"][0]["x"]) == 2  # Two tweets
    assert twitter_plot["data"][0]["name"] == "Twitter Mentions"
    assert "Twitter Activity" in twitter_plot["layout"]["title"]["text"]


def test_handle_export_json(component, compound):
    """Test JSON export."""
    # Export as JSON
    result = component.handle_export(
        compound,
        format="json",
    )

    # Check result
    assert result.success
    assert result.data["format"] == "json"
    
    # Parse data
    data = json.loads(result.data["export_data"])
    assert data["name"] == compound.name
    assert data["smiles"] == compound.smiles
    assert data["cas_number"] == compound.cas_number
    assert "binding_data" in data
    assert "social_data" in data


@patch("rdkit.Chem.MolFromSmiles")
@patch("rdkit.Chem.SDWriter")
def test_handle_export_sdf(mock_writer, mock_mol, component, compound):
    """Test SDF export."""
    # Mock RDKit
    mock_mol.return_value = Mock()
    mock_writer.return_value = "mock_sdf"

    # Export as SDF
    result = component.handle_export(
        compound,
        format="sdf",
    )

    # Check result
    assert result.success
    assert result.data["format"] == "sdf"
    assert result.data["export_data"] == "mock_sdf"


def test_metrics(component, compound):
    """Test metrics collection."""
    # Generate some activity
    component.render_detail(compound)
    component.handle_export(compound, format="json")

    # Get metrics
    metrics = component.get_metrics()

    # Check metrics
    assert metrics["view_stats"]["total_views"] == 1
    assert metrics["view_stats"]["structure_views"] == 1
    assert metrics["view_stats"]["binding_views"] == 1
    assert metrics["view_stats"]["social_views"] == 1
    assert metrics["view_stats"]["exported_views"] == 1


def test_error_handling(component, compound):
    """Test error handling."""
    # Test invalid export format
    result = component.handle_export(
        compound,
        format="invalid",
    )
    assert not result.success
    assert "Unsupported format" in result.error

    # Test template error
    with patch("flask.render_template", side_effect=Exception("Template error")):
        result = component.render_detail(compound)
        assert not result.success
        assert "Template error" in result.error
