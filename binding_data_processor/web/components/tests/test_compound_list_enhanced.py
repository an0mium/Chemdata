"""Tests for enhanced compound list view component."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import plotly.graph_objects as go

from ...components.compound_list_enhanced import CompoundListEnhanced
from ....models.compound import Compound


@pytest.fixture
def component(tmp_path):
    """Create test component."""
    return CompoundListEnhanced(
        template_dir=tmp_path / "templates",
        static_dir=tmp_path / "static",
    )


@pytest.fixture
def compounds():
    """Create test compounds."""
    compounds = []

    # Compound with binding data
    compound1 = Compound(
        name="Test Compound 1",
        smiles="CC1=CC=CC=C1",
        cas_number="100-00-0",
    )
    compound1.binding_data = [
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
    compounds.append(compound1)

    # Compound with social data
    compound2 = Compound(
        name="Test Compound 2",
        smiles="CC2=CC=CC=C2",
        cas_number="200-00-0",
    )
    compound2.social_data = {
        "reddit": {
            "posts": [{"id": "1", "title": "Test Post"}],
            "comments": [{"id": "1", "text": "Test Comment"}],
        },
        "twitter": {
            "tweets": [{"id": "1", "text": "Test Tweet"}],
        },
    }
    compounds.append(compound2)

    # Compound with both
    compound3 = Compound(
        name="Test Compound 3",
        smiles="CC3=CC=CC=C3",
        cas_number="300-00-0",
    )
    compound3.binding_data = [
        {
            "target": "5-HT2A",
            "affinity": "0.8",
            "confidence": 0.95,
        },
    ]
    compound3.social_data = {
        "reddit": {
            "posts": [{"id": "2", "title": "Another Post"}],
        },
    }
    compounds.append(compound3)

    return compounds


@patch("flask.render_template")
def test_render_list(mock_render, component, compounds):
    """Test rendering compound list."""
    # Mock template rendering
    mock_render.return_value = "<html>Test</html>"

    # Render list
    result = component.render_list(compounds)

    # Check result
    assert result.success
    assert result.data["html"] == "<html>Test</html>"
    assert len(result.data["compounds"]) == 3
    assert result.data["total_pages"] == 1
    assert result.data["current_page"] == 1
    assert "visualizations" in result.data

    # Check stats
    assert component.view_stats["total_views"] == 1
    assert component.view_stats["visualization_views"] == 1


def test_handle_filter(component, compounds):
    """Test filtering compounds."""
    # Filter by name
    result = component.handle_filter(
        compounds,
        filters={"name": "test compound 1"},
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 1
    assert result.data["compounds"][0].name == "Test Compound 1"
    assert result.data["filter_stats"]["total"] == 3
    assert result.data["filter_stats"]["filtered"] == 1

    # Filter by binding
    result = component.handle_filter(
        compounds,
        filters={"min_binding": 1.0},
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 2
    assert all(hasattr(c, "binding_data") for c in result.data["compounds"])

    # Filter by social data
    result = component.handle_filter(
        compounds,
        filters={"has_social_data": True},
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 2
    assert all(hasattr(c, "social_data") for c in result.data["compounds"])


def test_handle_sort(component, compounds):
    """Test sorting compounds."""
    # Sort by name ascending
    result = component.handle_sort(
        compounds,
        sort_by="name",
        ascending=True,
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 3
    assert [c.name for c in result.data["compounds"]] == [
        "Test Compound 1",
        "Test Compound 2",
        "Test Compound 3",
    ]

    # Sort by binding count descending
    result = component.handle_sort(
        compounds,
        sort_by="binding_count",
        ascending=False,
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 3
    assert result.data["compounds"][0].name == "Test Compound 1"  # 2 bindings

    # Sort by social count descending
    result = component.handle_sort(
        compounds,
        sort_by="social_count",
        ascending=False,
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 3
    assert result.data["compounds"][0].name == "Test Compound 2"  # 2 posts


def test_handle_export(component, compounds):
    """Test exporting compounds."""
    # Export as TSV
    result = component.handle_export(
        compounds,
        format="tsv",
        columns=["name", "cas_number"],
    )

    # Check result
    assert result.success
    assert result.data["format"] == "tsv"
    assert result.data["row_count"] == 3
    assert "name" in result.data["columns"]
    assert "cas_number" in result.data["columns"]

    # Export as JSON
    result = component.handle_export(
        compounds,
        format="json",
    )

    # Check result
    assert result.success
    assert result.data["format"] == "json"
    data = json.loads(result.data["export_data"])
    assert len(data) == 3
    assert all("binding_data" in c for c in data if "Test Compound 1" in c["name"])
    assert all("social_data" in c for c in data if "Test Compound 2" in c["name"])


def test_generate_visualizations(component, compounds):
    """Test visualization generation."""
    # Generate visualizations
    visualizations = component._generate_visualizations(compounds)

    # Check binding distribution
    assert "binding_distribution" in visualizations
    binding_plot = json.loads(visualizations["binding_distribution"])
    assert binding_plot["data"][0]["type"] == "box"
    assert "5-HT2A" in binding_plot["layout"]["title"]["text"]

    # Check social summary
    assert "social_summary" in visualizations
    social_plot = json.loads(visualizations["social_summary"])
    assert social_plot["data"][0]["type"] == "bar"
    assert "Social Media Mentions" in social_plot["layout"]["title"]["text"]


def test_metrics(component, compounds):
    """Test metrics collection."""
    # Generate some activity
    component.render_list(compounds)
    component.handle_filter(compounds, {"name": "test"})
    component.handle_export(compounds, format="tsv")

    # Get metrics
    metrics = component.get_metrics()

    # Check metrics
    assert metrics["view_stats"]["total_views"] == 1
    assert metrics["view_stats"]["filtered_views"] == 1
    assert metrics["view_stats"]["exported_views"] == 1
    assert metrics["view_stats"]["visualization_views"] == 1


def test_error_handling(component, compounds):
    """Test error handling."""
    # Test invalid sort field
    result = component.handle_sort(
        compounds,
        sort_by="invalid_field",
    )
    assert not result.success
    assert "Invalid sort field" in result.error

    # Test invalid export format
    result = component.handle_export(
        compounds,
        format="invalid",
    )
    assert not result.success
    assert "Unsupported format" in result.error

    # Test template error
    with patch("flask.render_template", side_effect=Exception("Template error")):
        result = component.render_list(compounds)
        assert not result.success
        assert "Template error" in result.error
