"""Tests for enhanced compound search component."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem

from ...components.compound_search_enhanced import CompoundSearchEnhanced
from ....models.compound import Compound


@pytest.fixture
def component(tmp_path):
    """Create test component."""
    return CompoundSearchEnhanced(
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
    now = datetime.now()
    compound2 = Compound(
        name="Test Compound 2",
        smiles="CC2=CC=CC=C2",
        cas_number="200-00-0",
    )
    compound2.social_data = {
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
            "posts": [
                {
                    "id": "3",
                    "title": "Another Post",
                    "created_utc": (now - timedelta(days=3)).isoformat(),
                },
            ],
        },
    }
    compounds.append(compound3)

    return compounds


@patch("flask.render_template")
def test_render_search(mock_render, component, compounds):
    """Test rendering search interface."""
    # Mock template rendering
    mock_render.return_value = "<html>Test</html>"

    # Render search
    result = component.render_search(compounds)

    # Check result
    assert result.success
    assert result.data["html"] == "<html>Test</html>"
    assert len(result.data["compounds"]) == 3
    assert "search_stats" in result.data

    # Check stats
    assert component.search_stats["total_searches"] == 1


def test_text_search(component, compounds):
    """Test text search functionality."""
    # Search by name
    result = component.render_search(
        compounds,
        query="test compound 1",
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 1
    assert result.data["compounds"][0].name == "Test Compound 1"
    assert component.search_stats["text_searches"] == 1

    # Search by target
    result = component.render_search(
        compounds,
        query="5-ht2a",
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 2
    assert all(
        any(b["target"] == "5-HT2A" for b in c.binding_data)
        for c in result.data["compounds"]
    )

    # Search social data
    result = component.render_search(
        compounds,
        query="test post",
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 1
    assert result.data["compounds"][0].name == "Test Compound 2"


def test_structure_search(component, compounds):
    """Test structure search functionality."""
    # Search by exact structure
    result = component.render_search(
        compounds,
        structure="CC1=CC=CC=C1",
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 1
    assert result.data["compounds"][0].name == "Test Compound 1"
    assert component.search_stats["structure_searches"] == 1

    # Check similarity sorting
    result = component.render_search(
        compounds,
        structure="CC1=CC=CC=C1",
        filters={"min_similarity": 0.5},
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 3
    similarities = [c.similarity for c in result.data["compounds"]]
    assert similarities == sorted(similarities, reverse=True)


def test_filter_search(component, compounds):
    """Test filter functionality."""
    # Filter by binding affinity
    result = component.render_search(
        compounds,
        filters={"min_binding": 1.0},
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 2
    assert all(
        any(float(b["affinity"]) >= 1.0 for b in c.binding_data)
        for c in result.data["compounds"]
    )
    assert component.search_stats["filter_searches"] == 1

    # Filter by target
    result = component.render_search(
        compounds,
        filters={"target": "d2"},
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 1
    assert result.data["compounds"][0].name == "Test Compound 1"

    # Filter by social data
    result = component.render_search(
        compounds,
        filters={"has_social_data": True},
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 2
    assert all(
        hasattr(c, "social_data")
        for c in result.data["compounds"]
    )

    # Filter by date range
    now = datetime.now()
    result = component.render_search(
        compounds,
        filters={
            "date_from": (now - timedelta(days=2)).isoformat(),
            "date_to": (now - timedelta(days=1)).isoformat(),
        },
    )

    # Check result
    assert result.success
    assert len(result.data["compounds"]) == 1
    assert result.data["compounds"][0].name == "Test Compound 2"


def test_search_history(component, compounds):
    """Test search history tracking."""
    # Perform searches
    component.render_search(compounds, query="test")
    component.render_search(compounds, structure="CC1=CC=CC=C1")
    component.render_search(compounds, filters={"target": "5-HT2A"})

    # Check history
    history = component.search_stats["search_history"]
    assert len(history) == 3
    assert history[0]["query"] == "test"
    assert history[1]["structure"] == "CC1=CC=CC=C1"
    assert history[2]["filters"]["target"] == "5-HT2A"


def test_error_handling(component, compounds):
    """Test error handling."""
    # Test invalid SMILES
    result = component.render_search(
        compounds,
        structure="invalid",
    )
    assert not result.success
    assert "Invalid SMILES" in result.error

    # Test template error
    with patch("flask.render_template", side_effect=Exception("Template error")):
        result = component.render_search(compounds)
        assert not result.success
        assert "Template error" in result.error
