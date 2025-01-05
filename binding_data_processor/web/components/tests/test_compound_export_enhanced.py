"""Tests for enhanced compound export component."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

import pandas as pd
from rdkit import Chem

from ...components.compound_export_enhanced import CompoundExportEnhanced
from ....models.compound import Compound


@pytest.fixture
def component(tmp_path):
    """Create test component."""
    return CompoundExportEnhanced(
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

    # Compound with predictions
    compound3 = Compound(
        name="Test Compound 3",
        smiles="CC3=CC=CC=C3",
        cas_number="300-00-0",
    )
    compound3.predictions = {
        "binding": {
            "5-HT2A": 0.8,
            "D2": 0.6,
        },
        "activity": {
            "stimulant": 0.7,
            "psychedelic": 0.3,
        },
        "safety": {
            "toxicity": 0.2,
            "addiction": 0.4,
        },
    }
    compounds.append(compound3)

    return compounds


@patch("flask.render_template")
def test_render_export(mock_render, component, compounds):
    """Test rendering export interface."""
    # Mock template rendering
    mock_render.return_value = "<html>Test</html>"

    # Render export
    result = component.render_export(compounds)

    # Check result
    assert result.success
    assert result.data["html"] == "<html>Test</html>"
    assert result.data["format"] == "tsv"
    assert "export_stats" in result.data

    # Check stats
    assert component.export_stats["total_exports"] == 1
    assert component.export_stats["tsv_exports"] == 1


def test_tsv_export(component, compounds):
    """Test TSV export functionality."""
    # Export all columns
    result = component.render_export(
        compounds,
        format="tsv",
    )

    # Check result
    assert result.success
    assert isinstance(result.data["export_data"], str)

    # Parse TSV
    df = pd.read_csv(pd.StringIO(result.data["export_data"]), sep="\t")
    assert len(df) == 3
    assert "name" in df.columns
    assert "smiles" in df.columns
    assert "cas_number" in df.columns

    # Export selected columns
    result = component.render_export(
        compounds,
        format="tsv",
        columns=["name", "cas_number"],
    )

    # Check result
    assert result.success
    df = pd.read_csv(pd.StringIO(result.data["export_data"]), sep="\t")
    assert list(df.columns) == ["name", "cas_number"]


def test_sdf_export(component, compounds):
    """Test SDF export functionality."""
    # Export compounds
    result = component.render_export(
        compounds,
        format="sdf",
    )

    # Check result
    assert result.success
    assert isinstance(result.data["export_data"], str)

    # Parse SDF
    suppl = Chem.SDMolSupplier()
    suppl.SetData(result.data["export_data"])
    mols = [mol for mol in suppl if mol]

    # Check molecules
    assert len(mols) == 3
    assert all(mol.HasProp("_Name") for mol in mols)
    assert all(mol.HasProp("CAS") for mol in mols)


def test_json_export(component, compounds):
    """Test JSON export functionality."""
    # Export all data
    result = component.render_export(
        compounds,
        format="json",
    )

    # Check result
    assert result.success
    data = json.loads(result.data["export_data"])
    assert len(data) == 3
    assert all(isinstance(d, dict) for d in data)

    # Export selected columns
    result = component.render_export(
        compounds,
        format="json",
        columns=["name", "cas_number"],
    )

    # Check result
    assert result.success
    data = json.loads(result.data["export_data"])
    assert all(set(d.keys()) == {"name", "cas_number"} for d in data)


def test_data_validation(component, compounds):
    """Test data validation."""
    # Test missing name
    compounds[0].name = ""
    result = component.render_export(compounds, validate=True)
    assert not result.success
    assert "validation_errors" in result.data
    assert any(e["field"] == "name" for e in result.data["validation_errors"])

    # Test invalid SMILES
    compounds[0].smiles = "invalid"
    result = component.render_export(compounds, validate=True)
    assert not result.success
    assert any(e["field"] == "smiles" for e in result.data["validation_errors"])

    # Test missing required column
    result = component.render_export(
        compounds,
        columns=["binding_data"],
        validate=True,
    )
    assert not result.success
    assert any(
        e["field"] == "binding_data" and "Missing binding data" in e["error"]
        for e in result.data["validation_errors"]
    )


def test_export_history(component, compounds):
    """Test export history tracking."""
    # Perform exports
    component.render_export(compounds, format="tsv")
    component.render_export(compounds, format="sdf")
    component.render_export(compounds, format="json")

    # Check history
    history = component.export_stats["export_history"]
    assert len(history) == 3
    assert history[0]["format"] == "tsv"
    assert history[1]["format"] == "sdf"
    assert history[2]["format"] == "json"


def test_error_handling(component, compounds):
    """Test error handling."""
    # Test invalid format
    result = component.render_export(
        compounds,
        format="invalid",
    )
    assert not result.success
    assert "Unsupported format" in result.error

    # Test template error
    with patch("flask.render_template", side_effect=Exception("Template error")):
        result = component.render_export(compounds)
        assert not result.success
        assert "Template error" in result.error
