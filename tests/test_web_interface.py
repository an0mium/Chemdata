"""Tests for web interface functionality."""

import json
import tempfile
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest
from flask import url_for

from binding_data_processor.config import Config
from web.app import ChemDataApp


@pytest.fixture
def test_config():
    """Create test configuration."""
    with tempfile.TemporaryDirectory() as temp_dir:
        config = Config()
        config.data_dir = Path(temp_dir) / "data"
        config.model_dir = Path(temp_dir) / "models"
        config.cache_dir = Path(temp_dir) / "cache"

        # Create directories
        for directory in [config.data_dir, config.model_dir, config.cache_dir]:
            directory.mkdir(parents=True)

        yield config


@pytest.fixture
def test_data():
    """Create test compound data."""
    compounds = pd.DataFrame(
        {
            "name": ["Test Compound 1", "Test Compound 2"],
            "smiles": [
                "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC",
                "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F",
            ],
            "molecular_weight": [395.47, 401.45],
            "logp": [3.2, 2.8],
            "target_type": ["serotonin", "dopamine"],
            "activity_type": ["Ki", "IC50"],
            "activity_value": [10.5, 25.3],
        }
    )
    return compounds


@pytest.fixture
def test_predictions():
    """Create test prediction data."""
    return {
        "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC": {
            "toxicity": {
                "hepatotoxicity": {"probability": 0.2, "confidence": 0.8},
            },
            "abuse": {
                "potential": {"probability": 0.3, "confidence": 0.7},
            },
        },
        "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F": {
            "toxicity": {
                "hepatotoxicity": {"probability": 0.4, "confidence": 0.6},
            },
            "abuse": {
                "potential": {"probability": 0.5, "confidence": 0.7},
            },
        },
    }


@pytest.fixture
def test_web_data():
    """Create test web data."""
    return {
        "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC": {
            "community": {
                "reports": [
                    {
                        "source": "erowid",
                        "text": "Test report 1",
                    }
                ],
            },
        },
        "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F": {
            "community": {
                "reports": [
                    {
                        "source": "psychonautwiki",
                        "text": "Test report 2",
                    }
                ],
            },
        },
    }


@pytest.fixture
def app(test_config, test_data, test_predictions, test_web_data):
    """Create test Flask application."""
    app = ChemDataApp(
        data_dir=str(test_config.data_dir),
        model_dir=str(test_config.model_dir),
        debug=True,
    )

    # Save test data
    test_data.to_csv(test_config.data_dir / "compounds.tsv", sep="\t", index=False)

    # Save test predictions
    with open(test_config.data_dir / "predictions.json", "w") as f:
        json.dump(test_predictions, f)

    # Save test web data
    with open(test_config.data_dir / "web_data.json", "w") as f:
        json.dump(test_web_data, f)

    return app.server


@pytest.fixture
def client(app):
    """Create test client."""
    return app.test_client()


def test_index_route(client):
    """Test index route."""
    response = client.get("/")
    assert response.status_code == 200
    assert b"ChemData" in response.data


def test_api_compounds(client, test_data):
    """Test compounds API endpoint."""
    response = client.get("/api/compounds")
    assert response.status_code == 200
    data = json.loads(response.data)
    assert len(data["data"]) == len(test_data)
    assert "draw" in data
    assert "recordsTotal" in data
    assert "recordsFiltered" in data


def test_api_compound_detail(client, test_data, test_predictions, test_web_data):
    """Test compound detail API endpoint."""
    smiles = test_data.iloc[0]["smiles"]
    response = client.get(f"/api/compounds/{smiles}")
    assert response.status_code == 200
    data = json.loads(response.data)
    assert data["name"] == test_data.iloc[0]["name"]
    assert "toxicity" in data["predictions"]
    assert "community" in data["web_data"]


def test_api_structure(client):
    """Test structure API endpoint."""
    smiles = "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
    response = client.get(f"/api/structure/{smiles}")
    assert response.status_code == 200
    assert response.mimetype == "image/svg+xml"


def test_api_filtering(client, test_data):
    """Test compound filtering."""
    # Test target type filter
    response = client.get("/api/compounds?target=serotonin")
    data = json.loads(response.data)
    assert len(data["data"]) == 1
    assert data["data"][0]["target_type"] == "serotonin"

    # Test activity type filter
    response = client.get("/api/compounds?activity=Ki")
    data = json.loads(response.data)
    assert len(data["data"]) == 1
    assert data["data"][0]["activity_type"] == "Ki"

    # Test text search
    response = client.get("/api/compounds?search=Test+Compound+1")
    data = json.loads(response.data)
    assert len(data["data"]) == 1
    assert data["data"][0]["name"] == "Test Compound 1"


def test_api_sorting(client):
    """Test compound sorting."""
    # Test name sorting
    response = client.get("/api/compounds?sort=name&order=asc")
    data = json.loads(response.data)
    names = [compound["name"] for compound in data["data"]]
    assert names == sorted(names)

    # Test molecular weight sorting
    response = client.get("/api/compounds?sort=molecular_weight&order=desc")
    data = json.loads(response.data)
    weights = [compound["molecular_weight"] for compound in data["data"]]
    assert weights == sorted(weights, reverse=True)


def test_api_pagination(client):
    """Test compound pagination."""
    # Test page size
    response = client.get("/api/compounds?length=1")
    data = json.loads(response.data)
    assert len(data["data"]) == 1

    # Test page number
    response = client.get("/api/compounds?start=1&length=1")
    data = json.loads(response.data)
    assert len(data["data"]) == 1
    assert data["data"][0]["name"] == "Test Compound 2"


def test_api_structure_similarity(client):
    """Test structure similarity search."""
    smiles = "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
    response = client.get(f"/api/compounds?structure={smiles}&similarity=0.7")
    data = json.loads(response.data)
    assert len(data["data"]) == 2  # Both test compounds are similar


def test_api_error_handling(client):
    """Test API error handling."""
    # Test invalid SMILES
    response = client.get("/api/structure/invalid")
    assert response.status_code == 400

    # Test missing compound
    response = client.get("/api/compounds/nonexistent")
    assert response.status_code == 404

    # Test invalid filter
    response = client.get("/api/compounds?invalid=filter")
    assert response.status_code == 400


def test_visualization_generation(client, test_data):
    """Test visualization generation."""
    # Test activity plot
    response = client.get("/api/plots/activity")
    assert response.status_code == 200
    assert response.mimetype == "image/png"

    # Test target distribution plot
    response = client.get("/api/plots/targets")
    assert response.status_code == 200
    assert response.mimetype == "image/png"

    # Test property distribution plot
    response = client.get("/api/plots/properties")
    assert response.status_code == 200
    assert response.mimetype == "image/png"


def test_data_export(client):
    """Test data export functionality."""
    # Test CSV export
    response = client.get("/api/export/csv")
    assert response.status_code == 200
    assert response.mimetype == "text/csv"

    # Test JSON export
    response = client.get("/api/export/json")
    assert response.status_code == 200
    assert response.mimetype == "application/json"

    # Test SDF export
    response = client.get("/api/export/sdf")
    assert response.status_code == 200
    assert response.mimetype == "chemical/x-mdl-sdfile"


def test_cache_control(client, test_data):
    """Test cache control headers."""
    # Static resources should be cached
    response = client.get("/static/css/base.css")
    assert response.status_code == 200
    assert "Cache-Control" in response.headers

    # Dynamic data should not be cached
    response = client.get("/api/compounds")
    assert response.status_code == 200
    assert "no-cache" in response.headers.get("Cache-Control", "")
