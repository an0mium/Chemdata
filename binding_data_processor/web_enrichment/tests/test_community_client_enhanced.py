"""Tests for enhanced community data source client."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime

import responses

from ..community_client_enhanced import CommunityClientEnhanced
from ...models.compound import Compound
from ...pipeline.infrastructure.circuit_breaker import CircuitConfig


@pytest.fixture
def client(tmp_path):
    """Create test client."""
    return CommunityClientEnhanced(
        cache_dir=tmp_path / "cache",
        circuit_config=CircuitConfig(
            failure_threshold=2,
            failure_timeout=1,
            reset_timeout=1,
        ),
    )


@pytest.fixture
def compound():
    """Create test compound."""
    return Compound(
        name="Test Compound",
        smiles="CC1=CC=CC=C1",
        cas_number="100-00-0",
    )


@pytest.fixture
def mock_psychonaut_response():
    """Create mock PsychonautWiki response."""
    return {
        "data": {
            "substances": [
                {
                    "name": "Test Compound",
                    "commonNames": ["Test"],
                    "class": {
                        "chemical": "Phenethylamine",
                        "psychoactive": "Stimulant",
                    },
                    "effects": [
                        {
                            "name": "Stimulation",
                            "url": "https://example.com",
                            "experience": {
                                "positive": ["Energy"],
                                "neutral": [],
                                "negative": ["Anxiety"],
                            },
                        }
                    ],
                    "roas": [
                        {
                            "name": "Oral",
                            "dose": {
                                "units": "mg",
                                "threshold": 10,
                                "light": {"min": 20, "max": 40},
                                "common": {"min": 40, "max": 80},
                                "strong": {"min": 80, "max": 120},
                                "heavy": 120,
                            },
                        }
                    ],
                }
            ]
        }
    }


@pytest.fixture
def mock_tripsit_response():
    """Create mock TripSit response."""
    return {
        "data": [
            {
                "name": "Test Compound",
                "properties": {
                    "dose": "20-80mg",
                    "duration": "4-6 hours",
                },
                "effects": {
                    "positive": ["Euphoria"],
                    "neutral": ["Stimulation"],
                    "negative": ["Anxiety"],
                },
            }
        ]
    }


@pytest.fixture
def mock_erowid_html():
    """Create mock Erowid HTML response."""
    return """
    <div class="experience-report">
        <h3>Test Report</h3>
        <div class="substance">Test Compound</div>
        <div class="author">Anonymous</div>
        <div class="date">2023-01-01</div>
        <div class="body">Test experience report text.</div>
    </div>
    """


@responses.activate
def test_get_psychonaut_data(client, compound, mock_psychonaut_response):
    """Test getting PsychonautWiki data."""
    # Mock request
    responses.add(
        responses.POST,
        client.PSYCHONAUT_API,
        json=mock_psychonaut_response,
        status=200,
    )

    # Get data
    data = client._get_psychonaut_data(compound.name)

    # Check data
    assert data is not None
    assert data["name"] == "Test Compound"
    assert data["class"]["psychoactive"] == "Stimulant"
    assert len(data["effects"]) == 1
    assert len(data["roas"]) == 1


@responses.activate
def test_get_tripsit_data(client, compound, mock_tripsit_response):
    """Test getting TripSit data."""
    # Mock request
    responses.add(
        responses.GET,
        f"{client.TRIPSIT_API}",
        json=mock_tripsit_response,
        status=200,
    )

    # Get data
    data = client._get_tripsit_data(compound.name)

    # Check data
    assert data is not None
    assert data["name"] == "Test Compound"
    assert "euphoria" in data["effects"]["positive"][0].lower()
    assert "anxiety" in data["effects"]["negative"][0].lower()


@responses.activate
def test_get_erowid_data(client, compound, mock_erowid_html):
    """Test getting Erowid data."""
    # Mock request
    responses.add(
        responses.GET,
        f"{client.EROWID_BASE}/search.php",
        body=mock_erowid_html,
        status=200,
    )

    # Get data
    data = client._get_erowid_data(compound.name)

    # Check data
    assert data is not None
    assert "reports" in data
    assert len(data["reports"]) == 1
    assert data["reports"][0]["title"] == "Test Report"
    assert data["reports"][0]["substance"] == "Test Compound"


def test_get_cached_psychonaut_data(client, compound, tmp_path):
    """Test getting cached PsychonautWiki data."""
    # Create cache file
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_file = cache_dir / f"psychonaut_{compound.name}.json"
    
    data = {
        "name": "Test Compound",
        "effects": ["Test Effect"],
        "roas": ["Test ROA"],
    }
    with cache_file.open("w") as f:
        json.dump(data, f)

    # Get cached data
    cached = client._get_cached_psychonaut_data(compound.name)

    # Check cached data
    assert cached is not None
    assert cached["name"] == "Test Compound"
    assert len(cached["effects"]) == 1
    assert len(cached["roas"]) == 1


def test_get_cached_tripsit_data(client, compound, tmp_path):
    """Test getting cached TripSit data."""
    # Create cache file
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_file = cache_dir / f"tripsit_{compound.name}.json"
    
    data = {
        "name": "Test Compound",
        "properties": {"test": "value"},
        "effects": ["Test Effect"],
    }
    with cache_file.open("w") as f:
        json.dump(data, f)

    # Get cached data
    cached = client._get_cached_tripsit_data(compound.name)

    # Check cached data
    assert cached is not None
    assert cached["name"] == "Test Compound"
    assert "test" in cached["properties"]
    assert len(cached["effects"]) == 1


def test_get_cached_erowid_data(client, compound, tmp_path):
    """Test getting cached Erowid data."""
    # Create cache file
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_file = cache_dir / f"erowid_{compound.name}.json"
    
    data = {
        "reports": [
            {
                "title": "Test Report",
                "substance": "Test Compound",
                "body_text": "Test text",
            }
        ]
    }
    with cache_file.open("w") as f:
        json.dump(data, f)

    # Get cached data
    cached = client._get_cached_erowid_data(compound.name)

    # Check cached data
    assert cached is not None
    assert len(cached["reports"]) == 1
    assert cached["reports"][0]["title"] == "Test Report"


@responses.activate
def test_process_compounds(
    client,
    compound,
    mock_psychonaut_response,
    mock_tripsit_response,
    mock_erowid_html,
):
    """Test processing compounds."""
    # Mock requests
    responses.add(
        responses.POST,
        client.PSYCHONAUT_API,
        json=mock_psychonaut_response,
        status=200,
    )
    responses.add(
        responses.GET,
        f"{client.TRIPSIT_API}",
        json=mock_tripsit_response,
        status=200,
    )
    responses.add(
        responses.GET,
        f"{client.EROWID_BASE}/search.php",
        body=mock_erowid_html,
        status=200,
    )

    # Process compound
    client.process_compounds([compound])

    # Check results
    assert compound.community_data is not None
    assert "psychonaut" in compound.community_data
    assert "tripsit" in compound.community_data
    assert "erowid" in compound.community_data
    assert len(client.processed_compounds) == 1
    assert len(client.failed_compounds) == 0


def test_metrics(client, compound):
    """Test metrics collection."""
    # Process successful compound
    client.processed_compounds.append(compound.name)
    client.source_stats["psychonaut"]["success"] += 1
    client.source_stats["tripsit"]["success"] += 1
    client.source_stats["erowid"]["failure"] += 1

    # Process failed compound
    client.failed_compounds.append("Failed Compound")

    # Get metrics
    metrics = client.get_metrics()

    # Check metrics
    assert metrics["processed_compounds"] == 1
    assert metrics["failed_compounds"] == 1
    assert metrics["success_rate"] == 0.5
    assert metrics["source_stats"]["psychonaut"]["success"] == 1
    assert metrics["source_stats"]["tripsit"]["success"] == 1
    assert metrics["source_stats"]["erowid"]["failure"] == 1


@responses.activate
def test_error_handling(client, compound):
    """Test error handling."""
    # Mock failing endpoints
    responses.add(
        responses.POST,
        client.PSYCHONAUT_API,
        status=500,
    )
    responses.add(
        responses.GET,
        f"{client.TRIPSIT_API}",
        status=500,
    )
    responses.add(
        responses.GET,
        f"{client.EROWID_BASE}/search.php",
        status=500,
    )

    # Process compound
    client.process_compounds([compound])

    # Check error handling
    assert len(client.processed_compounds) == 0
    assert len(client.failed_compounds) == 1
    assert client.failed_compounds[0] == compound.name
    assert client.source_stats["psychonaut"]["failure"] == 1
    assert client.source_stats["tripsit"]["failure"] == 1
    assert client.source_stats["erowid"]["failure"] == 1
