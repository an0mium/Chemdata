"""Tests for enhanced Swiss tools client."""

import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime

import pandas as pd
import responses

from ..swiss_client_enhanced import SwissClientEnhanced
from ...models.compound import Compound
from ...pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig


@pytest.fixture
def client(tmp_path):
    """Create test client."""
    return SwissClientEnhanced(
        cache_dir=tmp_path / "cache",
        circuit_config=CircuitBreakerConfig(
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
def mock_target_response():
    """Create mock target prediction response."""
    return {
        "results": [
            {
                "target": "5-HT2A",
                "uniprot": "P28223",
                "gene": "HTR2A",
                "probability": 0.95,
                "class": "GPCR",
                "known_actives": 100,
            }
        ]
    }


@pytest.fixture
def mock_adme_response():
    """Create mock ADME response."""
    return {
        "results": {
            "MW": "100.0",
            "LogP": "2.0",
            "HBD": "1",
            "HBA": "2",
            "TPSA": "20.0",
            "RotBonds": "2",
            "Lipinski": "Yes",
            "Ghose": "Yes",
            "Veber": "Yes",
            "Egan": "Yes",
            "Muegge": "Yes",
            "GI_absorption": "High",
            "BBB_permeant": "Yes",
            "Pgp_substrate": "No",
            "CYP1A2_inhibition": "No",
            "CYP2C19_inhibition": "No",
            "CYP2C9_inhibition": "No",
            "CYP2D6_inhibition": "No",
            "CYP3A4_inhibition": "No",
            "PAINS": "0",
            "Brenk": "0",
            "Leadlikeness": "Yes",
            "SA": "2.5",
        }
    }


@responses.activate
def test_get_target_predictions(client, compound, mock_target_response):
    """Test getting target predictions."""
    # Mock initial request
    responses.add(
        responses.POST,
        client.STP_URL,
        json={"job_id": "test_job"},
        status=200,
    )

    # Mock status check
    responses.add(
        responses.GET,
        f"{client.STP_URL}/status/test_job",
        json={"status": "completed", **mock_target_response},
        status=200,
    )

    # Get predictions
    targets = client._get_target_predictions(compound.smiles)

    # Check predictions
    assert targets is not None
    assert len(targets) == 1
    assert targets[0]["target"] == "5-HT2A"
    assert targets[0]["probability"] == 0.95


@responses.activate
def test_get_adme_properties(client, compound, mock_adme_response):
    """Test getting ADME properties."""
    # Mock initial request
    responses.add(
        responses.POST,
        client.ADME_URL,
        json={"job_id": "test_job"},
        status=200,
    )

    # Mock status check
    responses.add(
        responses.GET,
        f"{client.ADME_URL}/status/test_job",
        json={"status": "completed", **mock_adme_response},
        status=200,
    )

    # Get properties
    adme = client._get_adme_properties(compound.smiles)

    # Check properties
    assert adme is not None
    assert adme["molecular_weight"] == 100.0
    assert adme["logp"] == 2.0
    assert adme["bbb_permeant"] == "Yes"


@responses.activate
def test_process_compounds(client, compound, mock_target_response, mock_adme_response):
    """Test processing compounds."""
    # Mock target prediction
    responses.add(
        responses.POST,
        client.STP_URL,
        json={"job_id": "test_job"},
        status=200,
    )
    responses.add(
        responses.GET,
        f"{client.STP_URL}/status/test_job",
        json={"status": "completed", **mock_target_response},
        status=200,
    )

    # Mock ADME properties
    responses.add(
        responses.POST,
        client.ADME_URL,
        json={"job_id": "test_job"},
        status=200,
    )
    responses.add(
        responses.GET,
        f"{client.ADME_URL}/status/test_job",
        json={"status": "completed", **mock_adme_response},
        status=200,
    )

    # Process compound
    client.process_compounds([compound])

    # Check results
    assert compound.swiss_data is not None
    assert "targets" in compound.swiss_data
    assert "adme" in compound.swiss_data
    assert len(client.processed_compounds) == 1
    assert len(client.failed_compounds) == 0


def test_get_cached_targets(client, compound, tmp_path):
    """Test getting cached target predictions."""
    # Create cache file
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_file = cache_dir / f"targets_{compound.smiles}.json"

    targets = [
        {
            "target": "5-HT2A",
            "uniprot": "P28223",
            "gene": "HTR2A",
            "probability": 0.95,
            "class": "GPCR",
            "known_actives": 100,
        }
    ]
    pd.DataFrame(targets).to_json(cache_file)

    # Get cached targets
    cached = client._get_cached_targets(compound.smiles)

    # Check cached data
    assert cached is not None
    assert len(cached) == 1
    assert cached[0]["target"] == "5-HT2A"


def test_get_cached_adme(client, compound, tmp_path):
    """Test getting cached ADME properties."""
    # Create cache file
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_file = cache_dir / f"adme_{compound.smiles}.json"

    adme = {
        "molecular_weight": 100.0,
        "logp": 2.0,
        "bbb_permeant": "Yes",
    }
    pd.DataFrame([adme]).to_json(cache_file)

    # Get cached properties
    cached = client._get_cached_adme(compound.smiles)

    # Check cached data
    assert cached is not None
    assert cached["molecular_weight"] == 100.0
    assert cached["logp"] == 2.0
    assert cached["bbb_permeant"] == "Yes"


def test_metrics(client, compound):
    """Test metrics collection."""
    # Process successful compound
    client.processed_compounds.append(compound.name)

    # Process failed compound
    client.failed_compounds.append("Failed Compound")

    # Get metrics
    metrics = client.get_metrics()

    # Check metrics
    assert metrics["processed_compounds"] == 1
    assert metrics["failed_compounds"] == 1
    assert metrics["success_rate"] == 0.5


@responses.activate
def test_error_handling(client, compound):
    """Test error handling."""
    # Mock failing endpoint
    responses.add(
        responses.POST,
        client.STP_URL,
        status=500,
    )

    # Process compound
    client.process_compounds([compound])

    # Check error handling
    assert len(client.processed_compounds) == 0
    assert len(client.failed_compounds) == 1
    assert client.failed_compounds[0] == compound.name
