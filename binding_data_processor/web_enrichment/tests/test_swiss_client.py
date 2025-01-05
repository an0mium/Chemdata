"""Tests for Swiss web enrichment client."""

import logging
import pytest
from unittest.mock import Mock, patch

import requests

from ..clients.swiss import SwissClient, SwissValidationError
from ..validation.schema import ValidationLevel
from ...pipeline.infrastructure.circuit_breaker import CircuitConfig


@pytest.fixture
def client():
    """Create test client."""
    return SwissClient(
        base_url="http://test.swiss.com",
        requests_per_second=10,
        max_retries=1,
        timeout=1,
        cache_ttl=60,
        circuit_config=CircuitConfig(
            failure_threshold=3,
            recovery_timeout=60,
        ),
        validation_level=ValidationLevel.NORMAL,
        logger=logging.getLogger("test"),
    )


def test_client_initialization(client):
    """Test client initialization."""
    assert client.name == "swiss"
    assert client.base_url == "http://test.swiss.com"
    assert len(client.processed_compounds) == 0
    assert len(client.failed_compounds) == 0


def test_predict_targets_success(client):
    """Test successful target prediction."""
    mock_response = Mock()
    mock_response.json.return_value = {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "targets": [
            {
                "target": "COX-2",
                "probability": 0.95,
            },
        ],
    }

    with patch.object(client.http, "request", return_value=mock_response):
        data = client.predict_targets(
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
            threshold=0.5,
        )

    assert "targets" in data
    assert len(data["targets"]) == 1
    assert data["targets"][0]["target"] == "COX-2"
    assert "CC(=O)OC1=CC=CC=C1C(=O)O" in client.processed_compounds
    assert len(client.failed_compounds) == 0


def test_predict_targets_validation_error(client):
    """Test target prediction validation error."""
    mock_response = Mock()
    mock_response.json.return_value = {
        "smiles": "invalid",
        "error": "Invalid SMILES",
    }

    with patch.object(client.http, "request", return_value=mock_response):
        with pytest.raises(SwissValidationError) as exc:
            client.predict_targets(
                smiles="invalid",
                threshold=0.5,
            )

    assert exc.value.tool == "target"
    assert "invalid" in client.failed_compounds
    assert len(client.processed_compounds) == 0


def test_calculate_adme_success(client):
    """Test successful ADME calculation."""
    mock_response = Mock()
    mock_response.json.return_value = {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "properties": {
            "mw": 180.16,
            "logp": 1.43,
            "hbd": 1,
            "hba": 4,
            "psa": 63.6,
        },
    }

    with patch.object(client.http, "request", return_value=mock_response):
        data = client.calculate_adme(
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        )

    assert "properties" in data
    assert data["properties"]["mw"] == 180.16
    assert "CC(=O)OC1=CC=CC=C1C(=O)O" in client.processed_compounds
    assert len(client.failed_compounds) == 0


def test_calculate_adme_validation_error(client):
    """Test ADME calculation validation error."""
    mock_response = Mock()
    mock_response.json.return_value = {
        "smiles": "invalid",
        "error": "Invalid SMILES",
    }

    with patch.object(client.http, "request", return_value=mock_response):
        with pytest.raises(SwissValidationError) as exc:
            client.calculate_adme(
                smiles="invalid",
            )

    assert exc.value.tool == "adme"
    assert "invalid" in client.failed_compounds
    assert len(client.processed_compounds) == 0


def test_compare_structures_success(client):
    """Test successful structure comparison."""
    mock_response = Mock()
    mock_response.json.return_value = {
        "query": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "reference": "CC1=CC=CC=C1C(=O)O",
        "similarity": {
            "tanimoto": 0.85,
            "dice": 0.92,
        },
    }

    with patch.object(client.http, "request", return_value=mock_response):
        data = client.compare_structures(
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
            reference_smiles="CC1=CC=CC=C1C(=O)O",
        )

    assert "similarity" in data
    assert data["similarity"]["tanimoto"] == 0.85
    assert "CC(=O)OC1=CC=CC=C1C(=O)O" in client.processed_compounds
    assert len(client.failed_compounds) == 0


def test_compare_structures_validation_error(client):
    """Test structure comparison validation error."""
    mock_response = Mock()
    mock_response.json.return_value = {
        "query": "invalid",
        "reference": "CC1=CC=CC=C1C(=O)O",
        "error": "Invalid query SMILES",
    }

    with patch.object(client.http, "request", return_value=mock_response):
        with pytest.raises(SwissValidationError) as exc:
            client.compare_structures(
                smiles="invalid",
                reference_smiles="CC1=CC=CC=C1C(=O)O",
            )

    assert exc.value.tool == "similarity"
    assert "invalid" in client.failed_compounds
    assert len(client.processed_compounds) == 0


def test_http_error_handling(client):
    """Test HTTP error handling."""
    mock_error = requests.exceptions.RequestException()
    mock_error.response = Mock(status_code=500)

    with patch.object(client.http, "request", side_effect=mock_error):
        with pytest.raises(requests.exceptions.RequestException):
            client.predict_targets(
                smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
            )

    assert "CC(=O)OC1=CC=CC=C1C(=O)O" in client.failed_compounds
    assert len(client.processed_compounds) == 0


def test_metrics(client):
    """Test metrics collection."""
    # Successful request
    mock_response = Mock()
    mock_response.json.return_value = {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "targets": [],
    }

    with patch.object(client.http, "request", return_value=mock_response):
        client.predict_targets(
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        )

    # Failed request
    mock_error = requests.exceptions.RequestException()
    mock_error.response = Mock(status_code=500)

    with patch.object(client.http, "request", side_effect=mock_error):
        with pytest.raises(requests.exceptions.RequestException):
            client.predict_targets(
                smiles="invalid",
            )

    metrics = client.get_metrics()
    assert metrics["processed_compounds"] == 1
    assert metrics["failed_compounds"] == 1
    assert "http_client" in metrics
