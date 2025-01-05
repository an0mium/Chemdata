"""Tests for base web enrichment client."""

import logging
import pytest
from unittest.mock import Mock, patch

import requests

from ..clients.base import WebClient, WebClientError, ValidationError
from ..validation.schema import DataSource, ValidationLevel
from ..validation.data import ValidationConfig
from ...pipeline.infrastructure.circuit_breaker import CircuitConfig


class TestClient(WebClient):
    """Test client implementation."""

    def _get_validation_config(self) -> ValidationConfig:
        """Get validation configuration."""
        return ValidationConfig(
            rules={
                "name": [("required", None)],
                "value": [("range", (0, 100))],
            }
        )


@pytest.fixture
def client():
    """Create test client."""
    return TestClient(
        name="test",
        base_url="http://test.com",
        data_source=DataSource.SWISS,
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
    assert client.name == "test"
    assert client.base_url == "http://test.com"
    assert client.data_source == DataSource.SWISS
    assert client.processed_items == 0
    assert client.failed_items == 0
    assert client.validation_errors == 0
    assert client.http_errors == 0


def test_request_success(client):
    """Test successful request."""
    mock_response = Mock()
    mock_response.json.return_value = {"name": "test", "value": 50}

    with patch.object(client.http, "request", return_value=mock_response):
        data = client.request(
            method="GET",
            endpoint="/test",
            params={"id": 1},
            schema_name="test",
        )

    assert data == {"name": "test", "value": 50}
    assert client.processed_items == 1
    assert client.failed_items == 0
    assert client.validation_errors == 0
    assert client.http_errors == 0


def test_request_http_error(client):
    """Test HTTP error handling."""
    mock_error = requests.exceptions.RequestException()
    mock_error.response = Mock(status_code=500)

    with patch.object(client.http, "request", side_effect=mock_error):
        with pytest.raises(WebClientError) as exc:
            client.request(
                method="GET",
                endpoint="/test",
            )

    assert exc.value.status_code == 500
    assert client.processed_items == 0
    assert client.failed_items == 1
    assert client.validation_errors == 0
    assert client.http_errors == 1


def test_request_rate_limit_error(client):
    """Test rate limit error handling."""
    # Mock response with rate limit headers
    mock_response = Mock()
    mock_response.status_code = 429
    mock_response.headers = {"Retry-After": "60"}
    mock_response.json.return_value = {"error": "Rate limit exceeded"}

    # Mock request to raise rate limit error
    def mock_request(*args, **kwargs):
        error = requests.exceptions.RequestException("Rate limit exceeded")
        error.response = mock_response
        raise error

    with patch.object(client.http, "request", side_effect=mock_request):
        with pytest.raises(WebClientError) as exc:
            client.request(
                method="GET",
                endpoint="/test",
            )

    # Verify error details
    assert exc.value.status_code == 429
    assert "Rate limit exceeded" in str(exc.value)
    assert client.processed_items == 0
    assert client.failed_items == 1
    assert client.validation_errors == 0
    assert client.http_errors == 1


def test_validation_error(client):
    """Test validation error handling."""
    mock_response = Mock()
    mock_response.json.return_value = {"value": 200}  # Invalid value

    with patch.object(client.http, "request", return_value=mock_response):
        with pytest.raises(ValidationError) as exc:
            client.request(
                method="GET",
                endpoint="/test",
                schema_name="test",
            )

    assert "Data validation failed" in str(exc.value)
    assert client.processed_items == 0
    assert client.failed_items == 0
    assert client.validation_errors == 1
    assert client.http_errors == 0


def test_metrics(client):
    """Test metrics collection."""
    # Successful request
    mock_response = Mock()
    mock_response.json.return_value = {"name": "test", "value": 50}

    with patch.object(client.http, "request", return_value=mock_response):
        client.request(
            method="GET",
            endpoint="/test",
            schema_name="test",
        )

    # Failed request
    mock_error = requests.exceptions.RequestException()
    mock_error.response = Mock(status_code=500)

    with patch.object(client.http, "request", side_effect=mock_error):
        with pytest.raises(WebClientError):
            client.request(
                method="GET",
                endpoint="/test",
            )

    metrics = client.get_metrics()
    assert metrics["processed_items"] == 1
    assert metrics["failed_items"] == 1
    assert metrics["validation_errors"] == 0
    assert metrics["http_errors"] == 1
    assert "http_client" in metrics


def test_data_cleaning(client):
    """Test data cleaning utilities."""
    # Clean text
    text = "<p>Test  \n  text</p>"
    assert client.clean_text(text) == "Test text"

    # Clean number
    assert client.clean_number("123.45") == 123.45
    with pytest.raises(ValueError):
        client.clean_number("invalid")

    # Clean timestamp
    assert client.clean_timestamp("2023-01-01T00:00:00Z").year == 2023
    with pytest.raises(ValueError):
        client.clean_timestamp("invalid")

    # Clean duration
    assert client.clean_duration("1h30m") == "1h30m"
    assert client.clean_duration(90) == "1h30m"
    with pytest.raises(ValueError):
        client.clean_duration("invalid")


def test_cleanup(client):
    """Test client cleanup."""
    with patch.object(client.http, "close") as mock_close:
        client.close()
        mock_close.assert_called_once()
