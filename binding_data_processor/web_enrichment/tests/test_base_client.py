"""Tests for base web enrichment client."""

import logging
import pytest
from unittest.mock import Mock, patch

import requests

from ..clients.base import WebClient, WebClientError, ValidationError
from ..validation.schema import DataSource, ValidationLevel
from ..validation.data import ValidationConfig
from ...pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig


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
        circuit_config=CircuitBreakerConfig(
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
    """Test rate limit error handling with retries."""
    # Mock responses for retry sequence
    mock_responses = [
        # First attempt: Rate limit error
        Mock(
            status_code=429,
            headers={"Retry-After": "1"},
            json=lambda: {"error": "Rate limit exceeded"},
        ),
        # Second attempt: Success
        Mock(status_code=200, json=lambda: {"name": "test", "value": 50}),
    ]

    # Mock request to return sequence of responses
    with patch.object(client.http, "request", side_effect=mock_responses):
        data = client.request(
            method="GET",
            endpoint="/test",
            schema_name="test",
            retry_on_rate_limit=True,
        )

    # Verify successful retry
    assert data == {"name": "test", "value": 50}
    assert client.processed_items == 1
    assert client.failed_items == 0
    assert client.validation_errors == 0
    assert client.http_errors == 0


def test_circuit_breaker_state_transitions(client):
    """Test circuit breaker state transitions."""
    mock_responses = [
        # Initial failures to open circuit
        *[Mock(side_effect=requests.exceptions.RequestException("Error"))] * 3,
        # Recovery period
        Mock(side_effect=requests.exceptions.RequestException("Error")),
        # Half-open test request succeeds
        Mock(status_code=200, json=lambda: {"name": "test"}),
        # Circuit fully closed, normal operation
        Mock(status_code=200, json=lambda: {"name": "test"}),
    ]

    with patch.object(client.http, "request", side_effect=mock_responses):
        # Trigger circuit open
        for _ in range(3):
            with pytest.raises(WebClientError):
                client.request(method="GET", endpoint="/test")

        assert client.circuit.state == "open"

        # Wait for recovery timeout
        client.circuit._last_error_time -= 61  # Force timeout

        # Half-open test request
        data = client.request(method="GET", endpoint="/test")
        assert client.circuit.state == "closed"
        assert data == {"name": "test"}

        # Verify normal operation
        data = client.request(method="GET", endpoint="/test")
        assert data == {"name": "test"}


def test_request_batching(client):
    """Test request batching."""
    mock_response = Mock(status_code=200, json=lambda: [{"id": 1}, {"id": 2}, {"id": 3}])

    with patch.object(client.http, "request", return_value=mock_response):
        results = client.batch_request(
            method="GET",
            endpoint="/test",
            ids=[1, 2, 3],
            batch_size=2,
        )

    assert len(results) == 3
    assert all(r["id"] in [1, 2, 3] for r in results)
    # Verify batching
    assert client.http.request.call_count == 2


def test_enhanced_error_recovery(client):
    """Test enhanced error recovery with fallback strategies."""
    mock_responses = [
        # Primary endpoint fails
        Mock(side_effect=requests.exceptions.ConnectionError("Primary failed")),
        # Fallback endpoint succeeds
        Mock(status_code=200, json=lambda: {"name": "test", "source": "fallback"}),
    ]

    with patch.object(client.http, "request", side_effect=mock_responses):
        data = client.request(
            method="GET",
            endpoint="/test",
            fallback_endpoints=["/test-fallback"],
            retry_strategy="fallback_endpoints",
        )

    assert data["source"] == "fallback"
    assert client.processed_items == 1
    assert client.failed_items == 0


def test_response_parsing_with_retries(client):
    """Test response parsing with retry on parse error."""
    mock_responses = [
        # First attempt: Invalid JSON
        Mock(status_code=200, json=Mock(side_effect=ValueError("Invalid JSON"))),
        # Second attempt: Valid JSON
        Mock(status_code=200, json=lambda: {"name": "test", "value": 50}),
    ]

    with patch.object(client.http, "request", side_effect=mock_responses):
        data = client.request(
            method="GET",
            endpoint="/test",
            retry_on_parse_error=True,
        )

    assert data == {"name": "test", "value": 50}
    assert client.processed_items == 1
    assert client.validation_errors == 0


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
