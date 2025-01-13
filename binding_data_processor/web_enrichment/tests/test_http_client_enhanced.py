"""Tests for enhanced HTTP client."""

import pytest
import requests
import responses
from pathlib import Path
from datetime import datetime, timedelta
import json

from ..http_client_enhanced import HTTPClientEnhanced
from ...pipeline.infrastructure.circuit_breaker import (
    CircuitBreakerConfig,
    CircuitBreakerError,
)


@pytest.fixture
def client(tmp_path):
    """Create test client."""
    return HTTPClientEnhanced(
        name="test_client",
        cache_dir=tmp_path / "cache",
        rate_limit=0,  # Disable rate limiting for tests
        circuit_config=CircuitBreakerConfig(
            failure_threshold=2,
            failure_timeout=1,
            reset_timeout=1,
        ),
    )


@responses.activate
def test_get_success(client):
    """Test successful GET request."""
    # Mock response
    responses.add(
        responses.GET,
        "http://test.com/api",
        json={"status": "ok"},
        status=200,
    )

    # Make request
    response = client.get("http://test.com/api")
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}


@responses.activate
def test_get_with_circuit_breaker(client):
    """Test circuit breaker behavior."""
    # Mock failing endpoint
    responses.add(
        responses.GET,
        "http://test.com/api",
        status=500,
    )

    # Make requests until circuit opens
    with pytest.raises(requests.exceptions.HTTPError):
        client.get("http://test.com/api")
    with pytest.raises(requests.exceptions.HTTPError):
        client.get("http://test.com/api")

    # Circuit should be open
    with pytest.raises(CircuitBreakerError):
        client.get("http://test.com/api")


@responses.activate
def test_get_with_fallback(client):
    """Test fallback when circuit is open."""
    # Mock failing endpoint
    responses.add(
        responses.GET,
        "http://test.com/api",
        status=500,
    )

    # Create fallback response
    fallback_response = requests.Response()
    fallback_response.status_code = 200
    fallback_response._content = json.dumps({"status": "fallback"}).encode()

    # Make requests until circuit opens
    with pytest.raises(requests.exceptions.HTTPError):
        client.get("http://test.com/api")
    with pytest.raises(requests.exceptions.HTTPError):
        client.get("http://test.com/api")

    # Use fallback
    response = client.get(
        "http://test.com/api",
        fallback=lambda: fallback_response,
    )
    assert response.status_code == 200
    assert response.json() == {"status": "fallback"}


@responses.activate
def test_get_with_cache(client):
    """Test response caching."""
    # Mock response
    responses.add(
        responses.GET,
        "http://test.com/api",
        json={"status": "ok"},
        status=200,
    )

    # First request should hit network
    response1 = client.get("http://test.com/api")
    assert response1.status_code == 200

    # Second request should use cache
    response2 = client.get("http://test.com/api")
    assert response2.status_code == 200
    assert response2.json() == {"status": "ok"}

    # Only one actual request should have been made
    assert len(responses.calls) == 1


@responses.activate
def test_post_success(client):
    """Test successful POST request."""
    # Mock response
    responses.add(
        responses.POST,
        "http://test.com/api",
        json={"status": "created"},
        status=201,
    )

    # Make request
    response = client.post(
        "http://test.com/api",
        json_data={"key": "value"},
    )
    assert response.status_code == 201
    assert response.json() == {"status": "created"}


@responses.activate
def test_post_with_circuit_breaker(client):
    """Test circuit breaker with POST requests."""
    # Mock failing endpoint
    responses.add(
        responses.POST,
        "http://test.com/api",
        status=500,
    )

    # Make requests until circuit opens
    with pytest.raises(requests.exceptions.HTTPError):
        client.post("http://test.com/api")
    with pytest.raises(requests.exceptions.HTTPError):
        client.post("http://test.com/api")

    # Circuit should be open
    with pytest.raises(CircuitBreakerError):
        client.post("http://test.com/api")


def test_metrics(client):
    """Test metrics collection."""
    metrics = client.get_metrics()
    assert "circuit_breaker" in metrics
    assert "cache" in metrics
    assert "memory_size" in metrics["cache"]
    assert "disk_size" in metrics["cache"]


def test_cleanup(client):
    """Test cleanup on close."""
    # Add some cache entries
    client.memory_cache["key"] = "value"

    # Close client
    client.close()

    # Cache should be cleared
    assert len(client.memory_cache) == 0
