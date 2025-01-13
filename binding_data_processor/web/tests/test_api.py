"""Tests for FastAPI application."""

import pytest
from fastapi.testclient import TestClient
from ..api.app import create_app


@pytest.fixture
def app():
    """Create test application."""
    return create_app()


@pytest.fixture
def client(app):
    """Create test client."""
    return TestClient(app)


def test_health_check(client):
    """Test health check endpoint."""
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.json() == {"status": "healthy"}


def test_rate_limit_info(client):
    """Test rate limit info endpoint."""
    response = client.get("/api/rate-limit")
    assert response.status_code == 200
    data = response.json()
    assert "limit" in data
    assert "remaining" in data
    assert "reset" in data


def test_cache_info(client):
    """Test cache info endpoint."""
    response = client.get("/api/cache")
    assert response.status_code == 200
    data = response.json()
    assert "size" in data
    assert "hits" in data
    assert "misses" in data
    assert "ttl" in data


def test_openapi_docs(client):
    """Test OpenAPI documentation endpoints."""
    # Test Swagger UI
    response = client.get("/api/docs")
    assert response.status_code == 200
    assert "text/html" in response.headers["content-type"]

    # Test ReDoc
    response = client.get("/api/redoc")
    assert response.status_code == 200
    assert "text/html" in response.headers["content-type"]

    # Test OpenAPI JSON
    response = client.get("/api/openapi.json")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    data = response.json()
    assert "openapi" in data
    assert "paths" in data


def test_cors_headers(client):
    """Test CORS headers."""
    response = client.options(
        "/api/health",
        headers={
            "origin": "http://localhost:3000",
            "access-control-request-method": "GET",
        },
    )
    assert response.status_code == 200
    assert response.headers["access-control-allow-origin"] == "*"
    assert "GET" in response.headers["access-control-allow-methods"]


def test_rate_limiting(client):
    """Test rate limiting middleware."""
    # Make requests up to limit
    for _ in range(100):
        response = client.get("/api/health")
        assert response.status_code == 200
        assert "x-ratelimit-limit" in response.headers
        assert "x-ratelimit-remaining" in response.headers
        assert "x-ratelimit-reset" in response.headers

    # Next request should be rate limited
    response = client.get("/api/health")
    assert response.status_code == 429
    assert "retry-after" in response.headers


def test_caching(client):
    """Test caching middleware."""
    # First request should be uncached
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.headers["x-cache"] == "MISS"

    # Second request should be cached
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.headers["x-cache"] == "HIT"


def test_validation_error(client):
    """Test validation error handling."""
    response = client.post(
        "/api/search",
        json={"invalid": "query"},  # Invalid search query
    )
    assert response.status_code == 422
    data = response.json()
    assert "detail" in data


def test_error_with_rate_limit_headers(client):
    """Test rate limit headers on error responses."""
    response = client.get("/api/nonexistent")
    assert response.status_code == 404
    assert "x-ratelimit-limit" in response.headers
    assert "x-ratelimit-remaining" in response.headers
    assert "x-ratelimit-reset" in response.headers


def test_middleware_order(client):
    """Test middleware execution order."""
    # First request should be uncached and count against rate limit
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.headers["x-cache"] == "MISS"
    assert "x-ratelimit-remaining" in response.headers
    initial_remaining = int(response.headers["x-ratelimit-remaining"])

    # Second request should be cached but still count against rate limit
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.headers["x-cache"] == "HIT"
    assert int(response.headers["x-ratelimit-remaining"]) == initial_remaining - 1
