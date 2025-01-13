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
        "/api/compounds",
        headers={
            "origin": "http://localhost:3000",
            "access-control-request-method": "GET",
        },
    )
    assert response.status_code == 200
    assert "access-control-allow-origin" in response.headers
    assert "access-control-allow-methods" in response.headers
    assert "access-control-allow-headers" in response.headers


def test_get_compounds(client):
    """Test get compounds endpoint."""
    response = client.get("/api/compounds")
    assert response.status_code == 200
    data = response.json()
    assert "compounds" in data
    assert "total" in data
    assert "page" in data
    assert "per_page" in data
    assert "has_next" in data
    assert "has_prev" in data


def test_get_compound(client):
    """Test get single compound endpoint."""
    # Test non-existent compound
    response = client.get("/api/compounds/invalid-id")
    assert response.status_code == 404

    # TODO: Test existing compound once database is implemented
    # response = client.get("/api/compounds/valid-id")
    # assert response.status_code == 200
    # data = response.json()
    # assert data["id"] == "valid-id"


def test_search_compounds(client):
    """Test search compounds endpoint."""
    query = {
        "text": "test compound",
        "tags": ["validated"],
    }
    response = client.post("/api/search", json=query)
    assert response.status_code == 200
    data = response.json()
    assert "compounds" in data
    assert "total" in data
    assert "page" in data
    assert "per_page" in data


def test_export_compounds(client):
    """Test export compounds endpoint."""
    options = {
        "format": "json",
        "compound_ids": ["test-1", "test-2"],
        "include_fields": ["id", "name", "smiles"],
    }
    response = client.post("/api/export", json=options)
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


def test_validation_error(client):
    """Test validation error handling."""
    # Invalid compound ID format
    response = client.get("/api/compounds/123")  # Assuming IDs must be UUIDs
    assert response.status_code == 422
    data = response.json()
    assert "detail" in data

    # Invalid search query
    response = client.post("/api/search", json={"invalid": "query"})
    assert response.status_code == 422
    data = response.json()
    assert "detail" in data

    # Invalid export format
    response = client.post("/api/export", json={"format": "invalid"})
    assert response.status_code == 422
    data = response.json()
    assert "detail" in data


def test_rate_limiting(client):
    """Test rate limiting."""
    # Make requests up to limit
    for _ in range(100):
        response = client.get("/api/compounds")
        assert response.status_code == 200

    # Next request should be rate limited
    response = client.get("/api/compounds")
    assert response.status_code == 429
    assert "retry-after" in response.headers


def test_caching(client):
    """Test response caching."""
    # First request should be uncached
    response = client.get("/api/compounds/test-1")
    assert response.status_code == 404
    assert "x-cache" not in response.headers

    # Subsequent request should be cached
    response = client.get("/api/compounds/test-1")
    assert response.status_code == 404
    assert response.headers.get("x-cache") == "HIT"


def test_viewport_handling(client):
    """Test viewport-aware responses."""
    # Mobile viewport
    response = client.get(
        "/api/compounds",
        headers={"viewport-width": "375"},
    )
    assert response.status_code == 200
    data = response.json()
    for compound in data["compounds"]:
        assert "spectral_data" not in compound
        assert "crystal_data" not in compound

    # Desktop viewport
    response = client.get(
        "/api/compounds",
        headers={"viewport-width": "1920"},
    )
    assert response.status_code == 200
    data = response.json()
    for compound in data["compounds"]:
        if "spectral_data" in compound:
            assert isinstance(compound["spectral_data"], dict)
        if "crystal_data" in compound:
            assert isinstance(compound["crystal_data"], dict)
