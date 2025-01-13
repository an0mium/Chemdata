"""Tests for FastAPI middleware implementations."""

import pytest
from fastapi import FastAPI, Request
from fastapi.testclient import TestClient
from ..middleware import RateLimitMiddleware, CacheMiddleware


@pytest.fixture
def app():
    """Create test application."""
    app = FastAPI()

    @app.get("/test")
    async def test_endpoint():
        return {"message": "test"}

    @app.get("/test/{id}")
    async def test_with_param(id: str):
        return {"id": id}

    return app


@pytest.fixture
def rate_limited_app(app):
    """Create rate limited application."""
    app.add_middleware(RateLimitMiddleware, limit=2, window=1)
    return app


@pytest.fixture
def cached_app(app):
    """Create cached application."""
    app.add_middleware(CacheMiddleware, max_size=10, ttl=1)
    return app


def test_rate_limit_middleware(rate_limited_app):
    """Test rate limit middleware."""
    client = TestClient(rate_limited_app)

    # First request should succeed
    response = client.get("/test")
    assert response.status_code == 200
    assert "x-ratelimit-limit" in response.headers
    assert "x-ratelimit-remaining" in response.headers
    assert "x-ratelimit-reset" in response.headers

    # Second request should succeed
    response = client.get("/test")
    assert response.status_code == 200

    # Third request should be rate limited
    response = client.get("/test")
    assert response.status_code == 429
    assert "retry-after" in response.headers


def test_rate_limit_per_endpoint(rate_limited_app):
    """Test rate limiting applies per endpoint."""
    client = TestClient(rate_limited_app)

    # Use up limit on first endpoint
    client.get("/test")
    client.get("/test")
    response = client.get("/test")
    assert response.status_code == 429

    # Should still be able to access different endpoint
    response = client.get("/test/123")
    assert response.status_code == 200


def test_rate_limit_headers(rate_limited_app):
    """Test rate limit headers."""
    client = TestClient(rate_limited_app)

    response = client.get("/test")
    assert response.status_code == 200

    # Check header values
    assert int(response.headers["x-ratelimit-limit"]) == 2
    assert int(response.headers["x-ratelimit-remaining"]) == 1
    assert float(response.headers["x-ratelimit-reset"]) >= 0


def test_cache_middleware(cached_app):
    """Test cache middleware."""
    client = TestClient(cached_app)

    # First request should be a cache miss
    response = client.get("/test")
    assert response.status_code == 200
    assert response.headers["x-cache"] == "MISS"

    # Second request should be a cache hit
    response = client.get("/test")
    assert response.status_code == 200
    assert response.headers["x-cache"] == "HIT"


def test_cache_with_params(cached_app):
    """Test caching with URL parameters."""
    client = TestClient(cached_app)

    # Different parameters should cache separately
    response = client.get("/test/123")
    assert response.headers["x-cache"] == "MISS"

    response = client.get("/test/456")
    assert response.headers["x-cache"] == "MISS"

    response = client.get("/test/123")
    assert response.headers["x-cache"] == "HIT"


def test_cache_methods(cached_app):
    """Test caching different HTTP methods."""
    client = TestClient(cached_app)

    # GET requests should be cached
    response = client.get("/test")
    assert response.headers["x-cache"] == "MISS"
    response = client.get("/test")
    assert response.headers["x-cache"] == "HIT"

    # POST requests should not be cached
    @cached_app.post("/test")
    async def post_test():
        return {"message": "test"}

    response = client.post("/test")
    assert "x-cache" not in response.headers


def test_middleware_error_handling(rate_limited_app):
    """Test middleware error handling."""
    client = TestClient(rate_limited_app)

    # Add endpoint that raises error
    @rate_limited_app.get("/error")
    async def error_endpoint():
        raise ValueError("Test error")

    # Should still return rate limit headers on error
    response = client.get("/error")
    assert response.status_code == 500
    assert "x-ratelimit-limit" in response.headers
    assert "x-ratelimit-remaining" in response.headers
    assert "x-ratelimit-reset" in response.headers


def test_middleware_chaining(app):
    """Test chaining multiple middleware."""
    app.add_middleware(RateLimitMiddleware, limit=2, window=1)
    app.add_middleware(CacheMiddleware, max_size=10, ttl=1)
    client = TestClient(app)

    # First request
    response = client.get("/test")
    assert response.status_code == 200
    assert response.headers["x-cache"] == "MISS"
    assert "x-ratelimit-remaining" in response.headers

    # Second request should be cached but still rate limited
    response = client.get("/test")
    assert response.status_code == 200
    assert response.headers["x-cache"] == "HIT"
    assert "x-ratelimit-remaining" in response.headers

    # Third request should be rate limited even though cached
    response = client.get("/test")
    assert response.status_code == 429
