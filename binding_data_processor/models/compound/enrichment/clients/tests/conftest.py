"""Shared test fixtures for web enrichment clients."""

import json
from unittest.mock import Mock
import pytest

from ..base import WebClientError


@pytest.fixture
def mock_http_client():
    """Mock HTTP client."""
    mock = Mock()
    mock.get = Mock()
    mock.post = Mock()
    return mock


@pytest.fixture
def mock_circuit_breaker():
    """Mock circuit breaker."""
    mock = Mock()
    mock.is_open = False
    mock.record_success = Mock()
    mock.record_failure = Mock()
    return mock


@pytest.fixture
def mock_logger():
    """Mock logger."""
    mock = Mock()
    mock.info = Mock()
    mock.warning = Mock()
    mock.error = Mock()
    mock.debug = Mock()
    return mock


@pytest.fixture
def mock_response():
    """Mock HTTP response."""
    mock = Mock()
    mock.status_code = 200
    mock.json = Mock(return_value={"data": "test"})
    mock.text = "test response"
    return mock


@pytest.fixture
def mock_error_response():
    """Mock error HTTP response."""
    mock = Mock()
    mock.status_code = 429
    mock.json = Mock(return_value={"error": "Rate limit exceeded"})
    mock.text = "Rate limit exceeded"
    return mock


@pytest.fixture
def mock_cache_dir(tmp_path):
    """Create temporary cache directory."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    return cache_dir


@pytest.fixture
def mock_model_dir(tmp_path):
    """Create temporary model directory."""
    model_dir = tmp_path / "models"
    model_dir.mkdir()
    return model_dir


@pytest.fixture
def mock_cached_data(mock_cache_dir):
    """Create mock cached data."""
    data = {
        "test_data": {
            "id": "123",
            "name": "Test",
            "timestamp": "2023-01-01T00:00:00Z"
        }
    }
    
    cache_file = mock_cache_dir / "test_cache.json"
    cache_file.write_text(json.dumps(data))
    
    return data, cache_file


@pytest.fixture
def mock_web_client_error():
    """Create mock web client error."""
    return WebClientError(
        message="Test error",
        status_code=500,
        response={"error": "Internal server error"}
    )
