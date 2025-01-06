"""Tests for HTTP client.

This module tests the HTTPClient which handles network requests and responses
for all web enrichment clients, including retry logic, error handling, and
response validation.
"""

import json
import pytest
from unittest.mock import Mock, patch

import requests

from ..http import HTTPClient
from ..base import WebClientError


@pytest.fixture
def mock_session():
    """Mock requests session fixture."""
    return Mock()


@pytest.fixture
def mock_logger():
    """Mock logger fixture."""
    return Mock()


@pytest.fixture
def mock_success_response():
    """Mock successful response fixture."""
    response = Mock()
    response.status_code = 200
    response.json.return_value = {"data": "test"}
    return response


@pytest.fixture
def mock_error_response():
    """Mock error response fixture."""
    response = Mock()
    response.status_code = 429
    response.text = "Rate limit exceeded"
    response.json.return_value = {"error": "Rate limit exceeded"}
    return response


@pytest.fixture
def http_client(mock_session, mock_logger):
    """HTTP client fixture with mocked dependencies."""
    with patch("requests.Session") as mock_session_class:
        mock_session_class.return_value = mock_session
        return HTTPClient(
            base_url="https://api.example.com",
            logger=mock_logger,
            timeout=30,
            max_retries=3,
            retry_delay=1,
        )


def test_init(mock_logger):
    """Test client initialization."""
    client = HTTPClient(
        base_url="https://api.example.com",
        logger=mock_logger,
        timeout=30,
        max_retries=3,
        retry_delay=1,
    )
    
    assert client.base_url == "https://api.example.com"
    assert client.logger == mock_logger
    assert isinstance(client.session, requests.Session)
    assert client.timeout == 30
    assert client.max_retries == 3
    assert client.retry_delay == 1


def test_init_defaults():
    """Test initialization with default values."""
    client = HTTPClient()
    
    assert client.base_url == ""
    assert client.logger is not None
    assert isinstance(client.session, requests.Session)
    assert client.timeout == 30
    assert client.max_retries == 3
    assert client.retry_delay == 1


def test_init_custom_values():
    """Test initialization with custom values."""
    client = HTTPClient(
        base_url="https://api.example.com",
        timeout=60,
        max_retries=5,
        retry_delay=2,
    )
    
    assert client.base_url == "https://api.example.com"
    assert client.timeout == 60
    assert client.max_retries == 5
    assert client.retry_delay == 2


def test_get_success_base_url(http_client, mock_session, mock_success_response):
    """Test successful GET request with base URL."""
    mock_session.get.return_value = mock_success_response
    
    response = http_client.get(
        "/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
    )
    
    assert response == {"data": "test"}
    mock_session.get.assert_called_with(
        url="https://api.example.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
        timeout=30,
    )


def test_get_success_direct_url(http_client, mock_session, mock_success_response):
    """Test successful GET request with direct URL."""
    mock_session.get.return_value = mock_success_response
    
    response = http_client.get(
        url="https://other-api.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
    )
    
    assert response == {"data": "test"}
    mock_session.get.assert_called_with(
        url="https://other-api.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
        timeout=30,
    )


def test_get_retry_success(
    http_client,
    mock_session,
    mock_error_response,
    mock_success_response,
):
    """Test GET request with successful retry."""
    mock_session.get.side_effect = [
        mock_error_response,
        mock_success_response,
    ]
    
    response = http_client.get("/test")
    
    assert response == {"data": "test"}
    assert mock_session.get.call_count == 2


def test_get_retry_failure(http_client, mock_session, mock_error_response):
    """Test GET request with retry failure."""
    mock_session.get.return_value = mock_error_response
    
    with pytest.raises(WebClientError) as exc_info:
        http_client.get("/test")
    
    assert exc_info.value.status_code == 429
    assert "Rate limit exceeded" in str(exc_info.value)
    assert mock_session.get.call_count == 3


def test_get_connection_error(http_client, mock_session):
    """Test GET request with connection error."""
    mock_session.get.side_effect = requests.ConnectionError("Connection failed")
    
    with pytest.raises(WebClientError) as exc_info:
        http_client.get("/test")
    
    assert "Connection failed" in str(exc_info.value)
    assert exc_info.value.status_code == 503


def test_get_timeout_error(http_client, mock_session):
    """Test GET request with timeout error."""
    mock_session.get.side_effect = requests.Timeout("Request timed out")
    
    with pytest.raises(WebClientError) as exc_info:
        http_client.get("/test")
    
    assert "Request timed out" in str(exc_info.value)
    assert exc_info.value.status_code == 504


def test_post_success(http_client, mock_session, mock_success_response):
    """Test successful POST request."""
    mock_session.post.return_value = mock_success_response
    
    response = http_client.post(
        "/test",
        json={"key": "value"},
        headers={"Content-Type": "application/json"},
    )
    
    assert response == {"data": "test"}
    mock_session.post.assert_called_with(
        url="https://api.example.com/test",
        json={"key": "value"},
        headers={"Content-Type": "application/json"},
        timeout=30,
    )


def test_post_retry_success(
    http_client,
    mock_session,
    mock_error_response,
    mock_success_response,
):
    """Test POST request with successful retry."""
    mock_session.post.side_effect = [
        mock_error_response,
        mock_success_response,
    ]
    
    response = http_client.post("/test", json={})
    
    assert response == {"data": "test"}
    assert mock_session.post.call_count == 2


def test_post_retry_failure(http_client, mock_session, mock_error_response):
    """Test POST request with retry failure."""
    mock_session.post.return_value = mock_error_response
    
    with pytest.raises(WebClientError) as exc_info:
        http_client.post("/test", json={})
    
    assert exc_info.value.status_code == 429
    assert "Rate limit exceeded" in str(exc_info.value)
    assert mock_session.post.call_count == 3


def test_put_success(http_client, mock_session, mock_success_response):
    """Test successful PUT request."""
    mock_session.put.return_value = mock_success_response
    
    response = http_client.put(
        "/test",
        json={"key": "value"},
        headers={"Content-Type": "application/json"},
    )
    
    assert response == {"data": "test"}
    mock_session.put.assert_called_with(
        url="https://api.example.com/test",
        json={"key": "value"},
        headers={"Content-Type": "application/json"},
        timeout=30,
    )


def test_put_retry_success(
    http_client,
    mock_session,
    mock_error_response,
    mock_success_response,
):
    """Test PUT request with successful retry."""
    mock_session.put.side_effect = [
        mock_error_response,
        mock_success_response,
    ]
    
    response = http_client.put("/test", json={})
    
    assert response == {"data": "test"}
    assert mock_session.put.call_count == 2


def test_put_retry_failure(http_client, mock_session, mock_error_response):
    """Test PUT request with retry failure."""
    mock_session.put.return_value = mock_error_response
    
    with pytest.raises(WebClientError) as exc_info:
        http_client.put("/test", json={})
    
    assert exc_info.value.status_code == 429
    assert "Rate limit exceeded" in str(exc_info.value)
    assert mock_session.put.call_count == 3


def test_delete_success(http_client, mock_session, mock_success_response):
    """Test successful DELETE request."""
    mock_session.delete.return_value = mock_success_response
    
    response = http_client.delete(
        "/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
    )
    
    assert response == {"data": "test"}
    mock_session.delete.assert_called_with(
        url="https://api.example.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
        timeout=30,
    )


def test_delete_retry_success(
    http_client,
    mock_session,
    mock_error_response,
    mock_success_response,
):
    """Test DELETE request with successful retry."""
    mock_session.delete.side_effect = [
        mock_error_response,
        mock_success_response,
    ]
    
    response = http_client.delete("/test")
    
    assert response == {"data": "test"}
    assert mock_session.delete.call_count == 2


def test_delete_retry_failure(http_client, mock_session, mock_error_response):
    """Test DELETE request with retry failure."""
    mock_session.delete.return_value = mock_error_response
    
    with pytest.raises(WebClientError) as exc_info:
        http_client.delete("/test")
    
    assert exc_info.value.status_code == 429
    assert "Rate limit exceeded" in str(exc_info.value)
    assert mock_session.delete.call_count == 3


def test_build_url(http_client):
    """Test URL building."""
    # Test with leading slash
    url = http_client._build_url("/test")
    assert url == "https://api.example.com/test"
    
    # Test without leading slash
    url = http_client._build_url("test")
    assert url == "https://api.example.com/test"
    
    # Test with trailing slash in base_url
    http_client.base_url = "https://api.example.com/"
    url = http_client._build_url("/test")
    assert url == "https://api.example.com/test"
    
    # Test with full URL
    url = http_client._build_url("https://other-api.com/test")
    assert url == "https://other-api.com/test"


def test_handle_response_success(http_client, mock_success_response):
    """Test successful response handling."""
    response = http_client._handle_response(mock_success_response)
    assert response == {"data": "test"}


def test_handle_response_error(http_client, mock_error_response):
    """Test error response handling."""
    with pytest.raises(WebClientError) as exc_info:
        http_client._handle_response(mock_error_response)
    
    assert exc_info.value.status_code == 429
    assert "Rate limit exceeded" in str(exc_info.value)


def test_handle_response_invalid_json(http_client):
    """Test response handling with invalid JSON."""
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.json.side_effect = json.JSONDecodeError(
        "Invalid JSON",
        doc="not json",
        pos=0,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        http_client._handle_response(mock_response)
    
    assert "Invalid JSON response" in str(exc_info.value)
    assert exc_info.value.status_code == 502


def test_handle_response_no_json(http_client):
    """Test response handling with no JSON content."""
    mock_response = Mock()
    mock_response.status_code = 204
    mock_response.json.side_effect = ValueError("No JSON content")
    
    response = http_client._handle_response(mock_response)
    assert response is None


def test_custom_timeout(http_client, mock_session, mock_success_response):
    """Test custom timeout setting."""
    mock_session.get.return_value = mock_success_response
    
    http_client.get("/test", timeout=60)
    
    mock_session.get.assert_called_with(
        url="https://api.example.com/test",
        headers={},
        timeout=60,
    )


def test_custom_headers(http_client, mock_session, mock_success_response):
    """Test custom headers setting."""
    mock_session.get.return_value = mock_success_response
    
    headers = {
        "Authorization": "Bearer token",
        "Custom-Header": "value",
    }
    
    http_client.get("/test", headers=headers)
    
    mock_session.get.assert_called_with(
        url="https://api.example.com/test",
        headers=headers,
        timeout=30,
    )
