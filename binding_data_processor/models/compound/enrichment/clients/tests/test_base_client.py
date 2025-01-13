"""Tests for web client."""

import pytest
from unittest.mock import Mock

from ....web_enrichment.clients.base import WebClient, WebClientError


@pytest.fixture
def mock_http_client():
    """Mock HTTP client fixture."""
    return Mock()


@pytest.fixture
def mock_logger():
    """Mock logger fixture."""
    return Mock()


@pytest.fixture
def web_client(mock_http_client, mock_logger):
    """Web client fixture."""
    return WebClient(
        name="test_client",
        base_url="https://api.example.com",
        data_source="test",
        http_client=mock_http_client,
        logger=mock_logger,
    )


def test_init(mock_http_client, mock_logger):
    """Test client initialization."""
    client = WebClient(
        name="test_client",
        base_url="https://api.example.com",
        data_source="test",
        http_client=mock_http_client,
        logger=mock_logger,
    )

    assert client.name == "test_client"
    assert client.base_url == "https://api.example.com"
    assert client.data_source == "test"
    assert client.http == mock_http_client
    assert client.logger == mock_logger


def test_init_defaults():
    """Test initialization with default values."""
    client = WebClient(
        name="test_client",
        base_url="https://api.example.com",
        data_source="test",
    )

    assert client.name == "test_client"
    assert client.base_url == "https://api.example.com"
    assert client.data_source == "test"
    assert client.http is not None
    assert client.logger is not None


def test_get_success(web_client, mock_http_client):
    """Test successful GET request."""
    mock_http_client.get.return_value = {"data": "test"}

    response = web_client.get(
        url="test",
        params={"key": "value"},
    )

    assert response == {"data": "test"}
    mock_http_client.get.assert_called_with(
        url="https://api.example.com/test",
        params={"key": "value"},
        use_cache=True,
    )


def test_get_error(web_client, mock_http_client):
    """Test GET request error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "API error",
        status_code=500,
    )

    with pytest.raises(WebClientError) as exc_info:
        web_client.get(url="test")

    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_post_success(web_client, mock_http_client):
    """Test successful POST request."""
    mock_http_client.post.return_value = {"data": "test"}

    response = web_client.post(
        url="test",
        data={"key": "value"},
    )

    assert response == {"data": "test"}
    mock_http_client.post.assert_called_with(
        url="https://api.example.com/test",
        data={"key": "value"},
    )


def test_post_error(web_client, mock_http_client):
    """Test POST request error handling."""
    mock_http_client.post.side_effect = WebClientError(
        "API error",
        status_code=500,
    )

    with pytest.raises(WebClientError) as exc_info:
        web_client.post(
            url="test",
            data={"key": "value"},
        )

    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_validate_response_success(web_client):
    """Test successful response validation."""
    response = {
        "data": "test",
        "status": "success",
    }

    # Should not raise any errors
    web_client._validate_response(response)


def test_validate_response_error(web_client):
    """Test response validation error."""
    response = {
        "error": "Something went wrong",
        "status": "error",
    }

    with pytest.raises(WebClientError) as exc_info:
        web_client._validate_response(response)

    assert "API error response" in str(exc_info.value)


def test_validate_response_missing_required(web_client):
    """Test response validation with missing required fields."""
    response = {
        "status": "success",
        # Missing required data field
    }

    with pytest.raises(WebClientError) as exc_info:
        web_client._validate_response(
            response,
            required_fields=["data"],
        )

    assert "Missing required field" in str(exc_info.value)


def test_validate_response_invalid_type(web_client):
    """Test response validation with invalid type."""
    response = "not a dict"

    with pytest.raises(WebClientError) as exc_info:
        web_client._validate_response(response)

    assert "Invalid response type" in str(exc_info.value)


def test_log_request(web_client, mock_logger):
    """Test request logging."""
    web_client._log_request(
        method="GET",
        url="test",
        params={"key": "value"},
    )

    mock_logger.debug.assert_called_with("Making GET request to https://api.example.com/test with params {'key': 'value'}")


def test_log_response(web_client, mock_logger):
    """Test response logging."""
    web_client._log_response(
        method="GET",
        url="test",
        response={"data": "test"},
    )

    mock_logger.debug.assert_called_with("Received response from GET https://api.example.com/test: {'data': 'test'}")


def test_log_error(web_client, mock_logger):
    """Test error logging."""
    error = WebClientError("API error", status_code=500)

    web_client._log_error(
        method="GET",
        url="test",
        error=error,
    )

    mock_logger.error.assert_called_with("Error in GET request to https://api.example.com/test: API error (500)")
