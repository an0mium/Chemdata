"""Tests for base web client."""

import pytest
from unittest.mock import Mock

from ..base import BaseWebClient, WebClientError, RetryConfig


@pytest.fixture
def mock_http_client():
    """Mock HTTP client fixture."""
    return Mock()


@pytest.fixture
def mock_logger():
    """Mock logger fixture."""
    return Mock()


@pytest.fixture
def base_client(mock_http_client, mock_logger):
    """Base web client fixture."""
    return BaseWebClient(
        name="test_client",
        http_client=mock_http_client,
        logger=mock_logger,
    )


def test_init(mock_http_client, mock_logger):
    """Test client initialization."""
    client = BaseWebClient(
        name="test_client",
        http_client=mock_http_client,
        logger=mock_logger,
    )
    
    assert client.name == "test_client"
    assert client.http == mock_http_client
    assert client.logger == mock_logger


def test_init_defaults():
    """Test initialization with default values."""
    client = BaseWebClient(name="test_client")
    
    assert client.name == "test_client"
    assert client.http is not None
    assert client.logger is not None


def test_get_success(base_client, mock_http_client):
    """Test successful GET request."""
    mock_http_client.get.return_value = {"data": "test"}
    
    response = base_client.get(
        url="https://api.example.com/test",
        params={"key": "value"},
    )
    
    assert response == {"data": "test"}
    mock_http_client.get.assert_called_with(
        url="https://api.example.com/test",
        params={"key": "value"},
    )


def test_get_error(base_client, mock_http_client):
    """Test GET request error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "API error",
        status_code=500,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        base_client.get(url="https://api.example.com/test")
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_post_success(base_client, mock_http_client):
    """Test successful POST request."""
    mock_http_client.post.return_value = {"data": "test"}
    
    response = base_client.post(
        url="https://api.example.com/test",
        data={"key": "value"},
    )
    
    assert response == {"data": "test"}
    mock_http_client.post.assert_called_with(
        url="https://api.example.com/test",
        data={"key": "value"},
    )


def test_post_error(base_client, mock_http_client):
    """Test POST request error handling."""
    mock_http_client.post.side_effect = WebClientError(
        "API error",
        status_code=500,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        base_client.post(
            url="https://api.example.com/test",
            data={"key": "value"},
        )
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_retry_success(base_client, mock_http_client):
    """Test successful retry after temporary failure."""
    # Fail twice, succeed on third try
    mock_http_client.get.side_effect = [
        WebClientError("Rate limit", status_code=429),
        WebClientError("Rate limit", status_code=429),
        {"data": "test"},
    ]
    
    retry_config = RetryConfig(
        max_retries=3,
        retry_delay=0.1,
        retry_codes=[429],
    )
    
    response = base_client.get(
        url="https://api.example.com/test",
        retry_config=retry_config,
    )
    
    assert response == {"data": "test"}
    assert mock_http_client.get.call_count == 3


def test_retry_exhausted(base_client, mock_http_client):
    """Test retry exhaustion."""
    # Fail all attempts
    mock_http_client.get.side_effect = WebClientError(
        "Rate limit",
        status_code=429,
    )
    
    retry_config = RetryConfig(
        max_retries=3,
        retry_delay=0.1,
        retry_codes=[429],
    )
    
    with pytest.raises(WebClientError) as exc_info:
        base_client.get(
            url="https://api.example.com/test",
            retry_config=retry_config,
        )
    
    assert exc_info.value.status_code == 429
    assert "Rate limit" in str(exc_info.value)
    assert mock_http_client.get.call_count == 3


def test_retry_non_retryable(base_client, mock_http_client):
    """Test non-retryable error."""
    mock_http_client.get.side_effect = WebClientError(
        "Bad request",
        status_code=400,
    )
    
    retry_config = RetryConfig(
        max_retries=3,
        retry_delay=0.1,
        retry_codes=[429],
    )
    
    with pytest.raises(WebClientError) as exc_info:
        base_client.get(
            url="https://api.example.com/test",
            retry_config=retry_config,
        )
    
    assert exc_info.value.status_code == 400
    assert "Bad request" in str(exc_info.value)
    assert mock_http_client.get.call_count == 1


def test_validate_response_success(base_client):
    """Test successful response validation."""
    response = {
        "data": "test",
        "status": "success",
    }
    
    # Should not raise any errors
    base_client._validate_response(response)


def test_validate_response_error(base_client):
    """Test response validation error."""
    response = {
        "error": "Something went wrong",
        "status": "error",
    }
    
    with pytest.raises(WebClientError) as exc_info:
        base_client._validate_response(response)
    
    assert "API error response" in str(exc_info.value)


def test_validate_response_missing_required(base_client):
    """Test response validation with missing required fields."""
    response = {
        "status": "success",
        # Missing required data field
    }
    
    with pytest.raises(WebClientError) as exc_info:
        base_client._validate_response(
            response,
            required_fields=["data"],
        )
    
    assert "Missing required field" in str(exc_info.value)


def test_validate_response_invalid_type(base_client):
    """Test response validation with invalid type."""
    response = "not a dict"
    
    with pytest.raises(WebClientError) as exc_info:
        base_client._validate_response(response)
    
    assert "Invalid response type" in str(exc_info.value)


def test_log_request(base_client, mock_logger):
    """Test request logging."""
    base_client._log_request(
        method="GET",
        url="https://api.example.com/test",
        params={"key": "value"},
    )
    
    mock_logger.debug.assert_called_with(
        "Making GET request to https://api.example.com/test with params {'key': 'value'}"
    )


def test_log_response(base_client, mock_logger):
    """Test response logging."""
    base_client._log_response(
        method="GET",
        url="https://api.example.com/test",
        response={"data": "test"},
    )
    
    mock_logger.debug.assert_called_with(
        "Received response from GET https://api.example.com/test: {'data': 'test'}"
    )


def test_log_error(base_client, mock_logger):
    """Test error logging."""
    error = WebClientError("API error", status_code=500)
    
    base_client._log_error(
        method="GET",
        url="https://api.example.com/test",
        error=error,
    )
    
    mock_logger.error.assert_called_with(
        "Error in GET request to https://api.example.com/test: API error (500)"
    )
