"""Tests for base web enrichment client.

This module tests the BaseClient which provides common functionality for all
web enrichment clients.
"""

import pytest
from unittest.mock import Mock

from ..base import BaseClient, WebClientError


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
    """Base client fixture with mocked dependencies."""
    return BaseClient(
        name="test",
        http_client=mock_http_client,
        logger=mock_logger,
    )


def test_init(mock_http_client, mock_logger):
    """Test client initialization."""
    client = BaseClient(
        name="test",
        http_client=mock_http_client,
        logger=mock_logger,
    )
    
    assert client.name == "test"
    assert client.http == mock_http_client
    assert client.logger == mock_logger


def test_init_defaults():
    """Test initialization with default values."""
    client = BaseClient(name="test")
    
    assert client.name == "test"
    assert client.http is not None
    assert client.logger is not None


def test_init_invalid_name():
    """Test initialization with invalid name."""
    with pytest.raises(ValueError) as exc_info:
        BaseClient(name="")
    
    assert "Client name cannot be empty" in str(exc_info.value)


def test_get_success(base_client, mock_http_client):
    """Test successful GET request."""
    mock_response = {"data": "test"}
    mock_http_client.get.return_value = mock_response
    
    response = base_client.get(
        url="https://api.example.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
    )
    
    assert response == mock_response
    mock_http_client.get.assert_called_with(
        url="https://api.example.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
    )


def test_get_error(base_client, mock_http_client):
    """Test GET request error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "API error",
        status_code=500,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        base_client.get(
            url="https://api.example.com/test",
            params={"key": "value"},
        )
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_post_success(base_client, mock_http_client):
    """Test successful POST request."""
    mock_response = {"data": "test"}
    mock_http_client.post.return_value = mock_response
    
    response = base_client.post(
        url="https://api.example.com/test",
        json={"key": "value"},
        headers={"Content-Type": "application/json"},
    )
    
    assert response == mock_response
    mock_http_client.post.assert_called_with(
        url="https://api.example.com/test",
        json={"key": "value"},
        headers={"Content-Type": "application/json"},
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
            json={"key": "value"},
        )
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_put_success(base_client, mock_http_client):
    """Test successful PUT request."""
    mock_response = {"data": "test"}
    mock_http_client.put.return_value = mock_response
    
    response = base_client.put(
        url="https://api.example.com/test",
        json={"key": "value"},
        headers={"Content-Type": "application/json"},
    )
    
    assert response == mock_response
    mock_http_client.put.assert_called_with(
        url="https://api.example.com/test",
        json={"key": "value"},
        headers={"Content-Type": "application/json"},
    )


def test_put_error(base_client, mock_http_client):
    """Test PUT request error handling."""
    mock_http_client.put.side_effect = WebClientError(
        "API error",
        status_code=500,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        base_client.put(
            url="https://api.example.com/test",
            json={"key": "value"},
        )
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_delete_success(base_client, mock_http_client):
    """Test successful DELETE request."""
    mock_response = {"data": "test"}
    mock_http_client.delete.return_value = mock_response
    
    response = base_client.delete(
        url="https://api.example.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
    )
    
    assert response == mock_response
    mock_http_client.delete.assert_called_with(
        url="https://api.example.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
    )


def test_delete_error(base_client, mock_http_client):
    """Test DELETE request error handling."""
    mock_http_client.delete.side_effect = WebClientError(
        "API error",
        status_code=500,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        base_client.delete(
            url="https://api.example.com/test",
            params={"key": "value"},
        )
    
    assert exc_info.value.status_code == 500
    assert "API error" in str(exc_info.value)


def test_log_request(base_client, mock_logger):
    """Test request logging."""
    base_client._log_request(
        method="GET",
        url="https://api.example.com/test",
        params={"key": "value"},
        headers={"Accept": "application/json"},
    )
    
    mock_logger.debug.assert_called_with(
        "Making GET request to https://api.example.com/test",
        extra={
            "params": {"key": "value"},
            "headers": {"Accept": "application/json"},
        },
    )


def test_log_response(base_client, mock_logger):
    """Test response logging."""
    mock_response = {"data": "test"}
    
    base_client._log_response(
        method="GET",
        url="https://api.example.com/test",
        response=mock_response,
    )
    
    mock_logger.debug.assert_called_with(
        "Received response from GET https://api.example.com/test",
        extra={"response": mock_response},
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
        "Error making GET request to https://api.example.com/test: API error",
        extra={"status_code": 500},
    )


def test_validate_response_success(base_client):
    """Test successful response validation."""
    valid_response = {
        "data": {
            "field1": "value1",
            "field2": "value2",
        },
        "metadata": {
            "query_time": 0.15,
            "api_version": "2023.1",
        },
    }
    
    # Should not raise any errors
    base_client._validate_response(
        response=valid_response,
        required_fields=["data", "metadata"],
    )


def test_validate_response_missing_fields(base_client):
    """Test response validation with missing fields."""
    invalid_response = {
        "data": {
            "field1": "value1",
        },
    }
    
    with pytest.raises(WebClientError) as exc_info:
        base_client._validate_response(
            response=invalid_response,
            required_fields=["data", "metadata"],
        )
    
    assert "Missing required fields" in str(exc_info.value)


def test_validate_response_invalid_types(base_client):
    """Test response validation with invalid types."""
    invalid_response = {
        "data": "not an object",
        "metadata": "not an object",
    }
    
    with pytest.raises(WebClientError) as exc_info:
        base_client._validate_response(
            response=invalid_response,
            required_fields=["data", "metadata"],
            field_types={
                "data": dict,
                "metadata": dict,
            },
        )
    
    assert "Invalid field types" in str(exc_info.value)


def test_validate_response_nested_fields(base_client):
    """Test response validation with nested fields."""
    valid_response = {
        "data": {
            "nested": {
                "field1": "value1",
                "field2": "value2",
            },
        },
    }
    
    # Should not raise any errors
    base_client._validate_response(
        response=valid_response,
        required_fields=["data.nested.field1", "data.nested.field2"],
    )


def test_validate_response_nested_missing(base_client):
    """Test response validation with missing nested fields."""
    invalid_response = {
        "data": {
            "nested": {
                "field1": "value1",
            },
        },
    }
    
    with pytest.raises(WebClientError) as exc_info:
        base_client._validate_response(
            response=invalid_response,
            required_fields=["data.nested.field1", "data.nested.field2"],
        )
    
    assert "Missing required fields" in str(exc_info.value)


def test_validate_response_nested_types(base_client):
    """Test response validation with nested field types."""
    invalid_response = {
        "data": {
            "nested": {
                "field1": 123,  # Should be string
                "field2": "value2",
            },
        },
    }
    
    with pytest.raises(WebClientError) as exc_info:
        base_client._validate_response(
            response=invalid_response,
            required_fields=["data.nested.field1", "data.nested.field2"],
            field_types={
                "data.nested.field1": str,
                "data.nested.field2": str,
            },
        )
    
    assert "Invalid field types" in str(exc_info.value)
