"""Tests for compound analysis error handling."""

import pytest
from typing import Dict, Any

from ..errors import (
    AnalysisError,
    ValidationError,
    DataError,
    ConfigurationError,
    RecoveryError,
    format_error,
    handle_error,
    ErrorContext,
)


def test_analysis_error():
    """Test AnalysisError class."""
    # Test basic error
    error = AnalysisError("Test error")
    assert str(error) == "Test error"
    assert error.message == "Test error"
    assert error.context == {}

    # Test with context
    context = {"value": 1.2, "unit": "nM"}
    error = AnalysisError("Test error", context=context)
    assert error.context == context
    assert "value: 1.2" in str(error)
    assert "unit: nM" in str(error)

    # Test with cause
    cause = ValueError("Original error")
    error = AnalysisError("Test error", cause=cause)
    assert error.__cause__ == cause
    assert "Caused by ValueError: Original error" in str(error)


def test_validation_error():
    """Test ValidationError class."""
    # Test basic validation error
    error = ValidationError("Invalid value")
    assert isinstance(error, AnalysisError)
    assert str(error) == "Invalid value"

    # Test with field information
    error = ValidationError("Invalid value", field="concentration")
    assert "field 'concentration'" in str(error)
    assert error.field == "concentration"

    # Test with value information
    error = ValidationError("Invalid value", field="concentration", value="invalid")
    assert "field 'concentration'" in str(error)
    assert "value 'invalid'" in str(error)
    assert error.value == "invalid"


def test_data_error():
    """Test DataError class."""
    # Test basic data error
    error = DataError("Missing data")
    assert isinstance(error, AnalysisError)
    assert str(error) == "Missing data"

    # Test with data source
    error = DataError("Missing data", source="BindingDB")
    assert "source 'BindingDB'" in str(error)
    assert error.source == "BindingDB"

    # Test with record information
    error = DataError("Missing data", source="BindingDB", record_id="12345")
    assert "source 'BindingDB'" in str(error)
    assert "record '12345'" in str(error)
    assert error.record_id == "12345"


def test_configuration_error():
    """Test ConfigurationError class."""
    # Test basic configuration error
    error = ConfigurationError("Invalid configuration")
    assert isinstance(error, AnalysisError)
    assert str(error) == "Invalid configuration"

    # Test with parameter information
    error = ConfigurationError("Invalid configuration", parameter="cache_size")
    assert "parameter 'cache_size'" in str(error)
    assert error.parameter == "cache_size"

    # Test with value information
    error = ConfigurationError(
        "Invalid configuration",
        parameter="cache_size",
        value=-1
    )
    assert "parameter 'cache_size'" in str(error)
    assert "value '-1'" in str(error)
    assert error.value == -1


def test_recovery_error():
    """Test RecoveryError class."""
    # Test basic recovery error
    error = RecoveryError("Recovery failed")
    assert isinstance(error, AnalysisError)
    assert str(error) == "Recovery failed"

    # Test with operation information
    error = RecoveryError("Recovery failed", operation="cache_rebuild")
    assert "operation 'cache_rebuild'" in str(error)
    assert error.operation == "cache_rebuild"

    # Test with attempt information
    error = RecoveryError(
        "Recovery failed",
        operation="cache_rebuild",
        attempt=3
    )
    assert "operation 'cache_rebuild'" in str(error)
    assert "attempt 3" in str(error)
    assert error.attempt == 3


def test_error_formatting():
    """Test error message formatting."""
    # Test basic formatting
    msg = format_error("Test error")
    assert isinstance(msg, str)
    assert msg.startswith("[ERROR]")
    assert "Test error" in msg

    # Test with context
    context = {"value": 1.2, "unit": "nM"}
    msg = format_error("Test error", context)
    assert "value: 1.2" in msg
    assert "unit: nM" in msg

    # Test with exception
    try:
        raise ValueError("Original error")
    except ValueError as e:
        msg = format_error("Test error", exception=e)
        assert "Test error" in msg
        assert "Original error" in msg
        assert "ValueError" in msg


def test_error_handling():
    """Test error handling functionality."""
    def test_operation(context: Dict[str, Any]) -> None:
        if context.get("should_fail"):
            raise ValueError("Operation failed")
        return None

    # Test successful operation
    context = {"should_fail": False}
    result = handle_error(
        test_operation,
        context,
        error_msg="Test operation failed"
    )
    assert result is None

    # Test failed operation
    context = {"should_fail": True}
    with pytest.raises(AnalysisError) as exc:
        handle_error(
            test_operation,
            context,
            error_msg="Test operation failed"
        )
    assert "Test operation failed" in str(exc.value)
    assert "Operation failed" in str(exc.value)


def test_error_context():
    """Test error context management."""
    # Test context creation
    context = ErrorContext()
    assert context.data == {}

    # Test context update
    context.update({"value": 1.2})
    assert context.data["value"] == 1.2

    # Test context manager
    with ErrorContext() as ctx:
        ctx.update({"test": True})
        assert ctx.data["test"] is True

    # Test nested context
    with ErrorContext() as outer:
        outer.update({"outer": True})
        with ErrorContext(parent=outer) as inner:
            inner.update({"inner": True})
            assert "outer" in inner.data
            assert "inner" in inner.data


def test_error_recovery():
    """Test error recovery functionality."""
    def test_operation(should_fail: bool = False) -> str:
        if should_fail:
            raise ValueError("Operation failed")
        return "success"

    # Test successful recovery
    result = handle_error(
        test_operation,
        should_fail=False,
        error_msg="Operation failed",
        recovery_func=lambda: "recovered"
    )
    assert result == "success"

    # Test failed operation with recovery
    result = handle_error(
        test_operation,
        should_fail=True,
        error_msg="Operation failed",
        recovery_func=lambda: "recovered"
    )
    assert result == "recovered"

    # Test failed operation and recovery
    with pytest.raises(RecoveryError):
        handle_error(
            test_operation,
            should_fail=True,
            error_msg="Operation failed",
            recovery_func=lambda: exec('raise ValueError("Recovery failed")')
        )


def test_error_chaining():
    """Test error chaining functionality."""
    try:
        try:
            raise ValueError("Original error")
        except ValueError as e:
            raise DataError("Data processing failed", source="test") from e
    except DataError as e:
        assert isinstance(e.__cause__, ValueError)
        assert "Original error" in str(e.__cause__)
        assert "Data processing failed" in str(e)
        assert "source 'test'" in str(e)
