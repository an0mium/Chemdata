"""Tests for compound analysis utility functions."""

import pytest
from datetime import datetime
import json

from ..utils import (
    format_analysis_result,
    format_error_message,
    cache_result,
    convert_units,
    merge_data,
    validate_data_format,
)


def test_format_analysis_result():
    """Test analysis result formatting."""
    # Test basic result
    result = {
        "value": 1.2,
        "confidence": 0.9,
    }
    formatted = format_analysis_result(result)
    assert isinstance(formatted, dict)
    assert "timestamp" in formatted
    assert "version" in formatted
    assert formatted["data"] == result

    # Test nested result
    result = {
        "binding": {
            "value": 1.2,
            "confidence": 0.9,
        },
        "activity": {
            "value": 0.8,
            "confidence": 0.85,
        },
    }
    formatted = format_analysis_result(result)
    assert all(k in formatted["data"] for k in ["binding", "activity"])
    assert all(k in formatted for k in ["timestamp", "version"])

    # Test with metadata
    metadata = {
        "source": "test",
        "method": "test_method",
    }
    formatted = format_analysis_result(result, metadata=metadata)
    assert formatted["metadata"] == metadata


def test_format_error_message():
    """Test error message formatting."""
    # Test basic error
    msg = format_error_message("Test error")
    assert isinstance(msg, str)
    assert "Test error" in msg
    assert "[ERROR]" in msg

    # Test with context
    msg = format_error_message("Test error", context={"value": 1.2})
    assert "Test error" in msg
    assert "value: 1.2" in msg

    # Test with exception
    try:
        raise ValueError("Test exception")
    except ValueError as e:
        msg = format_error_message("Test error", exception=e)
        assert "Test error" in msg
        assert "Test exception" in msg
        assert "ValueError" in msg


def test_cache_result():
    """Test result caching functionality."""
    # Test basic caching
    @cache_result
    def test_func():
        return {"value": 1.2}

    result1 = test_func()
    result2 = test_func()
    assert result1 is result2  # Same object in memory

    # Test cache invalidation
    @cache_result(max_age=0)  # Immediate invalidation
    def test_func2():
        return {"value": datetime.now()}

    result1 = test_func2()
    result2 = test_func2()
    assert result1 is not result2  # Different objects

    # Test conditional caching
    @cache_result(condition=lambda x: x > 0)
    def test_func3(x):
        return {"value": x}

    result1 = test_func3(1)  # Should cache
    result2 = test_func3(1)  # Should use cache
    result3 = test_func3(-1)  # Should not cache
    result4 = test_func3(-1)  # Should not use cache
    assert result1 is result2
    assert result3 is not result4


def test_convert_units():
    """Test unit conversion functionality."""
    # Test concentration conversions
    assert convert_units(1000, "nM", "μM") == 1.0
    assert convert_units(1, "μM", "nM") == 1000.0
    assert convert_units(1000000, "nM", "mM") == 1.0

    # Test invalid conversions
    with pytest.raises(ValueError):
        convert_units(1, "invalid", "nM")
    with pytest.raises(ValueError):
        convert_units(1, "nM", "invalid")

    # Test same unit
    assert convert_units(1.2, "nM", "nM") == 1.2


def test_merge_data():
    """Test data merging functionality."""
    # Test basic merge
    data1 = {"a": 1, "b": 2}
    data2 = {"c": 3, "d": 4}
    merged = merge_data(data1, data2)
    assert all(k in merged for k in ["a", "b", "c", "d"])

    # Test overlapping keys
    data1 = {"a": 1, "b": 2}
    data2 = {"b": 3, "c": 4}
    merged = merge_data(data1, data2, strategy="keep_first")
    assert merged["b"] == 2
    merged = merge_data(data1, data2, strategy="keep_second")
    assert merged["b"] == 3

    # Test nested merge
    data1 = {"a": {"x": 1, "y": 2}}
    data2 = {"a": {"y": 3, "z": 4}}
    merged = merge_data(data1, data2, strategy="deep")
    assert merged["a"]["x"] == 1
    assert merged["a"]["y"] == 3
    assert merged["a"]["z"] == 4


def test_validate_data_format():
    """Test data format validation."""
    # Test valid formats
    valid_data = {
        "value": 1.2,
        "confidence": 0.9,
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "version": "1.0",
        },
    }
    assert validate_data_format(valid_data, ["value", "confidence"])
    assert validate_data_format(valid_data, ["value", "confidence", "metadata"])

    # Test missing required fields
    invalid_data = {
        "value": 1.2,
        # Missing confidence
    }
    with pytest.raises(ValueError):
        validate_data_format(invalid_data, ["value", "confidence"])

    # Test invalid types
    invalid_data = {
        "value": "not a number",
        "confidence": 0.9,
    }
    with pytest.raises(TypeError):
        validate_data_format(invalid_data, ["value", "confidence"])

    # Test nested validation
    nested_data = {
        "binding": {
            "value": 1.2,
            "confidence": 0.9,
        },
        "activity": {
            "value": 0.8,
            "confidence": 0.85,
        },
    }
    assert validate_data_format(
        nested_data,
        ["binding", "activity"],
        nested_fields={"binding": ["value", "confidence"]}
    )


def test_json_serialization():
    """Test JSON serialization of analysis results."""
    # Test basic serialization
    result = {
        "value": 1.2,
        "confidence": 0.9,
        "timestamp": datetime.now(),
    }
    formatted = format_analysis_result(result)
    
    # Should serialize without error
    json_str = json.dumps(formatted)
    assert isinstance(json_str, str)
    
    # Should deserialize correctly
    deserialized = json.loads(json_str)
    assert deserialized["data"]["value"] == 1.2
    assert deserialized["data"]["confidence"] == 0.9
    assert "timestamp" in deserialized


def test_error_handling():
    """Test error handling utilities."""
    # Test error collection
    errors = []
    
    def collect_error(msg):
        errors.append(msg)
    
    # Test validation with error collection
    invalid_data = {
        "value": "not a number",
        "confidence": 1.5,  # > 1
    }
    
    try:
        validate_data_format(
            invalid_data,
            ["value", "confidence"],
            on_error=collect_error
        )
    except (TypeError, ValueError):
        pass
    
    assert len(errors) == 2  # Should collect both errors
    assert any("not a number" in e for e in errors)
    assert any("confidence" in e for e in errors)


def test_unit_registry():
    """Test unit registry functionality."""
    from ..utils import UNIT_REGISTRY

    # Test basic units
    assert "nM" in UNIT_REGISTRY
    assert "μM" in UNIT_REGISTRY
    assert "mM" in UNIT_REGISTRY

    # Test conversions
    assert UNIT_REGISTRY["nM"]["to_base"] == 1
    assert UNIT_REGISTRY["μM"]["to_base"] == 1000
    assert UNIT_REGISTRY["mM"]["to_base"] == 1000000

    # Test unit validation
    assert "invalid" not in UNIT_REGISTRY
