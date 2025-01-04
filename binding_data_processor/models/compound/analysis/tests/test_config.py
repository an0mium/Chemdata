"""Tests for compound analysis configuration functionality."""

import pytest
import tempfile
import json
from pathlib import Path
from typing import Dict, Any

from ..config import (
    AnalysisConfig,
    ConfigError,
    load_config,
    validate_config,
    DEFAULT_CONFIG,
)


@pytest.fixture
def temp_config_file():
    """Create temporary config file."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        yield Path(f.name)
        Path(f.name).unlink()


@pytest.fixture
def valid_config() -> Dict[str, Any]:
    """Create valid configuration dictionary."""
    return {
        "analysis": {
            "binding": {
                "enabled": True,
                "min_confidence": 0.8,
                "max_targets": 10,
            },
            "activity": {
                "enabled": True,
                "min_confidence": 0.7,
                "effect_threshold": 0.6,
            },
            "safety": {
                "enabled": True,
                "risk_threshold": 0.8,
                "alert_level": "moderate",
            },
            "properties": {
                "enabled": True,
                "include_descriptors": True,
                "include_fingerprints": True,
            },
            "sar": {
                "enabled": True,
                "similarity_threshold": 0.7,
                "max_analogs": 5,
            },
        },
        "logging": {
            "level": "INFO",
            "file": "analysis.log",
            "format": "structured",
        },
        "caching": {
            "enabled": True,
            "max_size": 1000,
            "ttl": 3600,
        },
        "export": {
            "format": "tsv",
            "compress": False,
            "include_metadata": True,
        },
    }


def test_default_config():
    """Test default configuration."""
    config = AnalysisConfig()
    
    # Check default values
    assert config.analysis.binding.enabled is True
    assert config.analysis.activity.enabled is True
    assert config.analysis.safety.enabled is True
    assert config.analysis.properties.enabled is True
    assert config.analysis.sar.enabled is True

    # Check default thresholds
    assert 0 <= config.analysis.binding.min_confidence <= 1
    assert 0 <= config.analysis.activity.min_confidence <= 1
    assert 0 <= config.analysis.safety.risk_threshold <= 1
    assert 0 <= config.analysis.sar.similarity_threshold <= 1

    # Check default limits
    assert config.analysis.binding.max_targets > 0
    assert config.analysis.sar.max_analogs > 0
    assert config.caching.max_size > 0
    assert config.caching.ttl > 0


def test_config_loading(temp_config_file, valid_config):
    """Test configuration loading."""
    # Write config to file
    temp_config_file.write_text(json.dumps(valid_config))

    # Load config
    config = load_config(temp_config_file)
    
    # Check loaded values
    assert config.analysis.binding.enabled == valid_config["analysis"]["binding"]["enabled"]
    assert config.analysis.binding.min_confidence == valid_config["analysis"]["binding"]["min_confidence"]
    assert config.analysis.binding.max_targets == valid_config["analysis"]["binding"]["max_targets"]

    # Check nested values
    assert config.logging.level == valid_config["logging"]["level"]
    assert config.logging.file == valid_config["logging"]["file"]
    assert config.logging.format == valid_config["logging"]["format"]


def test_config_validation():
    """Test configuration validation."""
    # Test valid config
    valid = {
        "analysis": {
            "binding": {
                "enabled": True,
                "min_confidence": 0.8,
            }
        }
    }
    assert validate_config(valid)

    # Test invalid confidence
    invalid = {
        "analysis": {
            "binding": {
                "enabled": True,
                "min_confidence": 1.5,  # Should be <= 1
            }
        }
    }
    with pytest.raises(ConfigError):
        validate_config(invalid)

    # Test invalid type
    invalid = {
        "analysis": {
            "binding": {
                "enabled": "true",  # Should be bool
                "min_confidence": 0.8,
            }
        }
    }
    with pytest.raises(ConfigError):
        validate_config(invalid)


def test_config_overrides(valid_config):
    """Test configuration overrides."""
    # Create base config
    config = AnalysisConfig()
    
    # Override with valid config
    config.update(valid_config)
    
    # Check overridden values
    assert config.analysis.binding.min_confidence == valid_config["analysis"]["binding"]["min_confidence"]
    assert config.analysis.activity.min_confidence == valid_config["analysis"]["activity"]["min_confidence"]
    assert config.analysis.safety.risk_threshold == valid_config["analysis"]["safety"]["risk_threshold"]

    # Check unchanged defaults
    assert config.analysis.binding.enabled == DEFAULT_CONFIG["analysis"]["binding"]["enabled"]


def test_config_serialization(valid_config):
    """Test configuration serialization."""
    # Create config
    config = AnalysisConfig()
    config.update(valid_config)
    
    # Serialize to dict
    serialized = config.to_dict()
    
    # Check serialized values
    assert serialized["analysis"]["binding"]["min_confidence"] == valid_config["analysis"]["binding"]["min_confidence"]
    assert serialized["analysis"]["activity"]["min_confidence"] == valid_config["analysis"]["activity"]["min_confidence"]
    assert serialized["analysis"]["safety"]["risk_threshold"] == valid_config["analysis"]["safety"]["risk_threshold"]

    # Serialize to JSON
    json_str = config.to_json()
    deserialized = json.loads(json_str)
    assert deserialized == serialized


def test_config_validation_rules():
    """Test configuration validation rules."""
    # Test confidence values
    with pytest.raises(ConfigError):
        validate_config({
            "analysis": {
                "binding": {"min_confidence": -0.1}
            }
        })
    with pytest.raises(ConfigError):
        validate_config({
            "analysis": {
                "binding": {"min_confidence": 1.1}
            }
        })

    # Test threshold values
    with pytest.raises(ConfigError):
        validate_config({
            "analysis": {
                "sar": {"similarity_threshold": -0.1}
            }
        })
    with pytest.raises(ConfigError):
        validate_config({
            "analysis": {
                "sar": {"similarity_threshold": 1.1}
            }
        })

    # Test count values
    with pytest.raises(ConfigError):
        validate_config({
            "analysis": {
                "binding": {"max_targets": 0}
            }
        })
    with pytest.raises(ConfigError):
        validate_config({
            "analysis": {
                "binding": {"max_targets": -1}
            }
        })


def test_config_inheritance():
    """Test configuration inheritance."""
    # Create parent config
    parent = AnalysisConfig()
    parent.analysis.binding.min_confidence = 0.8
    
    # Create child config
    child = AnalysisConfig(parent=parent)
    
    # Check inherited values
    assert child.analysis.binding.min_confidence == parent.analysis.binding.min_confidence
    
    # Override value
    child.analysis.binding.min_confidence = 0.9
    
    # Check parent unchanged
    assert parent.analysis.binding.min_confidence == 0.8


def test_config_environment_override(monkeypatch):
    """Test configuration environment variable override."""
    # Set environment variables
    monkeypatch.setenv("ANALYSIS_BINDING_MIN_CONFIDENCE", "0.9")
    monkeypatch.setenv("ANALYSIS_ACTIVITY_MIN_CONFIDENCE", "0.8")
    
    # Create config
    config = AnalysisConfig()
    config.load_environment()
    
    # Check overridden values
    assert config.analysis.binding.min_confidence == 0.9
    assert config.analysis.activity.min_confidence == 0.8


def test_config_file_override(temp_config_file, valid_config):
    """Test configuration file override."""
    # Write override config
    override_config = {
        "analysis": {
            "binding": {
                "min_confidence": 0.9
            }
        }
    }
    temp_config_file.write_text(json.dumps(override_config))
    
    # Create config with base values
    config = AnalysisConfig()
    config.update(valid_config)
    
    # Apply override
    config.load_file(temp_config_file)
    
    # Check overridden value
    assert config.analysis.binding.min_confidence == 0.9
    
    # Check other values unchanged
    assert config.analysis.activity.min_confidence == valid_config["analysis"]["activity"]["min_confidence"]
