"""Tests for enrichment validation module."""

import pytest
from binding_data_processor.models.compound.enrichment import validation


@pytest.fixture
def mock_validator():
    """Create a mock enrichment validator for testing."""
    return validation.EnrichmentValidator()


def test_validator_initialization(mock_validator):
    """Test validator initialization."""
    assert isinstance(mock_validator, validation.EnrichmentValidator)
    assert hasattr(mock_validator, "validate")


def test_source_validation(mock_validator):
    """Test data source validation."""
    with pytest.raises(NotImplementedError):
        mock_validator.validate_source("pubchem")


def test_data_validation(mock_validator):
    """Test enrichment data validation."""
    with pytest.raises(NotImplementedError):
        mock_validator.validate_data({})


def test_schema_validation(mock_validator):
    """Test schema validation."""
    with pytest.raises(NotImplementedError):
        mock_validator.validate_schema({})


def test_field_validation(mock_validator):
    """Test field validation."""
    with pytest.raises(NotImplementedError):
        mock_validator.validate_fields([])


def test_source_compatibility(mock_validator):
    """Test source compatibility validation."""
    with pytest.raises(NotImplementedError):
        mock_validator.check_source_compatibility("pubchem", "chembl")


def test_data_consistency(mock_validator):
    """Test data consistency validation."""
    with pytest.raises(NotImplementedError):
        mock_validator.check_data_consistency({})


def test_enrichment_rules(mock_validator):
    """Test enrichment rules validation."""
    with pytest.raises(NotImplementedError):
        mock_validator.validate_enrichment_rules([])


def test_validation_report(mock_validator):
    """Test validation report generation."""
    with pytest.raises(NotImplementedError):
        mock_validator.generate_report({})
