"""Tests for data validation functionality."""

import pytest

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..validation import (
    DataValidator,
    ValidationRule,
    ValidationResult,
)


@pytest.fixture
def test_compound():
    """Create a test compound fixture."""
    compound = PsychoactiveCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    compound.psychoactive_class = PsychoactiveClass.STIMULANT
    compound.add_receptor_binding(
        "A2A",
        affinity=0.8,
        confidence=0.95,
        activity="antagonist"
    )
    compound.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
    }
    return compound


@pytest.fixture
def validator():
    """Create a test data validator fixture."""
    return DataValidator()


class TestDataValidator:
    """Tests for DataValidator class."""

    def test_initialization(self, validator):
        """Test initialization of DataValidator."""
        assert validator.rules is not None
        assert validator.stats == {}

    def test_basic_validation(self, test_compound, validator):
        """Test basic validation rules."""
        # Validate compound
        result = validator.validate_compound(test_compound)
        
        # Check validation result
        assert result.is_valid
        assert len(result.errors) == 0
        assert len(result.warnings) == 0

    def test_required_fields(self, test_compound, validator):
        """Test validation of required fields."""
        # Test missing name
        invalid_compound = PsychoactiveCompound(
            name="",  # Empty name
            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            cas_number="58-08-2",
        )
        
        result = validator.validate_compound(invalid_compound)
        assert not result.is_valid
        assert any("name" in error.message.lower() for error in result.errors)

    def test_smiles_validation(self, test_compound, validator):
        """Test SMILES structure validation."""
        # Test invalid SMILES
        invalid_compound = PsychoactiveCompound(
            name="Invalid",
            smiles="Invalid SMILES",
            cas_number="58-08-2",
        )
        
        result = validator.validate_compound(invalid_compound)
        assert not result.is_valid
        assert any("smiles" in error.message.lower() for error in result.errors)

    def test_cas_number_validation(self, test_compound, validator):
        """Test CAS number validation."""
        # Test invalid CAS format
        invalid_compound = PsychoactiveCompound(
            name="Caffeine",
            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            cas_number="invalid-cas",
        )
        
        result = validator.validate_compound(invalid_compound)
        assert not result.is_valid
        assert any("cas" in error.message.lower() for error in result.errors)

    def test_receptor_validation(self, test_compound, validator):
        """Test receptor binding data validation."""
        # Test invalid affinity value
        test_compound.add_receptor_binding(
            "5-HT2A",
            affinity=-1.0,  # Invalid negative affinity
            confidence=0.95,
            activity="agonist"
        )
        
        result = validator.validate_compound(test_compound)
        assert not result.is_valid
        assert any("affinity" in error.message.lower() for error in result.errors)

    def test_safety_validation(self, test_compound, validator):
        """Test safety data validation."""
        # Test invalid risk level
        test_compound.safety_alerts["test"] = "INVALID_LEVEL"  # Invalid enum value
        
        result = validator.validate_compound(test_compound)
        assert not result.is_valid
        assert any("risk level" in error.message.lower() for error in result.errors)

    def test_custom_validation_rules(self, test_compound, validator):
        """Test adding and using custom validation rules."""
        # Create custom rule
        custom_rule = ValidationRule(
            name="custom_rule",
            validate=lambda c: ValidationResult(
                is_valid=len(c.name) <= 50,
                errors=["Name too long"] if len(c.name) > 50 else []
            )
        )
        
        # Add custom rule
        validator.add_rule(custom_rule)
        
        # Test valid case
        result = validator.validate_compound(test_compound)
        assert result.is_valid
        
        # Test invalid case
        long_name_compound = PsychoactiveCompound(
            name="A" * 51,  # Name too long
            smiles=test_compound.smiles,
            cas_number=test_compound.cas_number,
        )
        result = validator.validate_compound(long_name_compound)
        assert not result.is_valid
        assert any("name too long" in error.message.lower() for error in result.errors)

    def test_warning_generation(self, test_compound, validator):
        """Test generation of validation warnings."""
        # Add warning rule for high risk compounds
        warning_rule = ValidationRule(
            name="high_risk_warning",
            validate=lambda c: ValidationResult(
                is_valid=True,
                warnings=["High risk compound"] if any(
                    risk == RiskLevel.HIGH for risk in c.safety_alerts.values()
                ) else []
            )
        )
        
        validator.add_rule(warning_rule)
        
        # Validate compound with warnings
        result = validator.validate_compound(test_compound)
        assert result.is_valid  # Warnings don't affect validity
        assert len(result.warnings) > 0
        assert any("high risk" in warning.lower() for warning in result.warnings)

    def test_batch_validation(self, test_compound, validator):
        """Test batch validation of multiple compounds."""
        compounds = [
            test_compound,
            PsychoactiveCompound(  # Invalid compound
                name="",
                smiles="Invalid",
                cas_number="invalid",
            )
        ]
        
        # Validate batch
        results = validator.validate_compounds(compounds)
        
        # Check results
        assert len(results) == 2
        assert results[0].is_valid  # First compound valid
        assert not results[1].is_valid  # Second compound invalid
        assert validator.stats["total_validated"] == 2
        assert validator.stats["invalid_compounds"] == 1

    def test_error_handling(self, test_compound, validator):
        """Test error handling during validation."""
        # Create failing rule
        failing_rule = ValidationRule(
            name="failing_rule",
            validate=lambda c: 1/0  # Force exception
        )
        
        validator.add_rule(failing_rule)
        
        # Validate with failing rule
        result = validator.validate_compound(test_compound)
        
        # Check error handling
        assert not result.is_valid
        assert any("validation error" in error.message.lower() for error in result.errors)
        assert "rule_errors" in validator.stats
        assert validator.stats["rule_errors"] > 0

    def test_validation_performance(self, test_compound, validator):
        """Test validation performance monitoring."""
        # Create slow rule
        slow_rule = ValidationRule(
            name="slow_rule",
            validate=lambda c: ValidationResult(is_valid=True)
        )
        
        validator.add_rule(slow_rule)
        
        # Validate with performance monitoring
        validator.validate_compound(test_compound)
        
        # Check performance stats
        assert "validation_time" in validator.stats
        assert isinstance(validator.stats["validation_time"], float)
        assert validator.stats["validation_time"] >= 0


if __name__ == "__main__":
    pytest.main([__file__])
