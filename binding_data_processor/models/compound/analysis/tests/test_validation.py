"""Tests for compound analysis validation functionality."""

import pytest
from rdkit import Chem

from ..base import CompoundAnalysis, AnalysisError
from ..validation import (
    validate_structure,
    validate_target_data,
    validate_activity_data,
    validate_property_data,
    validate_reference_data,
)


def test_structure_validation():
    """Test chemical structure validation."""
    # Valid SMILES
    assert validate_structure("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")  # Caffeine
    assert validate_structure("CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1")  # LSD

    # Invalid SMILES
    with pytest.raises(AnalysisError) as exc:
        validate_structure("invalid_smiles")
    assert "Invalid SMILES" in str(exc.value)

    # Empty SMILES
    with pytest.raises(AnalysisError) as exc:
        validate_structure("")
    assert "Empty SMILES" in str(exc.value)

    # None SMILES
    with pytest.raises(AnalysisError) as exc:
        validate_structure(None)
    assert "Missing SMILES" in str(exc.value)

    # Check molecule properties
    mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert mol.GetNumAtoms() > 0
    assert mol.GetNumBonds() > 0


def test_target_data_validation(lsd_data):
    """Test binding target data validation."""
    # Valid target data
    valid_targets = lsd_data["targets"]
    assert validate_target_data(valid_targets)

    # Missing required fields
    invalid_targets = [
        {
            "common_name": "5-HT2A",
            # Missing affinity_value
            "affinity_type": "Ki",
        }
    ]
    with pytest.raises(AnalysisError) as exc:
        validate_target_data(invalid_targets)
    assert "Missing required field" in str(exc.value)

    # Invalid affinity value
    invalid_targets = [
        {
            "common_name": "5-HT2A",
            "affinity_value": "invalid",  # Should be numeric
            "affinity_type": "Ki",
            "affinity_unit": "nM",
        }
    ]
    with pytest.raises(AnalysisError) as exc:
        validate_target_data(invalid_targets)
    assert "Invalid affinity value" in str(exc.value)

    # Invalid confidence value
    invalid_targets = [
        {
            "common_name": "5-HT2A",
            "affinity_value": 1.2,
            "affinity_type": "Ki",
            "affinity_unit": "nM",
            "confidence": 1.5,  # Should be between 0 and 1
        }
    ]
    with pytest.raises(AnalysisError) as exc:
        validate_target_data(invalid_targets)
    assert "Invalid confidence value" in str(exc.value)


def test_activity_data_validation(lsd_data):
    """Test activity data validation."""
    # Valid activity data
    valid_activity = {
        "primary_activity": lsd_data["primary_activity"],
        "mechanism_of_action": lsd_data["mechanism_of_action"],
        "effect_profile": lsd_data["effect_profile"],
    }
    assert validate_activity_data(valid_activity)

    # Invalid primary activity
    invalid_activity = {
        "primary_activity": "invalid",  # Should be numeric
        "mechanism_of_action": "test mechanism",
    }
    with pytest.raises(AnalysisError) as exc:
        validate_activity_data(invalid_activity)
    assert "Invalid primary activity" in str(exc.value)

    # Invalid effect profile
    invalid_activity = {
        "primary_activity": 1.2,
        "mechanism_of_action": "test mechanism",
        "effect_profile": {
            "effect1": (0.8, 1.2),  # Second value > 1
        },
    }
    with pytest.raises(AnalysisError) as exc:
        validate_activity_data(invalid_activity)
    assert "Invalid effect score" in str(exc.value)


def test_property_data_validation(drug_like_properties):
    """Test property data validation."""
    # Valid property data
    assert validate_property_data(drug_like_properties)

    # Invalid molecular weight
    invalid_properties = drug_like_properties.copy()
    invalid_properties["molecular_weight"] = "invalid"  # Should be numeric
    with pytest.raises(AnalysisError) as exc:
        validate_property_data(invalid_properties)
    assert "Invalid molecular weight" in str(exc.value)

    # Invalid LogP
    invalid_properties = drug_like_properties.copy()
    invalid_properties["logp"] = None  # Should be numeric
    with pytest.raises(AnalysisError) as exc:
        validate_property_data(invalid_properties)
    assert "Invalid LogP" in str(exc.value)

    # Invalid counts
    invalid_properties = drug_like_properties.copy()
    invalid_properties["hbd"] = -1  # Should be non-negative
    with pytest.raises(AnalysisError) as exc:
        validate_property_data(invalid_properties)
    assert "Invalid HBD count" in str(exc.value)


def test_reference_data_validation():
    """Test reference compound data validation."""
    # Valid reference data
    valid_references = [
        {
            "name": "Compound 1",
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "primary_activity": 1.2,
        },
        {
            "name": "Compound 2",
            "smiles": "CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1",
            "primary_activity": 0.8,
        },
    ]
    assert validate_reference_data(valid_references)

    # Invalid SMILES
    invalid_references = [
        {
            "name": "Compound 1",
            "smiles": "invalid_smiles",
            "primary_activity": 1.2,
        }
    ]
    with pytest.raises(AnalysisError) as exc:
        validate_reference_data(invalid_references)
    assert "Invalid SMILES" in str(exc.value)

    # Missing required fields
    invalid_references = [
        {
            "name": "Compound 1",
            # Missing SMILES
            "primary_activity": 1.2,
        }
    ]
    with pytest.raises(AnalysisError) as exc:
        validate_reference_data(invalid_references)
    assert "Missing required field" in str(exc.value)

    # Invalid activity value
    invalid_references = [
        {
            "name": "Compound 1",
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "primary_activity": "invalid",  # Should be numeric
        }
    ]
    with pytest.raises(AnalysisError) as exc:
        validate_reference_data(invalid_references)
    assert "Invalid activity value" in str(exc.value)


def test_compound_validation(lsd_data, drug_like_properties):
    """Test full compound validation."""
    compound = CompoundAnalysis()
    
    # Set valid data
    compound.smiles = lsd_data["smiles"]
    compound.targets = lsd_data["targets"]
    compound.primary_activity = lsd_data["primary_activity"]
    compound.mechanism_of_action = lsd_data["mechanism_of_action"]
    compound.effect_profile = lsd_data["effect_profile"]
    
    # Set property data
    for key, value in drug_like_properties.items():
        setattr(compound, key, value)
    
    # Validate all - should not raise
    compound.validate_all()

    # Test individual validations
    compound.validate_structure()
    compound.validate_targets()
    compound.validate_activity()
    compound.validate_properties()

    # Test validation order
    validation_order = [
        "validate_structure",
        "validate_targets",
        "validate_activity",
        "validate_properties",
    ]
    for method in validation_order:
        assert hasattr(compound, method)
        assert callable(getattr(compound, method))


def test_validation_error_messages():
    """Test validation error messages."""
    compound = CompoundAnalysis()

    # Structure validation
    with pytest.raises(AnalysisError) as exc:
        compound.validate_structure()
    assert "Missing SMILES" in str(exc.value)

    # Target validation
    compound.smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    with pytest.raises(AnalysisError) as exc:
        compound.validate_targets()
    assert "Missing target data" in str(exc.value)

    # Activity validation
    compound.targets = []
    with pytest.raises(AnalysisError) as exc:
        compound.validate_activity()
    assert "Missing activity data" in str(exc.value)

    # Property validation
    compound.primary_activity = 1.2
    with pytest.raises(AnalysisError) as exc:
        compound.validate_properties()
    assert "Missing property data" in str(exc.value)


def test_validation_warnings():
    """Test validation warnings."""
    compound = CompoundAnalysis()
    compound.smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    # Test with incomplete but valid data
    compound.targets = []
    compound.primary_activity = None
    compound.validate_all(strict=False)  # Should not raise

    # Test with strict validation
    with pytest.raises(AnalysisError):
        compound.validate_all(strict=True)
