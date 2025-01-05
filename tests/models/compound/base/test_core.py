"""Tests for core compound model functionality."""

import pytest
from typing import Dict, Any

from binding_data_processor.models.compound.base.core import CompoundBase
from binding_data_processor.models.compound.base.mixins import ValidationMixin, SerializationMixin
from binding_data_processor.models.compound.base.types import (
    CompoundIdentifiers,
    CompoundProperties,
    CompoundMetadata,
)

@pytest.fixture
def valid_identifiers() -> CompoundIdentifiers:
    """Create valid compound identifiers."""
    return {
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "cas_number": "58-08-2",
        "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
        "inchi_key": "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    }

@pytest.fixture
def valid_properties() -> CompoundProperties:
    """Create valid compound properties."""
    return {
        "molecular_weight": 194.19,
        "molecular_formula": "C8H10N4O2",
        "logp": -0.07,
        "hbd": 0,
        "hba": 6,
        "tpsa": 58.4,
        "rotatable_bonds": 0
    }

@pytest.fixture
def valid_metadata() -> CompoundMetadata:
    """Create valid compound metadata."""
    return {
        "source": "BindingDB",
        "source_id": "50105844",
        "last_updated": "2024-01-01",
        "data_quality": 0.95,
        "validation_status": "validated"
    }

@pytest.fixture
def valid_compound(
    valid_identifiers: CompoundIdentifiers,
    valid_properties: CompoundProperties,
    valid_metadata: CompoundMetadata
) -> CompoundBase:
    """Create valid compound instance."""
    compound = CompoundBase()
    compound.identifiers = valid_identifiers
    compound.properties = valid_properties
    compound.metadata = valid_metadata
    return compound

class TestCompoundBase:
    """Tests for CompoundBase class."""

    def test_initialization(self):
        """Test compound initialization."""
        compound = CompoundBase()
        assert isinstance(compound.identifiers, dict)
        assert isinstance(compound.properties, dict)
        assert isinstance(compound.metadata, dict)

    def test_valid_compound(self, valid_compound: CompoundBase):
        """Test valid compound data."""
        assert valid_compound.identifiers["name"] == "Caffeine"
        assert valid_compound.properties["molecular_weight"] == 194.19
        assert valid_compound.metadata["source"] == "BindingDB"

    def test_required_identifiers(self, valid_compound: CompoundBase):
        """Test required identifier validation."""
        # Remove required identifier
        del valid_compound.identifiers["smiles"]
        with pytest.raises(ValueError):
            valid_compound.validate()

    def test_property_types(self, valid_compound: CompoundBase):
        """Test property type validation."""
        # Set invalid property type
        valid_compound.properties["molecular_weight"] = "invalid"
        with pytest.raises(TypeError):
            valid_compound.validate()

    def test_metadata_validation(self, valid_compound: CompoundBase):
        """Test metadata validation."""
        # Set invalid metadata value
        valid_compound.metadata["data_quality"] = 2.0  # Should be 0-1
        with pytest.raises(ValueError):
            valid_compound.validate()

class TestValidationMixin:
    """Tests for ValidationMixin."""

    def test_validation_success(self, valid_compound: CompoundBase):
        """Test successful validation."""
        assert isinstance(valid_compound, ValidationMixin)
        assert valid_compound.validate() is True

    def test_validation_error_handling(self, valid_compound: CompoundBase):
        """Test validation error handling."""
        # Create invalid state
        valid_compound.identifiers = None
        with pytest.raises(ValueError) as exc_info:
            valid_compound.validate()
        assert "identifiers" in str(exc_info.value)

class TestSerializationMixin:
    """Tests for SerializationMixin."""

    def test_to_dict(self, valid_compound: CompoundBase):
        """Test dictionary serialization."""
        assert isinstance(valid_compound, SerializationMixin)
        data = valid_compound.to_dict()
        assert isinstance(data, dict)
        assert "identifiers" in data
        assert "properties" in data
        assert "metadata" in data

    def test_from_dict(self, valid_compound: CompoundBase):
        """Test dictionary deserialization."""
        data = valid_compound.to_dict()
        new_compound = CompoundBase.from_dict(data)
        assert new_compound.identifiers == valid_compound.identifiers
        assert new_compound.properties == valid_compound.properties
        assert new_compound.metadata == valid_compound.metadata

    def test_json_serialization(self, valid_compound: CompoundBase):
        """Test JSON serialization."""
        json_str = valid_compound.to_json()
        assert isinstance(json_str, str)
        new_compound = CompoundBase.from_json(json_str)
        assert new_compound.identifiers == valid_compound.identifiers

def test_compound_equality(valid_compound: CompoundBase):
    """Test compound equality comparison."""
    compound1 = valid_compound
    compound2 = CompoundBase.from_dict(valid_compound.to_dict())
    assert compound1 == compound2

    # Change a value
    compound2.identifiers["name"] = "Different"
    assert compound1 != compound2

def test_compound_hash(valid_compound: CompoundBase):
    """Test compound hash function."""
    compound1 = valid_compound
    compound2 = CompoundBase.from_dict(valid_compound.to_dict())
    
    # Same data should have same hash
    assert hash(compound1) == hash(compound2)
    
    # Different data should have different hash
    compound2.identifiers["name"] = "Different"
    assert hash(compound1) != hash(compound2)

def test_compound_copy(valid_compound: CompoundBase):
    """Test compound deep copy."""
    compound_copy = valid_compound.copy()
    assert compound_copy == valid_compound
    assert compound_copy is not valid_compound
    
    # Modify copy shouldn't affect original
    compound_copy.identifiers["name"] = "Different"
    assert compound_copy != valid_compound
    assert valid_compound.identifiers["name"] == "Caffeine"
