"""Tests for schema validation."""

import pytest
from datetime import datetime
from typing import Dict, Any

from ..validation.schema import (
    DataSource,
    ValidationLevel,
    ValidationResult,
    ResearchData,
    PatentData,
    CommunityData,
    WebDataSchema,
    SchemaRegistry,
    SchemaValidator,
)


@pytest.fixture
def valid_research_data() -> Dict[str, Any]:
    """Valid research paper data."""
    return {
        "title": "Test Paper",
        "abstract": "Test abstract",
        "authors": ["Author 1", "Author 2"],
        "date": datetime.now().isoformat(),
        "doi": "10.1234/test",
        "url": "https://example.com/paper",
        "source": "PubMed",
        "confidence": 0.95,
        "keywords": ["test", "research"],
        "full_text": "Full paper text",
    }


@pytest.fixture
def valid_patent_data() -> Dict[str, Any]:
    """Valid patent data."""
    return {
        "title": "Test Patent",
        "abstract": "Test abstract",
        "inventors": ["Inventor 1", "Inventor 2"],
        "assignee": "Test Company",
        "filing_date": datetime.now().isoformat(),
        "publication_date": datetime.now().isoformat(),
        "patent_number": "US123456",
        "url": "https://example.com/patent",
        "confidence": 0.9,
        "claims": ["claim1", "claim2"],
        "chemical_structures": [
            {"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "type": "SMILES"},
        ],
    }


@pytest.fixture
def valid_community_data() -> Dict[str, Any]:
    """Valid community data."""
    return {
        "title": "Test Post",
        "content": "Test content",
        "author": "user123",
        "date": datetime.now().isoformat(),
        "url": "https://example.com/post",
        "source": "Reddit",
        "subforum": "nootropics",
        "confidence": 0.85,
        "sentiment": 0.5,
        "effects": ["focus", "memory"],
        "dosage": "10mg",
        "roa": "oral",
        "safety_notes": ["note1", "note2"],
        "mentions": [{"text": "caffeine", "context": "positive"}],
        "tags": ["nootropic", "stimulant"],
    }


@pytest.fixture
def valid_web_data(
    valid_research_data,
    valid_patent_data,
    valid_community_data,
) -> Dict[str, Any]:
    """Valid combined web data."""
    return {
        "research": [valid_research_data],
        "patents": [valid_patent_data],
        "community": [valid_community_data],
        "metadata": {"timestamp": datetime.now().isoformat()},
        "confidence": 0.9,
    }


def test_research_data_model(valid_research_data):
    """Test ResearchData model validation."""
    data = ResearchData(**valid_research_data)
    assert data.title == valid_research_data["title"]
    assert data.abstract == valid_research_data["abstract"]
    assert data.confidence == valid_research_data["confidence"]


def test_patent_data_model(valid_patent_data):
    """Test PatentData model validation."""
    data = PatentData(**valid_patent_data)
    assert data.title == valid_patent_data["title"]
    assert data.abstract == valid_patent_data["abstract"]
    assert data.confidence == valid_patent_data["confidence"]


def test_community_data_model(valid_community_data):
    """Test CommunityData model validation."""
    data = CommunityData(**valid_community_data)
    assert data.content == valid_community_data["content"]
    assert data.source == valid_community_data["source"]
    assert data.confidence == valid_community_data["confidence"]


def test_web_data_schema(valid_web_data):
    """Test WebDataSchema model validation."""
    data = WebDataSchema(**valid_web_data)
    assert len(data.research) == 1
    assert len(data.patents) == 1
    assert len(data.community) == 1
    assert data.confidence == valid_web_data["confidence"]


def test_schema_registry():
    """Test SchemaRegistry schema retrieval."""
    # Test Swiss schemas
    swiss_target = SchemaRegistry.get_schema(DataSource.SWISS, "target")
    assert "targets" in swiss_target["required"]
    assert "scores" in swiss_target["required"]

    swiss_adme = SchemaRegistry.get_schema(DataSource.SWISS, "adme")
    assert "properties" in swiss_adme["required"]
    assert "predictions" in swiss_adme["required"]

    # Test community schemas
    web_schema = SchemaRegistry.get_schema(DataSource.COMMUNITY, "web")
    assert "research" in web_schema["properties"]
    assert "patents" in web_schema["properties"]

    # Test invalid schema
    with pytest.raises(ValueError):
        SchemaRegistry.get_schema(DataSource.SWISS, "invalid")


@pytest.mark.parametrize(
    "level,should_pass",
    [
        (ValidationLevel.STRICT, True),
        (ValidationLevel.NORMAL, True),
        (ValidationLevel.LENIENT, True),
    ],
)
def test_schema_validator_valid_data(valid_web_data, level, should_pass):
    """Test SchemaValidator with valid data."""
    validator = SchemaValidator(level=level)
    result = validator.validate(valid_web_data, DataSource.COMMUNITY, "web")
    assert result.is_valid == should_pass
    assert len(result.errors) == 0


@pytest.mark.parametrize(
    "level,field,value,should_pass",
    [
        (ValidationLevel.STRICT, "confidence", "0.9", False),  # Wrong type
        (ValidationLevel.NORMAL, "confidence", "0.9", True),  # Type coercion
        (ValidationLevel.LENIENT, "confidence", "invalid", True),  # Warning only
    ],
)
def test_schema_validator_type_coercion(
    valid_web_data,
    level,
    field,
    value,
    should_pass,
):
    """Test SchemaValidator type coercion."""
    validator = SchemaValidator(level=level)
    data = valid_web_data.copy()
    data[field] = value
    result = validator.validate(data, DataSource.COMMUNITY, "web")
    assert result.is_valid == should_pass


def test_schema_validator_missing_required():
    """Test SchemaValidator with missing required fields."""
    validator = SchemaValidator(level=ValidationLevel.STRICT)
    data = {"title": "Test"}  # Missing required fields
    result = validator.validate(data, DataSource.COMMUNITY, "research")
    assert not result.is_valid
    assert len(result.errors) > 0


def test_schema_validator_nested_validation():
    """Test SchemaValidator with nested objects."""
    validator = SchemaValidator(level=ValidationLevel.STRICT)
    data = {
        "properties": {
            "mw": 100,
            "logp": 2.5,
            # Missing required hbd, hba
        },
        "predictions": {
            "bbb": True,
            "pgp": False,
            "cyp": ["CYP3A4"],
        },
    }
    result = validator.validate(data, DataSource.SWISS, "adme")
    assert not result.is_valid
    assert any("hbd" in error for error in result.errors)
    assert any("hba" in error for error in result.errors)


def test_type_coercion_methods():
    """Test individual type coercion methods."""
    validator = SchemaValidator()

    # Test number coercion
    assert validator._coerce_number("1.5") == 1.5
    assert validator._coerce_number(1) == 1.0
    with pytest.raises(ValueError):
        validator._coerce_number("invalid")

    # Test integer coercion
    assert validator._coerce_integer("1") == 1
    assert validator._coerce_integer(1.0) == 1
    with pytest.raises(ValueError):
        validator._coerce_integer(1.5)
        validator._coerce_integer("invalid")

    # Test string coercion
    assert validator._coerce_string(123) == "123"
    assert validator._coerce_string(True) == "True"

    # Test boolean coercion
    assert validator._coerce_boolean("true")
    assert validator._coerce_boolean("1")
    assert not validator._coerce_boolean("false")
    assert not validator._coerce_boolean("0")
    with pytest.raises(ValueError):
        validator._coerce_boolean("invalid")


def test_validation_result_creation():
    """Test ValidationResult creation and properties."""
    result = ValidationResult(
        is_valid=True,
        errors=[],
        warnings=["warning1"],
        source=DataSource.SWISS,
        level=ValidationLevel.NORMAL,
    )
    assert result.is_valid
    assert len(result.errors) == 0
    assert len(result.warnings) == 1
    assert result.source == DataSource.SWISS
    assert result.level == ValidationLevel.NORMAL


def test_data_source_enum():
    """Test DataSource enum values."""
    assert DataSource.SWISS.value == "swiss"
    assert DataSource.COMMUNITY.value == "community"
    assert DataSource.SOCIAL.value == "social"


def test_validation_level_enum():
    """Test ValidationLevel enum values."""
    assert ValidationLevel.STRICT.value == "strict"
    assert ValidationLevel.NORMAL.value == "normal"
    assert ValidationLevel.LENIENT.value == "lenient"
