"""Schema validation for web enrichment data.

This module provides schema validation for:
- API responses
- Web scraping results
- Social media data
"""

import logging
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
from datetime import datetime
from enum import Enum

import jsonschema
from jsonschema import ValidationError
from pydantic import BaseModel, Field, ConfigDict


class BaseSchema:
    """Base class for all schema objects."""

    pass


class DataSource(Enum):
    """Supported data sources."""

    SWISS = "swiss"
    COMMUNITY = "community"
    SOCIAL = "social"


class ValidationLevel(Enum):
    """Validation strictness levels."""

    STRICT = "strict"  # All fields required, exact types
    NORMAL = "normal"  # Required fields only, type coercion
    LENIENT = "lenient"  # Best effort validation


@dataclass
class ValidationResult:
    """Result of schema validation."""

    is_valid: bool
    errors: List[str]
    warnings: List[str]
    source: DataSource
    level: ValidationLevel


# Pydantic Models for Data Validation


class ResearchData(BaseModel):
    """Research paper data model."""

    model_config = ConfigDict(extra="forbid")

    title: str = Field(..., description="Paper title")
    abstract: str = Field(..., description="Paper abstract")
    authors: List[str] = Field(default_factory=list, description="List of authors")
    date: Optional[datetime] = Field(None, description="Publication date")
    doi: Optional[str] = Field(None, description="Digital Object Identifier")
    url: str = Field(..., description="Source URL")
    source: str = Field(..., description="Source (e.g., PubMed, Google Scholar)")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Extraction confidence score")
    keywords: List[str] = Field(default_factory=list, description="Keywords")
    full_text: Optional[str] = Field(None, description="Full paper text if available")


class PatentData(BaseModel):
    """Patent document data model."""

    model_config = ConfigDict(extra="forbid")

    title: str = Field(..., description="Patent title")
    abstract: str = Field(..., description="Patent abstract")
    inventors: List[str] = Field(default_factory=list, description="List of inventors")
    assignee: Optional[str] = Field(None, description="Patent assignee")
    filing_date: Optional[datetime] = Field(None, description="Filing date")
    publication_date: Optional[datetime] = Field(None, description="Publication date")
    patent_number: str = Field(..., description="Patent number")
    url: str = Field(..., description="Source URL")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Extraction confidence score")
    claims: List[str] = Field(default_factory=list, description="Patent claims")
    chemical_structures: List[Dict[str, str]] = Field(default_factory=list, description="Chemical structures mentioned in patent")


class CommunityData(BaseModel):
    """Community forum data model."""

    model_config = ConfigDict(extra="forbid")

    title: Optional[str] = Field(None, description="Post title")
    content: str = Field(..., description="Post content")
    author: Optional[str] = Field(None, description="Post author")
    date: Optional[datetime] = Field(None, description="Post date")
    url: str = Field(..., description="Source URL")
    source: str = Field(..., description="Source (e.g., Reddit, Bluelight)")
    subforum: Optional[str] = Field(None, description="Subforum or category")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Extraction confidence score")
    sentiment: Optional[float] = Field(None, ge=-1.0, le=1.0, description="Sentiment score")
    effects: List[str] = Field(default_factory=list, description="Reported effects")
    dosage: Optional[str] = Field(None, description="Reported dosage")
    roa: Optional[str] = Field(None, description="Route of administration")
    safety_notes: List[str] = Field(default_factory=list, description="Safety-related notes")
    mentions: List[Dict[str, str]] = Field(
        default_factory=list,
        description="Entity mentions with context",
    )
    tags: List[str] = Field(default_factory=list, description="Post tags")


class WebDataSchema(BaseModel):
    """Combined web data schema."""

    model_config = ConfigDict(extra="forbid")

    research: List[ResearchData] = Field(default_factory=list, description="Research paper data")
    patents: List[PatentData] = Field(default_factory=list, description="Patent data")
    community: List[CommunityData] = Field(default_factory=list, description="Community data")
    metadata: dict = Field(default_factory=dict, description="Metadata about the extraction")
    confidence: Optional[float] = Field(
        None,
        ge=0.0,
        le=1.0,
        description="Overall confidence score",
    )


# Schema Registry for Flexible Validation


class SchemaRegistry:
    """Registry of JSON schemas for different data sources."""

    # Convert Pydantic models to JSON schemas
    RESEARCH_DATA_SCHEMA = ResearchData.model_json_schema()
    PATENT_DATA_SCHEMA = PatentData.model_json_schema()
    COMMUNITY_DATA_SCHEMA = CommunityData.model_json_schema()
    WEB_DATA_SCHEMA = WebDataSchema.model_json_schema()

    # Swiss tools schemas
    SWISS_TARGET_SCHEMA = {
        "type": "object",
        "required": ["targets", "scores"],
        "properties": {
            "targets": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": ["name", "uniprot_id", "score"],
                    "properties": {
                        "name": {"type": "string"},
                        "uniprot_id": {"type": "string"},
                        "score": {"type": "number"},
                    },
                },
            },
            "scores": {
                "type": "object",
                "required": ["probability", "reliability"],
                "properties": {
                    "probability": {"type": "number"},
                    "reliability": {"type": "number"},
                },
            },
        },
    }

    SWISS_ADME_SCHEMA = {
        "type": "object",
        "required": ["properties", "predictions"],
        "properties": {
            "properties": {
                "type": "object",
                "required": ["mw", "logp", "hbd", "hba"],
                "properties": {
                    "mw": {"type": "number"},
                    "logp": {"type": "number"},
                    "hbd": {"type": "integer"},
                    "hba": {"type": "integer"},
                },
            },
            "predictions": {
                "type": "object",
                "required": ["bbb", "pgp", "cyp"],
                "properties": {
                    "bbb": {"type": "boolean"},
                    "pgp": {"type": "boolean"},
                    "cyp": {"type": "array", "items": {"type": "string"}},
                },
            },
        },
    }

    @classmethod
    def get_schema(cls, source: DataSource, schema_name: str) -> Dict[str, Any]:
        """Get schema by source and name.

        Args:
            source: Data source
            schema_name: Schema name

        Returns:
            JSON schema

        Raises:
            ValueError: If schema not found
        """
        schema_map = {
            DataSource.SWISS: {
                "target": cls.SWISS_TARGET_SCHEMA,
                "adme": cls.SWISS_ADME_SCHEMA,
            },
            DataSource.COMMUNITY: {
                "web": cls.WEB_DATA_SCHEMA,
                "research": cls.RESEARCH_DATA_SCHEMA,
                "patent": cls.PATENT_DATA_SCHEMA,
                "forum": cls.COMMUNITY_DATA_SCHEMA,
            },
            DataSource.SOCIAL: {
                "post": cls.COMMUNITY_DATA_SCHEMA,  # Reuse community schema
            },
        }

        try:
            return schema_map[source][schema_name]
        except KeyError:
            raise ValueError(f"Schema not found: {source.value}/{schema_name}")


class SchemaValidator:
    """Schema validator with flexible validation levels."""

    def __init__(
        self,
        level: ValidationLevel = ValidationLevel.NORMAL,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize schema validator.

        Args:
            level: Validation strictness level
            logger: Optional logger instance
        """
        self.level = level
        self.logger = logger or logging.getLogger(self.__class__.__name__)

    def validate(
        self,
        data: Dict[str, Any],
        source: DataSource,
        schema_name: str,
    ) -> ValidationResult:
        """Validate data against schema.

        Args:
            data: Data to validate
            source: Data source
            schema_name: Schema name

        Returns:
            Validation result
        """
        schema = SchemaRegistry.get_schema(source, schema_name)
        errors = []
        warnings = []

        try:
            if self.level == ValidationLevel.STRICT:
                jsonschema.validate(data, schema)
            else:
                # Validate required fields
                self._validate_required(data, schema, errors, warnings)

                # Validate types with coercion
                if self.level == ValidationLevel.NORMAL:
                    self._validate_types(data, schema, errors, warnings)

        except ValidationError as e:
            errors.append(str(e))

        return ValidationResult(
            is_valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            source=source,
            level=self.level,
        )

    def _validate_required(
        self,
        data: Dict[str, Any],
        schema: Dict[str, Any],
        errors: List[str],
        warnings: List[str],
    ) -> None:
        """Validate required fields.

        Args:
            data: Data to validate
            schema: JSON schema
            errors: List to add errors to
            warnings: List to add warnings to
        """
        required = schema.get("required", [])
        properties = schema.get("properties", {})

        for field in required:
            if field not in data:
                if self.level == ValidationLevel.LENIENT:
                    warnings.append(f"Missing required field: {field}")
                else:
                    errors.append(f"Missing required field: {field}")
            elif field in properties:
                field_schema = properties[field]
                if field_schema.get("type") == "object":
                    self._validate_required(
                        data[field],
                        field_schema,
                        errors,
                        warnings,
                    )

    def _validate_types(
        self,
        data: Dict[str, Any],
        schema: Dict[str, Any],
        errors: List[str],
        warnings: List[str],
    ) -> None:
        """Validate and coerce types.

        Args:
            data: Data to validate
            schema: JSON schema
            errors: List to add errors to
            warnings: List to add warnings to
        """
        properties = schema.get("properties", {})

        for field, value in data.items():
            if field in properties:
                field_schema = properties[field]
                field_type = field_schema.get("type")

                if field_type == "object" and isinstance(value, dict):
                    self._validate_types(
                        value,
                        field_schema,
                        errors,
                        warnings,
                    )
                elif field_type == "array" and isinstance(value, list):
                    items_schema = field_schema.get("items", {})
                    for item in value:
                        if items_schema.get("type") == "object":
                            self._validate_types(
                                item,
                                items_schema,
                                errors,
                                warnings,
                            )
                else:
                    try:
                        self._coerce_type(data, field, field_type)
                    except (ValueError, TypeError) as e:
                        if self.level == ValidationLevel.LENIENT:
                            warnings.append(f"Type error for {field}: {str(e)}")
                        else:
                            errors.append(f"Type error for {field}: {str(e)}")

    def _coerce_type(
        self,
        data: Dict[str, Any],
        field: str,
        field_type: str,
    ) -> None:
        """Coerce value to expected type.

        Args:
            data: Data dictionary
            field: Field name
            field_type: Expected type

        Raises:
            ValueError: If coercion fails
        """
        value = data[field]
        coercion_map = {
            "number": self._coerce_number,
            "integer": self._coerce_integer,
            "string": self._coerce_string,
            "boolean": self._coerce_boolean,
        }

        if field_type in coercion_map:
            data[field] = coercion_map[field_type](value)

    def _coerce_number(self, value: Any) -> float:
        """Coerce value to number (float).

        Args:
            value: Value to coerce

        Returns:
            Coerced float value

        Raises:
            ValueError: If coercion fails
        """
        if isinstance(value, str):
            return float(value)
        elif isinstance(value, (int, float)):
            return float(value)
        raise ValueError(f"Expected number, got {type(value)}")

    def _coerce_integer(self, value: Any) -> int:
        """Coerce value to integer.

        Args:
            value: Value to coerce

        Returns:
            Coerced integer value

        Raises:
            ValueError: If coercion fails
        """
        if isinstance(value, str):
            return int(value)
        elif isinstance(value, float):
            if value.is_integer():
                return int(value)
            raise ValueError("Float value is not an integer")
        elif isinstance(value, int):
            return value
        raise ValueError(f"Expected integer, got {type(value)}")

    def _coerce_string(self, value: Any) -> str:
        """Coerce value to string.

        Args:
            value: Value to coerce

        Returns:
            Coerced string value
        """
        return str(value)

    def _coerce_boolean(self, value: Any) -> bool:
        """Coerce value to boolean.

        Args:
            value: Value to coerce

        Returns:
            Coerced boolean value

        Raises:
            ValueError: If coercion fails
        """
        if isinstance(value, str):
            value = value.lower()
            if value in ("true", "1", "yes"):
                return True
            elif value in ("false", "0", "no"):
                return False
            raise ValueError(f"Cannot convert to boolean: {value}")
        elif isinstance(value, bool):
            return value
        raise ValueError(f"Expected boolean, got {type(value)}")


def validate_community_data(data: Dict[str, Any], level: Optional[ValidationLevel] = None) -> ValidationResult:
    """Validate community data against schema.

    Args:
        data: Community data to validate
        level: Optional validation level (defaults to NORMAL)

    Returns:
        ValidationResult containing validation status and any errors/warnings
    """
    validator = SchemaValidator(level=level or ValidationLevel.NORMAL)
    return validator.validate(data, DataSource.COMMUNITY, "forum")
