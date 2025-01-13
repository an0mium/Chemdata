"""Validation module for structure analysis system.

This module provides:
1. Input validation
2. Configuration validation
3. State validation
4. Validation utilities
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Any, Type, Union, Callable
from pathlib import Path

from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from .errors import ValidationError
from .logging import get_logger

logger = get_logger(__name__)


@dataclass
class ValidationResult:
    """Result of validation operation."""

    valid: bool
    errors: List[str]
    warnings: List[str]
    details: Optional[Dict[str, Any]] = None

    @property
    def has_errors(self) -> bool:
        """Check if validation found errors."""
        return len(self.errors) > 0

    @property
    def has_warnings(self) -> bool:
        """Check if validation found warnings."""
        return len(self.warnings) > 0


class Validator:
    """Base class for validators."""

    def __init__(self, raise_on_error: bool = True):
        """Initialize validator.

        Args:
            raise_on_error: Whether to raise exception on validation error
        """
        self.raise_on_error = raise_on_error
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def validate(self, value: Any) -> ValidationResult:
        """Validate value.

        Args:
            value: Value to validate

        Returns:
            Validation result

        Raises:
            ValidationError: If validation fails and raise_on_error is True
        """
        self.errors.clear()
        self.warnings.clear()

        try:
            self._validate(value)
        except Exception as e:
            self.errors.append(str(e))

        result = ValidationResult(
            valid=len(self.errors) == 0,
            errors=self.errors.copy(),
            warnings=self.warnings.copy(),
        )

        if self.raise_on_error and result.has_errors:
            raise ValidationError("\n".join(result.errors))

        return result

    def _validate(self, value: Any):
        """Internal validation method to be implemented by subclasses.

        Args:
            value: Value to validate

        Raises:
            NotImplementedError: If not implemented by subclass
        """
        raise NotImplementedError


class StructureValidator(Validator):
    """Validator for protein structures."""

    def _validate(self, structure: Structure):
        """Validate protein structure.

        Args:
            structure: Structure to validate
        """
        if not isinstance(structure, Structure):
            self.errors.append(f"Expected Structure, got {type(structure)}")
            return

        # Check structure has models
        if len(structure) == 0:
            self.errors.append("Structure has no models")
            return

        # Check structure has chains
        model = structure[0]
        if len(model) == 0:
            self.errors.append("Structure has no chains")
            return

        # Check chains have residues
        for chain in model:
            if len(chain) == 0:
                self.warnings.append(f"Chain {chain.id} has no residues")

        # Check residues have atoms
        for chain in model:
            for residue in chain:
                if len(residue) == 0:
                    self.warnings.append(f"Residue {residue.id} in chain {chain.id} has no atoms")


class ConfigValidator(Validator):
    """Validator for configuration objects."""

    def __init__(
        self,
        required_fields: Optional[List[str]] = None,
        field_types: Optional[Dict[str, Type]] = None,
        field_validators: Optional[Dict[str, Callable]] = None,
        raise_on_error: bool = True,
    ):
        """Initialize config validator.

        Args:
            required_fields: List of required field names
            field_types: Dictionary mapping field names to expected types
            field_validators: Dictionary mapping field names to validator functions
            raise_on_error: Whether to raise exception on validation error
        """
        super().__init__(raise_on_error)
        self.required_fields = required_fields or []
        self.field_types = field_types or {}
        self.field_validators = field_validators or {}

    def _validate(self, config: Dict[str, Any]):
        """Validate configuration dictionary.

        Args:
            config: Configuration to validate
        """
        if not isinstance(config, dict):
            self.errors.append(f"Expected dict, got {type(config)}")
            return

        # Check required fields
        for field in self.required_fields:
            if field not in config:
                self.errors.append(f"Missing required field: {field}")

        # Check field types
        for field, expected_type in self.field_types.items():
            if field in config:
                value = config[field]
                if not isinstance(value, expected_type):
                    self.errors.append(f"Field {field} has wrong type: expected {expected_type}, got {type(value)}")

        # Run field validators
        for field, validator in self.field_validators.items():
            if field in config:
                try:
                    validator(config[field])
                except Exception as e:
                    self.errors.append(f"Validation failed for field {field}: {str(e)}")


class ComponentValidator(Validator):
    """Validator for component states."""

    def __init__(
        self,
        required_attributes: Optional[List[str]] = None,
        state_validators: Optional[Dict[str, Callable]] = None,
        raise_on_error: bool = True,
    ):
        """Initialize component validator.

        Args:
            required_attributes: List of required attribute names
            state_validators: Dictionary mapping state names to validator functions
            raise_on_error: Whether to raise exception on validation error
        """
        super().__init__(raise_on_error)
        self.required_attributes = required_attributes or []
        self.state_validators = state_validators or {}

    def _validate(self, component: Any):
        """Validate component state.

        Args:
            component: Component to validate
        """
        # Check required attributes
        for attr in self.required_attributes:
            if not hasattr(component, attr):
                self.errors.append(f"Missing required attribute: {attr}")

        # Run state validators
        for state, validator in self.state_validators.items():
            if hasattr(component, state):
                try:
                    validator(getattr(component, state))
                except Exception as e:
                    self.errors.append(f"Validation failed for state {state}: {str(e)}")


def validate_structure(
    structure: Structure,
    raise_on_error: bool = True,
) -> ValidationResult:
    """Validate protein structure.

    Args:
        structure: Structure to validate
        raise_on_error: Whether to raise exception on validation error

    Returns:
        Validation result
    """
    validator = StructureValidator(raise_on_error=raise_on_error)
    return validator.validate(structure)


def validate_config(
    config: Dict[str, Any],
    required_fields: Optional[List[str]] = None,
    field_types: Optional[Dict[str, Type]] = None,
    field_validators: Optional[Dict[str, Callable]] = None,
    raise_on_error: bool = True,
) -> ValidationResult:
    """Validate configuration dictionary.

    Args:
        config: Configuration to validate
        required_fields: List of required field names
        field_types: Dictionary mapping field names to expected types
        field_validators: Dictionary mapping field names to validator functions
        raise_on_error: Whether to raise exception on validation error

    Returns:
        Validation result
    """
    validator = ConfigValidator(
        required_fields=required_fields,
        field_types=field_types,
        field_validators=field_validators,
        raise_on_error=raise_on_error,
    )
    return validator.validate(config)


def validate_component(
    component: Any,
    required_attributes: Optional[List[str]] = None,
    state_validators: Optional[Dict[str, Callable]] = None,
    raise_on_error: bool = True,
) -> ValidationResult:
    """Validate component state.

    Args:
        component: Component to validate
        required_attributes: List of required attribute names
        state_validators: Dictionary mapping state names to validator functions
        raise_on_error: Whether to raise exception on validation error

    Returns:
        Validation result
    """
    validator = ComponentValidator(
        required_attributes=required_attributes,
        state_validators=state_validators,
        raise_on_error=raise_on_error,
    )
    return validator.validate(component)


def is_valid_structure(structure: Structure) -> bool:
    """Check if structure is valid.

    Args:
        structure: Structure to validate

    Returns:
        True if structure is valid
    """
    try:
        result = validate_structure(structure, raise_on_error=False)
        return result.valid
    except Exception:
        return False


def is_valid_config(
    config: Dict[str, Any],
    required_fields: Optional[List[str]] = None,
    field_types: Optional[Dict[str, Type]] = None,
) -> bool:
    """Check if configuration is valid.

    Args:
        config: Configuration to validate
        required_fields: List of required field names
        field_types: Dictionary mapping field names to expected types

    Returns:
        True if configuration is valid
    """
    try:
        result = validate_config(
            config,
            required_fields=required_fields,
            field_types=field_types,
            raise_on_error=False,
        )
        return result.valid
    except Exception:
        return False


def is_valid_component(
    component: Any,
    required_attributes: Optional[List[str]] = None,
) -> bool:
    """Check if component is valid.

    Args:
        component: Component to validate
        required_attributes: List of required attribute names

    Returns:
        True if component is valid
    """
    try:
        result = validate_component(
            component,
            required_attributes=required_attributes,
            raise_on_error=False,
        )
        return result.valid
    except Exception:
        return False
