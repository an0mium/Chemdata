"""Tests for compound analysis constants and enums."""

import pytest
from enum import Enum

from ..constants import (
    PsychoactiveClass,
    RiskLevel,
    ActivityType,
    DEFAULT_VALUES,
    UNIT_CONVERSIONS,
)


def test_psychoactive_classes():
    """Test psychoactive class enumeration."""
    # Test enum values
    assert isinstance(PsychoactiveClass, type(Enum))
    assert PsychoactiveClass.STIMULANT.value == "stimulant"
    assert PsychoactiveClass.DEPRESSANT.value == "depressant"
    assert PsychoactiveClass.PSYCHEDELIC.value == "psychedelic"
    assert PsychoactiveClass.DISSOCIATIVE.value == "dissociative"
    assert PsychoactiveClass.DELIRIANT.value == "deliriant"
    assert PsychoactiveClass.UNKNOWN.value == "unknown"

    # Test string conversion
    assert str(PsychoactiveClass.STIMULANT) == "STIMULANT"
    assert PsychoactiveClass.STIMULANT.name == "STIMULANT"

    # Test comparison
    assert PsychoactiveClass.STIMULANT != PsychoactiveClass.DEPRESSANT
    assert PsychoactiveClass.UNKNOWN == PsychoactiveClass.UNKNOWN

    # Test iteration
    classes = list(PsychoactiveClass)
    assert len(classes) > 0
    assert all(isinstance(c, PsychoactiveClass) for c in classes)


def test_risk_levels():
    """Test risk level enumeration."""
    # Test enum values
    assert isinstance(RiskLevel, type(Enum))
    assert RiskLevel.LOW.value == "low"
    assert RiskLevel.MODERATE.value == "moderate"
    assert RiskLevel.HIGH.value == "high"
    assert RiskLevel.SEVERE.value == "severe"
    assert RiskLevel.UNKNOWN.value == "unknown"

    # Test ordering
    assert RiskLevel.LOW < RiskLevel.MODERATE
    assert RiskLevel.MODERATE < RiskLevel.HIGH
    assert RiskLevel.HIGH < RiskLevel.SEVERE
    assert RiskLevel.UNKNOWN < RiskLevel.LOW

    # Test comparison methods
    assert RiskLevel.HIGH.is_higher_than(RiskLevel.MODERATE)
    assert RiskLevel.LOW.is_lower_than(RiskLevel.MODERATE)
    assert RiskLevel.MODERATE.is_between(RiskLevel.LOW, RiskLevel.HIGH)

    # Test string methods
    assert RiskLevel.HIGH.to_warning() == "HIGH RISK"
    assert RiskLevel.SEVERE.to_warning() == "SEVERE RISK - DANGER"


def test_activity_types():
    """Test activity type enumeration."""
    # Test enum values
    assert isinstance(ActivityType, type(Enum))
    assert ActivityType.AGONIST.value == "agonist"
    assert ActivityType.ANTAGONIST.value == "antagonist"
    assert ActivityType.PARTIAL_AGONIST.value == "partial_agonist"
    assert ActivityType.INVERSE_AGONIST.value == "inverse_agonist"
    assert ActivityType.UNKNOWN.value == "unknown"

    # Test classification methods
    assert ActivityType.AGONIST.is_activating()
    assert not ActivityType.ANTAGONIST.is_activating()
    assert ActivityType.PARTIAL_AGONIST.is_partial()
    assert ActivityType.INVERSE_AGONIST.is_inverse()

    # Test comparison
    assert ActivityType.AGONIST != ActivityType.ANTAGONIST
    assert ActivityType.UNKNOWN == ActivityType.UNKNOWN

    # Test string methods
    assert ActivityType.PARTIAL_AGONIST.to_display() == "Partial Agonist"
    assert ActivityType.INVERSE_AGONIST.to_display() == "Inverse Agonist"


def test_default_values():
    """Test default value constants."""
    # Test numeric defaults
    assert isinstance(DEFAULT_VALUES["confidence"], float)
    assert 0 <= DEFAULT_VALUES["confidence"] <= 1
    assert isinstance(DEFAULT_VALUES["activity"], float)
    assert DEFAULT_VALUES["activity"] > 0

    # Test string defaults
    assert isinstance(DEFAULT_VALUES["mechanism"], str)
    assert DEFAULT_VALUES["mechanism"] == "unknown"
    assert isinstance(DEFAULT_VALUES["activity_type"], str)
    assert DEFAULT_VALUES["activity_type"] == "unknown"

    # Test collection defaults
    assert isinstance(DEFAULT_VALUES["targets"], list)
    assert len(DEFAULT_VALUES["targets"]) == 0
    assert isinstance(DEFAULT_VALUES["effects"], dict)
    assert len(DEFAULT_VALUES["effects"]) == 0

    # Test enum defaults
    assert DEFAULT_VALUES["psychoactive_class"] == PsychoactiveClass.UNKNOWN
    assert DEFAULT_VALUES["risk_level"] == RiskLevel.UNKNOWN


def test_unit_conversions():
    """Test unit conversion constants."""
    # Test concentration conversions
    assert UNIT_CONVERSIONS["nM"]["to_base"] == 1
    assert UNIT_CONVERSIONS["μM"]["to_base"] == 1000
    assert UNIT_CONVERSIONS["mM"]["to_base"] == 1000000

    # Test conversion methods
    assert UNIT_CONVERSIONS["nM"]["convert_to"]["μM"](1000) == 1
    assert UNIT_CONVERSIONS["μM"]["convert_to"]["nM"](1) == 1000
    assert UNIT_CONVERSIONS["mM"]["convert_to"]["μM"](1) == 1000

    # Test validation
    with pytest.raises(KeyError):
        _ = UNIT_CONVERSIONS["invalid"]
    with pytest.raises(KeyError):
        _ = UNIT_CONVERSIONS["nM"]["convert_to"]["invalid"]

    # Test base unit identification
    assert UNIT_CONVERSIONS["nM"]["is_base"]
    assert not UNIT_CONVERSIONS["μM"]["is_base"]
    assert not UNIT_CONVERSIONS["mM"]["is_base"]


def test_enum_operations():
    """Test enum operation methods."""
    # Test PsychoactiveClass operations
    assert PsychoactiveClass.combine([
        PsychoactiveClass.STIMULANT,
        PsychoactiveClass.PSYCHEDELIC
    ]) == "stimulant/psychedelic"

    assert PsychoactiveClass.from_string("stimulant") == PsychoactiveClass.STIMULANT
    assert PsychoactiveClass.from_string("invalid") == PsychoactiveClass.UNKNOWN

    # Test RiskLevel operations
    assert RiskLevel.max([
        RiskLevel.LOW,
        RiskLevel.HIGH,
        RiskLevel.MODERATE
    ]) == RiskLevel.HIGH

    assert RiskLevel.min([
        RiskLevel.HIGH,
        RiskLevel.MODERATE,
        RiskLevel.SEVERE
    ]) == RiskLevel.MODERATE

    # Test ActivityType operations
    assert ActivityType.is_opposite(
        ActivityType.AGONIST,
        ActivityType.ANTAGONIST
    )
    assert not ActivityType.is_opposite(
        ActivityType.AGONIST,
        ActivityType.PARTIAL_AGONIST
    )


def test_constant_immutability():
    """Test that constants are immutable."""
    # Test enum immutability
    with pytest.raises(AttributeError):
        PsychoactiveClass.STIMULANT.value = "new_value"
    with pytest.raises(AttributeError):
        RiskLevel.HIGH.value = "new_value"
    with pytest.raises(AttributeError):
        ActivityType.AGONIST.value = "new_value"

    # Test default value immutability
    original = DEFAULT_VALUES.copy()
    with pytest.raises(TypeError):
        DEFAULT_VALUES["new_key"] = "new_value"
    assert DEFAULT_VALUES == original

    # Test unit conversion immutability
    original = UNIT_CONVERSIONS.copy()
    with pytest.raises(TypeError):
        UNIT_CONVERSIONS["new_unit"] = {"to_base": 1}
    assert UNIT_CONVERSIONS == original


def test_enum_validation():
    """Test enum validation methods."""
    # Test PsychoactiveClass validation
    assert PsychoactiveClass.is_valid("stimulant")
    assert not PsychoactiveClass.is_valid("invalid")
    assert PsychoactiveClass.validate("stimulant") == PsychoactiveClass.STIMULANT
    assert PsychoactiveClass.validate("invalid") == PsychoactiveClass.UNKNOWN

    # Test RiskLevel validation
    assert RiskLevel.is_valid("high")
    assert not RiskLevel.is_valid("invalid")
    assert RiskLevel.validate("high") == RiskLevel.HIGH
    assert RiskLevel.validate("invalid") == RiskLevel.UNKNOWN

    # Test ActivityType validation
    assert ActivityType.is_valid("agonist")
    assert not ActivityType.is_valid("invalid")
    assert ActivityType.validate("agonist") == ActivityType.AGONIST
    assert ActivityType.validate("invalid") == ActivityType.UNKNOWN
