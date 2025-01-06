"""Tests for compound type validation functionality."""

import pytest

from ..types import (
    CompoundType,
    LegalStatus,
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
    RiskLevel,
)
from ..validation import (
    validate_compound_type,
    validate_legal_status,
    validate_psychoactive_class,
    validate_nootropic_mechanism,
    validate_bbb_permeability,
)


def test_compound_type_validation():
    """Test compound type enum validation."""
    # Test all valid values
    for compound_type in CompoundType:
        is_valid, error = validate_compound_type(compound_type.value)
        assert is_valid
        assert error is None

    # Test specific new values
    valid_types = [
        "neurotransmitter",
        "psychoactive",
        "research_chemical",
        "novel_psychoactive_substance",
        "pharmaceutical",
        "natural_product",
        "small_molecule",
        "other",
        "unknown",
    ]
    for type_value in valid_types:
        is_valid, error = validate_compound_type(type_value)
        assert is_valid
        assert error is None

    # Test invalid values
    invalid_types = [
        "",
        "invalid_type",
        "INVALID",
        None,
        123,
    ]
    for type_value in invalid_types:
        is_valid, error = validate_compound_type(type_value)
        assert not is_valid
        assert error is not None
        assert "Invalid compound type" in error


def test_legal_status_validation():
    """Test legal status enum validation."""
    # Test all valid values
    for status in LegalStatus:
        is_valid, error = validate_legal_status(status.value)
        assert is_valid
        assert error is None

    # Test specific new values
    valid_statuses = [
        "legal",
        "controlled",
        "illegal",
        "research_only",
        "unscheduled",
        "prescription_only",
        "over_the_counter",
        "approved",
        "investigational",
        "withdrawn",
        "banned",
        "unknown",
    ]
    for status in valid_statuses:
        is_valid, error = validate_legal_status(status)
        assert is_valid
        assert error is None

    # Test invalid values
    invalid_statuses = [
        "",
        "invalid_status",
        "INVALID",
        None,
        123,
    ]
    for status in invalid_statuses:
        is_valid, error = validate_legal_status(status)
        assert not is_valid
        assert error is not None
        assert "Invalid legal status" in error


def test_psychoactive_class_validation():
    """Test psychoactive class enum validation."""
    # Test all valid values
    for pclass in PsychoactiveClass:
        is_valid, error = validate_psychoactive_class(pclass.value)
        assert is_valid
        assert error is None

    # Test specific new values
    valid_classes = [
        "psychedelic",
        "empathogen",
        "stimulant",
        "depressant",
        "dissociative",
        "deliriant",
        "nootropic",
        "anxiolytic",
        "antipsychotic",
        "antidepressant",
        "mood_stabilizer",
        "unknown",
    ]
    for pclass in valid_classes:
        is_valid, error = validate_psychoactive_class(pclass)
        assert is_valid
        assert error is None

    # Test invalid values
    invalid_classes = [
        "",
        "invalid_class",
        "INVALID",
        None,
        123,
    ]
    for pclass in invalid_classes:
        is_valid, error = validate_psychoactive_class(pclass)
        assert not is_valid
        assert error is not None
        assert "Invalid psychoactive class" in error


def test_nootropic_mechanism_validation():
    """Test nootropic mechanism enum validation."""
    # Test all valid values
    for mechanism in NootropicMechanism:
        is_valid, error = validate_nootropic_mechanism(mechanism.value)
        assert is_valid
        assert error is None

    # Test specific new values
    valid_mechanisms = [
        "cholinergic",
        "glutamatergic",
        "dopaminergic",
        "serotonergic",
        "gaba_modulation",
        "ampakine",
        "bdnf_modulation",
        "ngf_modulation",
        "neuroplasticity",
        "anti_inflammatory",
        "antioxidant",
        "memory_enhancement",
        "focus_improvement",
        "neuroprotection",
        "nootropic_synergy",
        "cognitive_modulation",
        "brain_metabolism",
        "unknown",
    ]
    for mechanism in valid_mechanisms:
        is_valid, error = validate_nootropic_mechanism(mechanism)
        assert is_valid
        assert error is None

    # Test invalid values
    invalid_mechanisms = [
        "",
        "invalid_mechanism",
        "INVALID",
        None,
        123,
    ]
    for mechanism in invalid_mechanisms:
        is_valid, error = validate_nootropic_mechanism(mechanism)
        assert not is_valid
        assert error is not None
        assert "Invalid nootropic mechanism" in error


def test_bbb_permeability_validation():
    """Test BBB permeability enum validation."""
    # Test all valid values
    for permeability in BBBPermeability:
        is_valid, error = validate_bbb_permeability(permeability.value)
        assert is_valid
        assert error is None

    # Test specific values
    valid_permeabilities = [
        "high",
        "moderate",
        "low",
        "negligible",
        "unknown",
    ]
    for permeability in valid_permeabilities:
        is_valid, error = validate_bbb_permeability(permeability)
        assert is_valid
        assert error is None

    # Test invalid values
    invalid_permeabilities = [
        "",
        "invalid_permeability",
        "INVALID",
        None,
        123,
    ]
    for permeability in invalid_permeabilities:
        is_valid, error = validate_bbb_permeability(permeability)
        assert not is_valid
        assert error is not None
        assert "Invalid BBB permeability" in error


def test_risk_level_validation():
    """Test risk level enum validation."""
    # Test all valid values
    for level in RiskLevel:
        assert level.value in {
            "severe",
            "high",
            "moderate",
            "low",
            "minimal",
            "unknown",
        }

    # Test specific values
    valid_levels = [
        RiskLevel.SEVERE,
        RiskLevel.HIGH,
        RiskLevel.MODERATE,
        RiskLevel.LOW,
        RiskLevel.MINIMAL,
        RiskLevel.UNKNOWN,
    ]
    for level in valid_levels:
        assert isinstance(level, RiskLevel)
        assert level.value in {
            "severe",
            "high",
            "moderate",
            "low",
            "minimal",
            "unknown",
        }

    # Test creation from string
    for value in ["severe", "high", "moderate", "low", "minimal", "unknown"]:
        level = RiskLevel(value)
        assert isinstance(level, RiskLevel)
        assert level.value == value

    # Test invalid values
    with pytest.raises(ValueError):
        RiskLevel("invalid_level")
