"""Compound data models.

This package provides a hierarchy of compound data models:

BaseCompound
    Basic chemical compound properties and validation

EnrichedCompound (extends BaseCompound)
    Web data enrichment capabilities

MLCompound (extends EnrichedCompound)
    Machine learning integration

AnalyzedCompound (extends MLCompound)
    Analysis capabilities

PsychopharmCompound (extends AnalyzedCompound)
    Psychopharmacology-specific capabilities
"""

from .base import BaseCompound, ValidationError
from .base.core import CompoundData  # Add CompoundData import
from .base.types import (
    CompoundType,
    LegalStatus,
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
    BindingType,
    ActivityType,
    RiskLevel,
    TargetData,
    TimeRange,
    DoseRange,
    EffectScore,
    RiskScore,
    ReceptorBinding,
    BindingData,
    ActivityData,
    PropertyData,
    SafetyData,
    PredictionResult,
)
from .enrichment import EnrichedCompound
from .ml import MLCompound
from .analysis import AnalyzedCompound
from .export import CompoundExporter

# Import psychopharm-specific functionality
from .psychopharm import PsychopharmCompound

# Main compound class is the fully featured version with psychopharm capabilities
Compound = PsychopharmCompound

# Export all public symbols
__all__ = [
    # Main classes
    "Compound",
    "CompoundData",  # Add CompoundData to exports
    # Base classes
    "BaseCompound",
    "EnrichedCompound",
    "MLCompound",
    "AnalyzedCompound",
    "PsychopharmCompound",
    "CompoundExporter",
    # Types
    "CompoundType",
    "LegalStatus",
    "PsychoactiveClass",
    "NootropicMechanism",
    "BBBPermeability",
    "BindingType",
    "ActivityType",
    "RiskLevel",
    "TargetData",
    "ValidationError",
    # Data structures
    "TimeRange",
    "DoseRange",
    "EffectScore",
    "RiskScore",
    "ReceptorBinding",
    "BindingData",
    "ActivityData",
    "PropertyData",
    "SafetyData",
    "PredictionResult",
]
