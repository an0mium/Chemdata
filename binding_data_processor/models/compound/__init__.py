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
"""

from .base import (
    BaseCompound,
    ValidationError,
    CompoundType,
    LegalStatus,
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
    BindingType,
    ActivityType,
    RiskLevel,
    TargetData,
)
from .enrichment import EnrichedCompound
from .ml import MLCompound
from .analysis import AnalyzedCompound

# Main compound class is the fully featured version
Compound = AnalyzedCompound

__all__ = [
    # Main class
    'Compound',
    
    # Base classes
    'BaseCompound',
    'EnrichedCompound',
    'MLCompound',
    'AnalyzedCompound',
    
    # Types
    'CompoundType',
    'LegalStatus',
    'PsychoactiveClass',
    'NootropicMechanism',
    'BBBPermeability',
    'BindingType',
    'ActivityType',
    'RiskLevel',
    'TargetData',
    
    # Exceptions
    'ValidationError',
]
