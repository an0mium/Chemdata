"""Base compound data models.

This module provides the core compound data model functionality:
- CompoundData class for basic compound properties
- Validation utilities
- Type definitions
"""

from .core import CompoundData as BaseCompound
from .validation import ValidationError
from .types import (
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

__all__ = [
    'BaseCompound',
    'ValidationError',
    'CompoundType',
    'LegalStatus',
    'PsychoactiveClass',
    'NootropicMechanism',
    'BBBPermeability',
    'BindingType',
    'ActivityType',
    'RiskLevel',
    'TargetData',
]
