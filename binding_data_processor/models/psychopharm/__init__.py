"""Psychopharmacological property functionality.

This package provides classes and utilities for handling psychopharmacological properties:
1. Receptor binding profiles
2. Blood-brain barrier penetration
3. Psychoactive effects classification
4. Nootropic activity metrics
5. Abuse potential assessment
6. Risk analysis
7. Community data integration

The functionality is split across several modules:
- base.py: Core data structures and base mixin class
- binding.py: Receptor binding analysis
- activity.py: Psychoactive and nootropic effects
- safety.py: Risk assessment and safety analysis
- community.py: Community data integration
- enrichment.py: Web data enrichment
- compound.py: Complete compound data model

Key Features:
- Comprehensive receptor binding analysis
- Blood-brain barrier penetration prediction
- Psychoactive effect classification
- Nootropic activity assessment
- Safety and risk analysis
- Community data integration
- Web data enrichment
- Structured data validation
"""

from ..compound.types import (
    CompoundType,
    LegalStatus,
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
    BindingType,
    RiskLevel,
    ToxicityClass,
    TargetData,
    StringSet,
    StringDict,
    ValidationErrors,
    OptionalStr,
    TargetDict,
    DoseRange,
    TimeRange,
    RiskScore,
    EffectScore,
    ReceptorBinding,
)
from .binding import ReceptorProfileMixin
from .activity import ActivityProfileMixin
from .safety import SafetyProfileMixin
from .community import CommunityDataMixin
from .enrichment import WebEnrichmentMixin
from .compound import PsychopharmCompound


class PsychopharmMixin(
    ReceptorProfileMixin,
    ActivityProfileMixin,
    SafetyProfileMixin,
    CommunityDataMixin,
    WebEnrichmentMixin,
):
    """Main mixin class combining all psychopharmacological functionality.

    This mixin combines all the specialized mixins into a single class that can be
    used to add complete psychopharmacological analysis capabilities to any class.

    For a complete compound data model that includes all this functionality, use
    PsychopharmCompound instead.
    """

    pass


__all__ = [
    # Core data types
    "PsychoactiveClass",
    "NootropicMechanism",
    "BBBPermeability",
    "RiskLevel",
    "ToxicityClass",
    "DoseRange",
    "TimeRange",
    "RiskScore",
    "EffectScore",
    "ReceptorBinding",
    # Individual mixins
    "ReceptorProfileMixin",
    "ActivityProfileMixin",
    "SafetyProfileMixin",
    "CommunityDataMixin",
    "WebEnrichmentMixin",
    # Combined mixin
    "PsychopharmMixin",
    # Complete compound model
    "PsychopharmCompound",
]

# Package metadata
__version__ = "0.1.0"
__author__ = "Armand"
__description__ = "Psychopharmacological analysis and data enrichment"
__license__ = "MIT"

# Module organization
__module_dependencies__ = {
    "base": [],
    "binding": ["base"],
    "activity": ["base", "binding"],
    "safety": ["base", "binding", "activity"],
    "community": ["base"],
    "enrichment": ["base", "community"],
    "compound": ["base", "binding", "activity", "safety", "community", "enrichment"],
}
