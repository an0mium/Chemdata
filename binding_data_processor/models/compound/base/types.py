"""Type definitions for compound data models.

This module defines:
1. Core enums for compound classification
2. Data structures for timing and scoring
3. Type aliases for common data types
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Set, Union, Any


class CompoundType(Enum):
    """Types of chemical compounds."""

    NEUROTRANSMITTER = "neurotransmitter"
    PSYCHOACTIVE = "psychoactive"
    RESEARCH_CHEMICAL = "research_chemical"
    NPS = "novel_psychoactive_substance"
    PHARMACEUTICAL = "pharmaceutical"
    NATURAL_PRODUCT = "natural_product"
    SMALL_MOLECULE = "small_molecule"
    OTHER = "other"
    UNKNOWN = "unknown"


class LegalStatus(Enum):
    """Legal status classifications."""

    LEGAL = "legal"
    CONTROLLED = "controlled"
    ILLEGAL = "illegal"
    RESEARCH_ONLY = "research_only"
    UNSCHEDULED = "unscheduled"
    PRESCRIPTION = "prescription_only"
    OTC = "over_the_counter"
    APPROVED = "approved"
    INVESTIGATIONAL = "investigational"
    WITHDRAWN = "withdrawn"
    BANNED = "banned"
    UNKNOWN = "unknown"


class PsychoactiveClass(Enum):
    """Classification of psychoactive effects."""

    PSYCHEDELIC = "psychedelic"
    EMPATHOGEN = "empathogen"
    STIMULANT = "stimulant"
    DEPRESSANT = "depressant"
    DISSOCIATIVE = "dissociative"
    DELIRIANT = "deliriant"
    NOOTROPIC = "nootropic"
    ANXIOLYTIC = "anxiolytic"
    ANTIPSYCHOTIC = "antipsychotic"
    ANTIDEPRESSANT = "antidepressant"
    MOOD_STABILIZER = "mood_stabilizer"
    UNKNOWN = "unknown"


class NootropicMechanism(Enum):
    """Mechanisms of nootropic activity."""

    CHOLINERGIC = "cholinergic"
    GLUTAMATERGIC = "glutamatergic"
    DOPAMINERGIC = "dopaminergic"
    SEROTONERGIC = "serotonergic"
    GABA = "gaba_modulation"
    AMPAKINE = "ampakine"
    BDNF = "bdnf_modulation"
    NGF = "ngf_modulation"
    NEUROPLASTICITY = "neuroplasticity"
    ANTI_INFLAMMATORY = "anti_inflammatory"
    ANTIOXIDANT = "antioxidant"
    MEMORY_ENHANCEMENT = "memory_enhancement"
    FOCUS_IMPROVEMENT = "focus_improvement"
    NEUROPROTECTION = "neuroprotection"
    NOOTROPIC_SYNERGY = "nootropic_synergy"
    COGNITIVE_MODULATION = "cognitive_modulation"
    BRAIN_METABOLISM = "brain_metabolism"
    UNKNOWN = "unknown"


class BBBPermeability(Enum):
    """Blood-brain barrier permeability classification."""

    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    NEGLIGIBLE = "negligible"
    UNKNOWN = "unknown"


class RiskLevel(Enum):
    """Risk assessment level."""

    SEVERE = "severe"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    MINIMAL = "minimal"
    UNKNOWN = "unknown"


@dataclass
class TimeRange:
    """Time range with min/max values."""

    min: float  # Minimum time in minutes
    max: float  # Maximum time in minutes
    typical: Optional[float] = None  # Typical time in minutes
    notes: Optional[str] = None  # Additional notes


@dataclass
class DoseRange:
    """Dose range with min/max values."""

    min: float  # Minimum dose in mg
    max: float  # Maximum dose in mg
    typical: Optional[float] = None  # Typical dose in mg
    unit: str = "mg"  # Dose unit
    route: Optional[str] = None  # Administration route
    notes: Optional[str] = None  # Additional notes


@dataclass
class EffectScore:
    """Effect intensity score."""

    value: float  # Score between 0-1
    confidence: Optional[float] = None  # Confidence between 0-1
    supporting_data: Optional[Dict[str, Any]] = None  # Evidence/sources
    notes: Optional[str] = None  # Additional notes


@dataclass
class RiskScore:
    """Risk assessment score."""

    value: float  # Score between 0-1
    confidence: Optional[float] = None  # Confidence between 0-1
    level: Optional[RiskLevel] = None  # Risk level category
    supporting_data: Optional[Dict[str, Any]] = None  # Evidence/sources
    notes: Optional[str] = None  # Additional notes


@dataclass
class ReceptorBinding:
    """Receptor binding data."""

    affinity: Optional[float] = None  # Binding affinity value
    affinity_type: Optional[str] = None  # Ki, IC50, EC50, etc.
    affinity_unit: Optional[str] = None  # nM, uM, etc.
    confidence: Optional[float] = None  # Confidence between 0-1
    supporting_data: Optional[Dict[str, Any]] = None  # Evidence/sources
    notes: Optional[str] = None  # Additional notes


@dataclass
class PredictionResult:
    """ML prediction result."""

    value: float  # Predicted value
    confidence: Optional[float] = None  # Confidence between 0-1
    probability: Optional[float] = None  # Probability between 0-1
    model_name: Optional[str] = None  # Name of model used
    model_version: Optional[str] = None  # Version of model
    features_used: Optional[List[str]] = None  # Features used in prediction
    metadata: Optional[Dict[str, Any]] = None  # Additional metadata


# Type aliases
CompoundID = str  # CAS number or other unique identifier
SMILES = str  # SMILES string
InChI = str  # InChI string
PubChemID = str  # PubChem CID
ChEMBLID = str  # ChEMBL ID

# Complex type aliases
Properties = Dict[str, Union[str, float, int, bool]]
Predictions = Dict[str, PredictionResult]
WebData = Dict[str, Any]
AnalysisData = Dict[str, Any]

# Social media data types
SocialPost = Dict[str, Any]  # Individual social media post
SocialData = Dict[str, List[SocialPost]]  # Platform -> Posts

# Literature data types
Paper = Dict[str, Any]  # Scientific paper data
Citation = Dict[str, Any]  # Citation data
LiteratureData = Dict[str, Union[Set[str], Dict[str, str], Dict[str, List[str]]]]

# Patent data types
Patent = Dict[str, Any]  # Patent data
PatentData = Dict[str, Union[Set[str], Dict[str, str], Dict[str, List[str]]]]

# Community data types
CommunityReport = Dict[str, Any]  # Community report data
CommunityData = Dict[str, Dict[str, Union[str, Dict[str, Any]]]]
