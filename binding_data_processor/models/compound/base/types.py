"""Type definitions for compound data models.

This module defines:
1. Core enums for compound classification and binding
2. Data structures for timing, scoring and binding data
3. Type aliases for common data types
4. Validation and data source types
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set, Union, Any


class CompoundType(str, Enum):
    """Types of chemical compounds."""

    NEUROTRANSMITTER = "neurotransmitter"
    PSYCHOACTIVE = "psychoactive"
    RESEARCH_CHEMICAL = "research_chemical"
    NPS = "novel_psychoactive_substance"
    PHARMACEUTICAL = "pharmaceutical"
    NATURAL_PRODUCT = "natural_product"
    SMALL_MOLECULE = "small_molecule"
    PEPTIDE = "peptide"
    PROTEIN = "protein"
    SYNTHETIC = "synthetic"
    NOOTROPIC = "nootropic"
    METABOLITE = "metabolite"
    OTHER = "other"
    UNKNOWN = "unknown"


class BindingType(str, Enum):
    """Type of binding interaction."""

    AGONIST = "agonist"
    ANTAGONIST = "antagonist"
    PARTIAL_AGONIST = "partial_agonist"
    INVERSE_AGONIST = "inverse_agonist"
    ALLOSTERIC = "allosteric"
    COMPETITIVE = "competitive"
    NONCOMPETITIVE = "noncompetitive"
    IRREVERSIBLE = "irreversible"
    UNKNOWN = "unknown"


class ActivityType(str, Enum):
    """Type of biological activity."""

    INHIBITOR = "inhibitor"
    ACTIVATOR = "activator"
    MODULATOR = "modulator"
    BLOCKER = "blocker"
    SUBSTRATE = "substrate"
    INDUCER = "inducer"
    UNKNOWN = "unknown"


class LegalStatus(str, Enum):
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


class PsychoactiveClass(str, Enum):
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


class NootropicMechanism(str, Enum):
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


class BBBPermeability(str, Enum):
    """Blood-brain barrier permeability classification."""

    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    NEGLIGIBLE = "negligible"
    UNKNOWN = "unknown"


class RiskLevel(str, Enum):
    """Risk assessment level."""

    SEVERE = "severe"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    MINIMAL = "minimal"
    UNKNOWN = "unknown"


class DataSource(str, Enum):
    """Source of compound data."""

    BINDINGDB = "bindingdb"
    PUBCHEM = "pubchem"
    CHEMBL = "chembl"
    CUSTOM = "custom"
    LITERATURE = "literature"
    PATENT = "patent"
    PREDICTED = "predicted"
    UNKNOWN = "unknown"


class ValidationStatus(str, Enum):
    """Data validation status."""

    VALIDATED = "validated"
    PARTIALLY_VALIDATED = "partially_validated"
    UNVALIDATED = "unvalidated"
    FAILED = "failed"
    PENDING = "pending"


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
class BindingData:
    """Container for binding interaction data."""

    target: str
    binding_type: BindingType
    affinity: Optional[float] = None
    units: Optional[str] = None
    confidence: Optional[float] = None
    method: Optional[str] = None
    conditions: Optional[Dict] = field(default_factory=dict)
    references: List[str] = field(default_factory=list)
    source: DataSource = DataSource.UNKNOWN
    validation_status: ValidationStatus = ValidationStatus.UNVALIDATED


@dataclass
class ActivityData:
    """Container for biological activity data."""

    activity_type: ActivityType
    target: str
    value: Optional[float] = None
    units: Optional[str] = None
    confidence: Optional[float] = None
    method: Optional[str] = None
    conditions: Optional[Dict] = field(default_factory=dict)
    references: List[str] = field(default_factory=list)
    source: DataSource = DataSource.UNKNOWN
    validation_status: ValidationStatus = ValidationStatus.UNVALIDATED


@dataclass
class PropertyData:
    """Container for physicochemical property data."""

    name: str
    value: Union[float, str, bool]
    units: Optional[str] = None
    confidence: Optional[float] = None
    method: Optional[str] = None
    conditions: Optional[Dict] = field(default_factory=dict)
    references: List[str] = field(default_factory=list)
    source: DataSource = DataSource.UNKNOWN
    validation_status: ValidationStatus = ValidationStatus.UNVALIDATED


@dataclass
class SafetyData:
    """Container for safety/toxicity data."""

    endpoint: str
    value: Union[float, str, bool]
    severity: Optional[str] = None
    confidence: Optional[float] = None
    method: Optional[str] = None
    conditions: Optional[Dict] = field(default_factory=dict)
    references: List[str] = field(default_factory=list)
    source: DataSource = DataSource.UNKNOWN
    validation_status: ValidationStatus = ValidationStatus.UNVALIDATED


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


@dataclass
class TargetData:
    """Container for target-specific data."""

    name: str  # Target name (e.g. receptor, enzyme)
    type: str  # Target type (e.g. GPCR, ion channel)
    organism: str = "human"  # Source organism
    uniprot_id: Optional[str] = None  # UniProt identifier
    gene_name: Optional[str] = None  # Gene name
    binding_data: Optional[BindingData] = None  # Binding interaction data
    activity_data: Optional[ActivityData] = None  # Activity measurements
    references: List[str] = field(default_factory=list)  # Literature references
    metadata: Dict[str, Any] = field(default_factory=dict)  # Additional metadata
    source: DataSource = DataSource.UNKNOWN
    validation_status: ValidationStatus = ValidationStatus.UNVALIDATED


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
