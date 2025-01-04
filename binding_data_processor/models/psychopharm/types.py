"""Consolidated type definitions for compound data models.

This module provides comprehensive type definitions used across the compound data models:

Core Types:
- Basic compound classifications (CompoundType)
- Legal and regulatory status (LegalStatus)
- Psychoactive classifications (PsychoactiveClass)
- Nootropic mechanisms (NootropicMechanism)
- BBB permeability levels (BBBPermeability)

Activity Types:
- Binding interactions (BindingType)
- Pharmacological activity (ActivityType)
- Risk classifications (RiskLevel)

Data Structures:
- Target binding data (TargetData)
- ML feature types
- Prediction result types
- Collection type aliases
- Export format types
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set, Union, Any, Tuple
import numpy as np


class CompoundType(Enum):
    """Types of chemical compounds."""
    # Core types
    NEUROTRANSMITTER = "neurotransmitter"
    PSYCHOACTIVE = "psychoactive"
    RESEARCH_CHEMICAL = "research_chemical"
    NPS = "novel_psychoactive_substance"
    PHARMACEUTICAL = "pharmaceutical"
    NATURAL_PRODUCT = "natural_product"
    PEPTIDE = "peptide"
    PROTEIN = "protein"
    ANTIBODY = "antibody"
    SMALL_MOLECULE = "small_molecule"
    METABOLITE = "metabolite"
    PRODRUG = "prodrug"
    OTHER = "other"


class LegalStatus(Enum):
    """Legal status classifications."""
    LEGAL = "legal"
    CONTROLLED = "controlled"
    ILLEGAL = "illegal"
    RESEARCH_ONLY = "research_only"
    UNSCHEDULED = "unscheduled"
    PRESCRIPTION = "prescription_only"
    OTC = "over_the_counter"
    INVESTIGATIONAL = "investigational"
    WITHDRAWN = "withdrawn"


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
    OPIOID = "opioid"
    CANNABINOID = "cannabinoid"
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
    MITOCHONDRIAL = "mitochondrial"
    UNKNOWN = "unknown"


class BBBPermeability(Enum):
    """Blood-brain barrier permeability classification."""
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    NEGLIGIBLE = "negligible"
    SUBSTRATE = "transporter_substrate"
    UNKNOWN = "unknown"


class BindingType(Enum):
    """Types of binding interactions."""
    AGONIST = "agonist"
    ANTAGONIST = "antagonist"
    PARTIAL_AGONIST = "partial_agonist"
    INVERSE_AGONIST = "inverse_agonist"
    ALLOSTERIC = "allosteric"
    POSITIVE_MODULATOR = "positive_modulator"
    NEGATIVE_MODULATOR = "negative_modulator"
    REUPTAKE_INHIBITOR = "reuptake_inhibitor"
    RELEASING_AGENT = "releasing_agent"
    ENZYME_INHIBITOR = "enzyme_inhibitor"
    UNKNOWN = "unknown"


class ActivityType(Enum):
    """Types of pharmacological activity."""
    STIMULANT = "stimulant"
    DEPRESSANT = "depressant"
    PSYCHEDELIC = "psychedelic"
    DISSOCIATIVE = "dissociative"
    NOOTROPIC = "nootropic"
    ANXIOLYTIC = "anxiolytic"
    ANTIPSYCHOTIC = "antipsychotic"
    ANTIDEPRESSANT = "antidepressant"
    MOOD_STABILIZER = "mood_stabilizer"
    UNKNOWN = "unknown"


class RiskLevel(Enum):
    """Risk level classifications."""
    NONE = "none"
    MINIMAL = "minimal"
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    SEVERE = "severe"
    UNKNOWN = "unknown"


@dataclass
class TargetData:
    """Structured data for a single target interaction."""
    # Target identification
    common_name: str = "N/A"
    protein_name: str = "N/A"
    gene_name: str = "N/A"
    organism: str = "human"
    
    # Binding data
    affinity_value: float = 0.0
    affinity_type: str = "N/A"  # Ki, IC50, Kd, EC50
    affinity_unit: str = "nM"
    activity_type: str = "N/A"  # agonist, antagonist, etc.
    confidence: float = 0.0
    
    # Metadata
    is_primary: bool = False
    pubmed_count: int = 0
    reference_dois: Set[str] = field(default_factory=set)
    assay_details: Dict = field(default_factory=dict)
    experimental_conditions: Dict = field(default_factory=dict)
    
    # Additional fields
    binding_type: BindingType = BindingType.UNKNOWN
    activity_level: float = 0.0  # 0.0-1.0 scale
    effect_confidence: float = 0.0  # 0.0-1.0 scale


# Value Types
AffinityValue = float  # Binding affinity measurements (Ki, IC50, etc.)
ConfidenceValue = float  # Confidence scores (0.0-1.0)
ProbabilityValue = float  # Probability values (0.0-1.0)
Features = np.ndarray  # ML feature vectors
ModelVersion = str  # Model version identifiers
Timestamp = str  # ISO format timestamps

# Structured Types
DoseRange = Tuple[float, float, float]  # min, max, recommended
TimeRange = Tuple[float, float]  # start, end
RiskScore = Tuple[float, float]  # severity, confidence
EffectScore = Tuple[float, float]  # magnitude, confidence
ReceptorBinding = Tuple[float, float, str]  # affinity, confidence, activity

# ML Types
PredictionResult = Dict[str, Union[float, str, Dict]]  # ML prediction results
ModelMetrics = Dict[str, float]  # Model performance metrics
FeatureImportances = Dict[str, float]  # Feature importance scores

# Collection Types
StringSet = Set[str]  # Set of strings
StringDict = Dict[str, str]  # String key-value pairs
MetricsDict = Dict[str, float]  # Metric key-value pairs
PredictionDict = Dict[str, PredictionResult]  # Prediction results by type
ValidationErrors = List[str]  # Validation error messages

# Optional Types
OptionalStr = Optional[str]  # Optional string values
OptionalFloat = Optional[float]  # Optional float values
OptionalDict = Optional[Dict]  # Optional dictionary values
OptionalList = Optional[List]  # Optional list values

# Export Types
ExportDict = Dict[str, Union[str, float, int, List, Dict]]  # Export data format

# Specialized Types
TargetDict = Dict[str, TargetData]  # Target data by name
ActivityDict = Dict[str, ActivityType]  # Activity types by target
RiskDict = Dict[str, RiskLevel]  # Risk levels by category
MechanismDict = Dict[str, NootropicMechanism]  # Mechanisms by type
