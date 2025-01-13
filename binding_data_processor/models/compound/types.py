"""Type definitions for compound data models.

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
from typing import Dict, List, Optional, Set, Union, Any, Tuple, Protocol
import numpy as np


class CompoundType(Enum):
    """Types of chemical compounds."""

    NEUROTRANSMITTER = "neurotransmitter"
    PSYCHOACTIVE = "psychoactive"
    RESEARCH_CHEMICAL = "research_chemical"
    NPS = "novel_psychoactive_substance"
    PHARMACEUTICAL = "pharmaceutical"
    NATURAL_PRODUCT = "natural_product"
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
    UNKNOWN = "unknown"


class BBBPermeability(Enum):
    """Blood-brain barrier permeability classification."""

    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    NEGLIGIBLE = "negligible"
    UNKNOWN = "unknown"


class BindingType(Enum):
    """Types of binding interactions."""

    AGONIST = "agonist"
    ANTAGONIST = "antagonist"
    PARTIAL_AGONIST = "partial_agonist"
    INVERSE_AGONIST = "inverse_agonist"
    ALLOSTERIC = "allosteric"
    UNKNOWN = "unknown"


class MechanismType(Enum):
    """Types of binding mechanisms."""

    AGONIST = "agonist"
    ANTAGONIST = "antagonist"
    PARTIAL_AGONIST = "partial_agonist"
    INVERSE_AGONIST = "inverse_agonist"
    ALLOSTERIC = "allosteric"
    UNKNOWN = "unknown"


class SelectivityType(Enum):
    """Types of target selectivity."""

    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    UNKNOWN = "unknown"


class AbusePotential(Enum):
    """Classification of abuse potential."""

    NONE = "none"
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    SEVERE = "severe"
    UNKNOWN = "unknown"


class TolerancePattern(Enum):
    """Classification of tolerance development patterns."""

    NONE = "none"
    ACUTE = "acute"
    CHRONIC = "chronic"
    BEHAVIORAL = "behavioral"
    MIXED = "mixed"
    UNKNOWN = "unknown"


class WithdrawalSeverity(Enum):
    """Classification of withdrawal severity."""

    NONE = "none"
    MILD = "mild"
    MODERATE = "moderate"
    SEVERE = "severe"
    LIFE_THREATENING = "life_threatening"
    UNKNOWN = "unknown"


class ToxicityClass(Enum):
    """Classification of toxicity mechanisms."""

    OXIDATIVE_STRESS = "oxidative_stress"
    MITOCHONDRIAL_TOXICITY = "mitochondrial_toxicity"
    DNA_DAMAGE = "dna_damage"
    PROTEIN_ADDUCTS = "protein_adducts"
    LIPID_PEROXIDATION = "lipid_peroxidation"
    MEMBRANE_DISRUPTION = "membrane_disruption"
    ENZYME_INHIBITION = "enzyme_inhibition"
    RECEPTOR_OVERSTIMULATION = "receptor_overstimulation"
    UNKNOWN = "unknown"


class ActivityType(Enum):
    """Types of pharmacological activity."""

    STIMULANT = "stimulant"
    DEPRESSANT = "depressant"
    PSYCHEDELIC = "psychedelic"
    DISSOCIATIVE = "dissociative"
    NOOTROPIC = "nootropic"
    UNKNOWN = "unknown"


class RiskLevel(Enum):
    """Risk level classifications."""

    NONE = "none"
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    SEVERE = "severe"
    UNKNOWN = "unknown"


class PropertyType(Enum):
    """Types of compound properties."""

    PHYSICOCHEMICAL = "physicochemical"
    STRUCTURAL = "structural"
    ELECTRONIC = "electronic"
    TOPOLOGICAL = "topological"
    GEOMETRIC = "geometric"
    QUANTUM = "quantum"
    UNKNOWN = "unknown"

    @classmethod
    def get_type(cls, property_name: str) -> "PropertyType":
        """Get property type from property name."""
        mapping = {
            "molecular_weight": cls.PHYSICOCHEMICAL,
            "logp": cls.PHYSICOCHEMICAL,
            "hbd": cls.PHYSICOCHEMICAL,
            "hba": cls.PHYSICOCHEMICAL,
            "tpsa": cls.PHYSICOCHEMICAL,
            "rotatable_bonds": cls.STRUCTURAL,
            "charge": cls.ELECTRONIC,
            "stereocenter_count": cls.STRUCTURAL,
            "ring_count": cls.STRUCTURAL,
            "atom_count": cls.STRUCTURAL,
        }
        return mapping.get(property_name, cls.UNKNOWN)


class DrugLikenessType(Enum):
    """Types of drug-likeness rules."""

    LIPINSKI = "lipinski"
    VEBER = "veber"
    GHOSE = "ghose"
    MUEGGE = "muegge"
    EGAN = "egan"
    QED = "qed"
    CUSTOM = "custom"
    UNKNOWN = "unknown"

    @classmethod
    def get_type(cls, rule_name: str) -> "DrugLikenessType":
        """Get drug-likeness type from rule name."""
        mapping = {
            "ro5": cls.LIPINSKI,
            "lipinski": cls.LIPINSKI,
            "veber": cls.VEBER,
            "ghose": cls.GHOSE,
            "muegge": cls.MUEGGE,
            "egan": cls.EGAN,
            "qed": cls.QED,
        }
        return mapping.get(rule_name.lower(), cls.UNKNOWN)


class ADMEType(Enum):
    """Types of ADME properties."""

    ABSORPTION = "absorption"
    DISTRIBUTION = "distribution"
    METABOLISM = "metabolism"
    EXCRETION = "excretion"
    BIOAVAILABILITY = "bioavailability"
    PERMEABILITY = "permeability"
    CLEARANCE = "clearance"
    HALF_LIFE = "half_life"
    PLASMA_BINDING = "plasma_binding"
    VOLUME_DISTRIBUTION = "volume_distribution"
    UNKNOWN = "unknown"


@dataclass
class TargetData:
    """Structured data for a single target interaction."""

    common_name: str = "N/A"
    protein_name: str = "N/A"
    gene_name: str = "N/A"
    organism: str = "human"
    affinity_value: float = 0.0
    affinity_type: str = "N/A"  # Ki, IC50, Kd, EC50
    affinity_unit: str = "nM"
    activity_type: str = "N/A"  # agonist, antagonist, etc.
    confidence: float = 0.0
    is_primary: bool = False
    pubmed_count: int = 0
    reference_dois: Set[str] = field(default_factory=set)
    assay_details: Dict = field(default_factory=dict)
    experimental_conditions: Dict = field(default_factory=dict)


class CompoundData(Protocol):
    """Protocol defining the interface for compound data."""

    smiles: str
    name: str
    primary_activity: Optional[float]
    reference_compounds: List["CompoundData"]


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
