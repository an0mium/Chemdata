"""Type definitions for compound data models.

This module provides enums and type definitions used across the compound data models:
- Basic compound classifications
- Legal status types
- Psychoactive classifications
- Binding and activity types
- Safety and risk levels
"""

from enum import Enum
from typing import Dict, List, Optional, Set, Union, Any
from dataclasses import dataclass


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
    reference_dois: Set[str] = None
    assay_details: Dict = None
    experimental_conditions: Dict = None

    def __post_init__(self):
        """Initialize collections."""
        if self.reference_dois is None:
            self.reference_dois = set()
        if self.assay_details is None:
            self.assay_details = {}
        if self.experimental_conditions is None:
            self.experimental_conditions = {}
