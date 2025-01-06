"""Type definitions for compound data models."""

from enum import Enum


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
