"""Base types and data structures for psychopharmacological properties."""

from dataclasses import field
from enum import Enum
from typing import Dict, List, Optional, Set, Tuple

# Type aliases for structured data
DoseRange = Tuple[float, float, float]  # min, max, recommended
TimeRange = Tuple[float, float]  # start, end
RiskScore = Tuple[float, float]  # severity, confidence
EffectScore = Tuple[float, float]  # magnitude, confidence
ReceptorBinding = Tuple[float, float, str]  # affinity, confidence, activity


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


class RiskLevel(Enum):
    """Risk level classification."""
    
    SEVERE = "severe"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    MINIMAL = "minimal"
    UNKNOWN = "unknown"


class BaseMixin:
    """Base mixin providing common functionality."""

    # Basic properties
    name: str = field(default="")
    smiles: Optional[str] = field(default=None)
    cas_number: Optional[str] = field(default=None)
    
    # Metadata
    data_sources: Set[str] = field(default_factory=set)
    last_updated: Optional[str] = field(default=None)
    confidence_scores: Dict[str, float] = field(default_factory=dict)
    # List of update records containing:
    # - timestamp: When the update occurred
    # - source: Where the update came from
    # - action: What type of update was performed
    update_history: List[Dict[str, str]] = field(default_factory=list)
    
    def get_base_dict(self) -> Dict:
        """Get dictionary of base properties."""
        return {
            "name": self.name,
            "smiles": self.smiles,
            "cas_number": self.cas_number,
            "data_sources": list(self.data_sources),
            "last_updated": self.last_updated,
            "confidence_scores": self.confidence_scores,
        }

    def merge_base_data(self, other: "BaseMixin", source: str = "merge") -> None:
        """Merge base data from another instance."""
        self.data_sources.update(other.data_sources)
        self.confidence_scores.update(other.confidence_scores)
        
        # Take most recent update
        if (
            other.last_updated and 
            (not self.last_updated or other.last_updated > self.last_updated)
        ):
            self.last_updated = other.last_updated
            self.update_history.append({
                "timestamp": self.last_updated,
                "source": source,
                "action": "merge",
            })
