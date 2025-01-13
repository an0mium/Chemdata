"""Psychopharmacology-specific compound model."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

from .analysis import AnalyzedCompound
from .base.types import (
    BBBPermeability,
    NootropicMechanism,
    PsychoactiveClass,
    EffectScore,
    RiskScore,
    ReceptorBinding,
    TimeRange,
    DoseRange,
)


@dataclass
class PsychopharmCompound(AnalyzedCompound):
    """Compound model with psychopharmacology-specific capabilities."""

    # Core psychoactive properties
    psychoactive_class: Optional[PsychoactiveClass] = None
    mechanism_of_action: Optional[NootropicMechanism] = None
    bbb_permeability: Optional[BBBPermeability] = None

    # Effects and timing
    onset: Optional[TimeRange] = None
    duration: Optional[TimeRange] = None
    dose_range: Optional[DoseRange] = None
    cognitive_effects: Dict[str, EffectScore] = field(default_factory=dict)
    side_effects: Dict[str, RiskScore] = field(default_factory=dict)

    # Receptor data
    receptor_affinities: Dict[str, ReceptorBinding] = field(default_factory=dict)
    receptor_activities: Dict[str, float] = field(default_factory=dict)

    # Clinical data
    clinical_notes: List[str] = field(default_factory=list)
    contraindications: Set[str] = field(default_factory=set)
    interactions: Dict[str, str] = field(default_factory=dict)

    # Literature and community data
    literature_references: List[str] = field(default_factory=list)
    community_reports: List[Dict] = field(default_factory=list)

    def __post_init__(self):
        """Validate required fields after initialization."""
        super().__post_init__()

        # Validate psychoactive class if provided
        if self.psychoactive_class and not isinstance(self.psychoactive_class, PsychoactiveClass):
            self.psychoactive_class = PsychoactiveClass(self.psychoactive_class)

        # Validate mechanism if provided
        if self.mechanism_of_action and not isinstance(self.mechanism_of_action, NootropicMechanism):
            self.mechanism_of_action = NootropicMechanism(self.mechanism_of_action)

        # Validate BBB permeability if provided
        if self.bbb_permeability and not isinstance(self.bbb_permeability, BBBPermeability):
            self.bbb_permeability = BBBPermeability(self.bbb_permeability)

    def get_receptor_profile(self) -> Dict[str, ReceptorBinding]:
        """Get complete receptor binding profile."""
        return self.receptor_affinities

    def get_cognitive_effects(self) -> Dict[str, EffectScore]:
        """Get cognitive effects with scores."""
        return self.cognitive_effects

    def get_safety_profile(self) -> Dict[str, RiskScore]:
        """Get safety profile with risk scores."""
        return self.side_effects

    def get_clinical_data(self) -> Dict:
        """Get clinical data including notes and contraindications."""
        return {"notes": self.clinical_notes, "contraindications": list(self.contraindications), "interactions": self.interactions}

    def get_literature_data(self) -> List[str]:
        """Get literature references."""
        return self.literature_references

    def get_community_data(self) -> List[Dict]:
        """Get community reports and experiences."""
        return self.community_reports
