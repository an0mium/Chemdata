"""Enhanced compound model combining all capabilities.

This module provides the EnhancedCompound class that combines:
1. Core compound data (CompoundData)
2. Common functionality (BaseMixin)
3. Web enrichment (WebEnrichmentMixin)
4. ML predictions (PredictionsMixin)
5. Analysis capabilities (AnalysisMixin)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Union

from .base import CompoundData
from .types import (
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
    RiskLevel,
    DoseRange,
    TimeRange,
    RiskScore,
    EffectScore,
    ReceptorBinding,
)
from ..mixins import (
    BaseMixin,
    WebEnrichmentMixin,
    PredictionsMixin,
    AnalysisMixin,
)


@dataclass
class EnhancedCompound(
    CompoundData,
    BaseMixin,
    WebEnrichmentMixin,
    PredictionsMixin,
    AnalysisMixin,
):
    """Enhanced compound with all capabilities."""

    # Psychopharmacology
    psychoactive_class: PsychoactiveClass = field(default=PsychoactiveClass.UNKNOWN)
    nootropic_mechanisms: Set[NootropicMechanism] = field(default_factory=set)
    bbb_permeability: BBBPermeability = field(default=BBBPermeability.UNKNOWN)

    # Dosage and timing
    dosage_ranges: Dict[str, DoseRange] = field(default_factory=dict)  # Route -> Range
    onset_times: Dict[str, TimeRange] = field(default_factory=dict)  # Route -> Time
    duration_times: Dict[str, TimeRange] = field(default_factory=dict)  # Route -> Time

    # Effects and risks
    effect_scores: Dict[str, EffectScore] = field(default_factory=dict)  # Effect -> Score
    risk_scores: Dict[str, RiskScore] = field(default_factory=dict)  # Risk -> Score
    receptor_bindings: Dict[str, ReceptorBinding] = field(default_factory=dict)

    def __post_init__(self):
        """Initialize all parent classes."""
        CompoundData.__post_init__(self)
        BaseMixin.__init__(self)
        WebEnrichmentMixin.__init__(self)
        PredictionsMixin.__init__(self)
        AnalysisMixin.__init__(self)

    def to_dict(self) -> Dict:
        """Convert to dictionary format."""
        data = super().to_dict()
        data.update({
            # Psychopharmacology
            "psychoactive_class": self.psychoactive_class.value,
            "nootropic_mechanisms": [m.value for m in self.nootropic_mechanisms],
            "bbb_permeability": self.bbb_permeability.value,
            # Dosage and timing
            "dosage_ranges": self.dosage_ranges,
            "onset_times": self.onset_times,
            "duration_times": self.duration_times,
            # Effects and risks
            "effect_scores": self.effect_scores,
            "risk_scores": self.risk_scores,
            "receptor_bindings": self.receptor_bindings,
        })
        return data

    def merge(self, other: "EnhancedCompound") -> None:
        """Merge data from another compound."""
        # Merge base data
        self.merge_base_data(other)

        # Merge web data
        self.merge_web_data(other.get_web_data())

        # Merge predictions
        for pred_type, pred in other._prediction_cache.items():
            if pred_type not in self._prediction_cache:
                self.cache_prediction(
                    pred_type,
                    pred,
                    other._prediction_history[pred_type]["confidence"],
                    other._prediction_history[pred_type]["supporting_data"],
                )

        # Merge analysis results
        if hasattr(other, "_binding_analysis"):
            if not hasattr(self, "_binding_analysis"):
                self._binding_analysis = {}
            self._binding_analysis.update(other._binding_analysis)

        if hasattr(other, "_activity_analysis"):
            if not hasattr(self, "_activity_analysis"):
                self._activity_analysis = {}
            self._activity_analysis.update(other._activity_analysis)

        if hasattr(other, "_safety_analysis"):
            if not hasattr(self, "_safety_analysis"):
                self._safety_analysis = {}
            self._safety_analysis.update(other._safety_analysis)

    def get_summary(self) -> Dict:
        """Get high-level summary of compound data."""
        return {
            # Basic info
            "name": self.name,
            "cas_number": self.cas_number,
            "psychoactive_class": self.psychoactive_class.value,
            # Key properties
            "bbb_permeability": self.bbb_permeability.value,
            "primary_mechanism": list(self.nootropic_mechanisms)[0].value
            if self.nootropic_mechanisms else None,
            # Analysis results
            "binding_summary": self.get_analysis_summary().get("binding", {}),
            "safety_summary": self.get_analysis_summary().get("safety", {}),
            # Predictions
            "predictions": {
                k: v for k, v in self._prediction_cache.items()
                if v is not None
            },
            # Web data
            "community_data": bool(self.community_data),
            "literature_data": bool(self.literature_data),
            # Metadata
            "data_sources": list(self.data_sources),
            "last_updated": self.last_updated,
        }

    def validate(self) -> List[str]:
        """Validate all compound data."""
        errors = []
        
        # Validate base data
        try:
            CompoundData._validate(self)
        except Exception as e:
            errors.append(f"Base validation error: {str(e)}")

        # Validate predictions
        if self._prediction_cache:
            for pred_type, pred in self._prediction_cache.items():
                if pred is None:
                    errors.append(f"Missing prediction for {pred_type}")

        # Validate analysis results
        if hasattr(self, "_binding_analysis"):
            if not self._binding_analysis.get("strongest_binding"):
                errors.append("Missing strongest binding analysis")

        if hasattr(self, "_safety_analysis"):
            if not self._safety_analysis.get("risk_assessment"):
                errors.append("Missing risk assessment")

        return errors
