"""Analysis functionality for compound data.

This module extends EnrichedCompoundData with analysis capabilities:
- Binding analysis
- Activity analysis
- Safety analysis
- SAR analysis
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from enum import Enum

from .compound_enrichment import EnrichedCompoundData


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


@dataclass
class AnalyzedCompoundData(EnrichedCompoundData):
    """CompoundData with analysis capabilities."""

    # Binding analysis
    binding_profiles: Dict[str, Dict] = field(default_factory=dict)
    binding_types: Dict[str, BindingType] = field(default_factory=dict)
    binding_affinities: Dict[str, float] = field(default_factory=dict)
    binding_confidences: Dict[str, float] = field(default_factory=dict)

    # Activity analysis
    activity_types: Set[ActivityType] = field(default_factory=set)
    activity_scores: Dict[str, float] = field(default_factory=dict)
    activity_confidences: Dict[str, float] = field(default_factory=dict)
    activity_mechanisms: Dict[str, List[str]] = field(default_factory=dict)

    # Safety analysis
    toxicity_alerts: List[Dict] = field(default_factory=list)
    safety_scores: Dict[str, float] = field(default_factory=dict)
    safety_warnings: List[str] = field(default_factory=list)
    contraindications: List[str] = field(default_factory=list)
    drug_interactions: List[Dict] = field(default_factory=list)

    # SAR analysis
    pharmacophores: List[Dict] = field(default_factory=list)
    structural_alerts: List[Dict] = field(default_factory=list)
    receptor_interactions: Dict[str, List[Dict]] = field(default_factory=dict)
    mechanism_predictions: List[Dict] = field(default_factory=list)
    sar_analysis: Dict = field(default_factory=dict)

    def analyze_binding(self) -> Dict:
        """Analyze binding data and generate summary."""
        summary = {
            "total_targets": len(self.binding_profiles),
            "strongest_binding": None,
            "primary_targets": [],
            "target_families": set(),
        }

        # Analyze each binding profile
        for target, profile in self.binding_profiles.items():
            # Track strongest binding
            affinity = self.binding_affinities.get(target)
            if affinity:
                if (not summary["strongest_binding"] or 
                    affinity < summary["strongest_binding"]["value"]):
                    summary["strongest_binding"] = {
                        "target": target,
                        "value": affinity,
                        "type": self.binding_types.get(target, BindingType.UNKNOWN).value,
                        "confidence": self.binding_confidences.get(target, 0.0),
                    }

            # Track primary targets
            if profile.get("is_primary"):
                summary["primary_targets"].append(target)

            # Track target families
            if "family" in profile:
                summary["target_families"].add(profile["family"])

        summary["target_families"] = list(summary["target_families"])
        return summary

    def analyze_activity(self) -> Dict:
        """Analyze activity data and generate summary."""
        summary = {
            "primary_type": None,
            "secondary_types": [],
            "mechanisms": [],
            "scores": {},
        }

        # Find primary activity type
        if self.activity_types:
            max_score = 0
            for activity_type in self.activity_types:
                score = self.activity_scores.get(activity_type.value, 0)
                if score > max_score:
                    max_score = score
                    summary["primary_type"] = activity_type.value
                elif score > 0:
                    summary["secondary_types"].append(activity_type.value)

        # Collect mechanisms
        for activity_type in self.activity_types:
            mechanisms = self.activity_mechanisms.get(activity_type.value, [])
            summary["mechanisms"].extend(mechanisms)

        # Collect scores
        for activity_type in self.activity_types:
            type_value = activity_type.value
            summary["scores"][type_value] = {
                "score": self.activity_scores.get(type_value, 0),
                "confidence": self.activity_confidences.get(type_value, 0),
            }

        return summary

    def analyze_safety(self) -> Dict:
        """Analyze safety data and generate summary."""
        summary = {
            "alerts": [],
            "warnings": [],
            "contraindications": [],
            "interactions": [],
            "risk_factors": [],
            "scores": {},
        }

        # Add structural alerts
        summary["alerts"].extend(
            alert["description"] for alert in self.structural_alerts
        )

        # Add toxicity alerts
        for alert in self.toxicity_alerts:
            if alert.get("severity", 0) > 0.7:  # High severity threshold
                summary["warnings"].append({
                    "type": alert["type"],
                    "severity": alert["severity"],
                    "description": alert["description"],
                })

        # Add contraindications
        summary["contraindications"].extend(self.contraindications)

        # Add drug interactions
        summary["interactions"].extend(self.drug_interactions)

        # Add risk factors
        summary["risk_factors"].extend(self.risk_factors)

        # Add safety scores
        summary["scores"] = self.safety_scores.copy()

        return summary

    def analyze_sar(self) -> Dict:
        """Analyze structure-activity relationships."""
        summary = {
            "pharmacophores": [],
            "key_fragments": [],
            "receptor_interactions": {},
            "mechanism_predictions": [],
        }

        # Add pharmacophores
        summary["pharmacophores"].extend(
            {
                "type": pharm["type"],
                "description": pharm["description"],
                "score": pharm["score"],
            }
            for pharm in self.pharmacophores
        )

        # Add key fragments from SAR analysis
        if self.sar_analysis:
            summary["key_fragments"].extend(
                self.sar_analysis.get("key_fragments", [])
            )

        # Add receptor interactions
        for receptor, interactions in self.receptor_interactions.items():
            summary["receptor_interactions"][receptor] = [
                {
                    "type": inter["type"],
                    "strength": inter["strength"],
                    "confidence": inter["confidence"],
                }
                for inter in interactions
            ]

        # Add mechanism predictions
        summary["mechanism_predictions"].extend(
            {
                "mechanism": pred["mechanism"],
                "probability": pred["probability"],
                "confidence": pred["confidence"],
            }
            for pred in self.mechanism_predictions
        )

        return summary

    def get_analysis_dict(self) -> Dict:
        """Get dictionary of all analysis results."""
        return {
            "binding_analysis": self.analyze_binding(),
            "activity_analysis": self.analyze_activity(),
            "safety_analysis": self.analyze_safety(),
            "sar_analysis": self.analyze_sar(),
        }

    def to_dict(self, include_predictions: bool = True) -> Dict:
        """Convert compound data to dictionary format."""
        data = super().to_dict(include_predictions=include_predictions)
        data.update(self.get_analysis_dict())
        return data
