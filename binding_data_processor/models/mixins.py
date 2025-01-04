"""Mixins for adding specific functionality to compound data models.

This module provides mixins that can be used to add specific functionality
to custom compound data implementations:

WebEnrichmentMixin:
    Web data enrichment capabilities

PredictionsMixin:
    ML prediction integration

AnalysisMixin:
    Analysis capabilities
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set
import json
import numpy as np
import pandas as pd

from .compound.types import (
    BindingType,
    ActivityType,
    TargetData,
)


@dataclass
class WebEnrichmentMixin:
    """Mixin for web data enrichment capabilities."""
    
    # Patent data
    patent_data: Dict = field(default_factory=dict)
    patent_count: int = 0

    # Web-enriched data
    community_data: Dict = field(default_factory=dict)
    literature_data: Dict = field(default_factory=dict)
    regulatory_data: Dict = field(default_factory=dict)
    experience_reports: List[Dict] = field(default_factory=list)
    safety_profile: Dict = field(default_factory=dict)
    dosage_info: Dict[str, Dict] = field(default_factory=dict)
    route_stats: Dict[str, int] = field(default_factory=dict)
    duration_stats: Dict[str, int] = field(default_factory=dict)

    def merge_web_data(self, other_data: Dict) -> None:
        """Merge web-enriched data."""
        if "patent_data" in other_data:
            self.patent_data.update(other_data["patent_data"])
            self.patent_count = max(self.patent_count, other_data.get("patent_count", 0))

        if "community_data" in other_data:
            self.community_data.update(other_data["community_data"])

        if "literature_data" in other_data:
            self.literature_data.update(other_data["literature_data"])

        if "regulatory_data" in other_data:
            self.regulatory_data.update(other_data["regulatory_data"])

        if "experience_reports" in other_data:
            self.experience_reports.extend(
                report for report in other_data["experience_reports"]
                if report not in self.experience_reports
            )

        if "safety_profile" in other_data:
            self.safety_profile.update(other_data["safety_profile"])

        if "dosage_info" in other_data:
            for route, stats in other_data["dosage_info"].items():
                if route not in self.dosage_info:
                    self.dosage_info[route] = stats
                else:
                    current = self.dosage_info[route]
                    current["min"] = min(current["min"], stats["min"])
                    current["max"] = max(current["max"], stats["max"])
                    current["avg"] = (current["avg"] * current["count"] + 
                                    stats["avg"] * stats["count"]) / (
                                        current["count"] + stats["count"]
                                    )
                    current["count"] += stats["count"]


@dataclass
class PredictionsMixin:
    """Mixin for ML prediction capabilities."""
    
    # Feature Management
    _feature_cache: Dict[str, np.ndarray] = field(default_factory=dict)
    _feature_importances: Dict[str, Dict[str, float]] = field(default_factory=dict)
    _feature_scalers: Dict[str, object] = field(default_factory=dict)

    # Prediction Integration
    _prediction_cache: Dict[str, object] = field(default_factory=dict)
    _prediction_history: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(
        columns=[
            'predictor_type',
            'prediction_value',
            'confidence',
            'timestamp',
            'supporting_data',
        ]
    ))

    def get_cached_prediction(
        self,
        predictor_type: str,
    ) -> Optional[object]:
        """Get cached prediction result if available."""
        return self._prediction_cache.get(predictor_type)

    def cache_prediction(
        self,
        predictor_type: str,
        result: object,
        confidence: float,
        supporting_data: Dict = None,
    ) -> None:
        """Cache prediction result."""
        self._prediction_cache[predictor_type] = result
        
        # Update history
        self._prediction_history = pd.concat([
            self._prediction_history,
            pd.DataFrame([{
                'predictor_type': predictor_type,
                'prediction_value': result,
                'confidence': confidence,
                'timestamp': pd.Timestamp.now(),
                'supporting_data': json.dumps(supporting_data or {}),
            }])
        ], ignore_index=True)

    def clear_prediction_cache(self) -> None:
        """Clear cached predictions."""
        self._prediction_cache.clear()


@dataclass
class AnalysisMixin:
    """Mixin for analysis capabilities."""
    
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
            "scores": {},
        }

        # Add toxicity alerts
        for alert in self.toxicity_alerts:
            if alert.get("severity", 0) > 0.7:  # High severity threshold
                summary["warnings"].append({
                    "type": alert["type"],
                    "severity": alert["severity"],
                    "description": alert["description"],
                })

        # Add contraindications and interactions
        summary["contraindications"].extend(self.contraindications)
        summary["interactions"].extend(self.drug_interactions)

        # Add safety scores
        summary["scores"] = self.safety_scores.copy()

        return summary
