"""Activity analysis functionality for compound data.

This module provides the ActivityAnalysisMixin class that adds activity analysis capabilities:
- Primary activity identification
- Activity pattern analysis
- Mechanism of action analysis
- Effect profile analysis
- Activity confidence assessment
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class ActivityAnalysisMixin:
    """Mixin class adding activity analysis capabilities."""

    _activity_analysis: Dict = field(default_factory=dict)

    def analyze_activity(self) -> Dict:
        """Analyze activity data to identify patterns and mechanisms."""
        analysis = {
            "primary_activity": self._find_primary_activity(),
            "activity_patterns": self._analyze_activity_patterns(),
            "mechanisms": self._analyze_mechanisms(),
            "effect_profile": self._analyze_effect_profile(),
            "activity_confidence": self._analyze_activity_confidence(),
        }
        self._activity_analysis = analysis
        return analysis

    def _find_primary_activity(self) -> Optional[Dict]:
        """Find primary activity type and mechanism."""
        if not hasattr(self, "primary_activity") or not self.primary_activity:
            return None

        return {
            "activity": self.primary_activity,
            "mechanism": self.mechanism_of_action,
            "confidence": self._calculate_activity_confidence(),
        }

    def _analyze_activity_patterns(self) -> List[Dict]:
        """Analyze patterns in activity data."""
        patterns = []

        # Analyze receptor-activity patterns
        for target in self.targets:
            if target.activity_type == "N/A":
                continue

            pattern = {
                "target": target.common_name,
                "activity": target.activity_type,
                "family": target.assay_details.get("family", "unknown"),
                "affinity": target.affinity_value,
                "confidence": target.confidence,
            }
            patterns.append(pattern)

        return patterns

    def _analyze_mechanisms(self) -> Dict:
        """Analyze mechanisms of action."""
        mechanisms = {}

        # Primary mechanism
        if self.mechanism_of_action:
            mechanisms["primary"] = {
                "mechanism": self.mechanism_of_action,
                "confidence": self._calculate_mechanism_confidence(),
            }

        # Secondary mechanisms from targets
        secondary = []
        for target in self.targets:
            if not target.activity_type or target.activity_type == "N/A":
                continue
            
            mechanism = {
                "target": target.common_name,
                "activity": target.activity_type,
                "affinity": target.affinity_value,
                "confidence": target.confidence,
            }
            secondary.append(mechanism)

        mechanisms["secondary"] = secondary
        return mechanisms

    def _analyze_effect_profile(self) -> Dict:
        """Analyze pharmacological effect profile."""
        if not hasattr(self, "effect_profile"):
            return {}

        effects = {}
        for effect, score in self.effect_profile.items():
            effects[effect] = {
                "score": score,
                "confidence": self._calculate_effect_confidence(effect),
            }

        return effects

    def _analyze_activity_confidence(self) -> float:
        """Calculate overall confidence in activity predictions."""
        confidences = []

        # Target data confidence
        if self.targets:
            confidences.append(
                sum(t.confidence for t in self.targets) / len(self.targets)
            )

        # Effect profile confidence
        if hasattr(self, "effect_profile") and self.effect_profile:
            effect_conf = self._calculate_effect_confidence(
                max(self.effect_profile, key=self.effect_profile.get)
            )
            confidences.append(effect_conf)

        # Mechanism confidence
        if self.mechanism_of_action:
            confidences.append(self._calculate_mechanism_confidence())

        return sum(confidences) / len(confidences) if confidences else 0.0

    def _calculate_effect_confidence(self, effect: str) -> float:
        """Calculate confidence for a specific effect."""
        confidence = 0.0
        factors = 0

        # Target evidence
        relevant_targets = [
            t for t in self.targets
            if effect.lower() in t.assay_details.get("effects", [])
        ]
        if relevant_targets:
            confidence += sum(t.confidence for t in relevant_targets) / len(relevant_targets)
            factors += 1

        # Literature evidence
        if hasattr(self, "literature_data"):
            effect_refs = len([
                ref for ref in self.literature_data.get("references", [])
                if effect.lower() in ref.get("effects", [])
            ])
            if effect_refs:
                confidence += min(effect_refs / 5, 1.0)  # Cap at 5 references
                factors += 1

        # Community reports
        if hasattr(self, "experience_reports"):
            effect_reports = len([
                report for report in self.experience_reports
                if effect.lower() in report.get("effects", [])
            ])
            if effect_reports:
                confidence += min(effect_reports / 10, 1.0)  # Cap at 10 reports
                factors += 1

        return confidence / factors if factors > 0 else 0.0

    def _calculate_mechanism_confidence(self) -> float:
        """Calculate confidence in mechanism of action."""
        confidence = 0.0
        factors = 0

        # Target evidence
        relevant_targets = [
            t for t in self.targets
            if self.mechanism_of_action.lower() in t.assay_details.get("mechanisms", [])
        ]
        if relevant_targets:
            confidence += sum(t.confidence for t in relevant_targets) / len(relevant_targets)
            factors += 1

        # Literature evidence
        if hasattr(self, "literature_data"):
            mech_refs = len([
                ref for ref in self.literature_data.get("references", [])
                if self.mechanism_of_action.lower() in ref.get("mechanisms", [])
            ])
            if mech_refs:
                confidence += min(mech_refs / 5, 1.0)  # Cap at 5 references
                factors += 1

        return confidence / factors if factors > 0 else 0.0
