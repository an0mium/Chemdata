"""Activity analysis functionality for compound data.

This module provides analysis capabilities for compound activity data:
- Pharmacological activity analysis
- Mechanism of action analysis
- Activity pattern analysis
- Functional selectivity analysis
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
import numpy as np

from ..types import (
    ActivityType,
    MechanismType,
    SelectivityType,
)


def default_activity_analysis() -> Dict:
    """Default empty activity analysis dictionary."""
    return {}


def default_activity_profiles() -> Dict[str, Tuple[float, float, str]]:
    """Default empty activity profiles dictionary."""
    return {}


def default_mechanisms() -> Dict[str, MechanismType]:
    """Default empty mechanisms dictionary."""
    return {}


def default_effect_mechanisms() -> Dict[str, List[str]]:
    """Default empty effect mechanisms dictionary."""
    return {"primary": [], "secondary": []}


def default_selectivity_data() -> Dict[str, SelectivityType]:
    """Default empty selectivity data dictionary."""
    return {}


@dataclass
class ActivityAnalyzer:
    """Analyzer for compound activity data."""

    # Core activity data
    _activity_analysis: Dict = field(default_factory=default_activity_analysis)
    activity_profiles: Dict[str, Tuple[float, float, str]] = field(default_factory=default_activity_profiles)
    mechanisms: Dict[str, MechanismType] = field(default_factory=default_mechanisms)
    effect_mechanisms: Dict[str, List[str]] = field(default_factory=default_effect_mechanisms)
    selectivity_data: Dict[str, SelectivityType] = field(default_factory=default_selectivity_data)

    def analyze_activity(self) -> Dict:
        """Analyze activity data to identify patterns and mechanisms."""
        analysis = {
            # Core activity analysis
            "activity_profiles": self._analyze_activity_profiles(),
            "mechanism_patterns": self._analyze_mechanism_patterns(),
            "selectivity_analysis": self._analyze_selectivity(),
            # Functional analysis
            "functional_activity": self._analyze_functional_activity(),
            "pathway_effects": self._analyze_pathway_effects(),
            # Overall assessment
            "activity_assessment": self._analyze_activity_assessment(),
        }
        self._activity_analysis = analysis
        return analysis

    def _analyze_activity_profiles(self) -> Dict:
        """Analyze pharmacological activity profiles."""
        profiles = {}

        for target, (potency, efficacy, mechanism) in self.activity_profiles.items():
            mechanism_type = self.mechanisms.get(target, MechanismType.UNKNOWN)

            if mechanism_type not in profiles:
                profiles[mechanism_type] = []

            profiles[mechanism_type].append(
                {
                    "target": target,
                    "potency": potency,
                    "efficacy": efficacy,
                    "mechanism": mechanism,
                }
            )

        return profiles

    def _analyze_mechanism_patterns(self) -> Dict:
        """Analyze patterns in mechanisms of action."""
        patterns = {}

        # Primary mechanisms
        patterns["primary"] = [
            {"mechanism": mechanism, "targets": [target for target, mech in self.mechanisms.items() if mech == mechanism and target in self.effect_mechanisms["primary"]]}
            for mechanism in set(self.mechanisms.values())
            if any(target in self.effect_mechanisms["primary"] for target, mech in self.mechanisms.items() if mech == mechanism)
        ]

        # Secondary mechanisms
        patterns["secondary"] = [
            {"mechanism": mechanism, "targets": [target for target, mech in self.mechanisms.items() if mech == mechanism and target in self.effect_mechanisms["secondary"]]}
            for mechanism in set(self.mechanisms.values())
            if any(target in self.effect_mechanisms["secondary"] for target, mech in self.mechanisms.items() if mech == mechanism)
        ]

        return patterns

    def _analyze_selectivity(self) -> Dict:
        """Analyze target selectivity patterns."""
        analysis = {}

        for target, selectivity in self.selectivity_data.items():
            if selectivity not in analysis:
                analysis[selectivity] = []

            if target in self.activity_profiles:
                potency, efficacy, mechanism = self.activity_profiles[target]
                analysis[selectivity].append(
                    {
                        "target": target,
                        "potency": potency,
                        "efficacy": efficacy,
                        "mechanism": mechanism,
                    }
                )

        return analysis

    def _analyze_functional_activity(self) -> Dict:
        """Analyze functional activity patterns."""
        functional = {}

        # Group by mechanism
        for target, (potency, efficacy, mechanism) in self.activity_profiles.items():
            if mechanism not in functional:
                functional[mechanism] = []

            functional[mechanism].append(
                {
                    "target": target,
                    "potency": potency,
                    "efficacy": efficacy,
                }
            )

        # Calculate mechanism statistics
        for mechanism, activities in functional.items():
            avg_potency = sum(a["potency"] for a in activities) / len(activities)
            avg_efficacy = sum(a["efficacy"] for a in activities) / len(activities)

            functional[mechanism] = {
                "activities": activities,
                "avg_potency": avg_potency,
                "avg_efficacy": avg_efficacy,
            }

        return functional

    def _analyze_pathway_effects(self) -> Dict:
        """Analyze pathway-level effects."""
        pathways = {}

        # Group by primary/secondary effects
        for effect_type, targets in self.effect_mechanisms.items():
            pathway_effects = []

            for target in targets:
                if target in self.activity_profiles:
                    potency, efficacy, mechanism = self.activity_profiles[target]
                    mechanism_type = self.mechanisms.get(target, MechanismType.UNKNOWN)

                    pathway_effects.append(
                        {
                            "target": target,
                            "potency": potency,
                            "efficacy": efficacy,
                            "mechanism": mechanism,
                            "mechanism_type": mechanism_type,
                        }
                    )

            if pathway_effects:
                pathways[effect_type] = pathway_effects

        return pathways

    def _analyze_activity_assessment(self) -> Dict:
        """Perform overall activity assessment."""
        assessment = {
            "primary_mechanisms": self._analyze_primary_mechanisms(),
            "secondary_mechanisms": self._analyze_secondary_mechanisms(),
            "selectivity_profile": self._analyze_selectivity_profile(),
        }

        # Add overall metrics
        if self.activity_profiles:
            potencies = [potency for potency, _, _ in self.activity_profiles.values()]
            efficacies = [efficacy for _, efficacy, _ in self.activity_profiles.values()]

            assessment["metrics"] = {
                "avg_potency": float(np.mean(potencies)),
                "avg_efficacy": float(np.mean(efficacies)),
                "potency_std": float(np.std(potencies)),
                "efficacy_std": float(np.std(efficacies)),
            }

        return assessment

    def _analyze_primary_mechanisms(self) -> List[Dict]:
        """Analyze primary mechanisms of action."""
        primary = []

        for target in self.effect_mechanisms["primary"]:
            if target in self.activity_profiles:
                potency, efficacy, mechanism = self.activity_profiles[target]
                mechanism_type = self.mechanisms.get(target, MechanismType.UNKNOWN)

                primary.append(
                    {
                        "target": target,
                        "potency": potency,
                        "efficacy": efficacy,
                        "mechanism": mechanism,
                        "mechanism_type": mechanism_type,
                    }
                )

        return sorted(primary, key=lambda x: x["potency"])

    def _analyze_secondary_mechanisms(self) -> List[Dict]:
        """Analyze secondary mechanisms of action."""
        secondary = []

        for target in self.effect_mechanisms["secondary"]:
            if target in self.activity_profiles:
                potency, efficacy, mechanism = self.activity_profiles[target]
                mechanism_type = self.mechanisms.get(target, MechanismType.UNKNOWN)

                secondary.append(
                    {
                        "target": target,
                        "potency": potency,
                        "efficacy": efficacy,
                        "mechanism": mechanism,
                        "mechanism_type": mechanism_type,
                    }
                )

        return sorted(secondary, key=lambda x: x["potency"])

    def _analyze_selectivity_profile(self) -> Dict:
        """Analyze overall selectivity profile."""
        profile = {}

        # Count selectivity types
        selectivity_counts = {}
        for selectivity in self.selectivity_data.values():
            if selectivity not in selectivity_counts:
                selectivity_counts[selectivity] = 0
            selectivity_counts[selectivity] += 1

        # Calculate selectivity ratios
        total = sum(selectivity_counts.values())
        if total > 0:
            profile["selectivity_ratios"] = {selectivity: count / total for selectivity, count in selectivity_counts.items()}

        # Identify dominant selectivity
        if selectivity_counts:
            profile["dominant_selectivity"] = max(selectivity_counts.items(), key=lambda x: x[1])[0]

        return profile

    def merge_activity_data(self, other: "ActivityAnalyzer") -> None:
        """Merge activity data from another instance."""
        # Merge activity profiles
        for target, (potency, efficacy, mechanism) in other.activity_profiles.items():
            if target not in self.activity_profiles:
                self.activity_profiles[target] = (potency, efficacy, mechanism)
            else:
                # Keep data with higher potency
                current_potency, current_efficacy, current_mechanism = self.activity_profiles[target]
                if potency < current_potency:  # Lower potency value = higher potency
                    self.activity_profiles[target] = (potency, efficacy, mechanism)

        # Merge mechanisms
        self.mechanisms.update(other.mechanisms)

        # Merge effect mechanisms
        for effect_type in ["primary", "secondary"]:
            self.effect_mechanisms[effect_type].extend(target for target in other.effect_mechanisms[effect_type] if target not in self.effect_mechanisms[effect_type])

        # Merge selectivity data
        self.selectivity_data.update(other.selectivity_data)
