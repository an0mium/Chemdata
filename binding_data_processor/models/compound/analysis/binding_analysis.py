"""Receptor binding analysis functionality.

This module provides comprehensive analysis of receptor binding profiles:
- Binding affinity analysis and classification
- Target selectivity analysis
- Binding pattern detection
- Risk assessment based on binding
- Effect prediction based on binding patterns
- Confidence assessment
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from collections import Counter
import numpy as np

from ..types import (
    BindingType,
    RiskLevel,
    ReceptorBinding,
    TargetData,
)


@dataclass
class BindingAnalysisMixin:
    """Mixin providing receptor binding analysis capabilities."""

    # Core binding data
    receptor_profiles: Dict[str, ReceptorBinding] = field(default_factory=dict)
    target_data: Dict[str, TargetData] = field(default_factory=dict)
    _binding_analysis: Dict = field(default_factory=dict)
    
    # Receptor families for analysis
    receptor_families: Dict[str, Set[str]] = {
        "serotonin": {
            "5-HT1A", "5-HT1B", "5-HT1D", "5-HT2A", "5-HT2B", "5-HT2C",
            "5-HT3", "5-HT4", "5-HT5A", "5-HT6", "5-HT7"
        },
        "dopamine": {"D1", "D2", "D3", "D4", "D5"},
        "norepinephrine": {"α1", "α2", "β1", "β2", "β3"},
        "histamine": {"H1", "H2", "H3", "H4"},
        "muscarinic": {"M1", "M2", "M3", "M4", "M5"},
        "nicotinic": {"α4β2", "α7"},
        "opioid": {"mu", "kappa", "delta", "NOP"},
        "glutamate": {"NMDA", "AMPA", "Kainate", "mGluR1", "mGluR5"},
        "gaba": {"GABA-A", "GABA-B"},
        "cannabinoid": {"CB1", "CB2"},
        "sigma": {"σ1", "σ2"},
        "transporter": {"SERT", "NET", "DAT"},
    }

    # High-risk binding patterns
    high_risk_combinations: Dict[str, Tuple[str, RiskLevel, str]] = {
        # Severe risks
        "5-HT2B": (
            "agonist",
            RiskLevel.SEVERE,
            "cardiotoxicity"
        ),
        "mu-opioid": (
            "agonist",
            RiskLevel.SEVERE,
            "respiratory_depression"
        ),
        # High risks
        "NMDA": (
            "antagonist",
            RiskLevel.HIGH,
            "neurotoxicity"
        ),
        "D2": (
            "agonist",
            RiskLevel.HIGH,
            "addiction"
        ),
        "5-HT2A": (
            "agonist",
            RiskLevel.HIGH,
            "psychosis_risk"
        ),
        "GABA-A": (
            "agonist",
            RiskLevel.HIGH,
            "dependence"
        ),
    }

    def analyze_binding(self) -> Dict:
        """Analyze binding data and generate comprehensive summary."""
        analysis = {
            "profiles": self._analyze_receptor_profiles(),
            "patterns": self._analyze_binding_patterns(),
            "selectivity": self._analyze_selectivity(),
            "risks": self._analyze_binding_risks(),
            "predictions": self._predict_binding_effects(),
            "confidence": self._analyze_binding_confidence(),
        }
        self._binding_analysis = analysis
        return analysis

    def _analyze_receptor_profiles(self) -> Dict:
        """Analyze individual receptor binding profiles."""
        profiles = {}
        
        for receptor, (affinity, confidence, activity) in self.receptor_profiles.items():
            # Get target data
            target = self.target_data.get(receptor, TargetData())
            
            # Analyze binding
            profiles[receptor] = {
                "affinity": affinity,
                "confidence": confidence,
                "activity": activity,
                "binding_type": self._classify_binding(affinity),
                "family": self._get_receptor_family(receptor),
                "is_primary": target.is_primary,
                "references": list(target.reference_dois),
                "risks": self._get_receptor_risks(receptor, activity),
                "experimental_data": bool(target.experimental_conditions),
            }
            
        return profiles

    def _analyze_binding_patterns(self) -> Dict:
        """Analyze binding patterns across receptor families."""
        patterns = {}
        
        for family, members in self.receptor_families.items():
            # Get binding data for family members
            family_data = [
                (receptor, data)
                for receptor, data in self.receptor_profiles.items()
                if any(member in receptor for member in members)
            ]
            
            if family_data:
                # Calculate family statistics
                affinities = [affinity for _, (affinity, _, _) in family_data]
                activities = [activity for _, (_, _, activity) in family_data]
                
                patterns[family] = {
                    "receptors": len(family_data),
                    "mean_affinity": float(np.mean(affinities)),
                    "max_affinity": float(min(affinities)),  # Lower is stronger
                    "activity_distribution": self._analyze_activities(activities),
                    "binding_profile": self._analyze_family_profile(family_data),
                }
                
        return patterns

    def _analyze_selectivity(self) -> Dict:
        """Analyze receptor selectivity patterns."""
        # Calculate individual receptor selectivity
        receptor_selectivity = {}
        for receptor, (affinity, confidence, activity) in self.receptor_profiles.items():
            # Calculate selectivity ratio
            other_affinities = [
                a for r, (a, _, _) in self.receptor_profiles.items()
                if r != receptor
            ]
            
            if other_affinities:
                ratio = float(np.mean(other_affinities) / affinity)
            else:
                ratio = 1.0
                
            receptor_selectivity[receptor] = {
                "ratio": ratio,
                "classification": self._classify_selectivity(ratio),
                "confidence": confidence,
            }
            
        # Calculate family selectivity
        family_affinities = self._get_family_affinities()
        family_selectivity = {}
        if family_affinities:
            selectivity_ratios = self._calculate_selectivity_ratios(family_affinities)
            if selectivity_ratios:
                avg_ratio = sum(selectivity_ratios) / len(selectivity_ratios)
                family_selectivity = {
                    "selectivity": self._get_selectivity_classification(avg_ratio),
                    "score": avg_ratio,
                    "ratios": selectivity_ratios,
                }
            
        return {
            "receptor_selectivity": receptor_selectivity,
            "family_selectivity": family_selectivity,
        }

    def _analyze_binding_risks(self) -> Dict:
        """Analyze binding-related risks."""
        risks = {}
        max_risk = RiskLevel.UNKNOWN
        
        # Check each receptor for risky binding patterns
        for receptor, (affinity, confidence, activity) in self.receptor_profiles.items():
            if affinity <= 100:  # Only consider strong binding
                receptor_risks = self._get_receptor_risks(receptor, activity)
                if receptor_risks:
                    risks[receptor] = {
                        "affinity": affinity,
                        "confidence": confidence,
                        "activity": activity,
                        "risks": receptor_risks,
                    }
                    
                    # Track highest risk level
                    for risk in receptor_risks:
                        risk_level = RiskLevel[risk["level"].upper()]
                        if risk_level.value > max_risk.value:
                            max_risk = risk_level
                            
        if risks:
            risks["overall_risk"] = max_risk.value
            risks["risk_distribution"] = self._analyze_risk_distribution(risks)
            
        return risks

    def _predict_binding_effects(self) -> Dict:
        """Predict effects based on binding patterns."""
        predictions = {}
        
        # Analyze each receptor family
        for family, members in self.receptor_families.items():
            family_data = [
                (receptor, data)
                for receptor, data in self.receptor_profiles.items()
                if any(member in receptor for member in members)
            ]
            
            if family_data:
                predictions[family] = self._predict_family_effects(family, family_data)
                
        return predictions

    def _analyze_binding_confidence(self) -> Dict:
        """Analyze confidence in binding data."""
        if not self.target_data:
            return {"overall": 0.0, "by_target": {}}

        confidences = {}
        for target_name, target in self.target_data.items():
            confidences[target_name] = {
                "confidence": target.confidence,
                "support": len(target.reference_dois),
                "experimental": bool(target.experimental_conditions),
                "assay_details": bool(target.assay_details),
            }

        return {
            "overall": sum(t.confidence for t in self.target_data.values()) / len(self.target_data),
            "by_target": confidences,
        }

    # Helper methods
    def _get_receptor_family(self, receptor: str) -> Optional[str]:
        """Get receptor family for a receptor."""
        for family, members in self.receptor_families.items():
            if any(member in receptor for member in members):
                return family
        return None

    def _classify_binding(self, affinity: float) -> str:
        """Classify binding strength from affinity value."""
        if affinity <= 1:  # Sub-nanomolar
            return "very_strong"
        elif affinity <= 10:  # Low nanomolar
            return "strong"
        elif affinity <= 100:  # High nanomolar
            return "moderate"
        elif affinity <= 1000:  # Micromolar
            return "weak"
        return "negligible"

    def _classify_selectivity(self, ratio: float) -> str:
        """Classify selectivity from ratio value."""
        if ratio >= 100:
            return "highly_selective"
        elif ratio >= 10:
            return "selective"
        elif ratio >= 3:
            return "moderately_selective"
        return "nonselective"

    def _get_receptor_risks(self, receptor: str, activity: str) -> List[Dict]:
        """Get risks associated with receptor binding."""
        risks = []
        
        # Check each known high-risk combination
        for risk_receptor, (risk_activity, risk_level, risk_type) in self.high_risk_combinations.items():
            if risk_receptor in receptor and activity.lower() == risk_activity.lower():
                risks.append({
                    "level": risk_level.value,
                    "type": risk_type,
                    "mechanism": f"{receptor}_{activity}",
                })
                
        return risks

    def _analyze_activities(self, activities: List[str]) -> Dict[str, float]:
        """Analyze distribution of activities."""
        total = len(activities)
        if not total:
            return {}
            
        counts = {
            "agonist": activities.count("agonist") / total,
            "antagonist": activities.count("antagonist") / total,
            "partial_agonist": activities.count("partial_agonist") / total,
            "inverse_agonist": activities.count("inverse_agonist") / total,
            "modulator": activities.count("modulator") / total,
        }
        
        return {k: v for k, v in counts.items() if v > 0}

    def _analyze_family_profile(
        self,
        family_data: List[Tuple[str, ReceptorBinding]]
    ) -> Dict:
        """Analyze binding profile for a receptor family."""
        if not family_data:
            return {}
            
        # Get strongest binding
        strongest = min(family_data, key=lambda x: x[1][0])
        receptor, (affinity, confidence, activity) = strongest
        
        return {
            "primary_activity": self._get_primary_activity(family_data),
            "strongest_binding": {
                "receptor": receptor,
                "affinity": affinity,
                "confidence": confidence,
                "activity": activity,
            },
            "selectivity": self._calculate_family_selectivity(family_data),
        }

    def _get_primary_activity(
        self,
        family_data: List[Tuple[str, ReceptorBinding]]
    ) -> Optional[str]:
        """Get primary activity type for a family."""
        activities = [activity for _, (_, _, activity) in family_data]
        if not activities:
            return None
            
        return Counter(activities).most_common(1)[0][0]

    def _calculate_family_selectivity(
        self,
        family_data: List[Tuple[str, ReceptorBinding]]
    ) -> float:
        """Calculate selectivity within a receptor family."""
        if not family_data:
            return 0.0
            
        # Get affinities
        family_affinities = [affinity for _, (affinity, _, _) in family_data]
        
        # Get affinities for other families
        other_affinities = [
            affinity
            for receptor, (affinity, _, _) in self.receptor_profiles.items()
            if not any(r[0] == receptor for r in family_data)
        ]
        
        if not other_affinities:
            return 1.0
            
        # Calculate selectivity ratio
        family_mean = np.mean(family_affinities)
        other_mean = np.mean(other_affinities)
        
        return float(other_mean / family_mean if family_mean > 0 else 0.0)

    def _get_family_affinities(self) -> Dict[str, List[float]]:
        """Get affinities grouped by receptor family."""
        family_affinities = {}
        for receptor, (affinity, _, _) in self.receptor_profiles.items():
            if affinity <= 0:
                continue
            family = self._get_receptor_family(receptor)
            if family:
                if family not in family_affinities:
                    family_affinities[family] = []
                family_affinities[family].append(affinity)
        return family_affinities

    def _calculate_selectivity_ratios(self, family_affinities: Dict[str, List[float]]) -> List[float]:
        """Calculate selectivity ratios between receptor families."""
        ratios = []
        for fam1, aff1 in family_affinities.items():
            min_aff1 = min(aff1)
            for fam2, aff2 in family_affinities.items():
                if fam1 == fam2:
                    continue
                min_aff2 = min(aff2)
                ratio = min_aff2 / min_aff1 if min_aff1 > 0 else 0
                ratios.append(ratio)
        return ratios

    def _get_selectivity_classification(self, avg_ratio: float) -> str:
        """Get selectivity classification based on average ratio."""
        if avg_ratio >= 100:
            return "highly selective"
        elif avg_ratio >= 10:
            return "moderately selective"
        return "non-selective"

    def _analyze_risk_distribution(self, risks: Dict) -> Dict[str, int]:
        """Analyze distribution of risk levels."""
        distribution = {level.value: 0 for level in RiskLevel}
        
        for receptor_data in risks.values():
            if isinstance(receptor_data, dict) and "risks" in receptor_data:
                for risk in receptor_data["risks"]:
                    distribution[risk["level"]] += 1
                    
        return {k: v for k, v in distribution.items() if v > 0}

    def _predict_family_effects(
        self,
        family: str,
        family_data: List[Tuple[str, ReceptorBinding]]
    ) -> Dict:
        """Predict effects for a receptor family."""
        if not family_data:
            return {}
            
        # Get primary activity and strongest binding
        primary_activity = self._get_primary_activity(family_data)
        strongest = min(family_data, key=lambda x: x[1][0])
        receptor, (affinity, confidence, activity) = strongest
        
        # Define effect patterns
        effect_patterns = {
            "serotonin": {
                "agonist": [
                    ("mood_elevation", 0.8),
                    ("anxiety_reduction", 0.6),
                ],
                "antagonist": [
                    ("mood_stabilization", 0.7),
                    ("antipsychotic", 0.6),
                ],
            },
            "dopamine": {
                "agonist": [
                    ("reward_enhancement", 0.8),
                    ("motivation_increase", 0.7),
                ],
                "antagonist": [
                    ("antipsychotic", 0.8),
                    ("mood_stabilization", 0.6),
                ],
            },
        }
        
        # Get predictions for family
        if family in effect_patterns and primary_activity in effect_patterns[family]:
            predictions = {
                effect: {
                    "probability": prob,
                    "confidence": confidence,
                    "mechanism": f"{family}_{primary_activity}",
                }
                for effect, prob in effect_patterns[family][primary_activity]
            }
            return predictions
            
        return {}

    def merge_binding_data(self, other: "BindingAnalysisMixin") -> None:
        """Merge binding data from another instance."""
        # Merge receptor profiles
        for receptor, (affinity, confidence, activity) in other.receptor_profiles.items():
            if receptor not in self.receptor_profiles:
                self.receptor_profiles[receptor] = (affinity, confidence, activity)
            else:
                # Keep data with higher confidence
                old_affinity, old_conf, old_activity = self.receptor_profiles[receptor]
                if confidence > old_conf:
                    self.receptor_profiles[receptor] = (affinity, confidence, activity)
                    
        # Merge target data
        for target, data in other.target_data.items():
            if target not in self.target_data:
                self.target_data[target] = data
            else:
                # Merge reference DOIs
                self.target_data[target].reference_dois.update(data.reference_dois)
