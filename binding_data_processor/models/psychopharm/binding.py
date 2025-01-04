"""Receptor binding analysis functionality."""

from dataclasses import field
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .base import ReceptorBinding, RiskLevel


class ReceptorProfileMixin:
    """Mixin class providing receptor binding profile analysis."""

    # Receptor binding profiles with structured data
    receptor_profiles: Dict[str, ReceptorBinding] = field(default_factory=dict)
    
    # Receptor families for grouping and analysis
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

    # High-risk receptor-activity combinations
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
        
        # Moderate risks
        "CB1": (
            "agonist",
            RiskLevel.MODERATE,
            "cognitive_impairment"
        ),
    }

    def get_binding_dict(self) -> Dict:
        """Get dictionary of binding data."""
        return {
            "receptor_profiles": self._format_receptor_profiles(),
            "binding_analysis": self._analyze_binding_patterns(),
            "selectivity_analysis": self._analyze_selectivity(),
            "predicted_effects": self._predict_binding_effects(),
            "risk_assessment": self._assess_binding_risks(),
        }

    def _format_receptor_profiles(self) -> Dict:
        """Format receptor binding profiles."""
        formatted = {}
        for receptor, (affinity, confidence, activity) in self.receptor_profiles.items():
            formatted[receptor] = {
                "affinity": affinity,
                "confidence": confidence,
                "activity": activity,
                "family": self._get_receptor_family(receptor),
                "binding_type": self._predict_binding_type(affinity),
                "risks": self._get_receptor_risks(receptor, activity),
            }
        return formatted

    def _get_receptor_family(self, receptor: str) -> Optional[str]:
        """Get receptor family for a receptor."""
        for family, members in self.receptor_families.items():
            if any(member in receptor for member in members):
                return family
        return None

    def _predict_binding_type(self, affinity: float) -> str:
        """Predict binding type from affinity value."""
        if affinity <= 1:  # Sub-nanomolar
            return "very_strong"
        elif affinity <= 10:  # Low nanomolar
            return "strong"
        elif affinity <= 100:  # High nanomolar
            return "moderate"
        elif affinity <= 1000:  # Micromolar
            return "weak"
        return "negligible"

    def _get_receptor_risks(self, receptor: str, activity: str) -> List[Dict]:
        """Get risks associated with receptor binding."""
        risks = []
        
        # Check each known high-risk combination
        for risk_receptor, risk_data in self.high_risk_combinations.items():
            risk_activity, risk_level, risk_type = risk_data
            if risk_receptor in receptor and activity.lower() == risk_activity.lower():
                risks.append({
                    "level": risk_level.value,
                    "type": risk_type,
                    "mechanism": f"{receptor}_{activity}",
                })
                
        return risks

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
                    "mean_affinity": np.mean(affinities),
                    "max_affinity": min(affinities),  # Lower is stronger
                    "activity_distribution": self._analyze_activities(activities),
                    "binding_profile": self._analyze_family_profile(family_data),
                    "risk_profile": self._analyze_family_risks(family_data),
                }
                
        return patterns

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
        return {
            "primary_activity": self._get_primary_activity(family_data),
            "strongest_binding": self._get_strongest_binding(family_data),
            "selectivity": self._calculate_family_selectivity(family_data),
        }

    def _analyze_family_risks(
        self,
        family_data: List[Tuple[str, ReceptorBinding]]
    ) -> Dict:
        """Analyze risks for a receptor family."""
        risks = []
        
        for receptor, (affinity, confidence, activity) in family_data:
            receptor_risks = self._get_receptor_risks(receptor, activity)
            if receptor_risks and affinity <= 100:  # Only include strong binding
                for risk in receptor_risks:
                    risks.append({
                        "receptor": receptor,
                        "affinity": affinity,
                        "confidence": confidence,
                        **risk,
                    })
                    
        if not risks:
            return {}
            
        # Get highest risk level
        risk_levels = [RiskLevel[risk["level"].upper()] for risk in risks]
        max_risk = max(risk_levels, key=lambda x: x.value)
        
        return {
            "overall_risk": max_risk.value,
            "risk_factors": risks,
            "risk_count": len(risks),
        }

    def _get_primary_activity(
        self,
        family_data: List[Tuple[str, ReceptorBinding]]
    ) -> Optional[str]:
        """Get primary activity type for a family."""
        activities = [activity for _, (_, _, activity) in family_data]
        if not activities:
            return None
            
        from collections import Counter
        return Counter(activities).most_common(1)[0][0]

    def _get_strongest_binding(
        self,
        family_data: List[Tuple[str, ReceptorBinding]]
    ) -> Optional[Dict]:
        """Get strongest binding receptor in family."""
        if not family_data:
            return None
            
        strongest = min(family_data, key=lambda x: x[1][0])  # Sort by affinity
        receptor, (affinity, confidence, activity) = strongest
        
        return {
            "receptor": receptor,
            "affinity": affinity,
            "confidence": confidence,
            "activity": activity,
            "risks": self._get_receptor_risks(receptor, activity),
        }

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
        
        return other_mean / family_mean if family_mean > 0 else 0.0

    def _analyze_selectivity(self) -> Dict:
        """Analyze receptor selectivity patterns."""
        selectivity = {}
        
        for receptor, (affinity, confidence, activity) in self.receptor_profiles.items():
            # Calculate selectivity ratio
            other_affinities = [
                a for r, (a, _, _) in self.receptor_profiles.items()
                if r != receptor
            ]
            
            if other_affinities:
                ratio = np.mean(other_affinities) / affinity
            else:
                ratio = 1.0
                
            selectivity[receptor] = {
                "ratio": ratio,
                "classification": self._classify_selectivity(ratio),
                "confidence": confidence,
                "risks": self._get_receptor_risks(receptor, activity),
            }
            
        return selectivity

    def _classify_selectivity(self, ratio: float) -> str:
        """Classify selectivity ratio."""
        if ratio >= 100:
            return "highly_selective"
        elif ratio >= 10:
            return "selective"
        elif ratio >= 3:
            return "moderately_selective"
        return "nonselective"

    def _predict_binding_effects(self) -> Dict:
        """Predict effects based on receptor binding patterns."""
        effects = {}
        
        # Analyze each receptor family
        for family, members in self.receptor_families.items():
            family_data = [
                (receptor, data)
                for receptor, data in self.receptor_profiles.items()
                if any(member in receptor for member in members)
            ]
            
            if family_data:
                effects[family] = self._predict_family_effects(family, family_data)
                
        return effects

    def _predict_family_effects(
        self,
        family: str,
        family_data: List[Tuple[str, ReceptorBinding]]
    ) -> Dict:
        """Predict effects for a receptor family."""
        # Get primary activity and strongest binding
        primary_activity = self._get_primary_activity(family_data)
        strongest = self._get_strongest_binding(family_data)
        
        if not primary_activity or not strongest:
            return {}
            
        # Predict effects based on family and activity
        predictions = []
        risks = []
        
        if family == "serotonin":
            if "5-HT2A" in strongest["receptor"]:
                if primary_activity == "agonist":
                    predictions.extend([
                        ("psychedelic", 0.9),
                        ("cognitive_enhancement", 0.7),
                    ])
                    risks.append({
                        "type": "psychosis_risk",
                        "level": RiskLevel.HIGH.value,
                        "probability": 0.8,
                    })
                elif primary_activity == "antagonist":
                    predictions.extend([
                        ("antipsychotic", 0.8),
                        ("neuroprotective", 0.6),
                    ])
                    
        elif family == "dopamine":
            if primary_activity == "agonist":
                predictions.extend([
                    ("stimulant", 0.8),
                    ("euphoric", 0.7),
                ])
                risks.append({
                    "type": "addiction_risk",
                    "level": RiskLevel.HIGH.value,
                    "probability": 0.7,
                })
            elif primary_activity == "antagonist":
                predictions.extend([
                    ("antipsychotic", 0.8),
                    ("sedating", 0.6),
                ])
                
        # Format predictions with confidence
        return {
            "effects": {
                effect: {
                    "probability": prob,
                    "confidence": strongest["confidence"],
                    "mechanism": f"{family}_{primary_activity}",
                }
                for effect, prob in predictions
            },
            "risks": risks,
        }

    def _assess_binding_risks(self) -> Dict:
        """Assess overall binding-related risks."""
        risks = {}
        max_risk = RiskLevel.UNKNOWN
        
        # Check each receptor for risks
        for receptor, (affinity, confidence, activity) in self.receptor_profiles.items():
            receptor_risks = self._get_receptor_risks(receptor, activity)
            if receptor_risks and affinity <= 100:  # Only include strong binding
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
                        
        return {
            "overall_risk": max_risk.value,
            "risk_factors": risks,
            "risk_count": len(risks),
            "risk_distribution": self._analyze_risk_distribution(risks),
        }

    def _analyze_risk_distribution(self, risks: Dict) -> Dict[str, int]:
        """Analyze distribution of risk levels."""
        distribution = {level.value: 0 for level in RiskLevel}
        
        for receptor_data in risks.values():
            for risk in receptor_data["risks"]:
                distribution[risk["level"]] += 1
                
        return {k: v for k, v in distribution.items() if v > 0}

    def merge_binding_data(self, other: "ReceptorProfileMixin") -> None:
        """Merge binding data from another instance."""
        for receptor, (affinity, confidence, activity) in other.receptor_profiles.items():
            if receptor not in self.receptor_profiles:
                self.receptor_profiles[receptor] = (affinity, confidence, activity)
            else:
                # Keep data with higher confidence
                old_affinity, old_conf, old_activity = self.receptor_profiles[receptor]
                if confidence > old_conf:
                    self.receptor_profiles[receptor] = (affinity, confidence, activity)
