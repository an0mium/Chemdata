"""Safety analysis functionality for compound data.

This module provides comprehensive safety analysis capabilities:
- Toxicity risk analysis
- Safety pattern analysis
- Drug interaction analysis
- Contraindication analysis
- Risk assessment and mitigation
- Psychopharmacological safety analysis
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import numpy as np

from ..types import (
    RiskLevel,
    PsychoactiveClass,
    NootropicMechanism,
)


def default_safety_analysis() -> Dict:
    """Default empty safety analysis dictionary."""
    return {}


def default_safety_alerts() -> Dict[str, RiskLevel]:
    """Default empty safety alerts dictionary."""
    return {}


def default_contraindications() -> Set[str]:
    """Default empty contraindications set."""
    return set()


def default_interaction_risks() -> Dict[str, Tuple[str, RiskLevel]]:
    """Default empty interaction risks dictionary."""
    return {}


def default_risk_thresholds() -> Dict[str, float]:
    """Default risk thresholds dictionary."""
    return {
        "cardiotoxicity": 0.7,
        "neurotoxicity": 0.8,
        "respiratory": 0.9,
        "addiction": 0.6,
    }


def default_interaction_matrix() -> Dict[str, Dict[str, RiskLevel]]:
    """Default interaction matrix with predefined drug interactions."""
    return {
        "SSRI": {
            "MAOI": RiskLevel.SEVERE,
            "MDMA": RiskLevel.SEVERE,
            "DRI": RiskLevel.HIGH,
            "TCA": RiskLevel.HIGH,
        },
        "MAOI": {
            "SSRI": RiskLevel.SEVERE,
            "MDMA": RiskLevel.SEVERE,
            "DRI": RiskLevel.HIGH,
            "TCA": RiskLevel.HIGH,
            "Tryptamine": RiskLevel.SEVERE,
        },
        "DRI": {
            "MAOI": RiskLevel.HIGH,
            "SSRI": RiskLevel.HIGH,
            "TCA": RiskLevel.MODERATE,
        },
        "TCA": {
            "MAOI": RiskLevel.HIGH,
            "SSRI": RiskLevel.HIGH,
            "DRI": RiskLevel.MODERATE,
        },
        "MDMA": {
            "MAOI": RiskLevel.SEVERE,
            "SSRI": RiskLevel.SEVERE,
            "DRI": RiskLevel.HIGH,
        },
        "Tryptamine": {
            "MAOI": RiskLevel.SEVERE,
            "SSRI": RiskLevel.HIGH,
        },
    }


@dataclass
class SafetyAnalyzer:
    """Analyzer for compound safety data."""

    # Core safety data
    _safety_analysis: Dict = field(default_factory=default_safety_analysis)
    safety_alerts: Dict[str, RiskLevel] = field(default_factory=default_safety_alerts)
    contraindications: Set[str] = field(default_factory=default_contraindications)
    interaction_risks: Dict[str, Tuple[str, RiskLevel]] = field(default_factory=default_interaction_risks)
    risk_thresholds: Dict[str, float] = field(
        default_factory=lambda: {
            "cardiotoxicity": 0.7,
            "neurotoxicity": 0.8,
            "respiratory": 0.9,
            "addiction": 0.6,
        }
    )
    interaction_matrix: Dict[str, Dict[str, RiskLevel]] = field(default_factory=default_interaction_matrix)

    def analyze_safety(self) -> Dict:
        """Analyze safety data to identify risks and interactions."""
        analysis = {
            # Core safety analysis
            "toxicity_risks": self._analyze_toxicity_risks(),
            "safety_patterns": self._analyze_safety_patterns(),
            "interactions": self._analyze_interactions(),
            "contraindications": self._analyze_contraindications(),
            # Psychopharm safety analysis
            "binding_risks": self._analyze_binding_risks(),
            "activity_risks": self._analyze_activity_risks(),
            "interaction_risks": self._analyze_interaction_risks(),
            # Overall assessment
            "risk_assessment": self._analyze_risk_assessment(),
            "recommendations": self._generate_safety_recommendations(),
        }
        self._safety_analysis = analysis
        return analysis

    def _analyze_toxicity_risks(self) -> List[Dict]:
        """Analyze toxicity risks."""
        risks = []

        # Structure-based alerts
        if hasattr(self, "structural_alerts"):
            for alert in self.structural_alerts:
                risks.append(
                    {
                        "type": "structural",
                        "alert": alert["description"],
                        "severity": alert.get("severity", "unknown"),
                        "confidence": alert.get("confidence", 0.5),
                    }
                )

        # Target-based risks
        for target in self.targets:
            if "toxicity" in target.assay_details:
                risks.append(
                    {
                        "type": "target",
                        "target": target.common_name,
                        "risk": target.assay_details["toxicity"],
                        "severity": target.assay_details.get("severity", "unknown"),
                        "confidence": target.confidence,
                    }
                )

        # Known risks
        if hasattr(self, "safety_profile"):
            for risk in self.safety_profile.get("known_risks", []):
                risks.append(
                    {
                        "type": "known",
                        "risk": risk["description"],
                        "severity": risk.get("severity", "unknown"),
                        "confidence": risk.get("confidence", 0.8),
                    }
                )

        return risks

    def _analyze_safety_patterns(self) -> Dict:
        """Analyze patterns in safety data."""
        return {
            "structural_patterns": self._analyze_structural_safety(),
            "target_patterns": self._analyze_target_safety(),
            "interaction_patterns": self._analyze_interaction_patterns(),
        }

    def _analyze_structural_safety(self) -> List[Dict]:
        """Analyze structural safety patterns."""
        patterns = []

        if hasattr(self, "structural_alerts"):
            # Group alerts by type
            alert_groups = {}
            for alert in self.structural_alerts:
                alert_type = alert.get("type", "unknown")
                if alert_type not in alert_groups:
                    alert_groups[alert_type] = []
                alert_groups[alert_type].append(alert)

            # Analyze patterns in each group
            for alert_type, alerts in alert_groups.items():
                pattern = {
                    "type": alert_type,
                    "count": len(alerts),
                    "severity": max((a.get("severity", "unknown") for a in alerts), key=lambda s: RiskLevel[s.upper()].value if s.upper() in RiskLevel.__members__ else 0),
                    "confidence": sum(a.get("confidence", 0.5) for a in alerts) / len(alerts),
                }
                patterns.append(pattern)

        return patterns

    def _analyze_target_safety(self) -> List[Dict]:
        """Analyze target-based safety patterns."""
        patterns = []

        # Group targets by safety implications
        safety_groups = {}
        for target in self.targets:
            if "safety" not in target.assay_details:
                continue

            safety_type = target.assay_details["safety"]
            if safety_type not in safety_groups:
                safety_groups[safety_type] = []
            safety_groups[safety_type].append(target)

        # Analyze patterns in each group
        for safety_type, targets in safety_groups.items():
            pattern = {
                "type": safety_type,
                "count": len(targets),
                "avg_affinity": sum(t.affinity_value for t in targets if t.affinity_value > 0) / len(targets),
                "confidence": sum(t.confidence for t in targets) / len(targets),
            }
            patterns.append(pattern)

        return patterns

    def _analyze_interaction_patterns(self) -> List[Dict]:
        """Analyze patterns in drug interactions."""
        patterns = []

        # Group interactions by type
        interaction_groups = {}
        for target in self.targets:
            if "interactions" not in target.assay_details:
                continue

            for interaction in target.assay_details["interactions"]:
                int_type = interaction.get("type", "unknown")
                if int_type not in interaction_groups:
                    interaction_groups[int_type] = []
                interaction_groups[int_type].append(interaction)

        # Analyze patterns in each group
        for int_type, interactions in interaction_groups.items():
            pattern = {
                "type": int_type,
                "count": len(interactions),
                "severity": max((i.get("severity", "unknown") for i in interactions), key=lambda s: RiskLevel[s.upper()].value if s.upper() in RiskLevel.__members__ else 0),
                "confidence": sum(i.get("confidence", 0.5) for i in interactions) / len(interactions),
            }
            patterns.append(pattern)

        return patterns

    def _analyze_binding_risks(self) -> Dict:
        """Analyze risks based on receptor binding patterns."""
        risks = {}
        max_risk = RiskLevel.UNKNOWN

        # Check each receptor for risky binding patterns
        for receptor, (affinity, confidence, activity) in self.receptor_profiles.items():
            if affinity <= 100:  # Only consider strong binding
                # Check for specific risk patterns
                if "5-HT2B" in receptor and activity == "agonist":
                    risks["cardiotoxicity"] = {
                        "risk_level": RiskLevel.SEVERE.value,
                        "mechanism": "5-HT2B agonism",
                        "affinity": affinity,
                        "confidence": confidence,
                    }
                    max_risk = max(max_risk, RiskLevel.SEVERE)

                elif "NMDA" in receptor and activity == "antagonist":
                    risks["neurotoxicity"] = {
                        "risk_level": RiskLevel.HIGH.value,
                        "mechanism": "NMDA antagonism",
                        "affinity": affinity,
                        "confidence": confidence,
                    }
                    max_risk = max(max_risk, RiskLevel.HIGH)

                elif "D2" in receptor and activity == "agonist":
                    risks["addiction"] = {
                        "risk_level": RiskLevel.HIGH.value,
                        "mechanism": "D2 agonism",
                        "affinity": affinity,
                        "confidence": confidence,
                    }
                    max_risk = max(max_risk, RiskLevel.HIGH)

        if risks:
            risks["risk_level"] = max_risk.value

        return risks

    def _analyze_activity_risks(self) -> Dict:
        """Analyze risks based on activity patterns."""
        risks = {}
        max_risk = RiskLevel.UNKNOWN

        # Check psychoactive class risks
        if self.psychoactive_class != PsychoactiveClass.UNKNOWN:
            class_risks = {
                PsychoactiveClass.PSYCHEDELIC: [
                    ("psychosis", RiskLevel.HIGH),
                    ("hppd", RiskLevel.MODERATE),
                ],
                PsychoactiveClass.STIMULANT: [
                    ("cardiovascular", RiskLevel.HIGH),
                    ("addiction", RiskLevel.HIGH),
                ],
                PsychoactiveClass.DEPRESSANT: [
                    ("respiratory", RiskLevel.SEVERE),
                    ("dependence", RiskLevel.HIGH),
                ],
            }

            if self.psychoactive_class in class_risks:
                for risk_type, risk_level in class_risks[self.psychoactive_class]:
                    risks[risk_type] = {
                        "risk_level": risk_level.value,
                        "mechanism": f"{self.psychoactive_class.value}_class",
                    }
                    max_risk = max(max_risk, risk_level)

        # Check mechanism-based risks
        for mechanism in self.nootropic_mechanisms:
            if mechanism == NootropicMechanism.GLUTAMATERGIC:
                risks["excitotoxicity"] = {
                    "risk_level": RiskLevel.MODERATE.value,
                    "mechanism": "glutamate_modulation",
                }
                max_risk = max(max_risk, RiskLevel.MODERATE)

        if risks:
            risks["risk_level"] = max_risk.value

        return risks

    def _analyze_interactions(self) -> List[Dict]:
        """Analyze potential drug interactions."""
        interactions = []

        # Target-based interactions
        for target in self.targets:
            if "interactions" not in target.assay_details:
                continue

            for interaction in target.assay_details["interactions"]:
                interactions.append(
                    {
                        "type": "target",
                        "target": target.common_name,
                        "interaction": interaction["description"],
                        "severity": interaction.get("severity", "unknown"),
                        "confidence": target.confidence,
                    }
                )

        # Known interactions
        if hasattr(self, "safety_profile"):
            for interaction in self.safety_profile.get("known_interactions", []):
                interactions.append(
                    {
                        "type": "known",
                        "interaction": interaction["description"],
                        "severity": interaction.get("severity", "unknown"),
                        "confidence": interaction.get("confidence", 0.8),
                    }
                )

        return interactions

    def _analyze_contraindications(self) -> List[Dict]:
        """Analyze contraindications."""
        contraindications = []

        # Target-based contraindications
        for target in self.targets:
            if "contraindications" not in target.assay_details:
                continue

            for ci in target.assay_details["contraindications"]:
                contraindications.append(
                    {
                        "type": "target",
                        "target": target.common_name,
                        "contraindication": ci["description"],
                        "severity": ci.get("severity", "unknown"),
                        "confidence": target.confidence,
                    }
                )

        # Known contraindications
        if hasattr(self, "safety_profile"):
            for ci in self.safety_profile.get("contraindications", []):
                contraindications.append(
                    {
                        "type": "known",
                        "contraindication": ci["description"],
                        "severity": ci.get("severity", "unknown"),
                        "confidence": ci.get("confidence", 0.8),
                    }
                )

        return contraindications

    def _get_drug_classes(self) -> Set[str]:
        """Determine compound's drug classes based on receptor profiles."""
        drug_classes = set()

        # Define class detection rules
        class_rules = {
            "SSRI": lambda r: r.startswith("SERT"),
            "MAOI": lambda r: r.startswith("MAO"),
            "DRI": lambda r: r.startswith("DAT"),
            "NMDA_antagonist": lambda r: r.startswith("NMDA"),
            "CNS_depressant": lambda r: r.startswith("GABA"),
        }

        # Apply rules to receptor profiles
        for class_name, rule in class_rules.items():
            if any(rule(r) for r in self.receptor_profiles):
                drug_classes.add(class_name)

        return drug_classes

    def _analyze_interaction_risks(self) -> Dict:
        """Analyze potential drug interaction risks."""
        risks = {}
        max_risk = RiskLevel.UNKNOWN

        # Get compound's drug classes
        drug_classes = self._get_drug_classes()

        # Analyze interactions for each drug class
        for drug_class in drug_classes:
            if drug_class in self.interaction_matrix:
                class_risks = {}
                for target, risk_level in self.interaction_matrix[drug_class].items():
                    class_risks[target] = {
                        "risk_level": risk_level.value,
                        "mechanism": f"{drug_class}_{target}_interaction",
                    }
                    max_risk = max(max_risk, risk_level)

                if class_risks:
                    risks[drug_class] = {
                        "interactions": class_risks,
                        "risk_level": max_risk.value,
                    }

        if risks:
            risks["risk_level"] = max_risk.value

        return risks

    def _analyze_risk_assessment(self) -> Dict:
        """Perform overall risk assessment."""
        return {
            "overall_risk": self._calculate_overall_risk(),
            "risk_factors": self._analyze_risk_factors(),
            "risk_mitigation": self._analyze_risk_mitigation(),
        }

    def _calculate_overall_risk(self) -> Dict:
        """Calculate overall risk level."""
        risks = []

        # Structural risks
        if hasattr(self, "structural_alerts"):
            risks.extend(RiskLevel[a.get("severity", "UNKNOWN").upper()].value for a in self.structural_alerts if a.get("severity", "").upper() in RiskLevel.__members__)

        # Target risks
        for target in self.targets:
            if "severity" in target.assay_details:
                severity = target.assay_details["severity"].upper()
                if severity in RiskLevel.__members__:
                    risks.append(RiskLevel[severity].value)

        # Known risks
        if hasattr(self, "safety_profile"):
            risks.extend(
                RiskLevel[r.get("severity", "UNKNOWN").upper()].value for r in self.safety_profile.get("known_risks", []) if r.get("severity", "").upper() in RiskLevel.__members__
            )

        if not risks:
            return {"level": "unknown", "score": 0.0}

        avg_risk = sum(risks) / len(risks)
        if avg_risk >= 4:
            level = "severe"
        elif avg_risk >= 3:
            level = "high"
        elif avg_risk >= 2:
            level = "moderate"
        else:
            level = "low"

        return {"level": level, "score": avg_risk}

    def _analyze_risk_factors(self) -> List[Dict]:
        """Analyze risk factors."""
        factors = []

        # Structural risk factors
        if hasattr(self, "structural_alerts"):
            for alert in self.structural_alerts:
                factors.append(
                    {
                        "type": "structural",
                        "factor": alert["description"],
                        "severity": alert.get("severity", "unknown"),
                        "confidence": alert.get("confidence", 0.5),
                    }
                )

        # Target risk factors
        for target in self.targets:
            if "risk_factors" in target.assay_details:
                for factor in target.assay_details["risk_factors"]:
                    factors.append(
                        {
                            "type": "target",
                            "target": target.common_name,
                            "factor": factor["description"],
                            "severity": factor.get("severity", "unknown"),
                            "confidence": target.confidence,
                        }
                    )

        # Known risk factors
        if hasattr(self, "safety_profile"):
            for factor in self.safety_profile.get("risk_factors", []):
                factors.append(
                    {
                        "type": "known",
                        "factor": factor["description"],
                        "severity": factor.get("severity", "unknown"),
                        "confidence": factor.get("confidence", 0.8),
                    }
                )

        return factors

    def _analyze_risk_mitigation(self) -> List[Dict]:
        """Analyze risk mitigation strategies."""
        strategies = []

        # Target-based strategies
        for target in self.targets:
            if "risk_mitigation" in target.assay_details:
                for strategy in target.assay_details["risk_mitigation"]:
                    strategies.append(
                        {
                            "type": "target",
                            "target": target.common_name,
                            "strategy": strategy["description"],
                            "effectiveness": strategy.get("effectiveness", "unknown"),
                            "confidence": target.confidence,
                        }
                    )

        # Known strategies
        if hasattr(self, "safety_profile"):
            for strategy in self.safety_profile.get("risk_mitigation", []):
                strategies.append(
                    {
                        "type": "known",
                        "strategy": strategy["description"],
                        "effectiveness": strategy.get("effectiveness", "unknown"),
                        "confidence": strategy.get("confidence", 0.8),
                    }
                )

        return strategies

    def _get_risk_recommendations(self, risk_type: str) -> Dict[str, List[str]]:
        """Get recommendations for a specific risk type."""
        risk_recommendations = {
            "cardiotoxicity": {
                "contraindications": ["Contraindicated in patients with cardiovascular conditions"],
                "monitoring": ["Regular cardiovascular monitoring required"],
            },
            "neurotoxicity": {
                "precautions": ["Use with caution due to potential neurotoxicity"],
                "monitoring": ["Monitor for cognitive and neurological effects"],
            },
            "respiratory": {
                "contraindications": ["Contraindicated with other respiratory depressants"],
                "monitoring": ["Monitor respiratory function"],
            },
            "addiction": {
                "precautions": ["High abuse potential - careful monitoring required"],
                "monitoring": ["Regular assessment of dependence risk"],
            },
        }

        return risk_recommendations.get(risk_type, {})

    def _analyze_risk_patterns(self, risks: Dict) -> List[str]:
        """Analyze patterns in risk data using numpy."""
        if not risks:
            return []

        # Convert risk levels to numerical scores
        risk_scores = [RiskLevel[risk["risk_level"].upper()].value for risk in risks.values() if "risk_level" in risk]

        if not risk_scores:
            return []

        # Calculate risk statistics
        mean_risk = np.mean(risk_scores)
        std_risk = np.std(risk_scores)
        max_risk = np.max(risk_scores)

        # Generate insights
        insights = []
        if max_risk >= RiskLevel.SEVERE.value:
            insights.append("Multiple severe risks detected - extreme caution required")
        if mean_risk > RiskLevel.MODERATE.value:
            insights.append("Overall high risk profile")
        if std_risk < 1.0 and mean_risk > RiskLevel.MODERATE.value:
            insights.append("Consistently high risk across multiple domains")

        return insights

    def _process_risk_recommendations(self, risks: Dict, valid_types: Set[str], recommendations: Dict[str, List[str]]) -> None:
        """Process recommendations for given risks."""
        if not risks:
            return

        for risk_type in risks:
            if risk_type in valid_types:
                risk_recs = self._get_risk_recommendations(risk_type)
                for category, recs in risk_recs.items():
                    recommendations[category].extend(recs)

    def _process_interaction_recommendations(self, interaction_risks: Dict, recommendations: Dict[str, List[str]]) -> None:
        """Process recommendations for interaction risks."""
        if not interaction_risks:
            return

        for drug_class, data in interaction_risks.items():
            if "interactions" in data:
                for target, risk in data["interactions"].items():
                    if risk["risk_level"] == RiskLevel.SEVERE.value:
                        recommendations["contraindications"].append(f"Contraindicated with {target}")
                    elif risk["risk_level"] == RiskLevel.HIGH.value:
                        recommendations["precautions"].append(f"Use with extreme caution with {target}")

    def _generate_safety_recommendations(self) -> Dict:
        """Generate safety recommendations based on risk analysis."""
        recommendations = defaultdict(list)

        # Process binding risks
        binding_risks = self._analyze_binding_risks()
        self._process_risk_recommendations(binding_risks, {"cardiotoxicity", "neurotoxicity"}, recommendations)

        # Process activity risks
        activity_risks = self._analyze_activity_risks()
        self._process_risk_recommendations(activity_risks, {"respiratory", "addiction"}, recommendations)

        # Process interaction risks
        interaction_risks = self._analyze_interaction_risks()
        self._process_interaction_recommendations(interaction_risks, recommendations)

        # Add risk pattern insights
        risk_patterns = self._analyze_risk_patterns({**binding_risks, **activity_risks, **(interaction_risks.get("interactions", {}))})
        if risk_patterns:
            recommendations["risk_patterns"] = risk_patterns

        return dict(recommendations)

    def merge_safety_data(self, other: "SafetyAnalyzer") -> None:
        """Merge safety data from another instance."""
        # Merge alerts with highest risk level
        for alert, level in other.safety_alerts.items():
            if alert not in self.safety_alerts:
                self.safety_alerts[alert] = level
            else:
                self.safety_alerts[alert] = max(self.safety_alerts[alert], level, key=lambda x: x.value)

        # Merge contraindications
        self.contraindications.update(other.contraindications)

        # Merge interaction risks with highest risk level
        for drug, (mechanism, level) in other.interaction_risks.items():
            if drug not in self.interaction_risks:
                self.interaction_risks[drug] = (mechanism, level)
            else:
                _, current_level = self.interaction_risks[drug]
                if level.value > current_level.value:
                    self.interaction_risks[drug] = (mechanism, level)
