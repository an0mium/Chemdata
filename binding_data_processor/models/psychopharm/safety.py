"""Safety assessment and risk analysis functionality."""

from dataclasses import field
from typing import Dict, List, Set, Tuple
import numpy as np
from collections import defaultdict

from .base import RiskLevel, PsychoactiveClass, NootropicMechanism


class SafetyProfileMixin:
    """Mixin class providing safety assessment and risk analysis."""

    # Safety data
    safety_alerts: Dict[str, RiskLevel] = field(default_factory=dict)
    contraindications: Set[str] = field(default_factory=set)
    interaction_risks: Dict[str, Tuple[str, RiskLevel]] = field(default_factory=dict)
    
    # Risk thresholds
    risk_thresholds: Dict[str, float] = {
        "cardiotoxicity": 0.7,
        "neurotoxicity": 0.6,
        "hepatotoxicity": 0.8,
        "nephrotoxicity": 0.8,
        "respiratory": 0.7,
        "addiction": 0.6,
        "psychosis": 0.7,
        "serotonin_syndrome": 0.8,
    }
    
    # Drug class interactions
    interaction_matrix: Dict[str, Dict[str, RiskLevel]] = {
        "SSRI": {
            "MAOI": RiskLevel.SEVERE,
            "MDMA": RiskLevel.SEVERE,
            "DRI": RiskLevel.HIGH,
            "TCA": RiskLevel.HIGH,
        },
        "MAOI": {
            "SSRI": RiskLevel.SEVERE,
            "DRI": RiskLevel.SEVERE,
            "TCA": RiskLevel.SEVERE,
            "SRA": RiskLevel.SEVERE,
        },
        "DRI": {
            "MAOI": RiskLevel.SEVERE,
            "SSRI": RiskLevel.HIGH,
            "TCA": RiskLevel.HIGH,
        },
        "NMDA_antagonist": {
            "CNS_depressant": RiskLevel.HIGH,
            "respiratory_depressant": RiskLevel.SEVERE,
        },
        "CNS_depressant": {
            "respiratory_depressant": RiskLevel.SEVERE,
            "NMDA_antagonist": RiskLevel.HIGH,
        },
    }

    def get_safety_dict(self) -> Dict:
        """Get dictionary of safety data."""
        return {
            "alerts": self._format_safety_alerts(),
            "contraindications": list(self.contraindications),
            "interactions": self._format_interactions(),
            "risk_assessment": self._assess_safety_risks(),
            "recommendations": self._generate_safety_recommendations(),
        }

    def _format_safety_alerts(self) -> Dict:
        """Format safety alerts."""
        return {
            alert: level.value
            for alert, level in self.safety_alerts.items()
        }

    def _format_interactions(self) -> Dict:
        """Format interaction risks."""
        return {
            drug: {
                "mechanism": mechanism,
                "risk_level": level.value,
            }
            for drug, (mechanism, level) in self.interaction_risks.items()
        }

    def _assess_safety_risks(self) -> Dict:
        """Assess comprehensive safety risks."""
        risks = {}
        
        # Analyze binding-based risks
        binding_risks = self._analyze_binding_risks()
        if binding_risks:
            risks["binding"] = binding_risks
            
        # Analyze activity-based risks
        activity_risks = self._analyze_activity_risks()
        if activity_risks:
            risks["activity"] = activity_risks
            
        # Analyze interaction risks
        interaction_risks = self._analyze_interaction_risks()
        if interaction_risks:
            risks["interactions"] = interaction_risks
            
        # Calculate overall risk level
        risk_levels = []
        for category in risks.values():
            if "risk_level" in category:
                risk_levels.append(RiskLevel[category["risk_level"].upper()])
        
        if risk_levels:
            max_risk = max(risk_levels, key=lambda x: x.value)
            risks["overall_risk"] = max_risk.value
            
        return risks

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

    def _analyze_class_interactions(
        self,
        drug_class: str,
        max_risk: RiskLevel
    ) -> Tuple[Dict, RiskLevel]:
        """Analyze interactions for a specific drug class."""
        class_risks = {}
        
        if drug_class in self.interaction_matrix:
            for target, risk_level in self.interaction_matrix[drug_class].items():
                class_risks[target] = {
                    "risk_level": risk_level.value,
                    "mechanism": f"{drug_class}_{target}_interaction",
                }
                max_risk = max(max_risk, risk_level)
                
        return class_risks, max_risk

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

    def _get_risk_recommendations(self, risk_type: str) -> List[Dict[str, str]]:
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
        risk_scores = [
            RiskLevel[risk["risk_level"].upper()].value
            for risk in risks.values()
            if "risk_level" in risk
        ]
        
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

    def _process_risk_recommendations(
        self,
        risks: Dict,
        valid_types: Set[str],
        recommendations: Dict[str, List[str]]
    ) -> None:
        """Process recommendations for given risks."""
        if not risks:
            return
            
        for risk_type in risks:
            if risk_type in valid_types:
                risk_recs = self._get_risk_recommendations(risk_type)
                for category, recs in risk_recs.items():
                    recommendations[category].extend(recs)

    def _process_interaction_recommendations(
        self,
        interaction_risks: Dict,
        recommendations: Dict[str, List[str]]
    ) -> None:
        """Process recommendations for interaction risks."""
        if not interaction_risks:
            return
            
        for drug_class, data in interaction_risks.items():
            if "interactions" in data:
                for target, risk in data["interactions"].items():
                    if risk["risk_level"] == RiskLevel.SEVERE.value:
                        recommendations["contraindications"].append(
                            f"Contraindicated with {target}"
                        )
                    elif risk["risk_level"] == RiskLevel.HIGH.value:
                        recommendations["precautions"].append(
                            f"Use with extreme caution with {target}"
                        )

    def _generate_safety_recommendations(self) -> Dict:
        """Generate safety recommendations based on risk analysis."""
        recommendations = defaultdict(list)
        
        # Process binding risks
        binding_risks = self._analyze_binding_risks()
        self._process_risk_recommendations(
            binding_risks,
            {"cardiotoxicity", "neurotoxicity"},
            recommendations
        )
        
        # Process activity risks
        activity_risks = self._analyze_activity_risks()
        self._process_risk_recommendations(
            activity_risks,
            {"respiratory", "addiction"},
            recommendations
        )
        
        # Process interaction risks
        interaction_risks = self._analyze_interaction_risks()
        self._process_interaction_recommendations(interaction_risks, recommendations)
        
        # Add risk pattern insights
        risk_patterns = self._analyze_risk_patterns({
            **binding_risks,
            **activity_risks,
            **(interaction_risks.get("interactions", {}))
        })
        if risk_patterns:
            recommendations["risk_patterns"] = risk_patterns
            
        return dict(recommendations)

    def merge_safety_data(self, other: "SafetyProfileMixin") -> None:
        """Merge safety data from another instance."""
        # Merge alerts with highest risk level
        for alert, level in other.safety_alerts.items():
            if alert not in self.safety_alerts:
                self.safety_alerts[alert] = level
            else:
                self.safety_alerts[alert] = max(
                    self.safety_alerts[alert],
                    level,
                    key=lambda x: x.value
                )
                
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
