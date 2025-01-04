"""Analysis mixin for compound data.

This module provides the AnalysisMixin class that adds analysis capabilities:
- Binding analysis (strongest binding, patterns, selectivity)
- Activity analysis (primary activity, patterns, mechanisms)
- Safety analysis (toxicity risks, patterns, interactions)
- SAR analysis (pharmacophores, similar compounds, activity cliffs)
- Property analysis (physicochemical properties, drug-likeness)

The mixin is designed to be used with CompoundData to add analysis functionality
while maintaining clean separation of concerns.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from collections import Counter

from .types import (
    TargetData,
    BindingType,
    ActivityType,
    RiskLevel,
    PsychoactiveClass,
    NootropicMechanism,
)


@dataclass
class AnalysisMixin:
    """Mixin class adding analysis capabilities to CompoundData."""

    # Analysis results
    _binding_analysis: Dict = field(default_factory=dict)
    _activity_analysis: Dict = field(default_factory=dict)
    _safety_analysis: Dict = field(default_factory=dict)
    _sar_analysis: Dict = field(default_factory=dict)
    _property_analysis: Dict = field(default_factory=dict)

    def analyze_binding(self) -> Dict:
        """Analyze binding data to identify patterns and key interactions."""
        analysis = {
            "strongest_binding": self._find_strongest_binding(),
            "binding_patterns": self._analyze_binding_patterns(),
            "target_selectivity": self._analyze_target_selectivity(),
            "receptor_families": self._analyze_receptor_families(),
            "binding_confidence": self._analyze_binding_confidence(),
        }
        self._binding_analysis = analysis
        return analysis

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

    def analyze_safety(self) -> Dict:
        """Analyze safety data to identify risks and interactions."""
        analysis = {
            "toxicity_risks": self._analyze_toxicity_risks(),
            "safety_patterns": self._analyze_safety_patterns(),
            "interactions": self._analyze_interactions(),
            "contraindications": self._analyze_contraindications(),
            "risk_assessment": self._analyze_risk_assessment(),
        }
        self._safety_analysis = analysis
        return analysis

    def analyze_sar(self) -> Dict:
        """Analyze structure-activity relationships."""
        analysis = {
            "pharmacophores": self._analyze_pharmacophores(),
            "similar_compounds": self._analyze_similar_compounds(),
            "activity_cliffs": self._analyze_activity_cliffs(),
            "structure_alerts": self._analyze_structure_alerts(),
            "scaffold_analysis": self._analyze_scaffolds(),
        }
        self._sar_analysis = analysis
        return analysis

    def analyze_properties(self) -> Dict:
        """Analyze physicochemical properties and drug-likeness."""
        analysis = {
            "physicochemical": self._analyze_physicochemical(),
            "drug_likeness": self._analyze_drug_likeness(),
            "admet_predictions": self._analyze_admet(),
            "property_alerts": self._analyze_property_alerts(),
            "bioavailability": self._analyze_bioavailability(),
        }
        self._property_analysis = analysis
        return analysis

    def _find_strongest_binding(self) -> Optional[Dict]:
        """Find strongest binding interaction."""
        if not self.targets:
            return None

        strongest = min(
            self.targets,
            key=lambda t: float('inf') if t.affinity_value <= 0 else t.affinity_value
        )

        if strongest.affinity_value <= 0:
            return None

        return {
            "target": strongest.common_name,
            "affinity": strongest.affinity_value,
            "type": strongest.affinity_type,
            "unit": strongest.affinity_unit,
            "activity": strongest.activity_type,
            "confidence": strongest.confidence,
        }

    def _analyze_binding_patterns(self) -> List[Dict]:
        """Analyze patterns in binding data."""
        patterns = []
        
        # Group by receptor family
        family_data = {}
        for target in self.targets:
            family = target.assay_details.get("family", "unknown")
            if family not in family_data:
                family_data[family] = []
            family_data[family].append(target)

        # Analyze patterns within families
        for family, targets in family_data.items():
            if len(targets) < 2:
                continue

            # Calculate average affinity
            affinities = [t.affinity_value for t in targets if t.affinity_value > 0]
            if not affinities:
                continue

            avg_affinity = sum(affinities) / len(affinities)
            
            # Find common activity types
            activities = [t.activity_type for t in targets if t.activity_type != "N/A"]
            if not activities:
                continue

            common_activity = Counter(activities).most_common(1)[0][0]

            patterns.append({
                "family": family,
                "target_count": len(targets),
                "avg_affinity": avg_affinity,
                "common_activity": common_activity,
                "confidence": sum(t.confidence for t in targets) / len(targets),
            })

        return patterns

    def _analyze_target_selectivity(self) -> Dict:
        """Analyze target selectivity."""
        if not self.targets:
            return {"selectivity": "unknown", "score": 0.0}

        # Group targets by family
        family_affinities = {}
        for target in self.targets:
            if target.affinity_value <= 0:
                continue
            family = target.assay_details.get("family", "unknown")
            if family not in family_affinities:
                family_affinities[family] = []
            family_affinities[family].append(target.affinity_value)

        if not family_affinities:
            return {"selectivity": "unknown", "score": 0.0}

        # Calculate selectivity metrics
        selectivity_ratios = []
        for fam1, aff1 in family_affinities.items():
            min_aff1 = min(aff1)
            for fam2, aff2 in family_affinities.items():
                if fam1 == fam2:
                    continue
                min_aff2 = min(aff2)
                ratio = min_aff2 / min_aff1 if min_aff1 > 0 else 0
                selectivity_ratios.append(ratio)

        if not selectivity_ratios:
            return {"selectivity": "unknown", "score": 0.0}

        # Calculate selectivity score
        avg_ratio = sum(selectivity_ratios) / len(selectivity_ratios)
        if avg_ratio >= 100:
            selectivity = "highly selective"
        elif avg_ratio >= 10:
            selectivity = "moderately selective"
        else:
            selectivity = "non-selective"

        return {
            "selectivity": selectivity,
            "score": avg_ratio,
            "ratios": selectivity_ratios,
        }

    def _analyze_receptor_families(self) -> Dict:
        """Analyze receptor family distribution."""
        families = {}
        for target in self.targets:
            family = target.assay_details.get("family", "unknown")
            if family not in families:
                families[family] = {
                    "count": 0,
                    "avg_affinity": 0.0,
                    "activities": Counter(),
                    "confidence": 0.0,
                }
            
            fam_data = families[family]
            fam_data["count"] += 1
            if target.affinity_value > 0:
                fam_data["avg_affinity"] += target.affinity_value
            if target.activity_type != "N/A":
                fam_data["activities"][target.activity_type] += 1
            fam_data["confidence"] += target.confidence

        # Calculate averages
        for fam_data in families.values():
            if fam_data["count"] > 0:
                fam_data["avg_affinity"] /= fam_data["count"]
                fam_data["confidence"] /= fam_data["count"]
                fam_data["activities"] = dict(fam_data["activities"].most_common())

        return families

    def _analyze_binding_confidence(self) -> Dict:
        """Analyze confidence in binding data."""
        if not self.targets:
            return {"overall": 0.0, "by_target": {}}

        confidences = {}
        for target in self.targets:
            confidences[target.common_name] = {
                "confidence": target.confidence,
                "support": len(target.reference_dois),
                "experimental": bool(target.experimental_conditions),
            }

        return {
            "overall": sum(t.confidence for t in self.targets) / len(self.targets),
            "by_target": confidences,
        }

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

    def _analyze_toxicity_risks(self) -> List[Dict]:
        """Analyze toxicity risks."""
        risks = []

        # Structure-based alerts
        if hasattr(self, "structural_alerts"):
            for alert in self.structural_alerts:
                risks.append({
                    "type": "structural",
                    "alert": alert["description"],
                    "severity": alert.get("severity", "unknown"),
                    "confidence": alert.get("confidence", 0.5),
                })

        # Target-based risks
        for target in self.targets:
            if "toxicity" in target.assay_details:
                risks.append({
                    "type": "target",
                    "target": target.common_name,
                    "risk": target.assay_details["toxicity"],
                    "severity": target.assay_details.get("severity", "unknown"),
                    "confidence": target.confidence,
                })

        # Known risks
        if hasattr(self, "safety_profile"):
            for risk in self.safety_profile.get("known_risks", []):
                risks.append({
                    "type": "known",
                    "risk": risk["description"],
                    "severity": risk.get("severity", "unknown"),
                    "confidence": risk.get("confidence", 0.8),
                })

        return risks

    def _analyze_safety_patterns(self) -> Dict:
        """Analyze patterns in safety data."""
        patterns = {
            "structural_patterns": self._analyze_structural_safety(),
            "target_patterns": self._analyze_target_safety(),
            "interaction_patterns": self._analyze_interaction_patterns(),
        }
        return patterns

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
                    "severity": max(
                        (a.get("severity", "unknown") for a in alerts),
                        key=lambda s: RiskLevel[s.upper()].value
                        if s.upper() in RiskLevel.__members__
                        else 0
                    ),
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
                "avg_affinity": sum(
                    t.affinity_value for t in targets
                    if t.affinity_value > 0
                ) / len(targets),
                "confidence": sum(t.confidence for t in targets) / len(targets),
            }
            patterns.append(pattern)

        return patterns

    def _analyze_interactions(self) -> List[Dict]:
        """Analyze potential drug interactions."""
        interactions = []

        # Target-based interactions
        for target in self.targets:
            if "interactions" not in target.assay_details:
                continue

            for interaction in target.assay_details["interactions"]:
                interactions.append({
                    "type": "target",
                    "target": target.common_name,
                    "interaction": interaction["description"],
                    "severity": interaction.get("severity", "unknown"),
                    "confidence": target.confidence,
                })

        # Known interactions
        if hasattr(self, "safety_profile"):
            for interaction in self.safety_profile.get("known_interactions", []):
                interactions.append({
                    "type": "known",
                    "interaction": interaction["description"],
                    "severity": interaction.get("severity", "unknown"),
                    "confidence": interaction.get("confidence", 0.8),
                })

        return interactions

    def _analyze_contraindications(self) -> List[Dict]:
        """Analyze contraindications."""
        contraindications = []

        # Target-based contraindications
        for target in self.targets:
            if "contraindications" not in target.assay_details:
                continue

            for ci in target.assay_details["contraindications"]:
                contraindications.append({
                    "type": "target",
                    "target": target.common_name,
                    "contraindication": ci["description"],
                    "severity": ci.get("severity", "unknown"),
                    "confidence": target.confidence,
                })

        # Known contraindications
        if hasattr(self, "safety_profile"):
            for ci in self.safety_profile.get("contraindications", []):
                contraindications.append({
                    "type": "known",
                    "contraindication": ci["description"],
                    "severity": ci.get("severity", "unknown"),
                    "confidence": ci.get("confidence", 0.8),
                })

        return contraindications

    def _analyze_risk_assessment(self) -> Dict:
        """Perform overall risk assessment."""
        assessment = {
            "overall_risk": self._calculate_overall_risk(),
            "risk_factors": self._analyze_risk_factors(),
            "risk_mitigation": self._analyze_risk_mitigation(),
        }
        return assessment

    def _calculate_overall_risk(self) -> Dict:
        """Calculate overall risk level."""
        risks = []

        # Structural risks
        if hasattr(self, "structural_alerts"):
            risks.extend(
                RiskLevel[a.get("severity", "UNKNOWN").upper()].value
                for a in self.structural_alerts
                if a.get("severity", "").upper() in RiskLevel.__members__
            )

        # Target risks
        for target in self.targets:
            if "severity" in target.assay_details:
                severity = target.assay_details["severity"].upper()
                if severity in RiskLevel.__members__:
                    risks.append(RiskLevel[severity].value)

        # Known risks
        if hasattr(self, "safety_profile"):
            risks.extend(
                RiskLevel[r.get("severity", "UNKNOWN").upper()].value
                for r in self.safety_profile.get("known_risks", [])
                if r.get("severity", "").upper() in RiskLevel.__members__
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
                factors.append({
                    "type": "structural",
                    "factor": alert["description"],
                    "severity": alert.get("severity", "unknown"),
                    "confidence": alert.get("confidence", 0.5),
                })

        # Target risk factors
        for target in self.targets:
            if "risk_factors" in target.assay_details:
                for factor in target.assay_details["risk_factors"]:
                    factors.append({
                        "type": "target",
                        "target": target.common_name,
                        "factor": factor["description"],
                        "severity": factor.get("severity", "unknown"),
                        "confidence": target.confidence,
                    })

        # Known risk factors
        if hasattr(self, "safety_profile"):
            for factor in self.safety_profile.get("risk_factors", []):
                factors.append({
                    "type": "known",
                    "factor": factor["description"],
                    "severity": factor.get("severity", "unknown"),
                    "confidence": factor.get("confidence", 0.8),
                })

        return factors

    def _analyze_risk_mitigation(self) -> List[Dict]:
        """Analyze risk mitigation strategies."""
        strategies = []

        # Target-based strategies
        for target in self.targets:
            if "risk_mitigation" in target.assay_details:
                for strategy in target.assay_details["risk_mitigation"]:
                    strategies.append({
                        "type": "target",
                        "target": target.common_name,
                        "strategy": strategy["description"],
                        "effectiveness": strategy.get("effectiveness", "unknown"),
                        "confidence": target.confidence,
                    })

        # Known strategies
        if hasattr(self, "safety_profile"):
            for strategy in self.safety_profile.get("risk_mitigation", []):
                strategies.append({
                    "type": "known",
                    "strategy": strategy["description"],
                    "effectiveness": strategy.get("effectiveness", "unknown"),
                    "confidence": strategy.get("confidence", 0.8),
                })

        return strategies

    def _analyze_pharmacophores(self) -> List[Dict]:
        """Analyze pharmacophore features."""
        if not hasattr(self, "pharmacophores"):
            return []

        return [
            {
                "type": pharm["type"],
                "description": pharm["description"],
                "score": pharm.get("score", 0.0),
                "confidence": pharm.get("confidence", 0.5),
            }
            for pharm in self.pharmacophores
        ]

    def _analyze_similar_compounds(self) -> List[Dict]:
        """Analyze similar compounds."""
        if not hasattr(self, "similar_compounds"):
            return []

        return [
            {
                "name": comp["name"],
                "similarity": comp["similarity"],
                "shared_features": comp.get("shared_features", []),
                "activity_difference": comp.get("activity_difference", 0.0),
            }
            for comp in self.similar_compounds
        ]

    def _analyze_activity_cliffs(self) -> List[Dict]:
        """Analyze activity cliffs."""
        cliffs = []

        if hasattr(self, "similar_compounds"):
            for comp in self.similar_compounds:
                if comp.get("activity_difference", 0.0) > 1.0:  # Significant difference
                    cliffs.append({
                        "compound": comp["name"],
                        "similarity": comp["similarity"],
                        "activity_difference": comp["activity_difference"],
                        "key_differences": comp.get("key_differences", []),
                    })

        return cliffs

    def _analyze_structure_alerts(self) -> List[Dict]:
        """Analyze structure-based alerts."""
        if not hasattr(self, "structural_alerts"):
            return []

        return [
            {
                "type": alert["type"],
                "description": alert["description"],
                "severity": alert.get("severity", "unknown"),
                "confidence": alert.get("confidence", 0.5),
            }
            for alert in self.structural_alerts
        ]

    def _analyze_scaffolds(self) -> Dict:
        """Analyze molecular scaffolds."""
        if not hasattr(self, "sar_analysis") or "scaffolds" not in self.sar_analysis:
            return {}

        scaffolds = self.sar_analysis["scaffolds"]
        return {
            "core_scaffold": scaffolds.get("core", "unknown"),
            "sub_scaffolds": scaffolds.get("sub_scaffolds", []),
            "scaffold_class": scaffolds.get("class", "unknown"),
            "similar_scaffolds": scaffolds.get("similar", []),
        }

    def _analyze_physicochemical(self) -> Dict:
        """Analyze physicochemical properties."""
        return {
            "molecular_weight": self.molecular_weight,
            "logp": self.logp,
            "hbd": self.hbd,
            "hba": self.hba,
            "tpsa": self.tpsa,
            "rotatable_bonds": self.rotatable_bonds,
            "charge": self.charge,
            "stereocenter_count": self.stereocenter_count,
            "ring_count": self.ring_count,
        }

    def _analyze_drug_likeness(self) -> Dict:
        """Analyze drug-likeness."""
        # Lipinski's Rule of 5
        lipinski = {
            "mw_ok": self.molecular_weight <= 500,
            "logp_ok": self.logp <= 5,
            "hbd_ok": self.hbd <= 5,
            "hba_ok": self.hba <= 10,
            "violations": 0,
        }
        
        lipinski["violations"] = sum(
            1 for ok in [
                lipinski["mw_ok"],
                lipinski["logp_ok"],
                lipinski["hbd_ok"],
                lipinski["hba_ok"],
            ]
            if not ok
        )

        # Veber's rules
        veber = {
            "rotatable_ok": self.rotatable_bonds <= 10,
            "tpsa_ok": self.tpsa <= 140,
        }

        return {
            "lipinski": lipinski,
            "veber": veber,
            "overall": "drug-like" if lipinski["violations"] <= 1 and all(veber.values())
                      else "non-drug-like",
        }

    def _analyze_admet(self) -> Dict:
        """Analyze ADMET predictions."""
        if not hasattr(self, "adme_properties"):
            return {}

        return {
            "absorption": self.adme_properties.get("absorption", {}),
            "distribution": self.adme_properties.get("distribution", {}),
            "metabolism": self.adme_properties.get("metabolism", {}),
            "excretion": self.adme_properties.get("excretion", {}),
            "toxicity": self.adme_properties.get("toxicity", {}),
        }

    def _analyze_property_alerts(self) -> List[Dict]:
        """Analyze property-based alerts."""
        alerts = []

        # Molecular weight alerts
        if self.molecular_weight > 500:
            alerts.append({
                "type": "molecular_weight",
                "description": "High molecular weight may reduce bioavailability",
                "value": self.molecular_weight,
                "threshold": 500,
                "severity": "moderate",
            })

        # LogP alerts
        if self.logp > 5:
            alerts.append({
                "type": "logp",
                "description": "High LogP may cause poor solubility",
                "value": self.logp,
                "threshold": 5,
                "severity": "moderate",
            })

        # TPSA alerts
        if self.tpsa > 140:
            alerts.append({
                "type": "tpsa",
                "description": "High TPSA may reduce membrane permeability",
                "value": self.tpsa,
                "threshold": 140,
                "severity": "moderate",
            })

        return alerts

    def _analyze_bioavailability(self) -> Dict:
        """Analyze predicted bioavailability."""
        # Basic bioavailability score based on properties
        score = 1.0

        # Molecular weight penalty
        if self.molecular_weight > 500:
            score *= 0.8

        # LogP penalty
        if self.logp > 5:
            score *= 0.8

        # TPSA penalty
        if self.tpsa > 140:
            score *= 0.8

        # Rotatable bonds penalty
        if self.rotatable_bonds > 10:
            score *= 0.9

        # H-bond penalties
        if self.hbd > 5:
            score *= 0.9
        if self.hba > 10:
            score *= 0.9

        # Classification
        if score >= 0.85:
            classification = "high"
        elif score >= 0.70:
            classification = "moderate"
        elif score >= 0.50:
            classification = "low"
        else:
            classification = "very low"

        # Absorption factors
        absorption_factors = []
        if self.molecular_weight > 500:
            absorption_factors.append("high molecular weight")
        if self.logp > 5:
            absorption_factors.append("high lipophilicity")
        if self.tpsa > 140:
            absorption_factors.append("high polar surface area")
        if self.rotatable_bonds > 10:
            absorption_factors.append("high flexibility")
        if self.hbd > 5:
            absorption_factors.append("many H-bond donors")
        if self.hba > 10:
            absorption_factors.append("many H-bond acceptors")

        # BBB permeability prediction
        bbb_score = 1.0
        if self.molecular_weight > 400:
            bbb_score *= 0.8
        if self.logp < 0 or self.logp > 6:
            bbb_score *= 0.7
        if self.tpsa > 90:
            bbb_score *= 0.6
        if self.rotatable_bonds > 8:
            bbb_score *= 0.9
        if self.hbd + self.hba > 8:
            bbb_score *= 0.8

        # BBB classification
        if bbb_score >= 0.80:
            bbb_class = "high"
        elif bbb_score >= 0.60:
            bbb_class = "moderate"
        else:
            bbb_class = "low"

        return {
            "score": score,
            "classification": classification,
            "limiting_factors": absorption_factors,
            "bbb_permeability": {
                "score": bbb_score,
                "classification": bbb_class,
            },
            "predictions": {
                "oral": score >= 0.70,
                "intestinal": score >= 0.60,
                "bbb": bbb_score >= 0.60,
            },
            "confidence": min(
                1.0,
                0.9 * (1.0 if hasattr(self, "experimental_data") else 0.7)
            ),
        }
