"""Safety analysis functionality for compound data.

This module provides the SafetyAnalysisMixin class that adds safety analysis capabilities:
- Toxicity risk analysis
- Safety pattern analysis
- Drug interaction analysis
- Contraindication analysis
- Risk assessment and mitigation
"""

from dataclasses import dataclass, field
from typing import Dict, List

from ..types import RiskLevel


@dataclass
class SafetyAnalysisMixin:
    """Mixin class adding safety analysis capabilities."""

    _safety_analysis: Dict = field(default_factory=dict)

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
                "severity": max(
                    (i.get("severity", "unknown") for i in interactions),
                    key=lambda s: RiskLevel[s.upper()].value
                    if s.upper() in RiskLevel.__members__
                    else 0
                ),
                "confidence": sum(i.get("confidence", 0.5) for i in interactions) / len(interactions),
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
