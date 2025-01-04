"""Base analysis functionality combining all analysis mixins.

This module provides the CompoundAnalysis class that combines all analysis capabilities:
- Binding analysis
- Activity analysis
- Safety analysis
- Property analysis
- SAR analysis
"""

from dataclasses import dataclass
from typing import Dict

from .activity_analysis import ActivityAnalysisMixin
from .binding_analysis import BindingAnalysisMixin
from .property_analysis import PropertyAnalysisMixin
from .safety_analysis import SafetyAnalysisMixin
from .sar_analysis import SARAnalysisMixin


@dataclass
class CompoundAnalysis(
    BindingAnalysisMixin,
    ActivityAnalysisMixin,
    SafetyAnalysisMixin,
    PropertyAnalysisMixin,
    SARAnalysisMixin,
):
    """Class combining all compound analysis capabilities."""

    def analyze_all(self) -> Dict:
        """Run all analysis methods and return combined results."""
        return {
            "binding": self.analyze_binding(),
            "activity": self.analyze_activity(),
            "safety": self.analyze_safety(),
            "properties": self.analyze_properties(),
            "sar": self.analyze_sar(),
        }

    def get_analysis_summary(self) -> Dict:
        """Get high-level summary of analysis results."""
        if not hasattr(self, "_binding_analysis"):
            self.analyze_all()

        return {
            "binding": {
                "strongest": self._binding_analysis.get("strongest_binding", {}),
                "selectivity": self._binding_analysis.get("target_selectivity", {}).get("selectivity"),
            },
            "activity": {
                "primary": self._activity_analysis.get("primary_activity", {}),
                "confidence": self._activity_analysis.get("activity_confidence"),
            },
            "safety": {
                "risk_level": self._safety_analysis.get("risk_assessment", {}).get("overall_risk", {}).get("level"),
                "alerts": len(self._safety_analysis.get("property_alerts", [])),
            },
            "properties": {
                "drug_likeness": self._property_analysis.get("drug_likeness", {}).get("overall", {}).get("classification"),
                "bioavailability": self._property_analysis.get("bioavailability", {}).get("classification"),
            },
            "sar": {
                "pharmacophores": len(self._sar_analysis.get("pharmacophores", [])),
                "activity_cliffs": len(self._sar_analysis.get("activity_cliffs", [])),
            },
        }

    def get_key_findings(self) -> Dict:
        """Get key findings and potential concerns."""
        if not hasattr(self, "_binding_analysis"):
            self.analyze_all()

        findings = {
            "highlights": [],
            "concerns": [],
            "opportunities": [],
        }

        # Analyze binding data
        strongest = self._binding_analysis.get("strongest_binding")
        if strongest and strongest.get("affinity", 0) < 100:  # High affinity (low nM)
            findings["highlights"].append(
                f"Strong binding to {strongest['target']} ({strongest['affinity']} {strongest['unit']})"
            )

        # Analyze safety
        risk_level = self._safety_analysis.get("risk_assessment", {}).get("overall_risk", {}).get("level")
        if risk_level in ["high", "severe"]:
            findings["concerns"].append(f"High safety risk level: {risk_level}")

        # Analyze properties
        bioavailability = self._property_analysis.get("bioavailability", {}).get("classification")
        if bioavailability in ["low", "very low"]:
            findings["concerns"].append(f"Poor bioavailability: {bioavailability}")

        # Analyze SAR
        activity_cliffs = self._sar_analysis.get("activity_cliffs", [])
        if activity_cliffs:
            cliff = activity_cliffs[0]  # Most significant cliff
            findings["opportunities"].append(
                f"Activity cliff identified vs {cliff['compound']} "
                f"(similarity: {cliff['similarity']:.2f}, activity ratio: {cliff['activity_ratio']:.1f})"
            )

        return findings

    def get_optimization_suggestions(self) -> Dict:
        """Get suggestions for compound optimization."""
        if not hasattr(self, "_binding_analysis"):
            self.analyze_all()

        suggestions = {
            "binding": [],
            "properties": [],
            "safety": [],
        }

        # Binding optimization
        selectivity = self._binding_analysis.get("target_selectivity", {}).get("selectivity")
        if selectivity in ["non-selective", "moderately selective"]:
            suggestions["binding"].append("Improve target selectivity")

        # Property optimization
        limiting_factors = self._property_analysis.get("bioavailability", {}).get("limiting_factors", [])
        for factor in limiting_factors:
            suggestions["properties"].append(f"Address {factor}")

        # Safety optimization
        alerts = self._safety_analysis.get("property_alerts", [])
        for alert in alerts:
            suggestions["safety"].append(
                f"Address {alert['type']}: {alert['description']}"
            )

        return suggestions
