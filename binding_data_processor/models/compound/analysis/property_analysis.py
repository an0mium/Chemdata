"""Property analysis functionality for compound data.

This module provides analysis capabilities for compound properties:
- Physicochemical property analysis
- Drug-likeness analysis
- ADME property analysis
- Structure-property relationships
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
import numpy as np

from ..types import (
    PropertyType,
    DrugLikenessType,
    ADMEType,
)


def default_property_analysis() -> Dict:
    """Default empty property analysis dictionary."""
    return {}


def default_physicochemical_properties() -> Dict[str, float]:
    """Default empty physicochemical properties dictionary."""
    return {}


def default_druglikeness_scores() -> Dict[str, float]:
    """Default empty druglikeness scores dictionary."""
    return {}


def default_adme_properties() -> Dict[str, ADMEType]:
    """Default empty ADME properties dictionary."""
    return {}


def default_property_thresholds() -> Dict[str, float]:
    """Default thresholds for physicochemical properties."""
    return {
        "molecular_weight": 500.0,
        "logp": 5.0,
        "hbd": 5.0,
        "hba": 10.0,
        "tpsa": 140.0,
        "rotatable_bonds": 10.0,
    }


@dataclass
class PropertyAnalyzer:
    """Analyzer for compound property data."""

    # Core property data
    _property_analysis: Dict = field(default_factory=default_property_analysis)
    physicochemical_properties: Dict[str, float] = field(default_factory=default_physicochemical_properties)
    druglikeness_scores: Dict[str, float] = field(default_factory=default_druglikeness_scores)
    adme_properties: Dict[str, ADMEType] = field(default_factory=default_adme_properties)
    property_thresholds: Dict[str, float] = field(default_factory=default_property_thresholds)

    def analyze_properties(self) -> Dict:
        """Analyze property data to identify patterns and relationships."""
        analysis = {
            # Core property analysis
            "physicochemical": self._analyze_physicochemical(),
            "druglikeness": self._analyze_druglikeness(),
            "adme": self._analyze_adme(),
            # Structure-property analysis
            "structure_property": self._analyze_structure_property(),
            "property_relationships": self._analyze_property_relationships(),
            # Overall assessment
            "property_assessment": self._analyze_property_assessment(),
        }
        self._property_analysis = analysis
        return analysis

    def _analyze_physicochemical(self) -> Dict:
        """Analyze physicochemical properties."""
        analysis = {}

        # Group properties by type
        for prop, value in self.physicochemical_properties.items():
            prop_type = PropertyType.get_type(prop)

            if prop_type not in analysis:
                analysis[prop_type] = []

            analysis[prop_type].append(
                {
                    "property": prop,
                    "value": value,
                    "threshold": self.property_thresholds.get(prop),
                }
            )

        return analysis

    def _analyze_druglikeness(self) -> Dict:
        """Analyze drug-likeness scores."""
        analysis = {}

        for rule, score in self.druglikeness_scores.items():
            rule_type = DrugLikenessType.get_type(rule)

            if rule_type not in analysis:
                analysis[rule_type] = []

            analysis[rule_type].append(
                {
                    "rule": rule,
                    "score": score,
                    "threshold": self.property_thresholds.get(rule),
                }
            )

        return analysis

    def _analyze_adme(self) -> Dict:
        """Analyze ADME properties."""
        analysis = {}

        for prop, adme_type in self.adme_properties.items():
            if adme_type not in analysis:
                analysis[adme_type] = []

            analysis[adme_type].append(
                {
                    "property": prop,
                    "type": adme_type.value,
                    "threshold": self.property_thresholds.get(prop),
                }
            )

        return analysis

    def _analyze_structure_property(self) -> Dict:
        """Analyze structure-property relationships."""
        relationships = {}

        # Group by property type
        for prop, value in self.physicochemical_properties.items():
            prop_type = PropertyType.get_type(prop)

            if prop_type not in relationships:
                relationships[prop_type] = []

            relationships[prop_type].append(
                {
                    "property": prop,
                    "value": value,
                    "threshold": self.property_thresholds.get(prop),
                }
            )

        return relationships

    def _analyze_property_relationships(self) -> Dict:
        """Analyze relationships between properties."""
        relationships = {}

        # Calculate correlations between properties
        properties = list(self.physicochemical_properties.keys())
        for i, prop1 in enumerate(properties):
            for prop2 in properties[i + 1 :]:
                val1 = self.physicochemical_properties[prop1]
                val2 = self.physicochemical_properties[prop2]

                if prop1 not in relationships:
                    relationships[prop1] = []

                relationships[prop1].append(
                    {
                        "property": prop2,
                        "correlation": np.corrcoef([val1], [val2])[0, 1],
                    }
                )

        return relationships

    def _analyze_property_assessment(self) -> Dict:
        """Perform overall property assessment."""
        assessment = {
            "violations": self._analyze_violations(),
            "compliance": self._analyze_compliance(),
            "optimization": self._analyze_optimization(),
        }

        # Add overall metrics
        if self.physicochemical_properties:
            assessment["metrics"] = {
                "violation_count": len(assessment["violations"]),
                "compliance_score": len(assessment["compliance"]) / len(self.property_thresholds),
                "optimization_count": len(assessment["optimization"]),
            }

        return assessment

    def _analyze_violations(self) -> List[Dict]:
        """Analyze property threshold violations."""
        violations = []

        for prop, value in self.physicochemical_properties.items():
            threshold = self.property_thresholds.get(prop)
            if threshold is not None and value > threshold:
                violations.append(
                    {
                        "property": prop,
                        "value": value,
                        "threshold": threshold,
                        "excess": value - threshold,
                    }
                )

        return sorted(violations, key=lambda x: x["excess"], reverse=True)

    def _analyze_compliance(self) -> List[Dict]:
        """Analyze property threshold compliance."""
        compliance = []

        for prop, value in self.physicochemical_properties.items():
            threshold = self.property_thresholds.get(prop)
            if threshold is not None and value <= threshold:
                compliance.append(
                    {
                        "property": prop,
                        "value": value,
                        "threshold": threshold,
                        "margin": threshold - value,
                    }
                )

        return sorted(compliance, key=lambda x: x["margin"])

    def _analyze_optimization(self) -> List[Dict]:
        """Analyze properties needing optimization."""
        optimization = []

        # Check physicochemical properties
        for prop, value in self.physicochemical_properties.items():
            threshold = self.property_thresholds.get(prop)
            if threshold is not None:
                if value > threshold * 1.2:  # 20% above threshold
                    optimization.append(
                        {
                            "property": prop,
                            "value": value,
                            "threshold": threshold,
                            "priority": "high",
                        }
                    )
                elif value > threshold:
                    optimization.append(
                        {
                            "property": prop,
                            "value": value,
                            "threshold": threshold,
                            "priority": "medium",
                        }
                    )

        return sorted(optimization, key=lambda x: x["priority"] == "high", reverse=True)

    def merge_property_data(self, other: "PropertyAnalyzer") -> None:
        """Merge property data from another instance."""
        # Merge physicochemical properties
        self.physicochemical_properties.update(other.physicochemical_properties)

        # Merge druglikeness scores
        self.druglikeness_scores.update(other.druglikeness_scores)

        # Merge ADME properties
        self.adme_properties.update(other.adme_properties)

        # Merge property thresholds, keeping stricter values
        for prop, threshold in other.property_thresholds.items():
            if prop not in self.property_thresholds:
                self.property_thresholds[prop] = threshold
            else:
                self.property_thresholds[prop] = min(self.property_thresholds[prop], threshold)
