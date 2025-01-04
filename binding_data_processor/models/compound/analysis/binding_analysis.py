"""Binding analysis functionality for compound data.

This module provides the BindingAnalysisMixin class that adds binding analysis capabilities:
- Strongest binding identification
- Binding pattern analysis
- Target selectivity analysis
- Receptor family analysis
- Binding confidence assessment
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
from collections import Counter


@dataclass
class BindingAnalysisMixin:
    """Mixin class adding binding analysis capabilities."""

    _binding_analysis: Dict = field(default_factory=dict)

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

    def _get_family_affinities(self) -> Dict[str, List[float]]:
        """Get affinities grouped by receptor family."""
        family_affinities = {}
        for target in self.targets:
            if target.affinity_value <= 0:
                continue
            family = target.assay_details.get("family", "unknown")
            if family not in family_affinities:
                family_affinities[family] = []
            family_affinities[family].append(target.affinity_value)
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

    def _analyze_target_selectivity(self) -> Dict:
        """Analyze target selectivity."""
        if not self.targets:
            return {"selectivity": "unknown", "score": 0.0}

        # Get affinities by family
        family_affinities = self._get_family_affinities()
        if not family_affinities:
            return {"selectivity": "unknown", "score": 0.0}

        # Calculate selectivity ratios
        selectivity_ratios = self._calculate_selectivity_ratios(family_affinities)
        if not selectivity_ratios:
            return {"selectivity": "unknown", "score": 0.0}

        # Calculate selectivity score and classification
        avg_ratio = sum(selectivity_ratios) / len(selectivity_ratios)
        selectivity = self._get_selectivity_classification(avg_ratio)

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
