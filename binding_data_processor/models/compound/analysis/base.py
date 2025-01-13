"""Analysis result functionality.

This module provides the AnalysisMixin class that implements methods for
handling various analysis results including:
- Pharmacophore analysis
- Structural alerts
- Receptor interactions
- SAR analysis
- Binding data analysis
"""

from dataclasses import dataclass, field
from typing import Dict, List

from ..enrichment.web import EnrichedCompound


class AnalysisMixin:
    """Mixin class providing analysis result methods."""

    # Analysis results
    pharmacophores: List[Dict] = field(default_factory=list)
    structural_alerts: List[Dict] = field(default_factory=list)
    receptor_interactions: Dict[str, List[Dict]] = field(default_factory=dict)
    mechanism_predictions: List[Dict] = field(default_factory=list)
    sar_analysis: Dict = field(default_factory=dict)

    def get_analysis_dict(self) -> Dict:
        """Get dictionary of analysis results."""
        return {
            "pharmacophores": self.pharmacophores,
            "structural_alerts": self.structural_alerts,
            "receptor_interactions": self.receptor_interactions,
            "mechanism_predictions": self.mechanism_predictions,
            "sar_analysis": self.sar_analysis,
            "binding_summary": self._summarize_binding_data(),
        }

    def _summarize_binding_data(self) -> Dict:
        """Create a summary of binding data."""
        summary = {
            "total_targets": len(self.targets),
            "strongest_binding": None,
            "primary_targets": [],
            "target_families": set(),
            "binding_distribution": {},
            "activity_types": {},
        }

        for target in self.targets:
            # Track strongest binding
            if target.affinity_value:
                if not summary["strongest_binding"] or target.affinity_value < summary["strongest_binding"]["value"]:
                    summary["strongest_binding"] = {
                        "target": target.common_name,
                        "value": target.affinity_value,
                        "type": target.affinity_type,
                        "unit": target.affinity_unit,
                    }

            # Track primary targets
            if target.is_primary:
                summary["primary_targets"].append(target.common_name)

            # Track target families
            if "family" in target.assay_details:
                summary["target_families"].add(target.assay_details["family"])

            # Track binding distribution
            affinity_range = self._get_affinity_range(target.affinity_value)
            if affinity_range:
                summary["binding_distribution"][affinity_range] = summary["binding_distribution"].get(affinity_range, 0) + 1

            # Track activity types
            if target.activity_type != "N/A":
                summary["activity_types"][target.activity_type] = summary["activity_types"].get(target.activity_type, 0) + 1

        # Convert sets to sorted lists
        summary["target_families"] = sorted(summary["target_families"])

        return summary

    def _get_affinity_range(self, value: float) -> str:
        """Get affinity range category for a value."""
        if not value or value <= 0:
            return None
        elif value < 1:
            return "sub-nM"
        elif value < 10:
            return "1-10 nM"
        elif value < 100:
            return "10-100 nM"
        elif value < 1000:
            return "100-1000 nM"
        elif value < 10000:
            return "1-10 µM"
        else:
            return ">10 µM"

    def merge_analysis_results(self, other: "AnalysisMixin") -> None:
        """Merge analysis results from another instance."""
        self._merge_pharmacophores(other)
        self._merge_structural_alerts(other)
        self._merge_receptor_interactions(other)
        self._merge_mechanism_predictions(other)
        self._merge_sar_analysis(other)

    def _merge_pharmacophores(self, other: "AnalysisMixin") -> None:
        """Merge pharmacophore data."""
        self.pharmacophores.extend(pharm for pharm in other.pharmacophores if pharm not in self.pharmacophores)

    def _merge_structural_alerts(self, other: "AnalysisMixin") -> None:
        """Merge structural alerts."""
        self.structural_alerts.extend(alert for alert in other.structural_alerts if alert not in self.structural_alerts)

    def _merge_receptor_interactions(self, other: "AnalysisMixin") -> None:
        """Merge receptor interaction data."""
        for receptor, interactions in other.receptor_interactions.items():
            if receptor not in self.receptor_interactions:
                self.receptor_interactions[receptor] = []
            self.receptor_interactions[receptor].extend(inter for inter in interactions if inter not in self.receptor_interactions[receptor])

    def _merge_mechanism_predictions(self, other: "AnalysisMixin") -> None:
        """Merge mechanism predictions."""
        self.mechanism_predictions.extend(pred for pred in other.mechanism_predictions if pred not in self.mechanism_predictions)

    def _merge_sar_analysis(self, other: "AnalysisMixin") -> None:
        """Merge SAR analysis data."""
        if not other.sar_analysis:
            return

        if not self.sar_analysis:
            self.sar_analysis = {}

        # Merge each analysis section
        for section, data in other.sar_analysis.items():
            if section not in self.sar_analysis:
                self.sar_analysis[section] = data
            else:
                if isinstance(data, dict):
                    self.sar_analysis[section].update(data)
                elif isinstance(data, list):
                    self.sar_analysis[section].extend(item for item in data if item not in self.sar_analysis[section])
                else:
                    # For scalar values, keep the most recent
                    self.sar_analysis[section] = data


@dataclass
class AnalyzedCompound(EnrichedCompound, AnalysisMixin):
    """Compound with analysis capabilities.

    Combines EnrichedCompound's web enrichment capabilities with
    AnalysisMixin's analysis result methods.
    """

    def __init__(self, smiles: str):
        """Initialize compound with analysis capabilities.

        Args:
            smiles: SMILES string representation
        """
        EnrichedCompound.__init__(self, smiles)
        # AnalysisMixin fields are initialized by dataclass

    def to_dict(self, include_predictions: bool = True) -> Dict:
        """Convert compound data to dictionary format."""
        data = super().to_dict(include_predictions=include_predictions)
        data.update(self.get_analysis_dict())
        return data
