"""Web data enrichment functionality.

This module provides the WebEnrichmentMixin class that implements methods for
handling web-enriched data from various sources including:
- Patent data
- Literature data
- Community reports
- Safety information
- Regulatory status
"""

from dataclasses import field
from typing import Dict, List


class WebEnrichmentMixin:
    """Mixin class providing web data enrichment methods."""

    # Patent data
    patent_data: Dict = field(default_factory=dict)
    patent_count: int = 0

    # Swiss tools data
    swiss_data: Dict = field(default_factory=dict)
    target_predictions: List[Dict] = field(default_factory=list)
    adme_properties: Dict = field(default_factory=dict)
    similar_compounds: List[Dict] = field(default_factory=list)

    # Web-enriched data
    community_data: Dict = field(default_factory=dict)
    literature_data: Dict = field(default_factory=dict)
    regulatory_data: Dict = field(default_factory=dict)
    experience_reports: List[Dict] = field(default_factory=list)
    safety_profile: Dict = field(default_factory=dict)

    def _summarize_community_data(self) -> Dict:
        """Create a summary of community data."""
        summary = self._init_community_summary()
        
        # Process each report
        for report in self.experience_reports:
            self._process_report_effects(report, summary)
            self._process_report_safety(report, summary)
            self._process_report_combinations(report, summary)
            self._process_report_stats(report, summary)

        # Finalize summary
        self._finalize_community_summary(summary)
        return summary

    def _init_community_summary(self) -> Dict:
        """Initialize community data summary structure."""
        return {
            "total_reports": len(self.experience_reports),
            "reported_effects": set(),
            "reported_mechanisms": set(),
            "safety_concerns": [],
            "common_combinations": [],
            "dosage_info": {},
            "route_stats": {},
            "duration_stats": {},
        }

    def _process_report_effects(self, report: Dict, summary: Dict) -> None:
        """Process effects and mechanisms from a report."""
        if "effects" in report:
            summary["reported_effects"].update(report["effects"])
        if "mechanisms" in report:
            summary["reported_mechanisms"].update(report["mechanisms"])

    def _process_report_safety(self, report: Dict, summary: Dict) -> None:
        """Process safety concerns from a report."""
        if "safety_concerns" in report:
            new_concerns = [
                concern for concern in report["safety_concerns"]
                if concern not in summary["safety_concerns"]
            ]
            summary["safety_concerns"].extend(new_concerns)

    def _process_report_combinations(self, report: Dict, summary: Dict) -> None:
        """Process drug combinations from a report."""
        if "combinations" in report:
            for combo in report["combinations"]:
                if combo not in summary["common_combinations"]:
                    summary["common_combinations"].append(combo)

    def _process_report_stats(self, report: Dict, summary: Dict) -> None:
        """Process dosage, route, and duration statistics from a report."""
        self._process_dosage_info(report, summary)
        self._process_route_stats(report, summary)
        self._process_duration_stats(report, summary)

    def _process_dosage_info(self, report: Dict, summary: Dict) -> None:
        """Process dosage information from a report."""
        if "dosage" in report:
            route = report["dosage"].get("route", "unknown")
            dose = report["dosage"].get("amount")
            if route and dose:
                if route not in summary["dosage_info"]:
                    summary["dosage_info"][route] = []
                summary["dosage_info"][route].append(dose)

    def _process_route_stats(self, report: Dict, summary: Dict) -> None:
        """Process administration route statistics from a report."""
        if "route" in report:
            route = report["route"]
            summary["route_stats"][route] = summary["route_stats"].get(route, 0) + 1

    def _process_duration_stats(self, report: Dict, summary: Dict) -> None:
        """Process duration statistics from a report."""
        if "duration" in report:
            duration = report["duration"]
            summary["duration_stats"][duration] = summary["duration_stats"].get(duration, 0) + 1

    def _finalize_community_summary(self, summary: Dict) -> None:
        """Finalize community data summary by sorting and calculating statistics."""
        # Convert sets to sorted lists
        summary["reported_effects"] = sorted(summary["reported_effects"])
        summary["reported_mechanisms"] = sorted(summary["reported_mechanisms"])
        
        # Sort combinations by frequency
        self._sort_combinations(summary)
        
        # Calculate dosage statistics
        self._calculate_dosage_stats(summary)

    def _sort_combinations(self, summary: Dict) -> None:
        """Sort combinations by frequency."""
        def get_combo_frequency(combo):
            return sum(
                1 for r in self.experience_reports
                if "combinations" in r and combo in r["combinations"]
            )
        
        summary["common_combinations"].sort(
            key=get_combo_frequency,
            reverse=True
        )

    def _calculate_dosage_stats(self, summary: Dict) -> None:
        """Calculate statistics for dosage information."""
        for route in summary["dosage_info"]:
            doses = summary["dosage_info"][route]
            summary["dosage_info"][route] = {
                "min": min(doses),
                "max": max(doses),
                "avg": sum(doses) / len(doses),
                "count": len(doses)
            }

    def _get_safety_summary(self) -> Dict:
        """Get combined safety summary."""
        summary = self._init_safety_summary()
        self._add_structural_alerts(summary)
        self._add_toxicity_warnings(summary)
        self._add_abuse_warnings(summary)
        self._add_safety_profile_data(summary)
        return summary

    def _init_safety_summary(self) -> Dict:
        """Initialize safety summary structure."""
        return {
            "alerts": [],
            "warnings": [],
            "contraindications": [],
            "interactions": [],
            "risk_factors": [],
            "overdose_risks": [],
            "long_term_risks": [],
        }

    def _add_structural_alerts(self, summary: Dict) -> None:
        """Add structural alerts to safety summary."""
        if hasattr(self, 'structural_alerts'):
            summary["alerts"].extend(
                alert["description"] for alert in self.structural_alerts
            )

    def _add_toxicity_warnings(self, summary: Dict) -> None:
        """Add toxicity warnings to safety summary."""
        if hasattr(self, 'toxicity_predictions'):
            for tox_type, tox_data in self.toxicity_predictions.items():
                if tox_data.get("probability", 0) > 0.7:
                    summary["warnings"].append({
                        "type": tox_type,
                        "probability": tox_data["probability"],
                        "severity": tox_data.get("severity", "unknown"),
                        "mechanisms": tox_data.get("mechanisms", [])
                    })

    def _add_abuse_warnings(self, summary: Dict) -> None:
        """Add abuse potential warnings to safety summary."""
        if hasattr(self, 'abuse_potential'):
            for abuse_type, abuse_data in self.abuse_potential.items():
                if abuse_data.get("probability", 0) > 0.7:
                    summary["warnings"].append({
                        "type": f"abuse_{abuse_type}",
                        "probability": abuse_data["probability"],
                        "risk_level": abuse_data.get("risk_level", "unknown"),
                        "mechanisms": abuse_data.get("mechanisms", [])
                    })

    def _add_safety_profile_data(self, summary: Dict) -> None:
        """Add safety profile data to safety summary."""
        if self.safety_profile:
            for key in ["contraindications", "interactions", "risk_factors",
                       "overdose_risks", "long_term_risks"]:
                summary[key].extend(self.safety_profile.get(key, []))

    def _get_regulatory_status(self) -> Dict:
        """Get combined regulatory status."""
        status = self._init_regulatory_status()
        self._add_regulatory_data(status)
        self._add_toxicity_status(status)
        self._add_abuse_status(status)
        return status

    def _init_regulatory_status(self) -> Dict:
        """Initialize regulatory status structure."""
        return {
            "scheduling": None,
            "control_status": None,
            "warnings": [],
            "restrictions": [],
            "approval_status": {},
            "clinical_status": {},
            "research_status": {},
        }

    def _add_regulatory_data(self, status: Dict) -> None:
        """Add regulatory data to status."""
        if self.regulatory_data:
            status.update(self.regulatory_data)

    def _add_toxicity_status(self, status: Dict) -> None:
        """Add toxicity warnings to regulatory status."""
        if hasattr(self, 'toxicity_predictions'):
            for tox_type, tox_data in self.toxicity_predictions.items():
                if tox_data.get("probability", 0) > 0.7:
                    status["warnings"].append({
                        "type": f"toxicity_{tox_type}",
                        "probability": tox_data["probability"],
                        "severity": tox_data.get("severity", "unknown")
                    })

    def _add_abuse_status(self, status: Dict) -> None:
        """Add abuse potential warnings to regulatory status."""
        if hasattr(self, 'abuse_potential'):
            for abuse_type, abuse_data in self.abuse_potential.items():
                if abuse_data.get("probability", 0) > 0.7:
                    status["warnings"].append({
                        "type": f"abuse_{abuse_type}",
                        "probability": abuse_data["probability"],
                        "risk_level": abuse_data.get("risk_level", "unknown")
                    })

    def get_enrichment_dict(self) -> Dict:
        """Get dictionary of web enrichment data."""
        return {
            # Patent data
            "patent_data": self.patent_data,
            "patent_count": self.patent_count,

            # Swiss data
            "swiss_data": self.swiss_data,
            "target_predictions": self.target_predictions,
            "adme_properties": self.adme_properties,
            "similar_compounds": self.similar_compounds,

            # Web data
            "community_data": self.community_data,
            "literature_data": self.literature_data,
            "regulatory_data": self.regulatory_data,
            "experience_reports": self.experience_reports,
            "safety_profile": self.safety_profile,
        }

    def merge_web_data(self, other: "WebEnrichmentMixin") -> None:
        """Merge web-enriched data from another instance."""
        self._merge_patent_data(other)
        self._merge_swiss_data(other)
        self._merge_community_data(other)
        self._merge_literature_data(other)
        self._merge_regulatory_data(other)
        self._merge_experience_reports(other)
        self._merge_safety_profile(other)

    def _merge_patent_data(self, other: "WebEnrichmentMixin") -> None:
        """Merge patent data."""
        if other.patent_data:
            if not self.patent_data:
                self.patent_data = {}
            self.patent_data.update(other.patent_data)
            self.patent_count = max(self.patent_count, other.patent_count)

    def _merge_swiss_data(self, other: "WebEnrichmentMixin") -> None:
        """Merge Swiss data."""
        if other.swiss_data:
            if not self.swiss_data:
                self.swiss_data = {}
            self.swiss_data.update(other.swiss_data)

        self.target_predictions.extend(
            pred for pred in other.target_predictions
            if pred not in self.target_predictions
        )

        if other.adme_properties:
            if not self.adme_properties:
                self.adme_properties = {}
            self.adme_properties.update(other.adme_properties)

        self.similar_compounds.extend(
            comp for comp in other.similar_compounds
            if comp not in self.similar_compounds
        )

    def _merge_community_data(self, other: "WebEnrichmentMixin") -> None:
        """Merge community data."""
        if not other.community_data:
            return

        if not self.community_data:
            self.community_data = {}

        for key, value in other.community_data.items():
            self._merge_data_field(self.community_data, key, value)

    def _merge_literature_data(self, other: "WebEnrichmentMixin") -> None:
        """Merge literature data."""
        if other.literature_data:
            if not self.literature_data:
                self.literature_data = {}
            self.literature_data.update(other.literature_data)

    def _merge_regulatory_data(self, other: "WebEnrichmentMixin") -> None:
        """Merge regulatory data."""
        if other.regulatory_data:
            if not self.regulatory_data:
                self.regulatory_data = {}
            self.regulatory_data.update(other.regulatory_data)

    def _merge_experience_reports(self, other: "WebEnrichmentMixin") -> None:
        """Merge experience reports."""
        self.experience_reports.extend(
            report for report in other.experience_reports
            if report not in self.experience_reports
        )

    def _merge_safety_profile(self, other: "WebEnrichmentMixin") -> None:
        """Merge safety profile data."""
        if not other.safety_profile:
            return

        if not self.safety_profile:
            self.safety_profile = {}

        for key, value in other.safety_profile.items():
            self._merge_data_field(self.safety_profile, key, value)

    def _merge_data_field(self, target: Dict, key: str, value) -> None:
        """Merge a single data field based on its type."""
        if isinstance(value, (list, set)):
            if key not in target:
                target[key] = type(value)()
            if isinstance(value, list):
                target[key].extend(
                    item for item in value
                    if item not in target[key]
                )
            else:  # set
                target[key].update(value)
        elif isinstance(value, dict):
            if key not in target:
                target[key] = {}
            target[key].update(value)
        else:
            target[key] = value
