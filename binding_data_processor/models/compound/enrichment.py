"""Web data enrichment for compound data.

This module extends BaseCompound with web enrichment capabilities:
- Patent data integration
- Literature data integration
- Community data integration
- Safety profile enrichment
"""

from dataclasses import dataclass, field
from typing import Dict, List, Set
from datetime import datetime

from .base import BaseCompound
from .types import TargetData


@dataclass
class EnrichedCompound(BaseCompound):
    """CompoundData with web enrichment capabilities."""

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
    dosage_info: Dict[str, Dict] = field(default_factory=dict)  # Route -> Stats
    route_stats: Dict[str, int] = field(default_factory=dict)  # Route -> Count
    duration_stats: Dict[str, int] = field(default_factory=dict)  # Duration -> Count
    common_combinations: List[Dict] = field(default_factory=list)
    risk_factors: List[str] = field(default_factory=list)
    overdose_risks: List[Dict] = field(default_factory=list)
    long_term_risks: List[Dict] = field(default_factory=list)
    approval_status: Dict[str, str] = field(default_factory=dict)  # Region -> Status
    clinical_status: Dict[str, str] = field(default_factory=dict)  # Region -> Status
    research_status: Dict[str, str] = field(default_factory=dict)  # Region -> Status

    # Reference URLs
    pubchem_url: str = "N/A"
    chembl_url: str = "N/A"
    psychonaut_url: str = "N/A"
    erowid_url: str = "N/A"
    wikipedia_url: str = "N/A"
    emcdda_url: str = "N/A"
    isomerdesign_url: str = "N/A"
    nida_url: str = "N/A"
    dea_url: str = "N/A"
    who_url: str = "N/A"

    def merge_web_data(self, other: 'EnrichedCompound') -> None:
        """Merge web-enriched data from another instance."""
        self._merge_patent_data(other)
        self._merge_swiss_data(other)
        self._merge_community_data(other)
        self._merge_literature_data(other)
        self._merge_regulatory_data(other)
        self._merge_experience_reports(other)
        self._merge_safety_profile(other)
        self._merge_status_data(other)

    def _merge_patent_data(self, other: 'EnrichedCompound') -> None:
        """Merge patent data."""
        if other.patent_data:
            if not self.patent_data:
                self.patent_data = {}
            self.patent_data.update(other.patent_data)
            self.patent_count = max(self.patent_count, other.patent_count)

    def _merge_swiss_data(self, other: 'EnrichedCompound') -> None:
        """Merge Swiss tools data."""
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

    def _merge_community_data(self, other: 'EnrichedCompound') -> None:
        """Merge community data."""
        if other.community_data:
            if not self.community_data:
                self.community_data = {}
            self.community_data.update(other.community_data)

        self.experience_reports.extend(
            report for report in other.experience_reports
            if report not in self.experience_reports
        )

        self._merge_dosage_info(other)
        self._merge_route_stats(other)
        self._merge_duration_stats(other)
        self._merge_combinations(other)

    def _merge_dosage_info(self, other: 'EnrichedCompound') -> None:
        """Merge dosage information."""
        for route, stats in other.dosage_info.items():
            if route not in self.dosage_info:
                self.dosage_info[route] = stats
            else:
                # Update stats
                current = self.dosage_info[route]
                current["min"] = min(current["min"], stats["min"])
                current["max"] = max(current["max"], stats["max"])
                current["avg"] = (current["avg"] * current["count"] + 
                                stats["avg"] * stats["count"]) / (
                                    current["count"] + stats["count"]
                                )
                current["count"] += stats["count"]

    def _merge_route_stats(self, other: 'EnrichedCompound') -> None:
        """Merge administration route statistics."""
        for route, count in other.route_stats.items():
            self.route_stats[route] = self.route_stats.get(route, 0) + count

    def _merge_duration_stats(self, other: 'EnrichedCompound') -> None:
        """Merge duration statistics."""
        for duration, count in other.duration_stats.items():
            self.duration_stats[duration] = self.duration_stats.get(duration, 0) + count

    def _merge_combinations(self, other: 'EnrichedCompound') -> None:
        """Merge drug combinations."""
        self.common_combinations.extend(
            combo for combo in other.common_combinations
            if combo not in self.common_combinations
        )

    def _merge_literature_data(self, other: 'EnrichedCompound') -> None:
        """Merge literature data."""
        if other.literature_data:
            if not self.literature_data:
                self.literature_data = {}
            self.literature_data.update(other.literature_data)

    def _merge_regulatory_data(self, other: 'EnrichedCompound') -> None:
        """Merge regulatory data."""
        if other.regulatory_data:
            if not self.regulatory_data:
                self.regulatory_data = {}
            self.regulatory_data.update(other.regulatory_data)

    def _merge_safety_profile(self, other: 'EnrichedCompound') -> None:
        """Merge safety profile data."""
        if other.safety_profile:
            if not self.safety_profile:
                self.safety_profile = {}
            self.safety_profile.update(other.safety_profile)

        self.risk_factors.extend(
            factor for factor in other.risk_factors
            if factor not in self.risk_factors
        )

        self.overdose_risks.extend(
            risk for risk in other.overdose_risks
            if risk not in self.overdose_risks
        )

        self.long_term_risks.extend(
            risk for risk in other.long_term_risks
            if risk not in self.long_term_risks
        )

    def _merge_status_data(self, other: 'EnrichedCompound') -> None:
        """Merge status data."""
        self.approval_status.update(other.approval_status)
        self.clinical_status.update(other.clinical_status)
        self.research_status.update(other.research_status)

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
            "dosage_info": self.dosage_info,
            "route_stats": self.route_stats,
            "duration_stats": self.duration_stats,
            "common_combinations": self.common_combinations,
            "risk_factors": self.risk_factors,
            "overdose_risks": self.overdose_risks,
            "long_term_risks": self.long_term_risks,
            "approval_status": self.approval_status,
            "clinical_status": self.clinical_status,
            "research_status": self.research_status,

            # Reference URLs
            "pubchem_url": self.pubchem_url,
            "chembl_url": self.chembl_url,
            "psychonaut_url": self.psychonaut_url,
            "erowid_url": self.erowid_url,
            "wikipedia_url": self.wikipedia_url,
            "emcdda_url": self.emcdda_url,
            "isomerdesign_url": self.isomerdesign_url,
            "nida_url": self.nida_url,
            "dea_url": self.dea_url,
            "who_url": self.who_url,
        }

    def to_dict(self) -> Dict:
        """Convert compound data to dictionary format."""
        data = super().to_dict()
        data.update(self.get_enrichment_dict())
        return data
