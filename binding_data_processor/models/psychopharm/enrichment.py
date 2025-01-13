"""Web data enrichment functionality for psychopharmacological compounds.

This module provides web enrichment capabilities through:
1. WebEnrichmentMixin - Core web data processing functionality

Key features:
- Patent data integration
- Literature data integration
- Community data integration
- Social media monitoring
- Web data validation
- Data merging and export
"""

from dataclasses import dataclass, field
from datetime import datetime
import json
from typing import Dict, List, Optional, Set

from .base import PsychopharmBase, RiskLevel
from .types import (
    StringSet,
    StringDict,
    ValidationErrors,
    OptionalStr,
    OptionalDict,
)


class WebEnrichmentMixin:
    """Mixin providing web data enrichment capabilities."""

    # Web data sources
    web_sources: Dict[str, Dict] = field(default_factory=dict)
    source_timestamps: Dict[str, datetime] = field(default_factory=dict)

    # Literature data
    pubmed_citations: List[Dict] = field(default_factory=list)
    patent_citations: List[Dict] = field(default_factory=list)

    # Community data
    psychonaut_data: Optional[Dict] = field(default=None)
    erowid_data: Optional[Dict] = field(default=None)
    tripsit_data: Optional[Dict] = field(default=None)

    # Social data
    reddit_mentions: List[Dict] = field(default_factory=list)
    twitter_mentions: List[Dict] = field(default_factory=list)

    # Regulatory data
    regulatory_status: Dict[str, str] = field(default_factory=dict)
    scheduling_info: Dict[str, str] = field(default_factory=dict)

    def get_enrichment_dict(self) -> Dict:
        """Get dictionary of enriched data."""
        return {
            "sources": {
                source: {
                    "data": data,
                    "timestamp": self.source_timestamps[source].isoformat() if source in self.source_timestamps else None,
                }
                for source, data in self.web_sources.items()
            },
            "literature": {
                "pubmed": self.pubmed_citations,
                "patents": self.patent_citations,
            },
            "community": {
                "psychonaut": self.psychonaut_data,
                "erowid": self.erowid_data,
                "tripsit": self.tripsit_data,
            },
            "social": {
                "reddit": self.reddit_mentions,
                "twitter": self.twitter_mentions,
            },
            "regulatory": {
                "status": self.regulatory_status,
                "scheduling": self.scheduling_info,
            },
        }

    def _update_source_data(self, source: str, data: Dict, timestamp: Optional[datetime] = None) -> None:
        """Update data from a specific source."""
        self.web_sources[source] = data
        if timestamp:
            self.source_timestamps[source] = timestamp
        else:
            self.source_timestamps[source] = datetime.now()

    def _extract_classification(self, web_data: Dict) -> Optional[Dict]:
        """Extract psychoactive classification from web data."""
        classification = {}
        confidence = 0.0

        # Check community sources
        if "psychonaut" in web_data.get("effects", {}):
            effects = web_data["effects"]["psychonaut"]
            if "class" in effects:
                classification["primary_class"] = effects["class"]
                classification["confidence"] = effects.get("confidence", 0.7)
                classification["source"] = "psychonaut"
                confidence = classification["confidence"]

        # Check literature sources
        if confidence < 0.8 and "pubmed" in web_data.get("literature", {}):
            for citation in web_data["literature"]["pubmed"]:
                if "classification" in citation:
                    new_conf = citation.get("confidence", 0.6)
                    if new_conf > confidence:
                        classification["primary_class"] = citation["classification"]
                        classification["confidence"] = new_conf
                        classification["source"] = "pubmed"
                        confidence = new_conf

        # Check regulatory sources
        if confidence < 0.9 and "scheduling" in web_data.get("regulatory", {}):
            scheduling = web_data["regulatory"]["scheduling"]
            if "class" in scheduling:
                classification["primary_class"] = scheduling["class"]
                classification["confidence"] = 0.9
                classification["source"] = "regulatory"

        return classification if classification else None

    def _extract_mechanisms(self, web_data: Dict) -> Dict[str, float]:
        """Extract mechanism information from web data."""
        mechanisms = {}

        # Check literature sources
        if "pubmed" in web_data.get("literature", {}):
            for citation in web_data["literature"]["pubmed"]:
                if "mechanisms" in citation:
                    for mech, conf in citation["mechanisms"].items():
                        if mech not in mechanisms or conf > mechanisms[mech]:
                            mechanisms[mech] = conf

        # Check community sources
        if "psychonaut" in web_data.get("effects", {}):
            effects = web_data["effects"]["psychonaut"]
            if "mechanisms" in effects:
                for mech, conf in effects["mechanisms"].items():
                    if mech not in mechanisms or conf > mechanisms[mech]:
                        mechanisms[mech] = conf

        return mechanisms

    def _extract_risks_from_source(self, source_data: Dict, existing_risks: Dict[str, str]) -> None:
        """Extract risk information from a data source."""
        if "risks" in source_data:
            for risk, level in source_data["risks"].items():
                try:
                    risk_level = RiskLevel[level.upper()]
                    if risk not in existing_risks or risk_level.value > RiskLevel[existing_risks[risk].upper()].value:
                        existing_risks[risk] = level
                except KeyError:
                    continue

    def _extract_safety_from_community(self, web_data: Dict, safety_data: Dict) -> None:
        """Extract safety information from community sources."""
        for source in ["psychonaut", "erowid", "tripsit"]:
            if source in web_data.get("community", {}):
                data = web_data["community"][source]

                # Extract risks
                self._extract_risks_from_source(data, safety_data["risks"])

                # Extract contraindications
                if "contraindications" in data:
                    safety_data["contraindications"].update(data["contraindications"])

                # Extract interactions
                if "interactions" in data:
                    for substance, interaction in data["interactions"].items():
                        if substance not in safety_data["interactions"]:
                            safety_data["interactions"][substance] = interaction

    def _extract_safety_from_literature(self, web_data: Dict, safety_data: Dict) -> None:
        """Extract safety information from literature sources."""
        if "pubmed" in web_data.get("literature", {}):
            for citation in web_data["literature"]["pubmed"]:
                if "safety" in citation:
                    safety = citation["safety"]

                    # Extract risks
                    self._extract_risks_from_source(safety, safety_data["risks"])

                    # Extract contraindications
                    if "contraindications" in safety:
                        safety_data["contraindications"].update(safety["contraindications"])

    def _extract_safety_data(self, web_data: Dict) -> Dict:
        """Extract safety information from web data."""
        safety_data = {
            "risks": {},
            "contraindications": set(),
            "interactions": {},
        }

        # Process community sources
        self._extract_safety_from_community(web_data, safety_data)

        # Process literature sources
        self._extract_safety_from_literature(web_data, safety_data)

        return safety_data

    def _update_literature_data(self, web_data: Dict) -> None:
        """Update literature-related data."""
        if "pubmed" in web_data.get("literature", {}):
            self.pubmed_citations.extend(web_data["literature"]["pubmed"])

        if "patents" in web_data.get("literature", {}):
            self.patent_citations.extend(web_data["literature"]["patents"])

    def _update_community_data(self, web_data: Dict) -> None:
        """Update community-related data."""
        if "psychonaut" in web_data.get("community", {}):
            self.psychonaut_data = web_data["community"]["psychonaut"]

        if "erowid" in web_data.get("community", {}):
            self.erowid_data = web_data["community"]["erowid"]

        if "tripsit" in web_data.get("community", {}):
            self.tripsit_data = web_data["community"]["tripsit"]

    def _update_social_data(self, web_data: Dict) -> None:
        """Update social media data."""
        if "reddit" in web_data.get("social", {}):
            self.reddit_mentions.extend(web_data["social"]["reddit"])

        if "twitter" in web_data.get("social", {}):
            self.twitter_mentions.extend(web_data["social"]["twitter"])

    def _update_regulatory_data(self, web_data: Dict) -> None:
        """Update regulatory data."""
        if "regulatory" in web_data:
            if "status" in web_data["regulatory"]:
                self.regulatory_status.update(web_data["regulatory"]["status"])
            if "scheduling" in web_data["regulatory"]:
                self.scheduling_info.update(web_data["regulatory"]["scheduling"])

    def update_classification(self, classification: Dict) -> None:
        """Update compound classification from web data."""
        if "primary_class" in classification:
            self.web_sources["classification"] = classification

    def update_mechanisms(self, mechanisms: Dict[str, float]) -> None:
        """Update mechanism data from web sources."""
        if mechanisms:
            self.web_sources["mechanisms"] = mechanisms

    def update_safety_data(self, safety_data: Dict) -> None:
        """Update safety data from web sources."""
        if safety_data:
            self.web_sources["safety"] = safety_data

    def update_from_web_sources(self, web_data: Dict, source: str, timestamp: Optional[datetime] = None) -> None:
        """Update compound data from web sources."""
        # Update source tracking
        self._update_source_data(source, web_data, timestamp)

        # Extract and update classification
        classification = self._extract_classification(web_data)
        if classification:
            self.update_classification(classification)

        # Extract and update mechanisms
        mechanisms = self._extract_mechanisms(web_data)
        if mechanisms:
            self.update_mechanisms(mechanisms)

        # Extract and update safety data
        safety_data = self._extract_safety_data(web_data)
        if safety_data:
            self.update_safety_data(safety_data)

        # Update specific data types
        self._update_literature_data(web_data)
        self._update_community_data(web_data)
        self._update_social_data(web_data)
        self._update_regulatory_data(web_data)

    def merge_enrichment_data(self, other: "WebEnrichmentMixin") -> None:
        """Merge enrichment data from another instance."""
        # Merge source data
        for source, data in other.web_sources.items():
            if source not in self.web_sources:
                self.web_sources[source] = data
                self.source_timestamps[source] = other.source_timestamps.get(source, datetime.now())
            elif source in other.source_timestamps and source in self.source_timestamps and other.source_timestamps[source] > self.source_timestamps[source]:
                self.web_sources[source] = data
                self.source_timestamps[source] = other.source_timestamps[source]

        # Merge literature data
        self.pubmed_citations.extend(other.pubmed_citations)
        self.patent_citations.extend(other.patent_citations)

        # Merge community data (take most recent)
        if other.psychonaut_data:
            self.psychonaut_data = other.psychonaut_data
        if other.erowid_data:
            self.erowid_data = other.erowid_data
        if other.tripsit_data:
            self.tripsit_data = other.tripsit_data

        # Merge social data
        self.reddit_mentions.extend(other.reddit_mentions)
        self.twitter_mentions.extend(other.twitter_mentions)

        # Merge regulatory data
        self.regulatory_status.update(other.regulatory_status)
        self.scheduling_info.update(other.scheduling_info)

    def to_json(self) -> str:
        """Convert enrichment data to JSON string."""
        data = self.get_enrichment_dict()

        # Convert datetime objects to ISO format strings
        for source in data["sources"]:
            if "timestamp" in data["sources"][source]:
                data["sources"][source]["timestamp"] = data["sources"][source]["timestamp"]

        return json.dumps(data, indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> "WebEnrichmentMixin":
        """Create instance from JSON string."""
        data = json.loads(json_str)
        instance = cls()

        # Convert ISO format strings back to datetime objects
        for source, source_data in data["sources"].items():
            if "timestamp" in source_data:
                instance.source_timestamps[source] = datetime.fromisoformat(source_data["timestamp"])
            instance.web_sources[source] = source_data["data"]

        # Load literature data
        instance.pubmed_citations = data["literature"]["pubmed"]
        instance.patent_citations = data["literature"]["patents"]

        # Load community data
        instance.psychonaut_data = data["community"]["psychonaut"]
        instance.erowid_data = data["community"]["erowid"]
        instance.tripsit_data = data["community"]["tripsit"]

        # Load social data
        instance.reddit_mentions = data["social"]["reddit"]
        instance.twitter_mentions = data["social"]["twitter"]

        # Load regulatory data
        instance.regulatory_status = data["regulatory"]["status"]
        instance.scheduling_info = data["regulatory"]["scheduling"]

        return instance


@dataclass
class EnrichedCompoundData(PsychopharmBase, WebEnrichmentMixin):
    """Compound data enriched with web information."""

    # Patent data
    patent_numbers: StringSet = field(default_factory=set)
    patent_titles: StringDict = field(default_factory=dict)  # number -> title
    patent_abstracts: StringDict = field(default_factory=dict)  # number -> abstract
    patent_claims: Dict[str, List[str]] = field(default_factory=dict)  # number -> claims
    patent_citations: Dict[str, List[str]] = field(default_factory=dict)  # number -> cited by

    # Literature data
    pubmed_ids: StringSet = field(default_factory=set)
    paper_titles: StringDict = field(default_factory=dict)  # pmid -> title
    paper_abstracts: StringDict = field(default_factory=dict)  # pmid -> abstract
    paper_citations: Dict[str, List[str]] = field(default_factory=dict)  # pmid -> cited by
    paper_keywords: Dict[str, List[str]] = field(default_factory=dict)  # pmid -> keywords

    # Community data
    psychonaut_url: OptionalStr = None
    psychonaut_data: OptionalDict = field(default_factory=dict)
    erowid_url: OptionalStr = None
    erowid_data: OptionalDict = field(default_factory=dict)
    tripsit_url: OptionalStr = None
    tripsit_data: OptionalDict = field(default_factory=dict)

    # Social media data
    reddit_mentions: List[Dict] = field(default_factory=list)  # [{subreddit, title, url, score, date}]
    twitter_mentions: List[Dict] = field(default_factory=list)  # [{user, text, url, date}]
    bluesky_mentions: List[Dict] = field(default_factory=list)  # [{user, text, url, date}]
    discord_mentions: List[Dict] = field(default_factory=list)  # [{server, channel, text, date}]

    # Web data metadata
    last_enriched: str = field(default_factory=lambda: datetime.now().isoformat())
    enrichment_sources: StringSet = field(default_factory=set)
    enrichment_stats: Dict[str, int] = field(default_factory=dict)  # source -> count

    def __post_init__(self):
        """Initialize and validate web enrichment data."""
        super().__post_init__()
        self._validate_web_data()

    def _validate_web_data(self) -> None:
        """Validate web enrichment data."""
        errors = []

        # Run web data validation checks
        errors.extend(self._validate_patents())
        errors.extend(self._validate_literature())
        errors.extend(self._validate_community())
        errors.extend(self._validate_social())

        if errors:
            raise ValidationError("\n".join(errors))

    def _validate_patents(self) -> ValidationErrors:
        """Validate patent data."""
        errors = []

        # Validate patent numbers
        for number in self.patent_numbers:
            if not self._validate_patent_number(number):
                errors.append(f"Invalid patent number format: {number}")

        # Validate patent data consistency
        for number in self.patent_numbers:
            if number not in self.patent_titles:
                errors.append(f"Missing title for patent: {number}")
            if number not in self.patent_abstracts:
                errors.append(f"Missing abstract for patent: {number}")

        return errors

    def _validate_literature(self) -> ValidationErrors:
        """Validate literature data."""
        errors = []

        # Validate PubMed IDs
        for pmid in self.pubmed_ids:
            if not pmid.isdigit():
                errors.append(f"Invalid PubMed ID format: {pmid}")

        # Validate literature data consistency
        for pmid in self.pubmed_ids:
            if pmid not in self.paper_titles:
                errors.append(f"Missing title for paper: {pmid}")
            if pmid not in self.paper_abstracts:
                errors.append(f"Missing abstract for paper: {pmid}")

        return errors

    def _validate_community(self) -> ValidationErrors:
        """Validate community data."""
        errors = []

        # Validate URLs
        if self.psychonaut_url and not self._validate_url(self.psychonaut_url):
            errors.append(f"Invalid PsychonautWiki URL: {self.psychonaut_url}")
        if self.erowid_url and not self._validate_url(self.erowid_url):
            errors.append(f"Invalid Erowid URL: {self.erowid_url}")
        if self.tripsit_url and not self._validate_url(self.tripsit_url):
            errors.append(f"Invalid TripSit URL: {self.tripsit_url}")

        return errors

    def _validate_social(self) -> ValidationErrors:
        """Validate social media data."""
        errors = []

        # Validate Reddit mentions
        for mention in self.reddit_mentions:
            if not all(k in mention for k in ["subreddit", "title", "url", "score", "date"]):
                errors.append(f"Invalid Reddit mention format: {mention}")

        # Validate Twitter mentions
        for mention in self.twitter_mentions:
            if not all(k in mention for k in ["user", "text", "url", "date"]):
                errors.append(f"Invalid Twitter mention format: {mention}")

        return errors

    def _validate_patent_number(self, number: str) -> bool:
        """Validate patent number format."""
        # Basic validation - can be enhanced
        return bool(number and len(number) >= 6)

    def _validate_url(self, url: str) -> bool:
        """Validate URL format."""
        # Basic validation - can be enhanced
        return url.startswith(("http://", "https://"))

    def merge_web_data(self, other: "EnrichedCompoundData", source: str = "merge") -> None:
        """Merge web enrichment data from another instance."""
        # Merge base web enrichment data
        super().merge_enrichment_data(other)

        # Merge patent data
        self.patent_numbers.update(other.patent_numbers)
        self.patent_titles.update(other.patent_titles)
        self.patent_abstracts.update(other.patent_abstracts)
        self.patent_claims.update(other.patent_claims)
        self.patent_citations.update(other.patent_citations)

        # Merge literature data
        self.pubmed_ids.update(other.pubmed_ids)
        self.paper_titles.update(other.paper_titles)
        self.paper_abstracts.update(other.paper_abstracts)
        self.paper_citations.update(other.paper_citations)
        self.paper_keywords.update(other.paper_keywords)

        # Merge community data
        if other.psychonaut_url:
            self.psychonaut_url = other.psychonaut_url
            self.psychonaut_data.update(other.psychonaut_data)
        if other.erowid_url:
            self.erowid_url = other.erowid_url
            self.erowid_data.update(other.erowid_data)
        if other.tripsit_url:
            self.tripsit_url = other.tripsit_url
            self.tripsit_data.update(other.tripsit_data)

        # Merge social data
        self.reddit_mentions.extend(other.reddit_mentions)
        self.twitter_mentions.extend(other.twitter_mentions)
        self.bluesky_mentions.extend(other.bluesky_mentions)
        self.discord_mentions.extend(other.discord_mentions)

        # Update metadata
        self.enrichment_sources.update(other.enrichment_sources)
        for k, v in other.enrichment_stats.items():
            self.enrichment_stats[k] = self.enrichment_stats.get(k, 0) + v

        self.last_enriched = datetime.now().isoformat()
        self.update_history.append(
            {
                "timestamp": self.last_enriched,
                "source": source,
                "action": "web_data_merge",
            }
        )

    def to_dict(self) -> Dict:
        """Convert enriched compound data to dictionary format."""
        base_dict = super().to_dict()
        enrichment_dict = self.get_enrichment_dict()

        web_dict = {
            # Patent data
            "patent_numbers": list(self.patent_numbers),
            "patent_titles": self.patent_titles,
            "patent_abstracts": self.patent_abstracts,
            "patent_claims": self.patent_claims,
            "patent_citations": self.patent_citations,
            # Literature data
            "pubmed_ids": list(self.pubmed_ids),
            "paper_titles": self.paper_titles,
            "paper_abstracts": self.paper_abstracts,
            "paper_citations": self.paper_citations,
            "paper_keywords": self.paper_keywords,
            # Community data
            "psychonaut_url": self.psychonaut_url,
            "psychonaut_data": self.psychonaut_data,
            "erowid_url": self.erowid_url,
            "erowid_data": self.erowid_data,
            "tripsit_url": self.tripsit_url,
            "tripsit_data": self.tripsit_data,
            # Social data
            "reddit_mentions": self.reddit_mentions,
            "twitter_mentions": self.twitter_mentions,
            "bluesky_mentions": self.bluesky_mentions,
            "discord_mentions": self.discord_mentions,
            # Metadata
            "last_enriched": self.last_enriched,
            "enrichment_sources": list(self.enrichment_sources),
            "enrichment_stats": self.enrichment_stats,
            # Enrichment data
            "enrichment": enrichment_dict,
        }

        return {**base_dict, **web_dict}


class ValidationError(Exception):
    """Raised when web data validation fails."""

    pass
