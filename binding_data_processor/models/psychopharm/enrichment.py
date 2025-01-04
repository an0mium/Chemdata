"""Web data enrichment functionality."""

from dataclasses import field
from typing import Dict, List, Optional
from datetime import datetime
import json

from .base import RiskLevel


class WebEnrichmentMixin:
    """Mixin class providing web data enrichment capabilities."""

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
                    "timestamp": self.source_timestamps[source].isoformat()
                    if source in self.source_timestamps else None,
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

    def _update_source_data(
        self,
        source: str,
        data: Dict,
        timestamp: Optional[datetime] = None
    ) -> None:
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

    def _extract_risks_from_source(
        self,
        source_data: Dict,
        existing_risks: Dict[str, str]
    ) -> None:
        """Extract risk information from a data source."""
        if "risks" in source_data:
            for risk, level in source_data["risks"].items():
                try:
                    risk_level = RiskLevel[level.upper()]
                    if (
                        risk not in existing_risks or
                        risk_level.value > RiskLevel[
                            existing_risks[risk].upper()
                        ].value
                    ):
                        existing_risks[risk] = level
                except KeyError:
                    continue

    def _extract_safety_from_community(
        self,
        web_data: Dict,
        safety_data: Dict
    ) -> None:
        """Extract safety information from community sources."""
        for source in ["psychonaut", "erowid", "tripsit"]:
            if source in web_data.get("community", {}):
                data = web_data["community"][source]
                
                # Extract risks
                self._extract_risks_from_source(data, safety_data["risks"])
                
                # Extract contraindications
                if "contraindications" in data:
                    safety_data["contraindications"].update(
                        data["contraindications"]
                    )
                    
                # Extract interactions
                if "interactions" in data:
                    for substance, interaction in data["interactions"].items():
                        if substance not in safety_data["interactions"]:
                            safety_data["interactions"][substance] = interaction

    def _extract_safety_from_literature(
        self,
        web_data: Dict,
        safety_data: Dict
    ) -> None:
        """Extract safety information from literature sources."""
        if "pubmed" in web_data.get("literature", {}):
            for citation in web_data["literature"]["pubmed"]:
                if "safety" in citation:
                    safety = citation["safety"]
                    
                    # Extract risks
                    self._extract_risks_from_source(safety, safety_data["risks"])
                    
                    # Extract contraindications
                    if "contraindications" in safety:
                        safety_data["contraindications"].update(
                            safety["contraindications"]
                        )

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

    def update_from_web_sources(
        self,
        web_data: Dict,
        source: str,
        timestamp: Optional[datetime] = None
    ) -> None:
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
                self.source_timestamps[source] = other.source_timestamps.get(
                    source, datetime.now()
                )
            elif (
                source in other.source_timestamps and
                source in self.source_timestamps and
                other.source_timestamps[source] > self.source_timestamps[source]
            ):
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
                data["sources"][source]["timestamp"] = (
                    data["sources"][source]["timestamp"]
                )
                
        return json.dumps(data, indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> "WebEnrichmentMixin":
        """Create instance from JSON string."""
        data = json.loads(json_str)
        instance = cls()
        
        # Convert ISO format strings back to datetime objects
        for source, source_data in data["sources"].items():
            if "timestamp" in source_data:
                instance.source_timestamps[source] = datetime.fromisoformat(
                    source_data["timestamp"]
                )
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
