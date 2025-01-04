"""Compound data model combining all psychopharmacological analysis capabilities."""

from dataclasses import dataclass, field
from typing import Dict, Optional, Set

from .base import BaseMixin, PsychoactiveClass, NootropicMechanism, RiskLevel
from .binding import ReceptorProfileMixin
from .activity import ActivityProfileMixin
from .safety import SafetyProfileMixin


@dataclass
class PsychopharmCompound(
    BaseMixin,
    ReceptorProfileMixin,
    ActivityProfileMixin,
    SafetyProfileMixin
):
    """Complete compound data model with psychopharmacological analysis."""

    # Core properties
    name: str = field(default="")
    smiles: Optional[str] = field(default=None)
    cas_number: Optional[str] = field(default=None)
    
    # Classification
    psychoactive_class: PsychoactiveClass = field(default=PsychoactiveClass.UNKNOWN)
    secondary_classes: Set[PsychoactiveClass] = field(default_factory=set)
    
    # Mechanisms
    nootropic_mechanisms: Set[NootropicMechanism] = field(default_factory=set)
    
    # Additional metadata
    source_urls: Set[str] = field(default_factory=set)
    literature_references: Set[str] = field(default_factory=set)
    community_reports: Set[str] = field(default_factory=set)

    def get_compound_dict(self) -> Dict:
        """Get complete dictionary of compound data."""
        return {
            # Base data
            **self.get_base_dict(),
            
            # Binding data
            "binding": self.get_binding_dict(),
            
            # Activity data
            "activity": self.get_activity_dict(),
            
            # Safety data
            "safety": self.get_safety_dict(),
            
            # Additional metadata
            "metadata": {
                "sources": list(self.source_urls),
                "references": list(self.literature_references),
                "reports": list(self.community_reports),
            },
        }

    def merge_compound_data(self, other: "PsychopharmCompound") -> None:
        """Merge all data from another compound instance."""
        # Merge base data
        self.merge_base_data(other)
        
        # Merge binding data
        self.merge_binding_data(other)
        
        # Merge activity data
        self.merge_activity_data(other)
        
        # Merge safety data
        self.merge_safety_data(other)
        
        # Merge additional metadata
        self.source_urls.update(other.source_urls)
        self.literature_references.update(other.literature_references)
        self.community_reports.update(other.community_reports)

    def _update_metadata(self, web_data: Dict) -> None:
        """Update metadata from web sources."""
        if "urls" in web_data:
            self.source_urls.update(web_data["urls"])
        if "references" in web_data:
            self.literature_references.update(web_data["references"])
        if "reports" in web_data:
            self.community_reports.update(web_data["reports"])

    def _update_classification(self, classification_data: Dict) -> None:
        """Update psychoactive classification from web data."""
        new_class = classification_data.get("primary_class")
        confidence = classification_data.get("confidence", 0.0)
        
        if not new_class:
            return
            
        if confidence > self.confidence_scores.get("classification", 0.0):
            # Store old class as secondary if it exists
            if self.psychoactive_class != PsychoactiveClass.UNKNOWN:
                self.secondary_classes.add(self.psychoactive_class)
                
            # Update to new class
            try:
                self.psychoactive_class = PsychoactiveClass[new_class.upper()]
                self.confidence_scores["classification"] = confidence
            except KeyError:
                return
                
            # Add any secondary classes
            for cls in classification_data.get("secondary_classes", []):
                try:
                    self.secondary_classes.add(PsychoactiveClass[cls.upper()])
                except KeyError:
                    continue

    def _update_mechanisms(self, mechanism_data: Dict) -> None:
        """Update nootropic mechanisms from web data."""
        for mech, conf in mechanism_data.items():
            try:
                mechanism = NootropicMechanism[mech.upper()]
                if conf > self.confidence_scores.get(f"mechanism_{mech}", 0.0):
                    self.nootropic_mechanisms.add(mechanism)
                    self.confidence_scores[f"mechanism_{mech}"] = conf
            except KeyError:
                continue

    def update_from_web_data(self, web_data: Dict) -> None:
        """Update compound data from web sources."""
        # Update metadata
        self._update_metadata(web_data)
        
        # Update classification
        if "classification" in web_data:
            self._update_classification(web_data["classification"])
            
        # Update mechanisms
        if "mechanisms" in web_data:
            self._update_mechanisms(web_data["mechanisms"])

    def validate(self) -> bool:
        """Validate compound data."""
        # Must have basic identification
        if not (self.name or self.smiles or self.cas_number):
            return False
            
        # Must have some binding data
        if not self.receptor_profiles:
            return False
            
        # Must have reasonable confidence scores
        if any(score > 1.0 for score in self.confidence_scores.values()):
            return False
            
        # Must have consistent risk levels
        safety_dict = self.get_safety_dict()
        if "risk_assessment" in safety_dict:
            risks = safety_dict["risk_assessment"]
            if "overall_risk" in risks:
                try:
                    RiskLevel[risks["overall_risk"].upper()]
                except KeyError:
                    return False
                    
        return True

    def __str__(self) -> str:
        """Get string representation."""
        return (
            f"PsychopharmCompound("
            f"name='{self.name}', "
            f"class={self.psychoactive_class.value}, "
            f"receptors={len(self.receptor_profiles)}, "
            f"mechanisms={len(self.nootropic_mechanisms)}"
            f")"
        )
