"""Analysis functionality for psychopharmacological compounds.

This module provides analysis capabilities that integrate with web enrichment:
1. Binding analysis using both experimental and predicted data
2. Activity analysis incorporating community reports
3. Safety analysis combining literature and community data
4. Property analysis with ML-enhanced predictions
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from .base import CompoundData, RiskLevel
from .types import (
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
    BindingType,
    ActivityType,
    TargetData,
    DoseRange,
    TimeRange,
    RiskScore,
    EffectScore,
    ReceptorBinding,
)


@dataclass
class AnalyzedCompoundData(CompoundData):
    """Compound data with comprehensive analysis results."""

    # Binding analysis
    binding_profiles: Dict[str, ReceptorBinding] = field(default_factory=dict)
    binding_confidence: Dict[str, float] = field(default_factory=dict)
    binding_sources: Dict[str, List[str]] = field(default_factory=dict)

    # Activity analysis
    activity_profiles: Dict[str, EffectScore] = field(default_factory=dict)
    activity_confidence: Dict[str, float] = field(default_factory=dict)
    activity_sources: Dict[str, List[str]] = field(default_factory=dict)

    # Safety analysis
    risk_profiles: Dict[str, RiskScore] = field(default_factory=dict)
    risk_confidence: Dict[str, float] = field(default_factory=dict)
    risk_sources: Dict[str, List[str]] = field(default_factory=dict)

    # Dosage analysis
    dose_ranges: Dict[str, DoseRange] = field(default_factory=dict)
    time_ranges: Dict[str, TimeRange] = field(default_factory=dict)
    dose_confidence: Dict[str, float] = field(default_factory=dict)

    # Analysis metadata
    analysis_version: str = "1.0.0"
    analysis_timestamp: Optional[str] = None
    analysis_sources: Set[str] = field(default_factory=set)

    def __post_init__(self):
        """Initialize and validate analysis data."""
        super().__post_init__()
        self._validate_analysis()

    def _validate_analysis(self) -> None:
        """Validate analysis data."""
        errors = []
        
        # Run analysis validation checks
        errors.extend(self._validate_binding_profiles())
        errors.extend(self._validate_activity_profiles())
        errors.extend(self._validate_risk_profiles())
        errors.extend(self._validate_dosage_data())
        
        if errors:
            raise ValidationError("\n".join(errors))

    def _validate_binding_profiles(self) -> List[str]:
        """Validate binding profile data."""
        errors = []
        
        for target, (affinity, confidence, activity) in self.binding_profiles.items():
            # Validate affinity value
            if affinity < 0:
                errors.append(f"Invalid binding affinity for {target}: {affinity}")
                
            # Validate confidence score
            if not 0 <= confidence <= 1:
                errors.append(f"Invalid confidence score for {target}: {confidence}")
                
            # Validate activity type
            try:
                BindingType(activity)
            except ValueError:
                errors.append(f"Invalid binding type for {target}: {activity}")
                
        return errors

    def _validate_activity_profiles(self) -> List[str]:
        """Validate activity profile data."""
        errors = []
        
        for effect, (magnitude, confidence) in self.activity_profiles.items():
            # Validate magnitude
            if not 0 <= magnitude <= 1:
                errors.append(f"Invalid effect magnitude for {effect}: {magnitude}")
                
            # Validate confidence
            if not 0 <= confidence <= 1:
                errors.append(f"Invalid confidence score for {effect}: {confidence}")
                
        return errors

    def _validate_risk_profiles(self) -> List[str]:
        """Validate risk profile data."""
        errors = []
        
        for risk, (severity, confidence) in self.risk_profiles.items():
            # Validate severity
            if not 0 <= severity <= 1:
                errors.append(f"Invalid risk severity for {risk}: {severity}")
                
            # Validate confidence
            if not 0 <= confidence <= 1:
                errors.append(f"Invalid confidence score for {risk}: {confidence}")
                
        return errors

    def _validate_dosage_data(self) -> List[str]:
        """Validate dosage data."""
        errors = []
        
        for route, (min_dose, max_dose, recommended) in self.dose_ranges.items():
            # Validate dose order
            if not (min_dose <= recommended <= max_dose):
                errors.append(
                    f"Invalid dose range for {route}: "
                    f"{min_dose}, {recommended}, {max_dose}"
                )
                
            # Validate confidence
            if route in self.dose_confidence:
                confidence = self.dose_confidence[route]
                if not 0 <= confidence <= 1:
                    errors.append(
                        f"Invalid confidence score for {route}: {confidence}"
                    )
                    
        return errors

    def analyze_binding(self, target_data: List[TargetData]) -> None:
        """Analyze binding data for targets."""
        for target in target_data:
            # Extract binding data
            affinity = target.affinity_value
            confidence = target.confidence
            activity = target.activity_type
            
            # Store binding profile
            self.binding_profiles[target.common_name] = (
                affinity, confidence, activity
            )
            self.binding_confidence[target.common_name] = confidence
            
            # Track sources
            if target.common_name not in self.binding_sources:
                self.binding_sources[target.common_name] = []
            self.binding_sources[target.common_name].extend(
                [doi for doi in target.reference_dois]
            )

    def _process_community_activity(self, source: str, data: Dict) -> None:
        """Process activity data from a community source."""
        if "effects" in data:
            for effect, details in data["effects"].items():
                magnitude = details.get("magnitude", 0.0)
                confidence = details.get("confidence", 0.5)
                
                self._update_activity_profile(effect, magnitude, confidence, source)

    def _process_literature_activity(self, citation: Dict) -> None:
        """Process activity data from a literature source."""
        if "activity" in citation:
            for effect, details in citation["activity"].items():
                magnitude = details.get("magnitude", 0.0)
                confidence = details.get("confidence", 0.7)
                
                self._update_activity_profile(
                    effect, magnitude, confidence, citation["pmid"]
                )

    def _update_activity_profile(
        self,
        effect: str,
        magnitude: float,
        confidence: float,
        source: str
    ) -> None:
        """Update activity profile with new data."""
        if effect not in self.activity_profiles:
            self.activity_profiles[effect] = (magnitude, confidence)
            self.activity_confidence[effect] = confidence
            self.activity_sources[effect] = [source]
        elif confidence > self.activity_confidence[effect]:
            self.activity_profiles[effect] = (magnitude, confidence)
            self.activity_confidence[effect] = confidence
            self.activity_sources[effect] = [source]

    def analyze_activity(self, web_data: Dict) -> None:
        """Analyze activity data from web sources."""
        # Process community data
        if "community" in web_data:
            for source, data in web_data["community"].items():
                self._process_community_activity(source, data)
                        
        # Process literature data
        if "literature" in web_data:
            for citation in web_data["literature"].get("pubmed", []):
                self._process_literature_activity(citation)

    def _process_community_safety(self, source: str, data: Dict) -> None:
        """Process safety data from a community source."""
        if "risks" in data:
            for risk, details in data["risks"].items():
                severity = details.get("severity", 0.0)
                confidence = details.get("confidence", 0.5)
                
                self._update_risk_profile(risk, severity, confidence, source)

    def _process_literature_safety(self, citation: Dict) -> None:
        """Process safety data from a literature source."""
        if "safety" in citation:
            for risk, details in citation["safety"].items():
                severity = details.get("severity", 0.0)
                confidence = details.get("confidence", 0.7)
                
                self._update_risk_profile(
                    risk, severity, confidence, citation["pmid"]
                )

    def _update_risk_profile(
        self,
        risk: str,
        severity: float,
        confidence: float,
        source: str
    ) -> None:
        """Update risk profile with new data."""
        if risk not in self.risk_profiles:
            self.risk_profiles[risk] = (severity, confidence)
            self.risk_confidence[risk] = confidence
            self.risk_sources[risk] = [source]
        elif confidence > self.risk_confidence[risk]:
            self.risk_profiles[risk] = (severity, confidence)
            self.risk_confidence[risk] = confidence
            self.risk_sources[risk] = [source]

    def analyze_safety(self, web_data: Dict) -> None:
        """Analyze safety data from web sources."""
        # Process community data
        if "community" in web_data:
            for source, data in web_data["community"].items():
                self._process_community_safety(source, data)
                            
        # Process literature data
        if "literature" in web_data:
            for citation in web_data["literature"].get("pubmed", []):
                self._process_literature_safety(citation)

    def analyze_dosage(self, web_data: Dict) -> None:
        """Analyze dosage data from web sources."""
        # Process community data
        if "community" in web_data:
            for source, data in web_data["community"].items():
                if "dosage" in data:
                    for route, details in data["dosage"].items():
                        min_dose = details.get("min", 0.0)
                        max_dose = details.get("max", 0.0)
                        recommended = details.get("recommended", 0.0)
                        confidence = details.get("confidence", 0.5)
                        
                        # Validate and store dose range
                        if min_dose <= recommended <= max_dose:
                            self.dose_ranges[route] = (
                                min_dose, max_dose, recommended
                            )
                            self.dose_confidence[route] = confidence
                            
                        # Store timing data if available
                        if "onset" in details and "duration" in details:
                            self.time_ranges[route] = (
                                details["onset"],
                                details["duration"]
                            )

    def get_analysis_summary(self) -> Dict:
        """Get summary of analysis results."""
        return {
            "binding": {
                target: {
                    "affinity": affinity,
                    "confidence": confidence,
                    "activity": activity,
                    "sources": self.binding_sources.get(target, [])
                }
                for target, (affinity, confidence, activity)
                in self.binding_profiles.items()
            },
            "activity": {
                effect: {
                    "magnitude": magnitude,
                    "confidence": confidence,
                    "sources": self.activity_sources.get(effect, [])
                }
                for effect, (magnitude, confidence)
                in self.activity_profiles.items()
            },
            "safety": {
                risk: {
                    "severity": severity,
                    "confidence": confidence,
                    "sources": self.risk_sources.get(risk, [])
                }
                for risk, (severity, confidence)
                in self.risk_profiles.items()
            },
            "dosage": {
                route: {
                    "range": self.dose_ranges[route],
                    "confidence": self.dose_confidence.get(route, 0.0),
                    "timing": self.time_ranges.get(route, None)
                }
                for route in self.dose_ranges
            },
            "metadata": {
                "version": self.analysis_version,
                "timestamp": self.analysis_timestamp,
                "sources": list(self.analysis_sources)
            }
        }


class ValidationError(Exception):
    """Raised when analysis data validation fails."""
    pass
