"""Base psychopharmacological compound data model.

This module provides the consolidated PsychopharmBase class that combines:
- Core compound data functionality (identifiers, properties, validation)
- Web enrichment capabilities (patents, literature, community data)
- ML prediction integration (features, predictions, uncertainty)
- Analysis features (binding, activity, safety)
- Enhanced validation and serialization

The base class serves as the foundation for all psychopharmacological compound
data models in the system, providing a comprehensive set of features while
maintaining tight integration with the rest of the codebase.
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Set, Any, Tuple
import json
import numpy as np
import pandas as pd
import re
from ..compound.types import (
    CompoundType,
    LegalStatus,
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
    BindingType,
    RiskLevel,
    TargetData,
    StringSet,
    StringDict,
    ValidationErrors,
    OptionalStr,
    TargetDict,
    DoseRange,
    TimeRange,
    RiskScore,
    EffectScore,
    ReceptorBinding,
)

class ValidationError(Exception):
    """Raised when compound data validation fails."""
    pass


@dataclass
class PsychopharmBase:
    """Base class for psychopharmacological compounds."""
    
    # Core identifiers (required)
    name: str
    smiles: str

    # Basic identifiers (optional)
    cas_number: OptionalStr = None
    inchi: OptionalStr = None
    inchi_key: OptionalStr = None
    iupac_name: OptionalStr = None
    compound_type: CompoundType = CompoundType.OTHER

    # Ranked common names with search results
    common_name_1: str = "N/A"
    common_name_2: str = "N/A"
    common_name_3: str = "N/A"
    common_name_1_results: int = 0
    common_name_2_results: int = 0
    common_name_3_results: int = 0
    other_names: StringSet = field(default_factory=set)

    # Chemical properties
    molecular_weight: float = 0.0
    logp: float = 0.0
    hbd: int = 0  # Hydrogen bond donors
    hba: int = 0  # Hydrogen bond acceptors
    tpsa: float = 0.0  # Topological polar surface area
    rotatable_bonds: int = 0
    charge: int = 0
    stereocenter_count: int = 0
    ring_count: int = 0

    # Database identifiers
    chembl_id: OptionalStr = None
    pubchem_cid: OptionalStr = None
    pubchem_sid: OptionalStr = None
    drugbank_id: OptionalStr = None
    bindingdb_id: OptionalStr = None

    # Target data
    targets: List[TargetData] = field(default_factory=list)
    target_data: TargetDict = field(default_factory=dict)
    primary_target: OptionalStr = None
    primary_activity: OptionalStr = None
    mechanism_of_action: OptionalStr = None
    pharmacology: str = "N/A"
    toxicity: str = "N/A"
    metabolism: str = "N/A"

    # Psychopharm-specific fields
    psychoactive_class: PsychoactiveClass = PsychoactiveClass.UNKNOWN
    nootropic_mechanisms: Set[NootropicMechanism] = field(default_factory=set)
    bbb_permeability: BBBPermeability = BBBPermeability.UNKNOWN
    risk_level: RiskLevel = RiskLevel.UNKNOWN

    # Web enrichment data
    patent_data: Dict = field(default_factory=dict)
    patent_count: int = 0
    community_data: Dict = field(default_factory=dict)
    literature_data: Dict = field(default_factory=dict)
    regulatory_data: Dict = field(default_factory=dict)
    experience_reports: List[Dict] = field(default_factory=list)
    safety_profile: Dict = field(default_factory=dict)
    dosage_info: Dict[str, Dict] = field(default_factory=dict)
    route_stats: Dict[str, int] = field(default_factory=dict)
    duration_stats: Dict[str, int] = field(default_factory=dict)

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
    psychonaut_data: Dict = field(default_factory=dict)
    erowid_url: OptionalStr = None
    erowid_data: Dict = field(default_factory=dict)
    tripsit_url: OptionalStr = None
    tripsit_data: Dict = field(default_factory=dict)

    # Social media data
    reddit_mentions: List[Dict] = field(default_factory=list)  # {subreddit, title, url,
                                                              # score, date}
    twitter_mentions: List[Dict] = field(default_factory=list)  # [{user, text, url, date}]
    bluesky_mentions: List[Dict] = field(default_factory=list)  # [{user, text, url, date}]
    discord_mentions: List[Dict] = field(default_factory=list)  # [{server, channel, text, date}]

    # ML prediction data
    _feature_cache: Dict[str, np.ndarray] = field(default_factory=dict)
    _feature_importances: Dict[str, Dict[str, float]] = field(default_factory=dict)
    _feature_scalers: Dict[str, object] = field(default_factory=dict)
    _prediction_cache: Dict[str, object] = field(default_factory=dict)
    _prediction_history: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(
            columns=[
                'predictor_type',
                'prediction_value',
                'confidence',
                'timestamp',
                'supporting_data',
            ]
        )
    )

    # Analysis data
    binding_profiles: Dict[str, ReceptorBinding] = field(default_factory=dict)
    binding_confidence: Dict[str, float] = field(default_factory=dict)
    binding_sources: Dict[str, List[str]] = field(default_factory=dict)
    activity_profiles: Dict[str, EffectScore] = field(default_factory=dict)
    activity_confidence: Dict[str, float] = field(default_factory=dict)
    activity_sources: Dict[str, List[str]] = field(default_factory=dict)
    risk_profiles: Dict[str, RiskScore] = field(default_factory=dict)
    risk_confidence: Dict[str, float] = field(default_factory=dict)
    risk_sources: Dict[str, List[str]] = field(default_factory=dict)
    dose_ranges: Dict[str, DoseRange] = field(default_factory=dict)
    time_ranges: Dict[str, TimeRange] = field(default_factory=dict)
    dose_confidence: Dict[str, float] = field(default_factory=dict)

    # Source information
    data_sources: StringSet = field(default_factory=set)
    reference_dois: StringSet = field(default_factory=set)
    reference_pmids: StringSet = field(default_factory=set)
    reference_urls: StringDict = field(default_factory=dict)

    # Legal & classification
    legal_status: Dict[str, LegalStatus] = field(default_factory=dict)  # Country -> Status
    scheduling: StringDict = field(default_factory=dict)  # Country -> Schedule

    # Metadata
    last_updated: str = field(default_factory=lambda: datetime.now().isoformat())
    version: str = "1.0.0"
    analysis_version: str = "1.0.0"
    analysis_timestamp: Optional[str] = None
    analysis_sources: Set[str] = field(default_factory=set)
    enrichment_sources: StringSet = field(default_factory=set)
    enrichment_stats: Dict[str, int] = field(default_factory=dict)  # source -> count
    last_enriched: str = field(default_factory=lambda: datetime.now().isoformat())

    def __post_init__(self):
        """Initialize and validate data."""
        self._validate()

    def _validate(self):
        """Validate compound data."""
        errors = []
        
        # Run all validation checks
        errors.extend(self._validate_identifiers())
        errors.extend(self._validate_properties())
        errors.extend(self._validate_targets())
        errors.extend(self._validate_psychopharm())
        errors.extend(self._validate_web_data())
        errors.extend(self._validate_analysis())
        
        if errors:
            raise ValidationError("\n".join(errors))
            
        # Initialize collections
        self._initialize_collections()

    def _validate_identifiers(self) -> ValidationErrors:
        """Validate chemical identifiers."""
        errors = []
        
        # Validate required fields
        if not self.name:
            errors.append("Compound name is required")
        if not self.smiles:
            errors.append("SMILES string is required")

        # Validate CAS number
        if self.cas_number and not self._validate_cas_format(self.cas_number):
            errors.append(f"Invalid CAS number format: {self.cas_number}")

        return errors

    def _validate_properties(self) -> ValidationErrors:
        """Validate chemical properties."""
        errors = []
        
        # Validate molecular weight
        if self.molecular_weight < 0:
            errors.append(f"Invalid molecular weight: {self.molecular_weight}")

        # Validate LogP
        if abs(self.logp) > 20:
            errors.append(f"Suspicious LogP value: {self.logp}")

        # Validate TPSA
        if self.tpsa < 0:
            errors.append(f"Invalid TPSA value: {self.tpsa}")

        return errors

    def _validate_targets(self) -> ValidationErrors:
        """Validate target data."""
        errors = []
        for i, target in enumerate(self.targets, 1):
            # Validate affinity value
            if target.affinity_value < 0:
                errors.append(
                    f"Invalid binding affinity value for target {i}: {target.affinity_value}"
                )
                
            # Validate confidence score
            if not 0 <= target.confidence <= 1:
                errors.append(
                    f"Invalid confidence score for target {i}: {target.confidence}"
                )
                
            # Validate affinity type
            valid_types = {'Ki', 'IC50', 'EC50', 'Kd'}
            if target.affinity_type not in valid_types and target.affinity_type != "N/A":
                errors.append(
                    f"Invalid affinity type for target {i}: {target.affinity_type}"
                )
                
            # Validate affinity unit
            valid_units = {'nM', 'uM', 'mM', 'pM'}
            if target.affinity_unit not in valid_units and target.affinity_unit != "N/A":
                errors.append(
                    f"Invalid affinity unit for target {i}: {target.affinity_unit}"
                )
        return errors

    def _validate_psychopharm(self) -> ValidationErrors:
        """Validate psychopharmacological properties."""
        errors = []
        
        # Validate binding type consistency
        if (
            any(bt != BindingType.UNKNOWN for bt in self.binding_types.values()) and
            not self.activity_types
        ):
            errors.append("Activity types must be specified when binding types are known")
            
        # Validate nootropic mechanism consistency
        if (
            self.nootropic_mechanisms and
            self.psychoactive_class != PsychoactiveClass.NOOTROPIC
        ):
            errors.append(
                "Nootropic mechanisms can only be specified for nootropic compounds"
            )
            
        return errors

    def _validate_web_data(self) -> ValidationErrors:
        """Validate web enrichment data."""
        errors = []
        
        # Validate patent data
        for number in self.patent_numbers:
            if not self._validate_patent_number(number):
                errors.append(f"Invalid patent number format: {number}")
                
        # Validate literature data
        for pmid in self.pubmed_ids:
            if not pmid.isdigit():
                errors.append(f"Invalid PubMed ID format: {pmid}")
                
        # Validate URLs
        if self.psychonaut_url and not self._validate_url(self.psychonaut_url):
            errors.append(f"Invalid PsychonautWiki URL: {self.psychonaut_url}")
        if self.erowid_url and not self._validate_url(self.erowid_url):
            errors.append(f"Invalid Erowid URL: {self.erowid_url}")
        if self.tripsit_url and not self._validate_url(self.tripsit_url):
            errors.append(f"Invalid TripSit URL: {self.tripsit_url}")
            
        return errors

    def _validate_analysis(self) -> ValidationErrors:
        """Validate analysis data."""
        errors = []
        errors.extend(self._validate_binding_profiles())
        errors.extend(self._validate_activity_profiles())
        errors.extend(self._validate_risk_profiles())
        errors.extend(self._validate_dosage_data())
        return errors

    def _validate_binding_profiles(self) -> ValidationErrors:
        """Validate binding profile data."""
        errors = []
        for target, (affinity, confidence, _) in self.binding_profiles.items():
            if affinity < 0:
                errors.append(f"Invalid binding affinity for {target}: {affinity}")
            if not 0 <= confidence <= 1:
                errors.append(f"Invalid confidence score for {target}: {confidence}")
        return errors

    def _validate_activity_profiles(self) -> ValidationErrors:
        """Validate activity profile data."""
        errors = []
        for effect, (magnitude, confidence) in self.activity_profiles.items():
            if not 0 <= magnitude <= 1:
                errors.append(f"Invalid effect magnitude for {effect}: {magnitude}")
            if not 0 <= confidence <= 1:
                errors.append(f"Invalid confidence score for {effect}: {confidence}")
        return errors

    def _validate_risk_profiles(self) -> ValidationErrors:
        """Validate risk profile data."""
        errors = []
        for risk, (severity, confidence) in self.risk_profiles.items():
            if not 0 <= severity <= 1:
                errors.append(f"Invalid risk severity for {risk}: {severity}")
            if not 0 <= confidence <= 1:
                errors.append(f"Invalid confidence score for {risk}: {confidence}")
        return errors

    def _validate_dosage_data(self) -> ValidationErrors:
        """Validate dosage data."""
        errors = []
        for route, (min_dose, max_dose, recommended) in self.dose_ranges.items():
            if not (min_dose <= recommended <= max_dose):
                msg = f"Invalid dose range for {route}: {min_dose}, {recommended}, {max_dose}"
                errors.append(msg)
        return errors

    def _initialize_collections(self) -> None:
        """Initialize collection fields."""
        for field_name, field_type in self.__annotations__.items():
            if hasattr(self, field_name):
                field_value = getattr(self, field_name)
                if field_value is None:
                    if "List" in str(field_type):
                        setattr(self, field_name, [])
                    elif "Dict" in str(field_type):
                        setattr(self, field_name, {})
                    elif "Set" in str(field_type):
                        setattr(self, field_name, set())

    def _validate_cas_format(self, cas: str) -> bool:
        """
        Validate CAS number format.
        
        Args:
            cas: CAS number to validate
            
        Returns:
            True if valid, False otherwise
        """
        pattern = r'^\d{1,7}-\d{2}-\d$'
        if not re.match(pattern, cas):
            return False

        # Validate checksum
        numbers = cas.replace('-', '')
        check_digit = int(numbers[-1])
        numbers = numbers[:-1]
        total = sum(
            int(num) * (i + 1)
            for i, num in enumerate(reversed(numbers))
        )
        return (total % 10) == check_digit

    def _validate_patent_number(self, number: str) -> bool:
        """Validate patent number format."""
        # Basic validation - can be enhanced
        return bool(number and len(number) >= 6)

    def _validate_url(self, url: str) -> bool:
        """Validate URL format."""
        # Basic validation - can be enhanced
        return url.startswith(("http://", "https://"))

    def format_numeric_values(self):
        """Format numeric values to specified precision."""
        self.logp = float(f"{self.logp:.5f}".rstrip('0').rstrip('.'))
        self.tpsa = float(f"{self.tpsa:.5f}".rstrip('0').rstrip('.'))
        self.molecular_weight = float(f"{self.molecular_weight:.5f}".rstrip('0').rstrip('.'))

    # Web Enrichment Methods
    def merge_web_data(self, other_data: Dict) -> None:
        """Merge web-enriched data."""
        if "patent_data" in other_data:
            self.patent_data.update(other_data["patent_data"])
            self.patent_count = max(self.patent_count, other_data.get("patent_count", 0))

        if "community_data" in other_data:
            self.community_data.update(other_data["community_data"])

        if "literature_data" in other_data:
            self.literature_data.update(other_data["literature_data"])

        if "regulatory_data" in other_data:
            self.regulatory_data.update(other_data["regulatory_data"])

        if "experience_reports" in other_data:
            self.experience_reports.extend(
                report for report in other_data["experience_reports"]
                if report not in self.experience_reports
            )

        if "safety_profile" in other_data:
            self.safety_profile.update(other_data["safety_profile"])

        if "dosage_info" in other_data:
            for route, stats in other_data["dosage_info"].items():
                if route not in self.dosage_info:
                    self.dosage_info[route] = stats
                else:
                    current = self.dosage_info[route]
                    current["min"] = min(current["min"], stats["min"])
                    current["max"] = max(current["max"], stats["max"])
                    total = current["avg"] * current["count"] + stats["avg"] * stats["count"]
                    count = current["count"] + stats["count"]
                    current["avg"] = total / count
                    current["count"] += stats["count"]

    def merge_enrichment_data(self, other: "PsychopharmBase") -> None:
        """Merge enrichment data from another instance."""
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

    # ML Prediction Methods
    def get_cached_prediction(
        self,
        predictor_type: str,
    ) -> Optional[object]:
        """Get cached prediction result if available."""
        return self._prediction_cache.get(predictor_type)

    def cache_prediction(
        self,
        predictor_type: str,
        result: object,
        confidence: float,
        supporting_data: Dict = None,
    ) -> None:
        """Cache prediction result."""
        self._prediction_cache[predictor_type] = result
        
        # Update history
        self._prediction_history = pd.concat([
            self._prediction_history,
            pd.DataFrame([{
                'predictor_type': predictor_type,
                'prediction_value': result,
                'confidence': confidence,
                'timestamp': pd.Timestamp.now(),
                'supporting_data': json.dumps(supporting_data or {}),
            }])
        ], ignore_index=True)

    def clear_prediction_cache(self) -> None:
        """Clear cached predictions."""
        self._prediction_cache.clear()

    def predict_binding(
        self,
        receptor: str,
        include_uncertainty: bool = True
    ) -> Tuple[float, float]:
        """
        Predict binding affinity with uncertainty.
        
        Args:
            receptor: Target receptor name
            include_uncertainty: Whether to include uncertainty estimate
            
        Returns:
            Tuple of (prediction, uncertainty)
        """
        cached = self.get_cached_prediction(f"binding_{receptor}")
        if cached is not None:
            return cached

        prediction = self._base_binding_prediction(receptor)
        uncertainty = self._estimate_binding_uncertainty(receptor) if include_uncertainty else 0.0
        
        self.cache_prediction(
            f"binding_{receptor}",
            (prediction, uncertainty),
            1.0 - uncertainty
        )
        
        return prediction, uncertainty

    def predict_ensemble(
        self,
        receptor: str,
        n_models: int = 5
    ) -> Dict[str, Any]:
        """
        Get ensemble predictions.
        
        Args:
            receptor: Target receptor name
            n_models: Number of models in ensemble
            
        Returns:
            Dictionary with prediction statistics
        """
        cached = self.get_cached_prediction(f"ensemble_{receptor}")
        if cached is not None:
            return cached

        predictions = []
        for _ in range(n_models):
            pred = self._single_model_prediction(receptor)
            predictions.append(pred)
            
        result = {
            "mean": float(np.mean(predictions)),
            "std": float(np.std(predictions)),
            "individual": predictions
        }
        
        self.cache_prediction(
            f"ensemble_{receptor}",
            result,
            1.0 / result["std"] if result["std"] > 0 else 1.0
        )
        
        return result

    # Analysis Methods
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
                    "range": dose_range,
                    "confidence": self.dose_confidence.get(route, 0.0),
                    "timing": self.time_ranges.get(route)
                }
                for route, dose_range in self.dose_ranges.items()
            },
            "metadata": {
                "version": self.analysis_version,
                "timestamp": self.analysis_timestamp,
                "sources": list(self.analysis_sources)
            }
        }
