"""Enhanced compound data model.

This module provides the CompoundData class for representing chemical compounds
and their associated data from various sources including:

Core Data:
- Basic identifiers (name, SMILES, InChI, CAS)
- Chemical properties (MW, LogP, etc.)
- Database IDs (ChEMBL, PubChem, DrugBank)
- Ranked common names with search result counts

Binding & Activity:
- BindingDB data
- Target predictions
- Activity classifications
- SAR analysis
- PubMed result counts

Predictions:
- ML-based binding predictions
- Toxicity predictions
- Abuse potential assessment
- Activity predictions
- Mechanism predictions

Psychopharmacology:
- Receptor binding profiles
- Blood-brain barrier properties
- Psychoactive classification
- Nootropic activity
- Duration metrics
- Tolerance data

Web-Enriched:
- Patent information
- Literature data
- Community reports
- Regulatory status
- Safety profiles
- Reference URLs

The model supports both flexible data structures for dynamic content and
structured fields for core data to ensure consistency and enable efficient
filtering and analysis.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Dict, List, Optional, Set
import json
import numpy as np
import pandas as pd
import re

from ..ml.predictors import PredictorBase, PredictionResult
from .validation import ValidationError, validate_smiles, validate_inchi


class CompoundType(Enum):
    """Types of chemical compounds."""

    NEUROTRANSMITTER = "neurotransmitter"
    PSYCHOACTIVE = "psychoactive"
    RESEARCH_CHEMICAL = "research_chemical"
    NPS = "novel_psychoactive_substance"
    PHARMACEUTICAL = "pharmaceutical"
    NATURAL_PRODUCT = "natural_product"
    SMALL_MOLECULE = "small_molecule"
    OTHER = "other"
    UNKNOWN = "unknown"


class LegalStatus(Enum):
    """Legal status classifications."""

    LEGAL = "legal"
    CONTROLLED = "controlled"
    ILLEGAL = "illegal"
    RESEARCH_ONLY = "research_only"
    UNSCHEDULED = "unscheduled"
    PRESCRIPTION = "prescription_only"
    OTC = "over_the_counter"
    APPROVED = "approved"
    INVESTIGATIONAL = "investigational"
    WITHDRAWN = "withdrawn"
    BANNED = "banned"
    UNKNOWN = "unknown"


class PsychoactiveClass(Enum):
    """Classification of psychoactive effects."""

    PSYCHEDELIC = "psychedelic"
    EMPATHOGEN = "empathogen"
    STIMULANT = "stimulant"
    DEPRESSANT = "depressant"
    DISSOCIATIVE = "dissociative"
    DELIRIANT = "deliriant"
    NOOTROPIC = "nootropic"
    ANXIOLYTIC = "anxiolytic"
    ANTIPSYCHOTIC = "antipsychotic"
    ANTIDEPRESSANT = "antidepressant"
    MOOD_STABILIZER = "mood_stabilizer"
    UNKNOWN = "unknown"


class NootropicMechanism(Enum):
    """Mechanisms of nootropic activity."""

    CHOLINERGIC = "cholinergic"
    GLUTAMATERGIC = "glutamatergic"
    DOPAMINERGIC = "dopaminergic"
    SEROTONERGIC = "serotonergic"
    GABA = "gaba_modulation"
    AMPAKINE = "ampakine"
    BDNF = "bdnf_modulation"
    NGF = "ngf_modulation"
    NEUROPLASTICITY = "neuroplasticity"
    ANTI_INFLAMMATORY = "anti_inflammatory"
    ANTIOXIDANT = "antioxidant"
    MEMORY_ENHANCEMENT = "memory_enhancement"
    FOCUS_IMPROVEMENT = "focus_improvement"
    NEUROPROTECTION = "neuroprotection"
    NOOTROPIC_SYNERGY = "nootropic_synergy"
    COGNITIVE_MODULATION = "cognitive_modulation"
    BRAIN_METABOLISM = "brain_metabolism"
    UNKNOWN = "unknown"


class BBBPermeability(Enum):
    """Blood-brain barrier permeability classification."""

    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    NEGLIGIBLE = "negligible"
    UNKNOWN = "unknown"


@dataclass
class TargetData:
    """Structured data for a single target interaction."""

    common_name: str = "N/A"
    protein_name: str = "N/A"
    gene_name: str = "N/A"
    organism: str = "human"
    affinity_value: float = 0.0
    affinity_type: str = "N/A"  # Ki, IC50, Kd, EC50
    affinity_unit: str = "nM"
    activity_type: str = "N/A"  # agonist, antagonist, etc.
    confidence: float = 0.0
    is_primary: bool = False
    pubmed_count: int = 0
    reference_dois: Set[str] = field(default_factory=set)
    assay_details: Dict = field(default_factory=dict)
    experimental_conditions: Dict = field(default_factory=dict)


@dataclass
class CompoundData:
    """Enhanced data class for chemical compound information."""

    # Core identifiers (required)
    name: str
    smiles: str

    # Basic identifiers (optional)
    cas_number: Optional[str] = None
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    iupac_name: Optional[str] = None
    compound_type: CompoundType = CompoundType.OTHER

    # Ranked common names with search results
    common_name_1: str = "N/A"
    common_name_2: str = "N/A"
    common_name_3: str = "N/A"
    common_name_1_results: int = 0
    common_name_2_results: int = 0
    common_name_3_results: int = 0
    other_names: Set[str] = field(default_factory=set)

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
    chembl_id: Optional[str] = None
    pubchem_cid: Optional[str] = None
    pubchem_sid: Optional[str] = None
    drugbank_id: Optional[str] = None
    bindingdb_id: Optional[str] = None

    # Binding and activity data
    targets: List[TargetData] = field(default_factory=list)
    binding_data: List[Dict] = field(default_factory=list)
    primary_target: Optional[str] = None
    primary_activity: Optional[str] = None
    mechanism_of_action: Optional[str] = None
    pharmacology: str = "N/A"
    toxicity: str = "N/A"
    metabolism: str = "N/A"

    # Source information
    data_sources: Set[str] = field(default_factory=set)
    reference_dois: Set[str] = field(default_factory=set)
    reference_pmids: Set[str] = field(default_factory=set)
    reference_urls: Dict[str, str] = field(default_factory=dict)

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

    # Patent data
    patent_data: Dict = field(default_factory=dict)
    patent_count: int = 0

    # Swiss tools data
    swiss_data: Dict = field(default_factory=dict)
    target_predictions: List[Dict] = field(default_factory=list)
    adme_properties: Dict = field(default_factory=dict)
    similar_compounds: List[Dict] = field(default_factory=list)

    # ML predictions
    toxicity_predictions: Dict[str, Dict] = field(default_factory=dict)
    abuse_potential: Dict[str, Dict] = field(default_factory=dict)
    binding_predictions: List[Dict] = field(default_factory=list)
    activity_predictions: List[Dict] = field(default_factory=list)

    # Psychopharmacological properties
    receptor_profiles: Dict[str, Dict] = field(default_factory=dict)
    bbb_permeability: BBBPermeability = BBBPermeability.UNKNOWN
    bbb_score: float = 0.0
    p_glycoprotein_substrate: bool = False
    psychoactive_class: PsychoactiveClass = PsychoactiveClass.UNKNOWN
    secondary_classes: Set[PsychoactiveClass] = field(default_factory=set)
    effect_profile: Dict[str, float] = field(default_factory=dict)
    nootropic_mechanisms: Set[NootropicMechanism] = field(default_factory=set)
    cognitive_effects: Dict[str, Dict] = field(default_factory=dict)
    side_effects: Dict[str, Dict] = field(default_factory=dict)
    onset_time: Optional[str] = None
    duration: Optional[str] = None
    half_life: Optional[str] = None
    tolerance_profile: Dict[str, str] = field(default_factory=dict)
    withdrawal_profile: Dict[str, str] = field(default_factory=dict)
    cross_tolerance: Set[str] = field(default_factory=set)

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

    # Legal & classification
    legal_status: Dict[str, LegalStatus] = field(default_factory=dict)  # Country -> Status
    scheduling: Dict[str, str] = field(default_factory=dict)  # Country -> Schedule

    # Analysis results
    pharmacophores: List[Dict] = field(default_factory=list)
    structural_alerts: List[Dict] = field(default_factory=list)
    receptor_interactions: Dict[str, List[Dict]] = field(default_factory=dict)
    mechanism_predictions: List[Dict] = field(default_factory=list)
    sar_analysis: Dict = field(default_factory=dict)

    # Metadata
    last_updated: str = field(default_factory=lambda: datetime.now().isoformat())
    version: str = "1.0.0"

    # ML Prediction Integration
    _predictors: Dict[str, PredictorBase] = field(default_factory=dict)
    _prediction_cache: Dict[str, PredictionResult] = field(default_factory=dict)
    _feature_cache: Dict[str, np.ndarray] = field(default_factory=dict)
    _prediction_history: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(
            columns=[
                "predictor_type",
                "prediction_value",
                "confidence",
                "timestamp",
                "supporting_data",
            ]
        )
    )

    def __post_init__(self):
        """Validate required fields and initialize collections."""
        self._validate()

    def predict(
        self,
        predictor_type: str,
        predictor: PredictorBase,
        use_cache: bool = True,
    ) -> PredictionResult:
        """Run prediction using specified predictor.

        Args:
            predictor_type: Type of predictor to use
            predictor: Predictor instance to use
            use_cache: Whether to use cached predictions

        Returns:
            PredictionResult containing prediction and confidence
        """
        # Check cache
        if use_cache and predictor_type in self._prediction_cache:
            return self._prediction_cache[predictor_type]

        # Run prediction
        result = predictor.predict(self)

        # Cache result
        self._prediction_cache[predictor_type] = result
        self._predictors[predictor_type] = predictor

        # Update history
        self._prediction_history = pd.concat(
            [
                self._prediction_history,
                pd.DataFrame(
                    [
                        {
                            "predictor_type": predictor_type,
                            "prediction_value": result.value,
                            "confidence": result.confidence,
                            "timestamp": pd.Timestamp.now(),
                            "supporting_data": json.dumps(result.supporting_data),
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        return result

    def get_cached_prediction(
        self,
        predictor_type: str,
    ) -> Optional[PredictionResult]:
        """Get cached prediction result if available."""
        return self._prediction_cache.get(predictor_type)

    def clear_prediction_cache(self) -> None:
        """Clear cached predictions."""
        self._prediction_cache.clear()

    def get_prediction_history(
        self,
        predictor_type: Optional[str] = None,
    ) -> pd.DataFrame:
        """Get prediction history, optionally filtered by type."""
        if predictor_type:
            return self._prediction_history[self._prediction_history["predictor_type"] == predictor_type]
        return self._prediction_history

    def get_cached_features(
        self,
        feature_type: str,
    ) -> Optional[np.ndarray]:
        """Get cached features if available."""
        return self._feature_cache.get(feature_type)

    def cache_features(
        self,
        feature_type: str,
        features: np.ndarray,
    ) -> None:
        """Cache features for reuse."""
        self._feature_cache[feature_type] = features

    def clear_feature_cache(self) -> None:
        """Clear cached features."""
        self._feature_cache.clear()

    def _validate(self):
        """Validate compound data."""
        errors = []

        # Run all validation checks
        errors.extend(self._validate_identifiers())
        errors.extend(self._validate_properties())
        errors.extend(self._validate_targets())
        errors.extend(self._validate_psychopharm_data())

        if errors:
            raise ValidationError("\n".join(errors))

        # Initialize collections
        self._initialize_collections()

    def _validate_identifiers(self) -> List[str]:
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

        # Validate SMILES
        if not validate_smiles(self.smiles):
            errors.append(f"Invalid SMILES format: {self.smiles}")

        # Validate InChI
        if self.inchi and not validate_inchi(self.inchi):
            errors.append(f"Invalid InChI format: {self.inchi}")

        return errors

    def _validate_properties(self) -> List[str]:
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

    def _validate_targets(self) -> List[str]:
        """Validate target data."""
        errors = []
        for i, target in enumerate(self.targets, 1):
            errors.extend(self._validate_single_target(i, target))
        return errors

    def _validate_single_target(self, index: int, target) -> List[str]:
        """Validate a single target entry."""
        errors = []

        # Validate affinity value
        if target.affinity_value < 0:
            errors.append(f"Invalid binding affinity value for target {index}: {target.affinity_value}")

        # Validate confidence score
        if not 0 <= target.confidence <= 1:
            errors.append(f"Invalid confidence score for target {index}: {target.confidence}")

        # Validate affinity type
        valid_types = {"Ki", "IC50", "EC50", "Kd"}
        if target.affinity_type not in valid_types and target.affinity_type != "N/A":
            errors.append(f"Invalid affinity type for target {index}: {target.affinity_type}")

        # Validate affinity unit
        valid_units = {"nM", "uM", "mM", "pM"}
        if target.affinity_unit not in valid_units and target.affinity_unit != "N/A":
            errors.append(f"Invalid affinity unit for target {index}: {target.affinity_unit}")

        return errors

    def _validate_psychopharm_data(self) -> List[str]:
        """Validate psychopharmacological data."""
        errors = []

        # Validate BBB score
        if not 0 <= self.bbb_score <= 1:
            errors.append(f"Invalid BBB score: {self.bbb_score}")

        # Validate effect profile scores
        for effect, score in self.effect_profile.items():
            if not 0 <= score <= 1:
                errors.append(f"Invalid effect score for {effect}: {score}")

        return errors

    def _validate_required_fields(self) -> List[str]:
        """Validate required fields."""
        errors = []
        if not self.name:
            errors.append("Compound name is required")
        if not self.smiles:
            errors.append("SMILES string is required")
        return errors

    def _validate_cas_number(self) -> List[str]:
        """Validate CAS number format."""
        errors = []
        if self.cas_number and not self._validate_cas_format(self.cas_number):
            errors.append(f"Invalid CAS number format: {self.cas_number}")
        return errors

    def _validate_numeric_values(self) -> List[str]:
        """Validate numeric values."""
        errors = []
        if self.molecular_weight < 0:
            errors.append(f"Invalid molecular weight: {self.molecular_weight}")
        return errors

    def _validate_binding_data(self) -> List[str]:
        """Validate binding data."""
        errors = []
        for i, target in enumerate(self.targets, 1):
            if target.affinity_value < 0:
                errors.append(f"Invalid binding affinity value for target {i}: {target.affinity_value}")
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
        pattern = r"^\d{1,7}-\d{2}-\d$"
        if not re.match(pattern, cas):
            return False

        # Validate checksum
        numbers = cas.replace("-", "")
        check_digit = int(numbers[-1])
        numbers = numbers[:-1]
        total = sum(int(num) * (i + 1) for i, num in enumerate(reversed(numbers)))
        return (total % 10) == check_digit

    def format_numeric_values(self):
        """Format numeric values to specified precision."""
        self.logp = float(f"{self.logp:.5f}".rstrip("0").rstrip("."))
        self.tpsa = float(f"{self.tpsa:.5f}".rstrip("0").rstrip("."))
        self.molecular_weight = float(f"{self.molecular_weight:.5f}".rstrip("0").rstrip("."))

    def get_psychopharm_dict(self) -> Dict:
        """Get dictionary of psychopharmacological properties."""
        return {
            "receptor_profiles": self._format_receptor_profiles(),
            "bbb_properties": self._format_bbb_properties(),
            "psychoactive_properties": self._format_psychoactive_properties(),
            "nootropic_properties": self._format_nootropic_properties(),
            "duration_metrics": self._format_duration_metrics(),
            "tolerance_data": self._format_tolerance_data(),
        }

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

    def merge(self, other: "CompoundData") -> None:
        """Merge data from another CompoundData object."""
        # Merge basic data
        self._merge_basic_data(other)

        # Merge analysis results
        self._merge_analysis_results(other)

        # Merge predictions
        self._merge_predictions(other)

        # Merge web data
        self._merge_web_data(other)

        # Merge psychopharmacological properties
        self._merge_psychopharm_data(other)

        # Update metadata
        self.last_updated = datetime.now().isoformat()

    def _merge_basic_data(self, other: "CompoundData") -> None:
        """Merge basic data fields."""
        # Merge sets
        self.other_names.update(other.other_names)
        self.data_sources.update(other.data_sources)
        self.reference_dois.update(other.reference_dois)
        self.reference_pmids.update(other.reference_pmids)

        # Merge dictionaries
        self.reference_urls.update(other.reference_urls)

        # Merge ranked names
        if other.common_name_1_results > self.common_name_1_results:
            self.common_name_1 = other.common_name_1
            self.common_name_1_results = other.common_name_1_results

        if other.common_name_2_results > self.common_name_2_results:
            self.common_name_2 = other.common_name_2
            self.common_name_2_results = other.common_name_2_results

        if other.common_name_3_results > self.common_name_3_results:
            self.common_name_3 = other.common_name_3
            self.common_name_3_results = other.common_name_3_results

        # Merge target data
        self._merge_target_data(other)

        # Merge all predictions
        self._merge_predictions(other)

        # Merge web data
        self._merge_web_data(other)

        # Merge analysis results
        self._merge_analysis_results(other)

        # Merge psychopharmacological properties
        self._merge_receptor_profiles(other)
        self._merge_bbb_properties(other)
        self._merge_psychoactive_properties(other)
        self._merge_nootropic_properties(other)
        self._merge_duration_metrics(other)
        self._merge_tolerance_data(other)

        # Update metadata
        self.last_updated = datetime.now().isoformat()

    def _merge_data_field(self, target: Dict, key: str, value) -> None:
        """Merge a single data field based on its type."""
        if isinstance(value, (list, set)):
            if key not in target:
                target[key] = type(value)()
            if isinstance(value, list):
                target[key].extend(item for item in value if item not in target[key])
            else:  # set
                target[key].update(value)
        elif isinstance(value, dict):
            if key not in target:
                target[key] = {}
            target[key].update(value)
        else:
            target[key] = value

    def _merge_target_data(self, other: "CompoundData") -> None:
        """Merge target data from another CompoundData object."""
        # Create lookup of existing targets by name
        existing_targets = {target.common_name: target for target in self.targets}

        # Merge target data
        for other_target in other.targets:
            if other_target.common_name in existing_targets:
                # Update existing target if new data has more references
                existing = existing_targets[other_target.common_name]
                if len(other_target.reference_dois) > len(existing.reference_dois):
                    existing_targets[other_target.common_name] = other_target
            else:
                # Add new target
                self.targets.append(other_target)

        # Update primary target if needed
        has_more_refs = len(other.reference_dois) > len(self.reference_dois)
        if other.primary_target and (not self.primary_target or has_more_refs):
            self.primary_target = other.primary_target

    def _merge_predictions(self, other: "CompoundData") -> None:
        """Merge all ML predictions from another instance."""
        self._merge_toxicity_predictions(other)
        self._merge_abuse_potential(other)
        self._merge_binding_predictions(other)
        self._merge_activity_predictions(other)
        self._merge_mechanism_predictions(other)

    def _merge_toxicity_predictions(self, other: "CompoundData") -> None:
        """Merge toxicity predictions."""
        for tox_type, tox_data in other.toxicity_predictions.items():
            if tox_type not in self.toxicity_predictions:
                self.toxicity_predictions[tox_type] = tox_data
            else:
                # Keep prediction with higher confidence
                old_conf = self.toxicity_predictions[tox_type].get("confidence", 0)
                if tox_data.get("confidence", 0) > old_conf:
                    self.toxicity_predictions[tox_type] = tox_data

    def _merge_abuse_potential(self, other: "CompoundData") -> None:
        """Merge abuse potential predictions."""
        for abuse_type, abuse_data in other.abuse_potential.items():
            if abuse_type not in self.abuse_potential:
                self.abuse_potential[abuse_type] = abuse_data
            else:
                # Keep prediction with higher confidence
                old_conf = self.abuse_potential[abuse_type].get("confidence", 0)
                if abuse_data.get("confidence", 0) > old_conf:
                    self.abuse_potential[abuse_type] = abuse_data

    def _merge_binding_predictions(self, other: "CompoundData") -> None:
        """Merge binding predictions."""
        self.binding_predictions.extend(
            pred for pred in other.binding_predictions if pred not in self.binding_predictions
        )

    def _merge_activity_predictions(self, other: "CompoundData") -> None:
        """Merge activity predictions."""
        self.activity_predictions.extend(
            pred for pred in other.activity_predictions if pred not in self.activity_predictions
        )

    def _merge_mechanism_predictions(self, other: "CompoundData") -> None:
        """Merge mechanism predictions."""
        self.mechanism_predictions.extend(
            pred for pred in other.mechanism_predictions if pred not in self.mechanism_predictions
        )

    def get_predictions_dict(self) -> Dict:
        """Get dictionary of all ML predictions."""
        return {
            "toxicity_predictions": self._format_toxicity_predictions(),
            "abuse_potential": self._format_abuse_potential(),
            "binding_predictions": self._format_binding_predictions(),
            "activity_predictions": self._format_activity_predictions(),
            "mechanism_predictions": self._format_mechanism_predictions(),
        }

    def _merge_web_data(self, other: "CompoundData") -> None:
        """Merge web-enriched data from another CompoundData object."""
        self._merge_community_data(other.community_data)
        self._merge_literature_data(other.literature_data)
        self._merge_regulatory_data(other.regulatory_data)
        self._merge_experience_reports(other.experience_reports)
        self._merge_safety_profile(other.safety_profile)

    def _merge_community_data(self, other_data: Dict) -> None:
        """Merge community data."""
        if not other_data:
            return

        if not self.community_data:
            self.community_data = {}

        for key, value in other_data.items():
            if isinstance(value, list):
                self._merge_list_field(self.community_data, key, value)
            elif isinstance(value, dict):
                self._merge_dict_field(self.community_data, key, value)
            else:
                self.community_data[key] = value

    def _merge_literature_data(self, other_data: Dict) -> None:
        """Merge literature data."""
        if other_data:
            if not self.literature_data:
                self.literature_data = {}
            self.literature_data.update(other_data)

    def _merge_regulatory_data(self, other_data: Dict) -> None:
        """Merge regulatory data."""
        if other_data:
            if not self.regulatory_data:
                self.regulatory_data = {}
            self.regulatory_data.update(other_data)

    def _merge_experience_reports(self, other_reports: List[Dict]) -> None:
        """Merge experience reports."""
        self.experience_reports.extend(report for report in other_reports if report not in self.experience_reports)

    def _merge_safety_profile(self, other_profile: Dict) -> None:
        """Merge safety profile data."""
        if other_profile:
            if not self.safety_profile:
                self.safety_profile = {}
            self.safety_profile.update(other_profile)

    def _merge_list_field(self, target_dict: Dict, key: str, value: List) -> None:
        """Merge a list field into a dictionary."""
        if key not in target_dict:
            target_dict[key] = []
        target_dict[key].extend(item for item in value if item not in target_dict[key])

    def _merge_dict_field(self, target_dict: Dict, key: str, value: Dict) -> None:
        """Merge a dictionary field into a dictionary."""
        if key not in target_dict:
            target_dict[key] = {}
        target_dict[key].update(value)

    def _merge_analysis_results(self, other: "CompoundData") -> None:
        """Merge analysis results from another CompoundData object."""
        # Merge pharmacophores
        self.pharmacophores.extend(pharm for pharm in other.pharmacophores if pharm not in self.pharmacophores)

        # Merge structural alerts
        self.structural_alerts.extend(alert for alert in other.structural_alerts if alert not in self.structural_alerts)

        # Merge receptor interactions
        for receptor, interactions in other.receptor_interactions.items():
            if receptor not in self.receptor_interactions:
                self.receptor_interactions[receptor] = []
            self.receptor_interactions[receptor].extend(
                inter for inter in interactions if inter not in self.receptor_interactions[receptor]
            )

        # Merge mechanism predictions
        self.mechanism_predictions.extend(
            mech for mech in other.mechanism_predictions if mech not in self.mechanism_predictions
        )

        # Merge SAR analysis
        if other.sar_analysis:
            if not self.sar_analysis:
                self.sar_analysis = {}
            self.sar_analysis.update(other.sar_analysis)

    def _merge_receptor_profiles(self, other: "CompoundData") -> None:
        """Merge receptor binding profiles."""
        for receptor, data in other.receptor_profiles.items():
            if receptor not in self.receptor_profiles:
                self.receptor_profiles[receptor] = data
            else:
                # Keep data with higher confidence
                old_conf = self.receptor_profiles[receptor].get("confidence", 0)
                if data.get("confidence", 0) > old_conf:
                    self.receptor_profiles[receptor] = data

    def _merge_bbb_properties(self, other: "CompoundData") -> None:
        """Merge blood-brain barrier properties."""
        if other.bbb_score > self.bbb_score:
            self.bbb_permeability = other.bbb_permeability
            self.bbb_score = other.bbb_score
            self.p_glycoprotein_substrate = other.p_glycoprotein_substrate

    def _merge_psychoactive_properties(self, other: "CompoundData") -> None:
        """Merge psychoactive classification data."""
        if other.psychoactive_class != PsychoactiveClass.UNKNOWN:
            if self.psychoactive_class == PsychoactiveClass.UNKNOWN:
                self.psychoactive_class = other.psychoactive_class
            else:
                self.secondary_classes.add(self.psychoactive_class)
                self.psychoactive_class = other.psychoactive_class
        self.secondary_classes.update(other.secondary_classes)
        self.effect_profile.update(other.effect_profile)

    def _merge_nootropic_properties(self, other: "CompoundData") -> None:
        """Merge nootropic properties."""
        self.nootropic_mechanisms.update(other.nootropic_mechanisms)
        self.cognitive_effects.update(other.cognitive_effects)
        self.side_effects.update(other.side_effects)

    def _merge_duration_metrics(self, other: "CompoundData") -> None:
        """Merge duration metrics."""
        if not self.onset_time:
            self.onset_time = other.onset_time
        if not self.duration:
            self.duration = other.duration
        if not self.half_life:
            self.half_life = other.half_life

    def _merge_tolerance_data(self, other: "CompoundData") -> None:
        """Merge tolerance and withdrawal data."""
        self.tolerance_profile.update(other.tolerance_profile)
        self.withdrawal_profile.update(other.withdrawal_profile)
        self.cross_tolerance.update(other.cross_tolerance)

    def _summarize_binding_data(self) -> Dict:
        """Create a summary of binding data."""
        summary = {
            "total_targets": len(self.targets),
            "strongest_binding": None,
            "primary_targets": [],
            "target_families": set(),
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

        summary["target_families"] = list(summary["target_families"])
        return summary

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
                concern for concern in report["safety_concerns"] if concern not in summary["safety_concerns"]
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
            return sum(1 for r in self.experience_reports if "combinations" in r and combo in r["combinations"])

        summary["common_combinations"].sort(key=get_combo_frequency, reverse=True)

    def _calculate_dosage_stats(self, summary: Dict) -> None:
        """Calculate statistics for dosage information."""
        for route in summary["dosage_info"]:
            doses = summary["dosage_info"][route]
            summary["dosage_info"][route] = {
                "min": min(doses),
                "max": max(doses),
                "avg": sum(doses) / len(doses),
                "count": len(doses),
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
        summary["alerts"].extend(alert["description"] for alert in self.structural_alerts)

    def _add_toxicity_warnings(self, summary: Dict) -> None:
        """Add toxicity warnings to safety summary."""
        for tox_type, tox_data in self.toxicity_predictions.items():
            if tox_data.get("probability", 0) > 0.7:
                summary["warnings"].append(
                    {
                        "type": tox_type,
                        "probability": tox_data["probability"],
                        "severity": tox_data.get("severity", "unknown"),
                        "mechanisms": tox_data.get("mechanisms", []),
                    }
                )

    def _add_abuse_warnings(self, summary: Dict) -> None:
        """Add abuse potential warnings to safety summary."""
        for abuse_type, abuse_data in self.abuse_potential.items():
            if abuse_data.get("probability", 0) > 0.7:
                summary["warnings"].append(
                    {
                        "type": f"abuse_{abuse_type}",
                        "probability": abuse_data["probability"],
                        "risk_level": abuse_data.get("risk_level", "unknown"),
                        "mechanisms": abuse_data.get("mechanisms", []),
                    }
                )

    def _add_safety_profile_data(self, summary: Dict) -> None:
        """Add safety profile data to safety summary."""
        if self.safety_profile:
            for key in ["contraindications", "interactions", "risk_factors", "overdose_risks", "long_term_risks"]:
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
            "approval_status": self.approval_status.copy(),
            "clinical_status": self.clinical_status.copy(),
            "research_status": self.research_status.copy(),
        }

    def _add_regulatory_data(self, status: Dict) -> None:
        """Add regulatory data to status."""
        if self.regulatory_data:
            status.update(self.regulatory_data)

    def _add_toxicity_status(self, status: Dict) -> None:
        """Add toxicity warnings to regulatory status."""
        for tox_type, tox_data in self.toxicity_predictions.items():
            if tox_data.get("probability", 0) > 0.7:
                status["warnings"].append(
                    {
                        "type": f"toxicity_{tox_type}",
                        "probability": tox_data["probability"],
                        "severity": tox_data.get("severity", "unknown"),
                    }
                )

    def _add_abuse_status(self, status: Dict) -> None:
        """Add abuse potential warnings to regulatory status."""
        for abuse_type, abuse_data in self.abuse_potential.items():
            if abuse_data.get("probability", 0) > 0.7:
                status["warnings"].append(
                    {
                        "type": f"abuse_{abuse_type}",
                        "probability": abuse_data["probability"],
                        "risk_level": abuse_data.get("risk_level", "unknown"),
                    }
                )

    def merge_web_data(self, other: "CompoundData") -> None:
        """Merge web-enriched data from another instance."""
        self._merge_patent_data(other)
        self._merge_swiss_data(other)
        self._merge_community_data(other)
        self._merge_literature_data(other)
        self._merge_regulatory_data(other)
        self._merge_experience_reports(other)
        self._merge_safety_profile(other)
        self._merge_status_data(other)

    def _merge_status_data(self, other: "CompoundData") -> None:
        """Merge status data."""
        self.approval_status.update(other.approval_status)
        self.clinical_status.update(other.clinical_status)
        self.research_status.update(other.research_status)
        self.risk_factors.extend(factor for factor in other.risk_factors if factor not in self.risk_factors)
        self.overdose_risks.extend(risk for risk in other.overdose_risks if risk not in self.overdose_risks)
        self.long_term_risks.extend(risk for risk in other.long_term_risks if risk not in self.long_term_risks)

    def _format_toxicity_predictions(self) -> Dict:
        """Format toxicity predictions for output."""
        formatted = {}
        for tox_type, prediction in self.toxicity_predictions.items():
            formatted[tox_type] = {
                "probability": prediction.get("probability", 0),
                "confidence": prediction.get("confidence", 0),
                "severity": prediction.get("severity", "unknown"),
                "mechanisms": prediction.get("mechanisms", []),
                "supporting_data": prediction.get("supporting_data", {}),
            }
        return formatted

    def _format_abuse_potential(self) -> Dict:
        """Format abuse potential data for output."""
        formatted = {}
        for abuse_type, prediction in self.abuse_potential.items():
            formatted[abuse_type] = {
                "probability": prediction.get("probability", 0),
                "confidence": prediction.get("confidence", 0),
                "risk_level": prediction.get("risk_level", "unknown"),
                "mechanisms": prediction.get("mechanisms", []),
                "receptor_involvement": prediction.get("receptor_involvement", {}),
                "supporting_data": prediction.get("supporting_data", {}),
            }
        return formatted

    def _format_binding_predictions(self) -> List[Dict]:
        """Format binding predictions for output."""
        return [
            {
                "target": pred.get("target", ""),
                "probability": pred.get("probability", 0),
                "affinity_type": pred.get("affinity_type", ""),
                "predicted_value": pred.get("predicted_value"),
                "confidence": pred.get("confidence", 0),
                "supporting_features": pred.get("supporting_features", []),
            }
            for pred in self.binding_predictions
        ]

    def _format_activity_predictions(self) -> List[Dict]:
        """Format activity predictions for output."""
        return [
            {
                "activity_type": pred.get("activity_type", ""),
                "probability": pred.get("probability", 0),
                "predicted_value": pred.get("predicted_value"),
                "confidence": pred.get("confidence", 0),
                "mechanism": pred.get("mechanism", ""),
                "supporting_data": pred.get("supporting_data", {}),
            }
            for pred in self.activity_predictions
        ]

    def _format_mechanism_predictions(self) -> List[Dict]:
        """Format mechanism predictions for output."""
        return [
            {
                "mechanism": pred.get("mechanism", ""),
                "probability": pred.get("probability", 0),
                "confidence": pred.get("confidence", 0),
                "supporting_evidence": pred.get("supporting_evidence", []),
                "receptor_systems": pred.get("receptor_systems", []),
            }
            for pred in self.mechanism_predictions
        ]

    def to_dict(self, include_predictions: bool = True) -> Dict:
        """Convert compound data to dictionary format."""
        data = {
            # Basic identifiers
            "name": self.name,
            "smiles": self.smiles,
            "inchi": self.inchi,
            "inchi_key": self.inchi_key,
            "cas_number": self.cas_number,
            "common_name_1": self.common_name_1,
            "common_name_2": self.common_name_2,
            "common_name_3": self.common_name_3,
            "common_name_1_results": self.common_name_1_results,
            "common_name_2_results": self.common_name_2_results,
            "common_name_3_results": self.common_name_3_results,
            "other_names": list(self.other_names),
            "iupac_name": self.iupac_name,
            "compound_type": self.compound_type.value,
            # Chemical properties
            "molecular_weight": self.molecular_weight,
            "logp": self.logp,
            "hbd": self.hbd,
            "hba": self.hba,
            "tpsa": self.tpsa,
            "rotatable_bonds": self.rotatable_bonds,
            "charge": self.charge,
            "stereocenter_count": self.stereocenter_count,
            "ring_count": self.ring_count,
            # Database IDs
            "chembl_id": self.chembl_id,
            "pubchem_cid": self.pubchem_cid,
            "pubchem_sid": self.pubchem_sid,
            "drugbank_id": self.drugbank_id,
            "bindingdb_id": self.bindingdb_id,
            # Target data
            "targets": [vars(target) for target in self.targets],
            "primary_target": self.primary_target,
            "primary_activity": self.primary_activity,
            "mechanism_of_action": self.mechanism_of_action,
            "pharmacology": self.pharmacology,
            "toxicity": self.toxicity,
            "metabolism": self.metabolism,
            # Source info
            "data_sources": list(self.data_sources),
            "reference_dois": list(self.reference_dois),
            "reference_pmids": list(self.reference_pmids),
            "reference_urls": self.reference_urls,
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
            # Patent data
            "patent_count": self.patent_count,
            "patent_numbers": self.patent_data.get("numbers", []),
            "patent_titles": self.patent_data.get("titles", []),
            # Legal status
            "legal_status": {country: status.value for country, status in self.legal_status.items()},
            "scheduling": self.scheduling,
            # Analysis summaries
            "binding_summary": self._summarize_binding_data(),
            "community_summary": self._summarize_community_data(),
            "safety_summary": self._get_safety_summary(),
            "regulatory_status": self._get_regulatory_status(),
            # Metadata
            "last_updated": self.last_updated,
            "version": self.version,
        }

        # Add ML predictions if requested
        if include_predictions:
            for pred_type, result in self._prediction_cache.items():
                data[f"{pred_type}_prediction"] = {
                    "value": result.value,
                    "confidence": result.confidence,
                    "supporting_data": result.supporting_data,
                }

        return data

    def to_json(self, include_predictions: bool = True) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(include_predictions=include_predictions))

    @classmethod
    def from_json(cls, json_str: str) -> "CompoundData":
        """Create CompoundData from JSON string."""
        data = json.loads(json_str)
        return cls(**data)
