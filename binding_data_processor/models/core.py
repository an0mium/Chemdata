"""Core compound data structures.

This module provides the base CompoundData class and related enums for
representing chemical compounds and their associated core data.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Dict, List, Optional, Set
import json

from .validation import ValidationMixin
from .enrichment import WebEnrichmentMixin
from .predictions import PredictionsMixin
from .analysis import AnalysisMixin
from .psychopharm import PsychopharmMixin


class CompoundType(Enum):
    """Types of chemical compounds."""

    NEUROTRANSMITTER = "neurotransmitter"
    PSYCHOACTIVE = "psychoactive"
    RESEARCH_CHEMICAL = "research_chemical"
    NPS = "novel_psychoactive_substance"
    PHARMACEUTICAL = "pharmaceutical"
    NATURAL_PRODUCT = "natural_product"
    NOOTROPIC = "nootropic"
    OTHER = "other"


class LegalStatus(Enum):
    """Legal status classifications."""

    LEGAL = "legal"
    CONTROLLED = "controlled"
    ILLEGAL = "illegal"
    RESEARCH_ONLY = "research_only"
    UNSCHEDULED = "unscheduled"
    PRESCRIPTION = "prescription_only"
    OTC = "over_the_counter"


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
class CompoundData(
    ValidationMixin,
    WebEnrichmentMixin,
    PredictionsMixin,
    AnalysisMixin,
    PsychopharmMixin,
):
    """Enhanced data class for chemical compound information."""

    # Core identifiers (required)
    name: str
    smiles: str

    # Basic identifiers (optional)
    cas_number: Optional[str] = None
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    common_names: List[str] = field(default_factory=list)
    iupac_name: Optional[str] = None
    compound_type: CompoundType = CompoundType.OTHER

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

    # Source information
    data_sources: Set[str] = field(default_factory=set)
    reference_dois: Set[str] = field(default_factory=set)
    reference_pmids: Set[str] = field(default_factory=set)
    reference_urls: Dict[str, str] = field(default_factory=dict)

    # Legal & classification
    legal_status: Dict[str, LegalStatus] = field(default_factory=dict)
    scheduling: Dict[str, str] = field(default_factory=dict)

    # Metadata
    last_updated: str = field(default_factory=lambda: datetime.now().isoformat())
    version: str = "2.0.0"

    def __post_init__(self):
        """Validate required fields and initialize collections."""
        if not self.name:
            raise ValueError("Compound name is required")
        if not self.smiles:
            raise ValueError("SMILES string is required")

        # Initialize all collection fields
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

        # Run validation
        self.validate()
        
        # Format numeric values
        self.format_numeric_values()

    def format_numeric_values(self):
        """Format numeric values to specified precision."""
        self.logp = float(f"{self.logp:.5f}".rstrip('0').rstrip('.'))
        self.tpsa = float(f"{self.tpsa:.5f}".rstrip('0').rstrip('.'))
        self.molecular_weight = float(f"{self.molecular_weight:.5f}".rstrip('0').rstrip('.'))

    def to_dict(self, include_predictions: bool = True, include_psychopharm: bool = True) -> Dict:
        """Convert compound data to dictionary format."""
        data = {
            # Basic identifiers
            "name": self.name,
            "smiles": self.smiles,
            "inchi": self.inchi,
            "inchi_key": self.inchi_key,
            "cas_number": self.cas_number,
            "common_names": self.common_names[:3],  # Top 3 common names
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

            # Source info
            "data_sources": list(self.data_sources),
            "reference_dois": list(self.reference_dois),
            "reference_pmids": list(self.reference_pmids),
            "reference_urls": self.reference_urls,

            # Legal status
            "legal_status": {
                country: status.value
                for country, status in self.legal_status.items()
            },
            "scheduling": self.scheduling,

            # Metadata
            "last_updated": self.last_updated,
            "version": self.version,
        }

        # Add summaries
        data.update({
            "binding_summary": self._summarize_binding_data(),
            "community_summary": self._summarize_community_data(),
            "safety_summary": self._get_safety_summary(),
            "regulatory_status": self._get_regulatory_status(),
        })

        # Add predictions if requested
        if include_predictions:
            data.update(self.get_predictions_dict())

        # Add web enrichment data
        data.update(self.get_enrichment_dict())

        # Add analysis results
        data.update(self.get_analysis_dict())

        # Add psychopharmacological properties
        if include_psychopharm:
            data.update(self.get_psychopharm_dict())

        return data

    def to_json(self, include_predictions: bool = True) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(include_predictions=include_predictions))

    @classmethod
    def from_json(cls, json_str: str) -> "CompoundData":
        """Create CompoundData from JSON string."""
        data = json.loads(json_str)
        return cls(**data)
