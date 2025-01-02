"""Data models and validation logic."""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Dict, List, Optional, Set


class CompoundType(Enum):
    """Types of chemical compounds."""
    NEUROTRANSMITTER = "neurotransmitter"
    PSYCHOACTIVE = "psychoactive"
    RESEARCH_CHEMICAL = "research_chemical"
    NPS = "novel_psychoactive_substance"
    OTHER = "other"


class LegalStatus(Enum):
    """Legal status classifications."""
    LEGAL = "legal"
    CONTROLLED = "controlled"
    ILLEGAL = "illegal"
    RESEARCH_ONLY = "research_only"
    UNSCHEDULED = "unscheduled"


@dataclass
class BindingData:
    """Represents binding affinity data for a compound-target interaction."""
    target: str
    affinity_type: str  # Ki, IC50, Kd, or EC50
    affinity_value: float
    assay_description: str
    reference: str


@dataclass
class CompoundData:
    """Represents comprehensive chemical compound data."""
    # Core identifiers
    cas: str = "N/A"
    name: str = "N/A"
    common_names: Set[str] = field(default_factory=set)  # Including street names
    iupac_name: str = "N/A"
    compound_type: CompoundType = CompoundType.OTHER
    
    # Chemical properties
    smiles: str = "N/A"
    inchi: str = "N/A"
    inchi_key: str = "N/A"
    molecular_weight: float = 0.0
    logp: float = 0.0
    hbd: int = 0  # hydrogen bond donors
    hba: int = 0  # hydrogen bond acceptors
    tpsa: float = 0.0  # topological polar surface area
    rotatable_bonds: int = 0
    
    # Biological & pharmacological data
    binding_data: List[BindingData] = field(default_factory=list)
    primary_activity: str = "N/A"
    mechanism_of_action: str = "N/A"
    pharmacology: str = "N/A"
    toxicity: str = "N/A"
    metabolism: str = "N/A"
    
    # Legal & classification
    legal_status: Dict[str, LegalStatus] = field(default_factory=dict)  # Country -> Status
    scheduling: Dict[str, str] = field(default_factory=dict)  # Country/Region -> Schedule
    
    # Database identifiers
    pubchem_cid: str = "N/A"
    chembl_id: str = "N/A"
    drugbank_id: str = "N/A"
    
    # Source links
    pubchem_url: str = "N/A"
    psychonaut_url: str = "N/A"
    erowid_url: str = "N/A"
    wikipedia_url: str = "N/A"
    
    # Additional metadata
    description: str = "N/A"
    primary_target: str = "N/A"
    data_sources: Set[str] = field(default_factory=set)  # Track where data came from
    last_updated: str = field(
        default_factory=lambda: datetime.now().isoformat()
    )
    
    def __post_init__(self):
        """Validate data after initialization."""
        self._validate()

    def _validate(self):
        """Validate compound data."""
        errors = []
        
        # Validate CAS number format if present
        if self.cas != "N/A" and not self._validate_cas_format(self.cas):
            errors.append(f"Invalid CAS number format: {self.cas}")
            
        # Validate molecular weight
        if self.molecular_weight < 0:
            errors.append(f"Invalid molecular weight: {self.molecular_weight}")
            
        # Validate binding data
        for binding in self.binding_data:
            if binding.affinity_value <= 0:
                errors.append(
                    f"Invalid binding affinity value: {binding.affinity_value}"
                )
                
        if errors:
            raise ValidationError("\n".join(errors))

    def _validate_cas_format(self, cas: str) -> bool:
        """
        Validate CAS number format.
        
        Args:
            cas: CAS number to validate
            
        Returns:
            True if valid, False otherwise
        """
        import re
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

    def merge(self, other: 'CompoundData') -> None:
        """
        Merge data from another compound instance.
        
        Args:
            other: CompoundData instance to merge from
        """
        # Merge sets
        self.common_names.update(other.common_names)
        self.data_sources.update(other.data_sources)
        
        # Merge dictionaries
        self.legal_status.update(other.legal_status)
        self.scheduling.update(other.scheduling)
        
        # Merge lists
        self.binding_data.extend(other.binding_data)
        
        # Update scalar fields if they have values
        for field in self.__dataclass_fields__:
            if field not in {'common_names', 'data_sources', 'legal_status',
                           'scheduling', 'binding_data'}:
                other_value = getattr(other, field)
                if other_value != "N/A" and other_value is not None:
                    setattr(self, field, other_value)
                    
        # Update timestamp
        self.last_updated = datetime.now().isoformat()


class ValidationError(Exception):
    """Raised when compound data validation fails."""
    pass
