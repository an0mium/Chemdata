"""Data models and validation logic."""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Dict, Set


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
    target_common_name: str = "N/A"
    target_protein_name: str = "N/A"
    target_gene_name: str = "N/A"
    affinity_value: float = 0.0
    affinity_type: str = "N/A"  # Ki, IC50, Kd, or EC50
    affinity_unit: str = "N/A"
    assay_description: str = "N/A"
    reference: str = "N/A"
    confidence_score: float = 0.0


@dataclass
class CompoundData:
    """Represents comprehensive chemical compound data."""
    # Core identifiers
    cas: str = "N/A"
    name: str = "N/A"
    common_name_1: str = "N/A"  # Most common name by web search results
    common_name_2: str = "N/A"  # Second most common name
    common_name_3: str = "N/A"  # Third most common name
    common_name_1_results: int = 0  # Number of search results for name 1
    common_name_2_results: int = 0  # Number of search results for name 2
    common_name_3_results: int = 0  # Number of search results for name 3
    other_names: Set[str] = field(default_factory=set)  # Additional names
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
    
    # Binding data (up to 12 targets, split into individual columns)
    target_1_common_name: str = "N/A"
    target_1_protein_name: str = "N/A"
    target_1_gene_name: str = "N/A"
    target_1_affinity: float = 0.0
    target_1_affinity_unit: str = "N/A"
    target_1_affinity_type: str = "N/A"
    target_1_pubmed_results: int = 0
    
    target_2_common_name: str = "N/A"
    target_2_protein_name: str = "N/A"
    target_2_gene_name: str = "N/A"
    target_2_affinity: float = 0.0
    target_2_affinity_unit: str = "N/A"
    target_2_affinity_type: str = "N/A"
    target_2_pubmed_results: int = 0
    
    target_3_common_name: str = "N/A"
    target_3_protein_name: str = "N/A"
    target_3_gene_name: str = "N/A"
    target_3_affinity: float = 0.0
    target_3_affinity_unit: str = "N/A"
    target_3_affinity_type: str = "N/A"
    target_3_pubmed_results: int = 0
    
    target_4_common_name: str = "N/A"
    target_4_protein_name: str = "N/A"
    target_4_gene_name: str = "N/A"
    target_4_affinity: float = 0.0
    target_4_affinity_unit: str = "N/A"
    target_4_affinity_type: str = "N/A"
    target_4_pubmed_results: int = 0
    
    target_5_common_name: str = "N/A"
    target_5_protein_name: str = "N/A"
    target_5_gene_name: str = "N/A"
    target_5_affinity: float = 0.0
    target_5_affinity_unit: str = "N/A"
    target_5_affinity_type: str = "N/A"
    target_5_pubmed_results: int = 0
    
    target_6_common_name: str = "N/A"
    target_6_protein_name: str = "N/A"
    target_6_gene_name: str = "N/A"
    target_6_affinity: float = 0.0
    target_6_affinity_unit: str = "N/A"
    target_6_affinity_type: str = "N/A"
    target_6_pubmed_results: int = 0
    
    target_7_common_name: str = "N/A"
    target_7_protein_name: str = "N/A"
    target_7_gene_name: str = "N/A"
    target_7_affinity: float = 0.0
    target_7_affinity_unit: str = "N/A"
    target_7_affinity_type: str = "N/A"
    target_7_pubmed_results: int = 0
    
    target_8_common_name: str = "N/A"
    target_8_protein_name: str = "N/A"
    target_8_gene_name: str = "N/A"
    target_8_affinity: float = 0.0
    target_8_affinity_unit: str = "N/A"
    target_8_affinity_type: str = "N/A"
    target_8_pubmed_results: int = 0
    
    target_9_common_name: str = "N/A"
    target_9_protein_name: str = "N/A"
    target_9_gene_name: str = "N/A"
    target_9_affinity: float = 0.0
    target_9_affinity_unit: str = "N/A"
    target_9_affinity_type: str = "N/A"
    target_9_pubmed_results: int = 0
    
    target_10_common_name: str = "N/A"
    target_10_protein_name: str = "N/A"
    target_10_gene_name: str = "N/A"
    target_10_affinity: float = 0.0
    target_10_affinity_unit: str = "N/A"
    target_10_affinity_type: str = "N/A"
    target_10_pubmed_results: int = 0
    
    target_11_common_name: str = "N/A"
    target_11_protein_name: str = "N/A"
    target_11_gene_name: str = "N/A"
    target_11_affinity: float = 0.0
    target_11_affinity_unit: str = "N/A"
    target_11_affinity_type: str = "N/A"
    target_11_pubmed_results: int = 0
    
    target_12_common_name: str = "N/A"
    target_12_protein_name: str = "N/A"
    target_12_gene_name: str = "N/A"
    target_12_affinity: float = 0.0
    target_12_affinity_unit: str = "N/A"
    target_12_affinity_type: str = "N/A"
    target_12_pubmed_results: int = 0
    
    # Pharmacology & activity
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
    
    # Source links and references
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
    
    # Additional metadata
    description: str = "N/A"
    primary_target: str = "N/A"
    data_sources: Set[str] = field(default_factory=set)  # Track where data came from
    last_updated: str = field(
        default_factory=lambda: datetime.now().isoformat()
    )
    
    def format_numeric_values(self):
        """Format numeric values to specified precision."""
        self.logp = float(f"{self.logp:.5f}".rstrip('0').rstrip('.'))
        self.tpsa = float(f"{self.tpsa:.5f}".rstrip('0').rstrip('.'))
        self.molecular_weight = float(f"{self.molecular_weight:.5f}".rstrip('0').rstrip('.'))
    
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
        for target_num in range(1, 13):
            affinity = getattr(self, f'target_{target_num}_affinity')
            if affinity < 0:
                errors.append(
                    f"Invalid binding affinity value for target {target_num}: {affinity}")
                
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
        # Merge names
        if other.common_name_1_results > self.common_name_1_results:
            self.common_name_1 = other.common_name_1
            self.common_name_1_results = other.common_name_1_results
            
        if other.common_name_2_results > self.common_name_2_results:
            self.common_name_2 = other.common_name_2
            self.common_name_2_results = other.common_name_2_results
            
        if other.common_name_3_results > self.common_name_3_results:
            self.common_name_3 = other.common_name_3
            self.common_name_3_results = other.common_name_3_results
            
        self.other_names.update(other.other_names)
        self.data_sources.update(other.data_sources)
        
        # Merge dictionaries
        self.legal_status.update(other.legal_status)
        self.scheduling.update(other.scheduling)
        
        # Merge binding data
        for i in range(1, 13):
            other_pubmed_results = getattr(other, f'target_{i}_pubmed_results')
            if other_pubmed_results > getattr(self, f'target_{i}_pubmed_results'):
                # Copy all target fields if the other has more PubMed results
                for field in ['common_name', 'protein_name', 'gene_name', 
                            'affinity', 'affinity_unit', 'affinity_type', 
                            'pubmed_results']:
                    setattr(self, f'target_{i}_{field}',
                           getattr(other, f'target_{i}_{field}'))
        
        # Update scalar fields if they have values
        excluded_fields = {'other_names', 'data_sources', 'legal_status', 'scheduling'}
        excluded_prefixes = ['common_name_', 'target_']
        
        def is_excluded_field(name: str) -> bool:
            """Check if field should be excluded from updates."""
            if name in excluded_fields:
                return True
            return any(name.startswith(prefix) for prefix in excluded_prefixes)
        
        # Update non-excluded fields
        for name, _ in self.__dataclass_fields__.items():
            if not is_excluded_field(name):
                other_value = getattr(other, name)
                if other_value != "N/A" and other_value is not None:
                    setattr(self, name, other_value)
                    
        # Update timestamp
        self.last_updated = datetime.now().isoformat()


class ValidationError(Exception):
    """Raised when compound data validation fails."""
    pass
