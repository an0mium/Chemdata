"""Base compound data model.

This module provides the base CompoundData class with core fields and methods for:
- Basic chemical identifiers
- Chemical properties
- Database IDs
- Common names
- Basic validation
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Set
import re

from .types import CompoundType, LegalStatus, TargetData


@dataclass
class BaseCompound:
    """Base data class for chemical compound information."""
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

    # Source information
    data_sources: Set[str] = field(default_factory=set)
    reference_dois: Set[str] = field(default_factory=set)
    reference_pmids: Set[str] = field(default_factory=set)
    reference_urls: Dict[str, str] = field(default_factory=dict)

    # Legal & classification
    legal_status: Dict[str, LegalStatus] = field(default_factory=dict)  # Country -> Status
    scheduling: Dict[str, str] = field(default_factory=dict)  # Country -> Schedule

    # Metadata
    last_updated: str = field(default_factory=lambda: datetime.now().isoformat())
    version: str = "1.0.0"

    def __post_init__(self):
        """Validate required fields and initialize collections."""
        self._validate()

    def _validate(self):
        """Validate compound data."""
        errors = []
        
        # Run all validation checks
        errors.extend(self._validate_identifiers())
        errors.extend(self._validate_properties())
        
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

    def format_numeric_values(self):
        """Format numeric values to specified precision."""
        self.logp = float(f"{self.logp:.5f}".rstrip('0').rstrip('.'))
        self.tpsa = float(f"{self.tpsa:.5f}".rstrip('0').rstrip('.'))
        self.molecular_weight = float(f"{self.molecular_weight:.5f}".rstrip('0').rstrip('.'))

    def to_dict(self) -> Dict:
        """Convert basic compound data to dictionary format."""
        return {
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
            # Source info
            "data_sources": list(self.data_sources),
            "reference_dois": list(self.reference_dois),
            "reference_pmids": list(self.reference_pmids),
            "reference_urls": self.reference_urls,
            # Legal status
            "legal_status": {
                country: status.value for country, status in self.legal_status.items()
            },
            "scheduling": self.scheduling,
            # Metadata
            "last_updated": self.last_updated,
            "version": self.version,
        }


class ValidationError(Exception):
    """Raised when compound data validation fails."""
    pass
