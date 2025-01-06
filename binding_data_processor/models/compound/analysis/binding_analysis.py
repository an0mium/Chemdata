"""Binding data model.

This module provides the BindingData class for representing binding affinity
and activity data for chemical compounds, including:
- Target information (name, gene, protein)
- Affinity measurements (Ki, IC50, Kd, EC50)
- Activity type classification
- Assay conditions
- Source information
"""

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class BindingData:
    """Data class for binding affinity information."""

    # Target information
    target_name: str
    target_gene: Optional[str] = None
    target_protein: Optional[str] = None
    target_organism: str = "human"
    target_type: Optional[str] = None

    # Affinity measurements
    affinity_type: str  # Ki, IC50, Kd, EC50
    affinity_value: float
    affinity_unit: str = "nM"
    affinity_modifier: Optional[str] = None  # >, <, ~
    affinity_error: Optional[float] = None
    affinity_temperature: Optional[float] = None
    affinity_ph: Optional[float] = None

    # Activity information
    activity_type: str = "unknown"
    activity_description: Optional[str] = None
    activity_concentration: Optional[float] = None
    activity_unit: str = "nM"

    # Kinetics
    kon: Optional[float] = None
    kon_unit: str = "M-1-s-1"
    koff: Optional[float] = None
    koff_unit: str = "s-1"

    # Source information
    source: str = "BindingDB"
    assay_type: Optional[str] = None
    assay_description: Optional[str] = None
    doi: Optional[str] = None
    pmid: Optional[str] = None
    patent_number: Optional[str] = None
    reference_type: Optional[str] = None
    reference_url: Optional[str] = None

    # Additional data
    conditions: Dict = field(default_factory=dict)
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        """Validate required fields and initialize collections."""
        if not self.target_name:
            raise ValueError("Target name is required")
        if not self.affinity_type or self.affinity_value is None:
            raise ValueError("Affinity type and value are required")

        # Initialize collections if None
        if self.conditions is None:
            self.conditions = {}
        if self.metadata is None:
            self.metadata = {}

    def to_dict(self) -> Dict:
        """
        Convert binding data to dictionary format.

        Returns:
            Dictionary representation of binding data
        """
        data = {
            # Target info
            "target_name": self.target_name,
            "target_gene": self.target_gene,
            "target_protein": self.target_protein,
            "target_organism": self.target_organism,
            "target_type": self.target_type,
            # Affinity data
            "affinity_type": self.affinity_type,
            "affinity_value": self.affinity_value,
            "affinity_unit": self.affinity_unit,
            "affinity_modifier": self.affinity_modifier,
            "affinity_error": self.affinity_error,
            "affinity_temperature": self.affinity_temperature,
            "affinity_ph": self.affinity_ph,
            # Activity data
            "activity_type": self.activity_type,
            "activity_description": self.activity_description,
            "activity_concentration": self.activity_concentration,
            "activity_unit": self.activity_unit,
            # Kinetics
            "kon": self.kon,
            "kon_unit": self.kon_unit,
            "koff": self.koff,
            "koff_unit": self.koff_unit,
            # Source info
            "source": self.source,
            "assay_type": self.assay_type,
            "assay_description": self.assay_description,
            "doi": self.doi,
            "pmid": self.pmid,
            "patent_number": self.patent_number,
            "reference_type": self.reference_type,
            "reference_url": self.reference_url,
        }

        # Add conditions and metadata if present
        if self.conditions:
            data["conditions"] = self.conditions
        if self.metadata:
            data["metadata"] = self.metadata

        return data

    def merge(self, other: "BindingData") -> None:
        """
        Merge data from another BindingData object.

        Args:
            other: BindingData object to merge from
        """
        # Only merge if target names match
        if self.target_name != other.target_name:
            return

        # Update target info if missing
        if not self.target_gene and other.target_gene:
            self.target_gene = other.target_gene
        if not self.target_protein and other.target_protein:
            self.target_protein = other.target_protein
        if not self.target_type and other.target_type:
            self.target_type = other.target_type

        # Update affinity data if other has stronger affinity
        if other.affinity_value is not None and (
            self.affinity_value is None or other.affinity_value < self.affinity_value
        ):
            self.affinity_type = other.affinity_type
            self.affinity_value = other.affinity_value
            self.affinity_unit = other.affinity_unit
            self.affinity_modifier = other.affinity_modifier
            self.affinity_error = other.affinity_error
            self.affinity_temperature = other.affinity_temperature
            self.affinity_ph = other.affinity_ph

        # Update activity info if missing
        if not self.activity_type or self.activity_type == "unknown":
            self.activity_type = other.activity_type
        if not self.activity_description and other.activity_description:
            self.activity_description = other.activity_description
        if not self.activity_concentration and other.activity_concentration:
            self.activity_concentration = other.activity_concentration
            self.activity_unit = other.activity_unit

        # Update kinetics if missing
        if not self.kon and other.kon:
            self.kon = other.kon
            self.kon_unit = other.kon_unit
        if not self.koff and other.koff:
            self.koff = other.koff
            self.koff_unit = other.koff_unit

        # Merge conditions and metadata
        if other.conditions:
            self.conditions.update(other.conditions)
        if other.metadata:
            self.metadata.update(other.metadata)

    @property
    def kd_calculated(self) -> Optional[float]:
        """Calculate Kd from kon and koff if available."""
        if self.kon is not None and self.koff is not None:
            try:
                return self.koff / self.kon
            except ZeroDivisionError:
                return None
        return None
