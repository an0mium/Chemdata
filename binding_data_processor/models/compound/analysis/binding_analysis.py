"""Binding data model and analyzer.

This module provides:
1. BindingData class for representing binding affinity and activity data
2. BindingAnalyzer class for analyzing binding data, including:
   - Affinity calculations and normalization
   - Activity classification
   - Statistical analysis
   - Data validation and quality assessment
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np
from scipy import stats


def default_conditions() -> Dict:
    """Default empty conditions dictionary."""
    return {}


def default_metadata() -> Dict:
    """Default empty metadata dictionary."""
    return {}


@dataclass
class BindingData:
    """Data class for binding affinity information."""

    # Required fields - must come first
    target_name: str
    affinity_type: str  # Ki, IC50, Kd, EC50
    affinity_value: float
    affinity_unit: str
    activity_type: str
    activity_unit: str
    source: str
    kon_unit: str
    koff_unit: str
    target_organism: str

    # Optional target information
    target_gene: Optional[str] = None
    target_protein: Optional[str] = None
    target_type: Optional[str] = None

    # Optional affinity measurements
    affinity_modifier: Optional[str] = None  # >, <, ~
    affinity_error: Optional[float] = None
    affinity_temperature: Optional[float] = None
    affinity_ph: Optional[float] = None

    # Activity information
    activity_description: Optional[str] = None
    activity_concentration: Optional[float] = None

    # Kinetics
    kon: Optional[float] = None
    koff: Optional[float] = None

    # Source information
    assay_type: Optional[str] = None
    assay_description: Optional[str] = None
    doi: Optional[str] = None
    pmid: Optional[str] = None
    patent_number: Optional[str] = None
    reference_type: Optional[str] = None
    reference_url: Optional[str] = None

    # Additional data
    conditions: Dict = field(default_factory=default_conditions)
    metadata: Dict = field(default_factory=default_metadata)

    # Default values for required fields
    def __init__(
        self,
        target_name: str,
        affinity_type: str,
        affinity_value: float,
        affinity_unit: str = "nM",
        activity_type: str = "unknown",
        activity_unit: str = "nM",
        source: str = "BindingDB",
        kon_unit: str = "M-1-s-1",
        koff_unit: str = "s-1",
        target_organism: str = "human",
        **kwargs,
    ):
        """Initialize binding data with required and optional fields."""
        # Initialize required fields
        self.target_name = target_name
        self.affinity_type = affinity_type
        self.affinity_value = affinity_value
        self.affinity_unit = affinity_unit
        self.activity_type = activity_type
        self.activity_unit = activity_unit
        self.source = source
        self.kon_unit = kon_unit
        self.koff_unit = koff_unit
        self.target_organism = target_organism

        # Initialize optional fields with defaults
        self.target_gene = kwargs.get("target_gene")
        self.target_protein = kwargs.get("target_protein")
        self.target_type = kwargs.get("target_type")
        self.affinity_modifier = kwargs.get("affinity_modifier")
        self.affinity_error = kwargs.get("affinity_error")
        self.affinity_temperature = kwargs.get("affinity_temperature")
        self.affinity_ph = kwargs.get("affinity_ph")
        self.activity_description = kwargs.get("activity_description")
        self.activity_concentration = kwargs.get("activity_concentration")
        self.kon = kwargs.get("kon")
        self.koff = kwargs.get("koff")
        self.assay_type = kwargs.get("assay_type")
        self.assay_description = kwargs.get("assay_description")
        self.doi = kwargs.get("doi")
        self.pmid = kwargs.get("pmid")
        self.patent_number = kwargs.get("patent_number")
        self.reference_type = kwargs.get("reference_type")
        self.reference_url = kwargs.get("reference_url")

        # Initialize collections
        self.conditions = kwargs.get("conditions", {})
        self.metadata = kwargs.get("metadata", {})

    def __post_init__(self):
        """Validate required fields and their values."""
        # Validate target information
        if not self.target_name:
            raise ValueError("Target name is required")
        if not self.target_organism:
            raise ValueError("Target organism is required")

        # Validate affinity information
        if not self.affinity_type:
            raise ValueError("Affinity type is required")
        if self.affinity_type not in ["Ki", "IC50", "Kd", "EC50"]:
            raise ValueError("Invalid affinity type. Must be one of: Ki, IC50, Kd, EC50")
        if self.affinity_value is None:
            raise ValueError("Affinity value is required")
        if not isinstance(self.affinity_value, (int, float)):
            raise ValueError("Affinity value must be a number")
        if not self.affinity_unit:
            raise ValueError("Affinity unit is required")

        # Validate activity information
        if not self.activity_type:
            raise ValueError("Activity type is required")
        if not self.activity_unit:
            raise ValueError("Activity unit is required")

        # Validate source information
        if not self.source:
            raise ValueError("Source is required")

        # Validate kinetics units
        if not self.kon_unit:
            raise ValueError("Kon unit is required")
        if not self.koff_unit:
            raise ValueError("Koff unit is required")

        # Validate optional numeric fields if present
        if self.affinity_error is not None and not isinstance(self.affinity_error, (int, float)):
            raise ValueError("Affinity error must be a number")
        if self.affinity_temperature is not None and not isinstance(self.affinity_temperature, (int, float)):
            raise ValueError("Affinity temperature must be a number")
        if self.affinity_ph is not None and not isinstance(self.affinity_ph, (int, float)):
            raise ValueError("Affinity pH must be a number")
        if self.activity_concentration is not None and not isinstance(self.activity_concentration, (int, float)):
            raise ValueError("Activity concentration must be a number")
        if self.kon is not None and not isinstance(self.kon, (int, float)):
            raise ValueError("Kon must be a number")
        if self.koff is not None and not isinstance(self.koff, (int, float)):
            raise ValueError("Koff must be a number")

        # Validate affinity modifier if present
        if self.affinity_modifier is not None and self.affinity_modifier not in [">", "<", "~"]:
            raise ValueError("Invalid affinity modifier. Must be one of: >, <, ~")

        # Initialize collections
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
        # Required fields first
        data = {
            # Required fields
            "target_name": self.target_name,
            "affinity_type": self.affinity_type,
            "affinity_value": self.affinity_value,
            "affinity_unit": self.affinity_unit,
            "activity_type": self.activity_type,
            "activity_unit": self.activity_unit,
            "source": self.source,
            "kon_unit": self.kon_unit,
            "koff_unit": self.koff_unit,
            "target_organism": self.target_organism,
        }

        # Optional target information
        target_info = {
            "target_gene": self.target_gene,
            "target_protein": self.target_protein,
            "target_type": self.target_type,
        }
        data.update({k: v for k, v in target_info.items() if v is not None})

        # Optional affinity measurements
        affinity_info = {
            "affinity_modifier": self.affinity_modifier,
            "affinity_error": self.affinity_error,
            "affinity_temperature": self.affinity_temperature,
            "affinity_ph": self.affinity_ph,
        }
        data.update({k: v for k, v in affinity_info.items() if v is not None})

        # Optional activity information
        activity_info = {
            "activity_description": self.activity_description,
            "activity_concentration": self.activity_concentration,
        }
        data.update({k: v for k, v in activity_info.items() if v is not None})

        # Optional kinetics
        kinetics_info = {
            "kon": self.kon,
            "koff": self.koff,
        }
        data.update({k: v for k, v in kinetics_info.items() if v is not None})

        # Optional source information
        source_info = {
            "assay_type": self.assay_type,
            "assay_description": self.assay_description,
            "doi": self.doi,
            "pmid": self.pmid,
            "patent_number": self.patent_number,
            "reference_type": self.reference_type,
            "reference_url": self.reference_url,
        }
        data.update({k: v for k, v in source_info.items() if v is not None})

        # Additional data
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

        # Update affinity data if other has stronger affinity
        if other.affinity_value is not None and (self.affinity_value is None or other.affinity_value < self.affinity_value):
            # Update all affinity-related fields together
            self.affinity_type = other.affinity_type
            self.affinity_value = other.affinity_value
            self.affinity_unit = other.affinity_unit
            self.affinity_modifier = other.affinity_modifier
            self.affinity_error = other.affinity_error
            self.affinity_temperature = other.affinity_temperature
            self.affinity_ph = other.affinity_ph

        # Update activity fields if current data is unknown/missing
        if self.activity_type == "unknown" and other.activity_type != "unknown":
            self.activity_type = other.activity_type
            self.activity_unit = other.activity_unit
            self.activity_description = other.activity_description
            self.activity_concentration = other.activity_concentration

        # Update optional target info if missing
        if not self.target_gene and other.target_gene:
            self.target_gene = other.target_gene
        if not self.target_protein and other.target_protein:
            self.target_protein = other.target_protein
        if not self.target_type and other.target_type:
            self.target_type = other.target_type

        # Update kinetics if missing
        if not self.kon and other.kon:
            self.kon = other.kon
            self.kon_unit = other.kon_unit
        if not self.koff and other.koff:
            self.koff = other.koff
            self.koff_unit = other.koff_unit

        # Update source info if more detailed
        if other.source != "BindingDB" and self.source == "BindingDB" or len(other.source) > len(self.source):
            self.source = other.source
            self.assay_type = other.assay_type
            self.assay_description = other.assay_description
            self.doi = other.doi
            self.pmid = other.pmid
            self.patent_number = other.patent_number
            self.reference_type = other.reference_type
            self.reference_url = other.reference_url

        # Merge collections
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


class BindingAnalyzer:
    """Analyzer for compound binding data."""

    def __init__(self):
        """Initialize binding analyzer."""
        self.binding_data: List[BindingData] = []
        self.statistics: Dict = {}
        self.quality_scores: Dict = {}

    def add_binding_data(self, data: BindingData) -> None:
        """Add binding data for analysis."""
        self.binding_data.append(data)

    def add_multiple_binding_data(self, data_list: List[BindingData]) -> None:
        """Add multiple binding data entries."""
        self.binding_data.extend(data_list)

    def normalize_affinity_values(self) -> List[float]:
        """Normalize affinity values to a common unit and type."""
        normalized = []
        for data in self.binding_data:
            value = data.affinity_value

            # Convert to nM if needed
            if data.affinity_unit == "μM":
                value *= 1000
            elif data.affinity_unit == "mM":
                value *= 1000000
            elif data.affinity_unit == "pM":
                value /= 1000

            # Adjust based on affinity type
            if data.affinity_type == "IC50":
                value *= 0.5  # Approximate Ki from IC50
            elif data.affinity_type == "EC50":
                value *= 0.5  # Similar approximation

            normalized.append(value)
        return normalized

    def calculate_statistics(self) -> Dict:
        """Calculate statistical measures for binding data."""
        values = self.normalize_affinity_values()
        if not values:
            return {}

        self.statistics = {
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
            "std": float(np.std(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "n": len(values),
            "geometric_mean": float(stats.gmean(values)) if all(v > 0 for v in values) else None,
        }
        return self.statistics

    def assess_data_quality(self) -> Dict:
        """Assess quality of binding data."""
        quality_scores = {}
        for i, data in enumerate(self.binding_data):
            score = 100  # Start with perfect score

            # Deduct points for missing important data
            if data.affinity_error is None:
                score -= 10
            if data.affinity_temperature is None:
                score -= 5
            if data.affinity_ph is None:
                score -= 5

            # Deduct for less reliable sources
            if data.source == "BindingDB":
                score -= 0
            elif "patent" in data.source.lower():
                score -= 15
            elif "predicted" in data.source.lower():
                score -= 25

            # Deduct for missing references
            if not any([data.doi, data.pmid, data.patent_number]):
                score -= 20

            # Deduct for missing experimental details
            if not data.assay_description:
                score -= 10

            # Deduct for approximated values
            if data.affinity_modifier in ["~", ">", "<"]:
                score -= 15

            quality_scores[f"binding_data_{i}"] = max(0, score)  # Ensure non-negative

        self.quality_scores = quality_scores
        return quality_scores

    def classify_activity(self, threshold_high: float = 100, threshold_low: float = 1000) -> Dict[str, List[BindingData]]:
        """
        Classify binding data into activity categories.

        Args:
            threshold_high: Threshold for high activity in nM (default 100)
            threshold_low: Threshold for low activity in nM (default 1000)

        Returns:
            Dictionary with activity classifications
        """
        classifications = {"high_activity": [], "moderate_activity": [], "low_activity": []}

        for data in self.binding_data:
            value = data.affinity_value
            if data.affinity_unit == "μM":
                value *= 1000
            elif data.affinity_unit == "mM":
                value *= 1000000
            elif data.affinity_unit == "pM":
                value /= 1000

            if value <= threshold_high:
                classifications["high_activity"].append(data)
            elif value <= threshold_low:
                classifications["moderate_activity"].append(data)
            else:
                classifications["low_activity"].append(data)

        return classifications

    def get_summary(self) -> Dict:
        """Get comprehensive summary of binding data analysis."""
        return {
            "statistics": self.calculate_statistics(),
            "quality_scores": self.assess_data_quality(),
            "activity_classification": self.classify_activity(),
            "data_count": len(self.binding_data),
            "unique_targets": len(set(data.target_name for data in self.binding_data)),
            "unique_sources": len(set(data.source for data in self.binding_data)),
        }

    def get_confidence_interval(self, confidence: float = 0.95) -> Tuple[float, float]:
        """
        Calculate confidence interval for binding affinities.

        Args:
            confidence: Confidence level (default 0.95 for 95% CI)

        Returns:
            Tuple of (lower bound, upper bound)
        """
        values = self.normalize_affinity_values()
        if not values:
            return (0.0, 0.0)

        mean = np.mean(values)
        std_err = stats.sem(values)
        ci = stats.t.interval(confidence, len(values) - 1, mean, std_err)
        return (float(ci[0]), float(ci[1]))

    def get_outliers(self, threshold: float = 2.0) -> List[BindingData]:
        """
        Identify outliers in binding data.

        Args:
            threshold: Z-score threshold for outlier detection

        Returns:
            List of outlier BindingData objects
        """
        values = self.normalize_affinity_values()
        if not values:
            return []

        z_scores = stats.zscore(values)
        outlier_indices = np.where(np.abs(z_scores) > threshold)[0]
        return [self.binding_data[i] for i in outlier_indices]
