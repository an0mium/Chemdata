"""Validation functionality for data export.

This module provides functionality to:
1. Validate export fields and formats
2. Validate compound data before export
3. Ensure data consistency
4. Type check exported data
"""

import logging
from typing import Dict, List, Optional, Any, Set
from dataclasses import dataclass

from ....models.validation import ValidationResult, Validator
from .data_enrichment import EnrichedData


@dataclass
class ExportValidationResult(ValidationResult):
    """Result of export data validation."""
    
    valid_compounds: List[EnrichedData]
    invalid_compounds: List[EnrichedData]
    missing_fields: Set[str]
    invalid_fields: Dict[str, List[str]]
    field_coverage: Dict[str, float]


class ExportValidator(Validator):
    """Validator for export data."""

    def __init__(
        self,
        log_level: int = logging.INFO,
    ):
        """Initialize export validator."""
        super().__init__(log_level)

    def validate_export_data(
        self,
        compounds: List[EnrichedData],
        fields: Set[str],
        required_fields: Optional[Set[str]] = None,
    ) -> ExportValidationResult:
        """Validate compounds for export."""
        self.logger.debug(f"Validating {len(compounds)} compounds for export")
        
        valid_compounds = []
        invalid_compounds = []
        missing_fields = set()
        invalid_fields = {}
        field_coverage = {}
        
        # Check required fields
        if required_fields:
            missing_fields = required_fields - fields
        
        # Validate each compound
        for compound in compounds:
            field_issues = self._validate_compound_fields(
                compound, fields
            )
            
            if field_issues:
                invalid_compounds.append(compound)
                invalid_fields[compound.compound.name] = field_issues
            else:
                valid_compounds.append(compound)
        
        # Calculate field coverage
        field_coverage = self._calculate_field_coverage(
            compounds, fields
        )
        
        return ExportValidationResult(
            is_valid=len(invalid_compounds) == 0 and len(missing_fields) == 0,
            valid_compounds=valid_compounds,
            invalid_compounds=invalid_compounds,
            missing_fields=missing_fields,
            invalid_fields=invalid_fields,
            field_coverage=field_coverage,
            issues=self._format_validation_issues(
                missing_fields,
                invalid_compounds,
                field_coverage,
            ),
        )

    def _validate_compound_fields(
        self,
        compound: EnrichedData,
        fields: Set[str],
    ) -> List[str]:
        """Validate compound fields."""
        issues = []
        
        # Validate basic fields
        if "name" in fields and not self._is_valid_string(
            compound.compound.name
        ):
            issues.append("Invalid name")
        
        if "smiles" in fields and not self._is_valid_string(
            compound.compound.smiles
        ):
            issues.append("Invalid SMILES")
        
        if "cas" in fields and not self._is_valid_cas(
            compound.standardized.cas
        ):
            issues.append("Invalid CAS")
        
        # Validate numeric fields
        if "molecular_weight" in fields and not self._is_valid_numeric(
            compound.properties.get("molecular_weight")
        ):
            issues.append("Invalid molecular weight")
        
        if "logp" in fields and not self._is_valid_numeric(
            compound.properties.get("logp")
        ):
            issues.append("Invalid LogP")
        
        # Validate prediction fields
        if "bbb_predictions" in fields and not self._is_valid_predictions(
            compound.bbb_predictions
        ):
            issues.append("Invalid BBB predictions")
        
        if "activity_predictions" in fields and not self._is_valid_predictions(
            compound.activity_predictions
        ):
            issues.append("Invalid activity predictions")
        
        if "toxicity_predictions" in fields and not self._is_valid_predictions(
            compound.toxicity_predictions
        ):
            issues.append("Invalid toxicity predictions")
        
        if "abuse_predictions" in fields and not self._is_valid_predictions(
            compound.abuse_predictions
        ):
            issues.append("Invalid abuse predictions")
        
        return issues

    def _is_valid_string(self, value: Any) -> bool:
        """Check if value is valid string."""
        return isinstance(value, str) and len(value.strip()) > 0

    def _is_valid_cas(self, value: Any) -> bool:
        """Check if value is valid CAS number."""
        if not self._is_valid_string(value):
            return False
        
        # Basic CAS format validation
        parts = value.split("-")
        if len(parts) != 3:
            return False
        
        try:
            # Check each part is numeric
            for part in parts[:-1]:
                int(part)
            
            # Validate check digit
            digits = "".join(parts[:-1])
            check = int(parts[-1])
            
            total = sum(
                int(d) * (i + 1)
                for i, d in enumerate(reversed(digits))
            )
            return (total % 10) == check
            
        except ValueError:
            return False

    def _is_valid_numeric(self, value: Any) -> bool:
        """Check if value is valid number."""
        if value is None:
            return True
        
        try:
            float(value)
            return True
        except (ValueError, TypeError):
            return False

    def _is_valid_predictions(self, predictions: Any) -> bool:
        """Check if predictions are valid."""
        if not isinstance(predictions, dict):
            return False
        
        # Check required prediction fields
        required = {"value", "confidence"}
        if not all(
            isinstance(pred, dict) and required.issubset(pred.keys())
            for pred in predictions.values()
            if isinstance(pred, dict)
        ):
            return False
        
        # Validate confidence values
        if not all(
            isinstance(pred.get("confidence"), (int, float)) and
            0 <= pred["confidence"] <= 1
            for pred in predictions.values()
            if isinstance(pred, dict)
        ):
            return False
        
        return True

    def _format_validation_issues(
        self,
        missing_fields: Set[str],
        invalid_compounds: List[EnrichedData],
        field_coverage: Dict[str, float],
    ) -> List[str]:
        """Format validation issues into messages."""
        issues = []
        
        if missing_fields:
            issues.append(f"Missing required fields: {missing_fields}")
        
        if invalid_compounds:
            issues.append(f"Invalid compounds: {len(invalid_compounds)}")
        
        if field_coverage:
            low_coverage = [
                f for f, c in field_coverage.items() if c < 0.5
            ]
            if low_coverage:
                issues.append(f"Low field coverage: {low_coverage}")
        
        return issues

    def _calculate_field_coverage(
        self,
        compounds: List[EnrichedData],
        fields: Set[str],
    ) -> Dict[str, float]:
        """Calculate coverage for each field."""
        if not compounds:
            return {}
        
        coverage = {}
        total = len(compounds)
        
        for field in fields:
            valid = self._count_valid_values(compounds, field)
            coverage[field] = valid / total
        
        return coverage

    def _count_valid_values(
        self,
        compounds: List[EnrichedData],
        field: str,
    ) -> int:
        """Count compounds with valid values for a field."""
        valid = 0
        
        for compound in compounds:
            value = self._get_field_value(compound, field)
            if self._is_valid_value(value, field):
                valid += 1
        
        return valid

    def _get_field_value(
        self,
        compound: EnrichedData,
        field: str,
    ) -> Any:
        """Get value for a field from compound."""
        if field in {"name", "smiles"}:
            return getattr(compound.compound, field, None)
        elif field == "cas":
            return compound.standardized.cas
        elif field in {"molecular_weight", "logp"}:
            return compound.properties.get(field)
        elif field.endswith("_predictions"):
            predictor = field.replace("_predictions", "")
            return getattr(compound, f"{predictor}_predictions", None)
        return None

    def _is_valid_value(
        self,
        value: Any,
        field: str,
    ) -> bool:
        """Check if value is valid for field type."""
        if field in {"name", "smiles", "cas"}:
            return self._is_valid_string(value)
        elif field in {"molecular_weight", "logp"}:
            return self._is_valid_numeric(value)
        elif field.endswith("_predictions"):
            return self._is_valid_predictions(value)
        return False
