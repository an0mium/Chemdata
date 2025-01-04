"""Validation manager for pipeline.

This module provides the ValidationManager class that handles:
1. Structure validation
2. Data validation
3. Format validation
4. Consistency checks
5. Quality assessment

The manager supports multiple validation types:
- Chemical structure validation
- Property validation
- Data completeness checks
- Format consistency
- Quality metrics
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Set
from dataclasses import dataclass, field
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw

from ..models import CompoundData
from ..processors.structure import (
    StructureValidator,
    PharmacophoreDetector,
    SimilaritySearcher,
)


@dataclass
class ValidationConfig:
    """Validation configuration."""
    
    # Structure validation
    validate_structure: bool = True
    require_valid_structure: bool = True
    check_stereochemistry: bool = True
    check_aromaticity: bool = True
    
    # Property validation
    validate_properties: bool = True
    require_properties: bool = False
    property_ranges: Dict[str, Dict[str, float]] = field(default_factory=lambda: {
        "molecular_weight": {"min": 0, "max": 2000},
        "logp": {"min": -10, "max": 10},
        "tpsa": {"min": 0, "max": 500},
        "hbd": {"min": 0, "max": 20},
        "hba": {"min": 0, "max": 20},
        "rotatable_bonds": {"min": 0, "max": 50},
    })
    
    # Data validation
    validate_data: bool = True
    required_fields: Set[str] = field(default_factory=lambda: {
        "name",
        "smiles",
    })
    optional_fields: Set[str] = field(default_factory=lambda: {
        "cas_number",
        "inchi",
        "inchi_key",
        "molecular_weight",
        "logp",
    })
    
    # Format validation
    validate_format: bool = True
    cas_pattern: str = r'^\d{1,7}-\d{2}-\d$'
    name_pattern: str = r'^[A-Za-z0-9\-\(\)\[\] ]+$'
    
    # Quality thresholds
    min_completeness: float = 0.5
    min_consistency: float = 0.8
    min_quality: float = 0.7


@dataclass
class ValidationStats:
    """Validation statistics."""
    
    # Validation counts
    total_validations: int = 0
    successful_validations: int = 0
    failed_validations: int = 0
    
    # Issue counts
    structure_issues: int = 0
    property_issues: int = 0
    data_issues: int = 0
    format_issues: int = 0
    quality_issues: int = 0
    
    # Quality metrics
    completeness_scores: List[float] = field(default_factory=list)
    consistency_scores: List[float] = field(default_factory=list)
    quality_scores: List[float] = field(default_factory=list)
    
    # Error tracking
    errors: List[Dict[str, Any]] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "validations": {
                "total": self.total_validations,
                "successful": self.successful_validations,
                "failed": self.failed_validations,
                "success_rate": self._get_success_rate(),
            },
            "issues": {
                "structure": self.structure_issues,
                "property": self.property_issues,
                "data": self.data_issues,
                "format": self.format_issues,
                "quality": self.quality_issues,
            },
            "quality": {
                "completeness": self._get_average(self.completeness_scores),
                "consistency": self._get_average(self.consistency_scores),
                "quality": self._get_average(self.quality_scores),
            },
            "errors": self.errors,
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get validation success rate."""
        if not self.total_validations:
            return None
        return self.successful_validations / self.total_validations
    
    def _get_average(self, scores: List[float]) -> Optional[float]:
        """Get average score."""
        if not scores:
            return None
        return sum(scores) / len(scores)


class ValidationManager:
    """Manager for data validation."""

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[ValidationConfig] = None,
    ):
        """Initialize validation manager.
        
        Args:
            cache_dir: Optional directory for caching
            config: Optional validation configuration
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.config = config or ValidationConfig()
        
        # Initialize validators
        self._init_validators()
        
        # Initialize stats
        self.stats = ValidationStats()

    def _init_validators(self) -> None:
        """Initialize validators."""
        try:
            # Structure validators
            self.structure_validator = StructureValidator()
            self.pharmacophore_detector = PharmacophoreDetector()
            self.similarity_searcher = SimilaritySearcher()
            
            self.logger.info("Successfully initialized validators")
            
        except Exception as e:
            self.logger.error(f"Failed to initialize validators: {str(e)}")
            raise

    def validate_compound(
        self,
        compound: CompoundData,
    ) -> CompoundData:
        """Validate compound data.
        
        Args:
            compound: CompoundData instance to validate
            
        Returns:
            Validated CompoundData instance
            
        Raises:
            ValidationError: If validation fails and require_valid is True
        """
        try:
            self.stats.total_validations += 1
            issues = []
            
            # Structure validation
            if self.config.validate_structure:
                structure_issues = self._validate_structure(compound)
                if structure_issues:
                    issues.extend(structure_issues)
                    self.stats.structure_issues += len(structure_issues)
            
            # Property validation
            if self.config.validate_properties:
                property_issues = self._validate_properties(compound)
                if property_issues:
                    issues.extend(property_issues)
                    self.stats.property_issues += len(property_issues)
            
            # Data validation
            if self.config.validate_data:
                data_issues = self._validate_data(compound)
                if data_issues:
                    issues.extend(data_issues)
                    self.stats.data_issues += len(data_issues)
            
            # Format validation
            if self.config.validate_format:
                format_issues = self._validate_format(compound)
                if format_issues:
                    issues.extend(format_issues)
                    self.stats.format_issues += len(format_issues)
            
            # Quality assessment
            quality_issues = self._assess_quality(compound)
            if quality_issues:
                issues.extend(quality_issues)
                self.stats.quality_issues += len(quality_issues)
            
            # Update stats
            if issues:
                self.stats.failed_validations += 1
                self.stats.errors.append({
                    "type": "validation_error",
                    "compound": compound.name,
                    "issues": issues,
                    "timestamp": datetime.now().isoformat(),
                })
                if self.config.require_valid_structure and any(
                    i["type"] == "structure" for i in issues
                ):
                    raise ValidationError(
                        f"Structure validation failed for {compound.name}: {issues}"
                    )
            else:
                self.stats.successful_validations += 1
            
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to validate compound {compound.name}: {str(e)}"
            )
            self.stats.failed_validations += 1
            self.stats.errors.append({
                "type": "validation_error",
                "compound": compound.name,
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            })
            raise

    def _validate_structure(
        self,
        compound: CompoundData,
    ) -> List[Dict[str, Any]]:
        """Validate chemical structure."""
        issues = []
        
        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                issues.append({
                    "type": "structure",
                    "severity": "error",
                    "message": "Invalid SMILES string",
                })
                return issues
            
            # Check stereochemistry
            if self.config.check_stereochemistry:
                if not self.structure_validator.check_stereochemistry(mol):
                    issues.append({
                        "type": "structure",
                        "severity": "warning",
                        "message": "Incomplete stereochemistry",
                    })
            
            # Check aromaticity
            if self.config.check_aromaticity:
                if not self.structure_validator.check_aromaticity(mol):
                    issues.append({
                        "type": "structure",
                        "severity": "warning",
                        "message": "Invalid aromaticity",
                    })
            
            return issues
            
        except Exception as e:
            self.logger.error(f"Structure validation error: {str(e)}")
            issues.append({
                "type": "structure",
                "severity": "error",
                "message": str(e),
            })
            return issues

    def _validate_properties(
        self,
        compound: CompoundData,
    ) -> List[Dict[str, Any]]:
        """Validate chemical properties."""
        issues = []
        
        try:
            for prop, ranges in self.config.property_ranges.items():
                value = getattr(compound, prop, None)
                if value is None:
                    if self.config.require_properties:
                        issues.append({
                            "type": "property",
                            "severity": "error",
                            "message": f"Missing property: {prop}",
                        })
                else:
                    if not ranges["min"] <= value <= ranges["max"]:
                        issues.append({
                            "type": "property",
                            "severity": "warning",
                            "message": (
                                f"Property {prop} ({value}) outside "
                                f"range [{ranges['min']}, {ranges['max']}]"
                            ),
                        })
            
            return issues
            
        except Exception as e:
            self.logger.error(f"Property validation error: {str(e)}")
            issues.append({
                "type": "property",
                "severity": "error",
                "message": str(e),
            })
            return issues

    def _validate_data(
        self,
        compound: CompoundData,
    ) -> List[Dict[str, Any]]:
        """Validate compound data."""
        issues = []
        
        try:
            # Check required fields
            for field in self.config.required_fields:
                if not hasattr(compound, field) or not getattr(compound, field):
                    issues.append({
                        "type": "data",
                        "severity": "error",
                        "message": f"Missing required field: {field}",
                    })
            
            # Check optional fields
            for field in self.config.optional_fields:
                if not hasattr(compound, field):
                    issues.append({
                        "type": "data",
                        "severity": "warning",
                        "message": f"Missing optional field: {field}",
                    })
            
            return issues
            
        except Exception as e:
            self.logger.error(f"Data validation error: {str(e)}")
            issues.append({
                "type": "data",
                "severity": "error",
                "message": str(e),
            })
            return issues

    def _validate_format(
        self,
        compound: CompoundData,
    ) -> List[Dict[str, Any]]:
        """Validate data formats."""
        issues = []
        
        try:
            # Validate CAS number
            if compound.cas_number:
                if not self.structure_validator.validate_cas(
                    compound.cas_number,
                    pattern=self.config.cas_pattern,
                ):
                    issues.append({
                        "type": "format",
                        "severity": "error",
                        "message": "Invalid CAS number format",
                    })
            
            # Validate name
            if compound.name:
                if not self.structure_validator.validate_name(
                    compound.name,
                    pattern=self.config.name_pattern,
                ):
                    issues.append({
                        "type": "format",
                        "severity": "warning",
                        "message": "Invalid name format",
                    })
            
            return issues
            
        except Exception as e:
            self.logger.error(f"Format validation error: {str(e)}")
            issues.append({
                "type": "format",
                "severity": "error",
                "message": str(e),
            })
            return issues

    def _assess_quality(
        self,
        compound: CompoundData,
    ) -> List[Dict[str, Any]]:
        """Assess data quality."""
        issues = []
        
        try:
            # Calculate completeness
            completeness = self._calculate_completeness(compound)
            self.stats.completeness_scores.append(completeness)
            if completeness < self.config.min_completeness:
                issues.append({
                    "type": "quality",
                    "severity": "warning",
                    "message": f"Low completeness score: {completeness:.2f}",
                })
            
            # Calculate consistency
            consistency = self._calculate_consistency(compound)
            self.stats.consistency_scores.append(consistency)
            if consistency < self.config.min_consistency:
                issues.append({
                    "type": "quality",
                    "severity": "warning",
                    "message": f"Low consistency score: {consistency:.2f}",
                })
            
            # Calculate quality
            quality = self._calculate_quality(compound)
            self.stats.quality_scores.append(quality)
            if quality < self.config.min_quality:
                issues.append({
                    "type": "quality",
                    "severity": "warning",
                    "message": f"Low quality score: {quality:.2f}",
                })
            
            return issues
            
        except Exception as e:
            self.logger.error(f"Quality assessment error: {str(e)}")
            issues.append({
                "type": "quality",
                "severity": "error",
                "message": str(e),
            })
            return issues

    def _calculate_completeness(self, compound: CompoundData) -> float:
        """Calculate data completeness score."""
        total_fields = len(self.config.required_fields) + len(self.config.optional_fields)
        present_fields = sum(
            1 for f in self.config.required_fields | self.config.optional_fields
            if hasattr(compound, f) and getattr(compound, f)
        )
        return present_fields / total_fields

    def _calculate_consistency(self, compound: CompoundData) -> float:
        """Calculate data consistency score."""
        # TODO: Implement consistency checks
        return 1.0

    def _calculate_quality(self, compound: CompoundData) -> float:
        """Calculate overall quality score."""
        completeness = self._calculate_completeness(compound)
        consistency = self._calculate_consistency(compound)
        return (completeness + consistency) / 2


class ValidationError(Exception):
    """Raised when validation fails."""
    pass
