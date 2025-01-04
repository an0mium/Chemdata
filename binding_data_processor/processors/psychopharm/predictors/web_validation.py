"""Web data validation for BBB permeability prediction.

This module provides functionality to:
1. Validate web-sourced data
2. Standardize data formats
3. Cross-reference data sources
4. Detect and resolve conflicts
5. Score data reliability
"""

import logging
from typing import Dict, List, Optional, Any
import pandas as pd
from dataclasses import dataclass


@dataclass
class WebDataValidationResult:
    """Result of web data validation."""
    
    is_valid: bool
    confidence: float
    issues: List[str]
    source_scores: Dict[str, float]
    cross_references: Dict[str, List[str]]
    supporting_data: Dict[str, Any]


class WebDataValidator:
    """Validator for web-sourced compound data."""

    def __init__(
        self,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize web data validator."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.cache_dir = cache_dir

        # Source reliability scores (0-1)
        self.source_weights = {
            "chembl": 1.0,      # Highly reliable
            "pubchem": 0.9,     # Very reliable
            "swiss": 0.9,       # Very reliable
            "community": 0.6,   # Moderately reliable
            "social": 0.3,      # Less reliable
            "web_search": 0.4,  # Less reliable
        }

        # Required fields by source
        self.required_fields = {
            "chembl": [
                "compound_id",
                "activity_data",
                "target_data",
            ],
            "pubchem": [
                "cid",
                "properties",
                "bioactivity",
            ],
            "swiss": [
                "predictions",
                "properties",
                "targets",
            ],
            "community": [
                "reports",
                "effects",
                "safety",
            ],
            "social": [
                "mentions",
                "discussions",
                "reports",
            ],
            "web_search": [
                "references",
                "mentions",
                "context",
            ],
        }

        # Value range validations
        self.value_ranges = {
            "activity_value": (0, 1e6),
            "confidence": (0, 1),
            "probability": (0, 1),
            "score": (0, 100),
        }

    def validate_web_data(
        self, web_data: Dict[str, Any], compound_name: str
    ) -> WebDataValidationResult:
        """Validate web-sourced data for a compound."""
        self.logger.debug(f"Validating web data for {compound_name}")
        
        issues = []
        source_scores = {}
        cross_references = {}
        supporting_data = {}
        
        # Validate each source
        for source, data in web_data.items():
            if source not in self.source_weights:
                issues.append(f"Unknown source: {source}")
                continue
            
            # Check required fields
            source_issues = self._validate_source_data(source, data)
            issues.extend(source_issues)
            
            # Calculate source score
            base_score = self.source_weights[source]
            completeness = self._calculate_completeness(source, data)
            quality = self._calculate_quality(source, data)
            source_scores[source] = base_score * completeness * quality
            
            # Find cross-references
            refs = self._find_cross_references(source, data, web_data)
            if refs:
                cross_references[source] = refs
            
            # Collect supporting data
            if source_scores[source] > 0.5:  # Only include reliable sources
                supporting_data[source] = self._extract_supporting_data(data)
        
        # Calculate overall confidence
        confidence = self._calculate_overall_confidence(source_scores)
        
        # Determine overall validity
        is_valid = (
            confidence > 0.5 and
            len(issues) < 5 and
            any(score > 0.7 for score in source_scores.values())
        )
        
        return WebDataValidationResult(
            is_valid=is_valid,
            confidence=confidence,
            issues=issues,
            source_scores=source_scores,
            cross_references=cross_references,
            supporting_data=supporting_data,
        )

    def _validate_source_data(
        self, source: str, data: Dict[str, Any]
    ) -> List[str]:
        """Validate data from a specific source."""
        issues = []
        
        # Check required fields
        for field in self.required_fields[source]:
            if field not in data:
                issues.append(f"Missing required field '{field}' for {source}")
                continue
            
            # Check for empty/null values
            if data[field] is None or data[field] == "":
                issues.append(f"Empty required field '{field}' for {source}")
            
            # Check value ranges
            if field in self.value_ranges:
                min_val, max_val = self.value_ranges[field]
                try:
                    val = float(data[field])
                    if not min_val <= val <= max_val:
                        issues.append(
                            f"Value {val} for '{field}' in {source} "
                            f"outside range [{min_val}, {max_val}]"
                        )
                except (ValueError, TypeError):
                    issues.append(
                        f"Invalid numeric value for '{field}' in {source}"
                    )
        
        return issues

    def _calculate_completeness(
        self, source: str, data: Dict[str, Any]
    ) -> float:
        """Calculate data completeness score."""
        required = set(self.required_fields[source])
        present = set(data.keys())
        
        # Calculate basic completeness
        completeness = len(present & required) / len(required)
        
        # Bonus for additional useful fields
        bonus = len(present - required) * 0.1
        
        return min(1.0, completeness + bonus)

    def _calculate_quality(
        self, source: str, data: Dict[str, Any]
    ) -> float:
        """Calculate data quality score."""
        quality_scores = []
        
        # Check data format quality
        format_score = self._check_format_quality(data)
        quality_scores.append(format_score)
        
        # Check data consistency
        consistency_score = self._check_data_consistency(data)
        quality_scores.append(consistency_score)
        
        # Check data recency
        recency_score = self._check_data_recency(data)
        quality_scores.append(recency_score)
        
        # Calculate weighted average
        weights = [0.4, 0.4, 0.2]  # Format, consistency, recency
        return sum(s * w for s, w in zip(quality_scores, weights))

    def _check_format_quality(self, data: Dict[str, Any]) -> float:
        """Check quality of data formats."""
        format_issues = 0
        checked_fields = 0
        
        for field, value in data.items():
            checked_fields += 1
            
            # Check JSON structure
            if isinstance(value, dict):
                if not all(isinstance(k, str) for k in value.keys()):
                    format_issues += 1
            
            # Check list structure
            elif isinstance(value, list):
                if not all(isinstance(x, (str, dict)) for x in value):
                    format_issues += 1
            
            # Check string values
            elif isinstance(value, str):
                if len(value.strip()) == 0:
                    format_issues += 1
            
            # Check numeric values
            elif isinstance(value, (int, float)):
                if value < 0 and field not in ["delta", "change"]:
                    format_issues += 1
        
        return max(0.0, 1.0 - (format_issues / checked_fields))

    def _check_data_consistency(self, data: Dict[str, Any]) -> float:
        """Check internal consistency of data."""
        consistency_issues = 0
        checked_pairs = 0
        
        # Check numeric relationships
        numeric_fields = {
            k: v for k, v in data.items()
            if isinstance(v, (int, float))
        }
        for field1, value1 in numeric_fields.items():
            for field2, value2 in numeric_fields.items():
                if field1 < field2:  # Avoid checking same pair twice
                    checked_pairs += 1
                    # Check for impossible relationships
                    if "probability" in field1 and "probability" in field2:
                        if value1 + value2 > 1.1:  # Allow small error
                            consistency_issues += 1
                    elif "confidence" in field1 and "confidence" in field2:
                        if abs(value1 - value2) > 0.5:  # Large confidence diff
                            consistency_issues += 1
        
        # Check categorical consistency
        if "category" in data and "subcategory" in data:
            checked_pairs += 1
            if not str(data["subcategory"]).startswith(str(data["category"])):
                consistency_issues += 1
        
        if checked_pairs == 0:
            return 1.0
        
        return max(0.0, 1.0 - (consistency_issues / checked_pairs))

    def _check_data_recency(self, data: Dict[str, Any]) -> float:
        """Check recency of data."""
        if "timestamp" not in data:
            return 0.5  # Neutral score if no timestamp
        
        try:
            timestamp = pd.to_datetime(data["timestamp"])
            now = pd.Timestamp.now()
            age_days = (now - timestamp).days
            
            # Score decreases with age
            if age_days < 30:
                return 1.0
            elif age_days < 90:
                return 0.8
            elif age_days < 180:
                return 0.6
            elif age_days < 365:
                return 0.4
            else:
                return 0.2
                
        except (ValueError, TypeError):
            return 0.5  # Neutral score if invalid timestamp

    def _find_cross_references(
        self,
        source: str,
        data: Dict[str, Any],
        all_data: Dict[str, Dict[str, Any]],
    ) -> List[str]:
        """Find cross-references between sources."""
        refs = []
        
        # Look for matching IDs
        for other_source, other_data in all_data.items():
            if other_source == source:
                continue
            
            # Check for matching compound IDs
            if (
                "compound_id" in data and
                "compound_id" in other_data and
                data["compound_id"] == other_data["compound_id"]
            ):
                refs.append(other_source)
            
            # Check for matching CAS numbers
            elif (
                "cas" in data and
                "cas" in other_data and
                data["cas"] == other_data["cas"]
            ):
                refs.append(other_source)
            
            # Check for matching SMILES
            elif (
                "smiles" in data and
                "smiles" in other_data and
                data["smiles"] == other_data["smiles"]
            ):
                refs.append(other_source)
        
        return refs

    def _extract_supporting_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Extract relevant supporting data."""
        supporting = {}
        
        # Extract activity data
        if "activity_data" in data:
            supporting["activity"] = data["activity_data"]
        
        # Extract target data
        if "target_data" in data:
            supporting["targets"] = data["target_data"]
        
        # Extract property data
        if "properties" in data:
            supporting["properties"] = data["properties"]
        
        # Extract safety data
        if "safety" in data:
            supporting["safety"] = data["safety"]
        
        # Extract reference data
        if "references" in data:
            supporting["references"] = data["references"]
        
        return supporting

    def _calculate_overall_confidence(
        self, source_scores: Dict[str, float]
    ) -> float:
        """Calculate overall confidence score."""
        if not source_scores:
            return 0.0
        
        # Weight scores by source reliability
        weighted_scores = [
            score * self.source_weights[source]
            for source, score in source_scores.items()
        ]
        
        # Calculate weighted average
        total_weight = sum(self.source_weights[s] for s in source_scores)
        if total_weight == 0:
            return 0.0
            
        return sum(weighted_scores) / total_weight
