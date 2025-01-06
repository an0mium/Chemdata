"""Data validation for web enrichment data.

This module provides data validation and cleaning for:
- Text data (HTML, markdown, etc.)
- Numeric data (ranges, units, etc.)
- Temporal data (timestamps, durations, etc.)
"""

import logging
import re
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass
from enum import Enum

import bs4
import markdown
import dateutil.parser


class DataType(Enum):
    """Supported data types."""
    TEXT = "text"
    NUMBER = "number"
    TIMESTAMP = "timestamp"
    DURATION = "duration"
    IDENTIFIER = "identifier"


class ValidationRule(Enum):
    """Data validation rules."""
    REQUIRED = "required"
    NON_EMPTY = "non_empty"
    RANGE = "range"
    PATTERN = "pattern"
    UNIQUE = "unique"
    FORMAT = "format"


@dataclass
class ValidationConfig:
    """Configuration for data validation."""
    rules: Dict[str, List[Tuple[ValidationRule, Any]]]
    strip_html: bool = True
    normalize_whitespace: bool = True
    convert_markdown: bool = True
    timezone: str = "UTC"


@dataclass
class ValidationIssue:
    """Data validation issue."""
    field: str
    rule: ValidationRule
    message: str
    severity: str


class DataValidator:
    """Data validator and cleaner."""

    def __init__(
        self,
        config: ValidationConfig,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize data validator.
        
        Args:
            config: Validation configuration
            logger: Optional logger instance
        """
        self.config = config
        self.logger = logger or logging.getLogger(self.__class__.__name__)

    def validate(
        self,
        data: Dict[str, Any],
    ) -> Tuple[Dict[str, Any], List[ValidationIssue]]:
        """Validate and clean data.
        
        Args:
            data: Data to validate
            
        Returns:
            Tuple of (cleaned data, validation issues)
        """
        cleaned = {}
        issues = []

        for field, value in data.items():
            if field in self.config.rules:
                cleaned_value = self._clean_value(field, value)
                validation_issues = self._validate_value(
                    field,
                    cleaned_value,
                    self.config.rules[field],
                )
                cleaned[field] = cleaned_value
                issues.extend(validation_issues)
            else:
                cleaned[field] = value

        return cleaned, issues

    def _clean_value(self, field: str, value: Any) -> Any:
        """Clean value based on type and configuration.
        
        Args:
            field: Field name
            value: Value to clean
            
        Returns:
            Cleaned value
        """
        if value is None:
            return None

        if isinstance(value, str):
            # Strip HTML
            if self.config.strip_html:
                value = self._strip_html(value)

            # Convert markdown
            if self.config.convert_markdown:
                value = self._convert_markdown(value)

            # Normalize whitespace
            if self.config.normalize_whitespace:
                value = self._normalize_whitespace(value)

        return value

    def _validate_value(
        self,
        field: str,
        value: Any,
        rules: List[Tuple[ValidationRule, Any]],
    ) -> List[ValidationIssue]:
        """Validate value against rules.
        
        Args:
            field: Field name
            value: Value to validate
            rules: List of (rule, rule_value) tuples
            
        Returns:
            List of validation issues
        """
        issues = []
        validation_map = {
            ValidationRule.REQUIRED: self._validate_required,
            ValidationRule.NON_EMPTY: self._validate_non_empty,
            ValidationRule.RANGE: self._validate_range,
            ValidationRule.PATTERN: self._validate_pattern,
            ValidationRule.UNIQUE: self._validate_unique,
            ValidationRule.FORMAT: self._validate_format_rule,
        }

        for rule, rule_value in rules:
            if rule in validation_map:
                issue = validation_map[rule](field, value, rule_value)
                if issue:
                    issues.append(issue)

        return issues

    def _validate_required(
        self,
        field: str,
        value: Any,
        _: Any,
    ) -> Optional[ValidationIssue]:
        """Validate required field."""
        if value is None:
            return ValidationIssue(
                field=field,
                rule=ValidationRule.REQUIRED,
                message="Field is required",
                severity="error",
            )
        return None

    def _validate_non_empty(
        self,
        field: str,
        value: Any,
        _: Any,
    ) -> Optional[ValidationIssue]:
        """Validate non-empty field."""
        if not value:
            return ValidationIssue(
                field=field,
                rule=ValidationRule.NON_EMPTY,
                message="Field cannot be empty",
                severity="error",
            )
        return None

    def _validate_range(
        self,
        field: str,
        value: Any,
        range_value: Tuple[float, float],
    ) -> Optional[ValidationIssue]:
        """Validate value range."""
        if value is not None:
            min_val, max_val = range_value
            if not min_val <= value <= max_val:
                return ValidationIssue(
                    field=field,
                    rule=ValidationRule.RANGE,
                    message=f"Value must be between {min_val} and {max_val}",
                    severity="error",
                )
        return None

    def _validate_pattern(
        self,
        field: str,
        value: Any,
        pattern: str,
    ) -> Optional[ValidationIssue]:
        """Validate value pattern."""
        if value is not None and not re.match(pattern, str(value)):
            return ValidationIssue(
                field=field,
                rule=ValidationRule.PATTERN,
                message=f"Value must match pattern: {pattern}",
                severity="error",
            )
        return None

    def _validate_unique(
        self,
        field: str,
        value: Any,
        seen_values: set,
    ) -> Optional[ValidationIssue]:
        """Validate unique value."""
        if value is not None:
            if value in seen_values:
                return ValidationIssue(
                    field=field,
                    rule=ValidationRule.UNIQUE,
                    message="Value must be unique",
                    severity="error",
                )
            seen_values.add(value)
        return None

    def _validate_format_rule(
        self,
        field: str,
        value: Any,
        format_type: str,
    ) -> Optional[ValidationIssue]:
        """Validate value format."""
        if value is not None and not self._validate_format(value, format_type):
            return ValidationIssue(
                field=field,
                rule=ValidationRule.FORMAT,
                message=f"Invalid {format_type} format",
                severity="error",
            )
        return None

    def _strip_html(self, text: str) -> str:
        """Strip HTML tags from text.
        
        Args:
            text: Text to clean
            
        Returns:
            Text with HTML removed
        """
        soup = bs4.BeautifulSoup(text, "html.parser")
        return soup.get_text()

    def _convert_markdown(self, text: str) -> str:
        """Convert markdown to plain text.
        
        Args:
            text: Text to convert
            
        Returns:
            Plain text version
        """
        html = markdown.markdown(text)
        soup = bs4.BeautifulSoup(html, "html.parser")
        return soup.get_text()

    def _normalize_whitespace(self, text: str) -> str:
        """Normalize whitespace in text.
        
        Args:
            text: Text to normalize
            
        Returns:
            Text with normalized whitespace
        """
        return " ".join(text.split())

    def _validate_format(self, value: Any, format_type: str) -> bool:
        """Validate value format.
        
        Args:
            value: Value to validate
            format_type: Expected format type
            
        Returns:
            True if valid, False otherwise
        """
        try:
            if format_type == "timestamp":
                dateutil.parser.parse(str(value))
                return True

            elif format_type == "duration":
                # Format: 1h30m, 45m, 2h, etc.
                pattern = r"^(\d+h)?(\d+m)?$"
                return bool(re.match(pattern, str(value)))

            elif format_type == "identifier":
                # Format: letters, numbers, dashes, underscores
                pattern = r"^[a-zA-Z0-9_-]+$"
                return bool(re.match(pattern, str(value)))

            return False

        except (ValueError, TypeError):
            return False


class DataCleaner:
    """Data cleaning utilities."""

    @staticmethod
    def clean_text(
        text: str,
        strip_html: bool = True,
        convert_markdown: bool = True,
        normalize_whitespace: bool = True,
    ) -> str:
        """Clean text data.
        
        Args:
            text: Text to clean
            strip_html: Whether to strip HTML
            convert_markdown: Whether to convert markdown
            normalize_whitespace: Whether to normalize whitespace
            
        Returns:
            Cleaned text
        """
        if strip_html:
            soup = bs4.BeautifulSoup(text, "html.parser")
            text = soup.get_text()

        if convert_markdown:
            html = markdown.markdown(text)
            soup = bs4.BeautifulSoup(html, "html.parser")
            text = soup.get_text()

        if normalize_whitespace:
            text = " ".join(text.split())

        return text

    @staticmethod
    def clean_number(
        value: Any,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None,
        round_digits: Optional[int] = None,
    ) -> float:
        """Clean numeric data.
        
        Args:
            value: Value to clean
            min_value: Optional minimum value
            max_value: Optional maximum value
            round_digits: Optional number of decimal places
            
        Returns:
            Cleaned number
            
        Raises:
            ValueError: If value is invalid
        """
        if isinstance(value, str):
            value = float(value)
        elif not isinstance(value, (int, float)):
            raise ValueError(f"Invalid numeric value: {value}")

        if min_value is not None and value < min_value:
            raise ValueError(f"Value below minimum: {value} < {min_value}")

        if max_value is not None and value > max_value:
            raise ValueError(f"Value above maximum: {value} > {max_value}")

        if round_digits is not None:
            value = round(float(value), round_digits)

        return float(value)

    @staticmethod
    def clean_timestamp(
        value: Any,
        timezone: str = "UTC",
    ) -> datetime:
        """Clean timestamp data.
        
        Args:
            value: Value to clean
            timezone: Timezone name
            
        Returns:
            Cleaned datetime
            
        Raises:
            ValueError: If value is invalid
        """
        if isinstance(value, datetime):
            return value

        try:
            dt = dateutil.parser.parse(str(value))
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=dateutil.tz.gettz(timezone))
            return dt
        except (ValueError, TypeError) as e:
            raise ValueError(f"Invalid timestamp: {value}") from e

    @staticmethod
    def clean_duration(value: Any) -> str:
        """Clean duration data.
        
        Args:
            value: Value to clean
            
        Returns:
            Cleaned duration string
            
        Raises:
            ValueError: If value is invalid
        """
        if isinstance(value, (int, float)):
            # Convert minutes to format
            minutes = int(value)
            hours = minutes // 60
            minutes = minutes % 60
            if hours > 0:
                return f"{hours}h{minutes}m" if minutes > 0 else f"{hours}h"
            return f"{minutes}m"

        value = str(value).lower().replace(" ", "")
        pattern = r"^(\d+h)?(\d+m)?$"
        if not re.match(pattern, value):
            raise ValueError(f"Invalid duration format: {value}")

        return value
