"""Data validation for web enrichment data.

This module provides comprehensive data validation and cleaning for:
- Text data (HTML, markdown, etc.)
- Numeric data (ranges, units, etc.)
- Temporal data (timestamps, durations, etc.)
- Validation rules and configurations
- Data cleaning utilities
"""

import logging
import re
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple, Union, Callable
from dataclasses import dataclass, field
from enum import Enum

import bs4
import markdown
import dateutil.parser
from dateutil import tz


class DataType(Enum):
    """Supported data types."""

    TEXT = "text"
    NUMBER = "number"
    TIMESTAMP = "timestamp"
    DURATION = "duration"
    IDENTIFIER = "identifier"
    LIST = "list"
    DICT = "dict"


class ValidationLevel(Enum):
    """Validation strictness levels."""

    STRICT = "strict"  # All rules enforced
    NORMAL = "normal"  # Most rules enforced, some warnings
    LENIENT = "lenient"  # Basic validation only


class ValidationRule(Enum):
    """Core validation rules."""

    REQUIRED = "required"
    NON_EMPTY = "non_empty"
    RANGE = "range"
    PATTERN = "pattern"
    UNIQUE = "unique"
    FORMAT = "format"
    TYPE = "type"
    LENGTH = "length"
    ITEMS = "items"


@dataclass
class ValidationConfig:
    """Configuration for validation."""

    level: ValidationLevel = ValidationLevel.NORMAL
    rules: Dict[str, List[Tuple[ValidationRule, Any]]] = field(default_factory=dict)
    custom_validators: List[Callable] = field(default_factory=list)
    error_messages: Dict[str, str] = field(default_factory=dict)
    strip_html: bool = True
    normalize_whitespace: bool = True
    convert_markdown: bool = True
    timezone: str = "UTC"


@dataclass
class FieldRule:
    """Rule for validating a single field."""

    field: str
    required: bool = False
    data_type: Optional[DataType] = None
    min_length: Optional[int] = None
    max_length: Optional[int] = None
    min_value: Optional[float] = None
    max_value: Optional[float] = None
    min_items: Optional[int] = None
    max_items: Optional[int] = None
    pattern: Optional[str] = None
    allowed_values: Optional[List[Any]] = None
    custom_validator: Optional[Callable] = None
    error_message: Optional[str] = None


@dataclass
class ValidationIssue:
    """Data validation issue."""

    field: str
    rule: Union[ValidationRule, str]
    message: str
    severity: str = "error"


@dataclass
class ValidationResult:
    """Result of validation."""

    is_valid: bool
    errors: List[ValidationIssue] = field(default_factory=list)
    warnings: List[ValidationIssue] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


def clean_community_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """Clean and validate community data.

    Args:
        data: Raw community data dictionary

    Returns:
        Cleaned and validated data dictionary
    """
    if not data:
        return {}

    cleaned = {}
    cleaner = DataCleaner()

    # Clean text fields
    text_fields = ["title", "content", "summary", "source", "platform", "author"]
    for field in text_fields:
        if field in data:
            cleaned[field] = cleaner.clean_text(data[field], strip_html=True, convert_markdown=True, normalize_whitespace=True)

    # Clean timestamps
    time_fields = ["timestamp", "created_at", "updated_at", "last_activity"]
    for field in time_fields:
        if field in data:
            try:
                cleaned[field] = cleaner.clean_timestamp(data[field])
            except ValueError:
                cleaned[field] = None

    # Clean numeric fields
    numeric_fields = ["score", "views", "replies", "mentions"]
    for field in numeric_fields:
        if field in data:
            try:
                cleaned[field] = cleaner.clean_number(data[field])
            except ValueError:
                cleaned[field] = 0

    # Clean lists
    list_fields = ["tags", "categories", "references"]
    for field in list_fields:
        if field in data:
            if isinstance(data[field], list):
                cleaned[field] = [str(item).strip() for item in data[field] if item]
            else:
                cleaned[field] = []

    # Clean URLs
    url_fields = ["url", "permalink", "image_url"]
    for field in url_fields:
        if field in data:
            cleaned[field] = str(data[field]).strip()

    # Copy through any other fields
    for key, value in data.items():
        if key not in cleaned:
            cleaned[key] = value

    return cleaned


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
        if not text:
            return ""

        result = text

        if strip_html:
            soup = bs4.BeautifulSoup(text, "html.parser")
            result = soup.get_text()

        if convert_markdown:
            html = markdown.markdown(result)
            soup = bs4.BeautifulSoup(html, "html.parser")
            result = soup.get_text()

        if normalize_whitespace:
            result = " ".join(result.split())

        return result

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
                dt = dt.replace(tzinfo=tz.gettz(timezone))
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


class DataValidator:
    """Data validator."""

    def __init__(
        self,
        config: ValidationConfig,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize validator.

        Args:
            config: Validation configuration
            logger: Optional logger instance
        """
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        self.cleaner = DataCleaner()

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

        # Apply field rules
        for field, rules in self.config.rules.items():
            try:
                value = data.get(field)
                cleaned_value = self._clean_value(field, value)
                validation_issues = self._validate_value(field, cleaned_value, rules)

                if validation_issues:
                    issues.extend(validation_issues)
                else:
                    cleaned[field] = cleaned_value

            except Exception as e:
                self.logger.error(f"Error validating {field}: {str(e)}")
                issues.append(ValidationIssue(field=field, rule="error", message=str(e), severity="error"))

        # Apply custom validators
        for validator in self.config.custom_validators:
            try:
                result = validator(cleaned)
                if isinstance(result, ValidationResult) and not result.is_valid:
                    issues.extend(result.errors)
            except Exception as e:
                self.logger.error(f"Error in custom validator: {str(e)}")
                issues.append(ValidationIssue(field="custom", rule="error", message=str(e), severity="error"))

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

        # Get field rules
        rules = self.config.rules.get(field, [])
        data_type = None

        # Find data type rule
        for rule, rule_value in rules:
            if rule == ValidationRule.TYPE:
                data_type = rule_value
                break

        if data_type == DataType.TEXT:
            return self.cleaner.clean_text(
                str(value), strip_html=self.config.strip_html, convert_markdown=self.config.convert_markdown, normalize_whitespace=self.config.normalize_whitespace
            )
        elif data_type == DataType.NUMBER:
            return self.cleaner.clean_number(value)
        elif data_type == DataType.TIMESTAMP:
            return self.cleaner.clean_timestamp(value, self.config.timezone)
        elif data_type == DataType.DURATION:
            return self.cleaner.clean_duration(value)

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

        for rule, rule_value in rules:
            if rule == ValidationRule.REQUIRED and value is None:
                issues.append(ValidationIssue(field=field, rule=rule, message=f"Field {field} is required", severity="error"))

            elif value is not None:
                if rule == ValidationRule.NON_EMPTY and not value:
                    issues.append(ValidationIssue(field=field, rule=rule, message=f"Field {field} cannot be empty", severity="error"))

                elif rule == ValidationRule.RANGE:
                    min_val, max_val = rule_value
                    if not min_val <= value <= max_val:
                        issues.append(ValidationIssue(field=field, rule=rule, message=f"Value must be between {min_val} and {max_val}", severity="error"))

                elif rule == ValidationRule.PATTERN and not re.match(rule_value, str(value)):
                    issues.append(ValidationIssue(field=field, rule=rule, message=f"Value must match pattern: {rule_value}", severity="error"))

                elif rule == ValidationRule.LENGTH:
                    min_len, max_len = rule_value
                    length = len(str(value))
                    if length < min_len or length > max_len:
                        issues.append(ValidationIssue(field=field, rule=rule, message=f"Length must be between {min_len} and {max_len}", severity="error"))

                elif rule == ValidationRule.ITEMS and isinstance(value, (list, tuple)):
                    min_items, max_items = rule_value
                    if len(value) < min_items or len(value) > max_items:
                        issues.append(ValidationIssue(field=field, rule=rule, message=f"Number of items must be between {min_items} and {max_items}", severity="error"))

        return issues
