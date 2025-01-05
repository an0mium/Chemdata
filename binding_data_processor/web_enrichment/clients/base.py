"""Base client for web enrichment data sources.

This module provides a base client with:
- HTTP client with circuit breaker
- Rate limiting and retries
- Response caching
- Schema validation
- Data validation
- Error handling
- Metrics collection
"""

import logging
from abc import ABC, abstractmethod
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import requests

from ..http_client_enhanced import HTTPClientEnhanced
from ..validation.schema import (
    DataSource,
    SchemaValidator,
    ValidationLevel,
)
from ..validation.data import (
    DataCleaner,
    DataValidator,
    ValidationConfig,
    ValidationIssue,
)
from ...pipeline.infrastructure.circuit_breaker import CircuitConfig


class WebClientError(Exception):
    """Base error for web clients."""

    def __init__(
        self,
        message: str,
        source: str,
        status_code: Optional[int] = None,
        details: Optional[Dict[str, Any]] = None,
    ):
        """Initialize error.
        
        Args:
            message: Error message
            source: Error source
            status_code: Optional HTTP status code
            details: Optional error details
        """
        super().__init__(message)
        self.source = source
        self.status_code = status_code
        self.details = details or {}
        self.timestamp = datetime.now().isoformat()


class RateLimitError(WebClientError):
    """Rate limit exceeded error."""

    def __init__(self, source: str, retry_after: Optional[int] = None):
        """Initialize error.
        
        Args:
            source: Error source
            retry_after: Optional seconds to wait before retry
        """
        super().__init__(
            message="Rate limit exceeded",
            source=source,
            status_code=429,
            details={"retry_after": retry_after} if retry_after else None,
        )
        self.retry_after = retry_after


class ValidationError(WebClientError):
    """Validation error."""

    def __init__(
        self,
        message: str,
        source: str,
        issues: List[ValidationIssue],
    ):
        """Initialize error.
        
        Args:
            message: Error message
            source: Error source
            issues: Validation issues
        """
        super().__init__(
            message=message,
            source=source,
            details={"issues": [vars(issue) for issue in issues]},
        )
        self.issues = issues


class WebClient(ABC):
    """Base client for web enrichment data sources."""

    def __init__(
        self,
        name: str,
        base_url: str,
        data_source: DataSource,
        requests_per_second: float = 1.0,
        max_retries: int = 3,
        timeout: int = 30,
        cache_ttl: int = 3600,
        circuit_config: Optional[CircuitConfig] = None,
        validation_level: ValidationLevel = ValidationLevel.NORMAL,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize client.
        
        Args:
            name: Client name
            base_url: Base URL for all requests
            data_source: Data source type
            requests_per_second: Maximum requests per second
            max_retries: Maximum number of retries
            timeout: Request timeout in seconds
            cache_ttl: Cache TTL in seconds
            circuit_config: Optional circuit breaker config
            validation_level: Validation strictness level
            logger: Optional logger instance
        """
        self.name = name
        self.base_url = base_url.rstrip("/")
        self.data_source = data_source
        self.logger = logger or logging.getLogger(self.__class__.__name__)

        # Initialize HTTP client
        self.http = HTTPClientEnhanced(
            name=name,
            requests_per_second=requests_per_second,
            max_retries=max_retries,
            timeout=timeout,
            cache_ttl=cache_ttl,
            circuit_config=circuit_config,
            logger=self.logger,
        )

        # Initialize validators
        self.schema_validator = SchemaValidator()
        self.data_validator = DataValidator(
            config=self._get_validation_config(),
            logger=self.logger,
        )

        # Initialize metrics
        self.processed_items = 0
        self.failed_items = 0
        self.validation_errors = 0
        self.http_errors = 0

    @abstractmethod
    def _get_validation_config(self) -> ValidationConfig:
        """Get validation configuration.
        
        Returns:
            Validation configuration
        """
        pass

    def request(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        json: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        use_cache: bool = True,
        schema_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Make HTTP request with validation and error handling.
        
        Args:
            method: HTTP method (GET, POST, etc)
            endpoint: API endpoint (will be joined with base_url)
            params: Query parameters
            json: JSON body for POST/PUT
            headers: Additional headers
            use_cache: Whether to use cache for GET requests
            schema_name: Optional schema name for validation
            
        Returns:
            Response data
            
        Raises:
            WebClientError: For client errors
            ValidationError: For validation errors
        """
        try:
            # Make request
            url = f"{self.base_url}/{endpoint.lstrip('/')}"
            response = self.http.request(
                method=method,
                url=url,
                params=params,
                json=json,
                headers=headers,
                use_cache=use_cache,
            )

            # Parse response
            data = response.json()

            # Validate response if schema provided
            if schema_name:
                data, issues = self.validate_response(data, schema_name)
                if issues:
                    self.logger.warning(
                        f"Validation issues for {endpoint}: {issues}"
                    )

            self.processed_items += 1
            return data

        except requests.exceptions.RequestException as e:
            self.http_errors += 1
            self.failed_items += 1
            raise WebClientError(
                message=str(e),
                source=self.name,
                status_code=e.response.status_code if hasattr(e, "response") else None,
            )

    def validate_response(
        self,
        data: Dict[str, Any],
        schema_name: str,
    ) -> Tuple[Dict[str, Any], List[ValidationIssue]]:
        """Validate response data.
        
        Args:
            data: Response data
            schema_name: Schema name
            
        Returns:
            Tuple of (cleaned data, validation issues)
            
        Raises:
            ValidationError: If validation fails
        """
        # Validate schema
        schema_result = self.schema_validator.validate(
            data=data,
            source=self.data_source,
            schema_name=schema_name,
        )
        if not schema_result.is_valid:
            self.validation_errors += 1
            raise ValidationError(
                message="Schema validation failed",
                source=self.name,
                issues=[
                    ValidationIssue(
                        field="schema",
                        rule=issue.rule,
                        message=issue.message,
                        severity=issue.severity,
                    )
                    for issue in schema_result.errors
                ],
            )

        # Validate and clean data
        cleaned, issues = self.data_validator.validate(data)
        if issues:
            self.validation_errors += 1
            raise ValidationError(
                message="Data validation failed",
                source=self.name,
                issues=issues,
            )

        return cleaned, issues

    def clean_text(
        self,
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
        return DataCleaner.clean_text(
            text=text,
            strip_html=strip_html,
            convert_markdown=convert_markdown,
            normalize_whitespace=normalize_whitespace,
        )

    def clean_number(
        self,
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
        return DataCleaner.clean_number(
            value=value,
            min_value=min_value,
            max_value=max_value,
            round_digits=round_digits,
        )

    def clean_timestamp(
        self,
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
        return DataCleaner.clean_timestamp(
            value=value,
            timezone=timezone,
        )

    def clean_duration(self, value: Any) -> str:
        """Clean duration data.
        
        Args:
            value: Value to clean
            
        Returns:
            Cleaned duration string
            
        Raises:
            ValueError: If value is invalid
        """
        return DataCleaner.clean_duration(value)

    def get_metrics(self) -> Dict[str, Any]:
        """Get client metrics.
        
        Returns:
            Client metrics
        """
        return {
            "processed_items": self.processed_items,
            "failed_items": self.failed_items,
            "validation_errors": self.validation_errors,
            "http_errors": self.http_errors,
            "http_client": self.http.get_metrics(),
        }

    def close(self) -> None:
        """Close client and cleanup."""
        self.http.close()
