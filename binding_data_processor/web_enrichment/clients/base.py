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

    def _make_request(
        self,
        method: str,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        json: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        use_cache: bool = True,
    ) -> requests.Response:
        """Make HTTP request.

        Args:
            method: HTTP method
            url: Request URL
            params: Query parameters
            json: JSON body
            headers: Request headers
            use_cache: Whether to use cache

        Returns:
            Response object

        Raises:
            WebClientError: For request errors
        """
        try:
            return self.http.request(
                method=method,
                url=url,
                params=params,
                json=json,
                headers=headers,
                use_cache=use_cache,
            )
        except requests.exceptions.RequestException as e:
            self.http_errors += 1
            raise WebClientError(
                message=str(e),
                source=self.name,
                status_code=e.response.status_code if hasattr(e, "response") else None,
            )

    def _handle_rate_limit(
        self,
        response: requests.Response,
        retry_on_rate_limit: bool = True,
    ) -> bool:
        """Handle rate limiting.

        Args:
            response: Response object
            retry_on_rate_limit: Whether to retry

        Returns:
            Whether to retry request
        """
        if response.status_code == 429 and retry_on_rate_limit:
            retry_after = int(response.headers.get("Retry-After", "60"))
            self.logger.warning(f"Rate limited, waiting {retry_after}s before retry")
            import time

            time.sleep(retry_after)
            return True
        return False

    def _parse_response(
        self,
        response: requests.Response,
        endpoint: str,
        retry_on_parse_error: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Parse response data.

        Args:
            response: Response object
            endpoint: API endpoint
            retry_on_parse_error: Whether to retry on parse error

        Returns:
            Parsed data or None to retry

        Raises:
            ValueError: For parse errors
        """
        try:
            return response.json()
        except ValueError as e:
            if retry_on_parse_error:
                self.logger.warning(f"Parse error for {endpoint}, retrying: {e}")
                return None
            raise

    def _validate_data(
        self,
        data: Dict[str, Any],
        schema_name: str,
        endpoint: str,
    ) -> Dict[str, Any]:
        """Validate response data.

        Args:
            data: Response data
            schema_name: Schema name
            endpoint: API endpoint

        Returns:
            Validated data
        """
        validated_data, issues = self.validate_response(data, schema_name)
        if issues:
            self.logger.warning(f"Validation issues for {endpoint}: {issues}")
        return validated_data

    def _try_endpoint(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        json: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        use_cache: bool = True,
        schema_name: Optional[str] = None,
        retry_on_rate_limit: bool = True,
        retry_on_parse_error: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Try a single endpoint with retries.

        Args:
            method: HTTP method
            endpoint: API endpoint
            params: Query parameters
            json: JSON body
            headers: Request headers
            use_cache: Whether to use cache
            schema_name: Schema name for validation
            retry_on_rate_limit: Whether to retry on rate limit
            retry_on_parse_error: Whether to retry on parse error

        Returns:
            Response data or None if should retry with next endpoint

        Raises:
            WebClientError: For request errors
        """
        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        response = self._make_request(
            method=method,
            url=url,
            params=params,
            json=json,
            headers=headers,
            use_cache=use_cache,
        )

        # Handle rate limiting
        if self._handle_rate_limit(response, retry_on_rate_limit):
            return None

        # Parse response
        data = self._parse_response(response, endpoint, retry_on_parse_error)
        if data is None:  # Retry on parse error
            return None

        # Validate response
        if schema_name:
            data = self._validate_data(data, schema_name, endpoint)

        self.processed_items += 1
        return data

    def request(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        json: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        use_cache: bool = True,
        schema_name: Optional[str] = None,
        retry_on_rate_limit: bool = True,
        retry_on_parse_error: bool = True,
        fallback_endpoints: Optional[List[str]] = None,
        retry_strategy: str = "exponential_backoff",
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
            retry_on_rate_limit: Whether to retry on rate limit errors
            retry_on_parse_error: Whether to retry on parse errors
            fallback_endpoints: List of fallback endpoints to try
            retry_strategy: Retry strategy (exponential_backoff, fallback_endpoints)

        Returns:
            Response data

        Raises:
            WebClientError: For client errors
            ValidationError: For validation errors
        """
        endpoints = [endpoint] + (fallback_endpoints or [])
        last_error = None

        for current_endpoint in endpoints:
            try:
                data = self._try_endpoint(
                    method=method,
                    endpoint=current_endpoint,
                    params=params,
                    json=json,
                    headers=headers,
                    use_cache=use_cache,
                    schema_name=schema_name,
                    retry_on_rate_limit=retry_on_rate_limit,
                    retry_on_parse_error=retry_on_parse_error,
                )
                if data is not None:
                    return data
            except WebClientError as e:
                last_error = e
                self.logger.warning(f"Request failed for {current_endpoint}: {e}")
                if current_endpoint == endpoints[-1]:
                    self.failed_items += 1
                    raise

        # Should never reach here
        assert last_error is not None
        raise last_error

    def batch_request(
        self,
        method: str,
        endpoint: str,
        ids: List[Any],
        batch_size: int = 50,
        **kwargs: Any,
    ) -> List[Dict[str, Any]]:
        """Make batched requests.

        Args:
            method: HTTP method
            endpoint: API endpoint
            ids: List of IDs to batch
            batch_size: Batch size
            **kwargs: Additional arguments for request()

        Returns:
            List of response data
        """
        results = []
        for i in range(0, len(ids), batch_size):
            batch = ids[i : i + batch_size]
            if "params" in kwargs:
                kwargs["params"]["ids"] = batch
            else:
                kwargs["params"] = {"ids": batch}

            batch_data = self.request(
                method=method,
                endpoint=endpoint,
                **kwargs,
            )
            if isinstance(batch_data, list):
                results.extend(batch_data)
            else:
                results.append(batch_data)

        return results

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
