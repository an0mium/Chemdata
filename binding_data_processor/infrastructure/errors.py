"""Error handling for infrastructure components.

This module defines custom exceptions for:
1. Cache errors
2. Rate limit errors
3. Circuit breaker errors
4. Monitoring errors

These exceptions provide detailed error information and consistent error handling.
"""

from typing import Any, Dict, Optional


class InfrastructureError(Exception):
    """Base class for all infrastructure errors."""

    def __init__(self, message: str, code: Optional[str] = None, details: Optional[Dict[str, Any]] = None):
        """Initialize error.

        Args:
            message: Error message
            code: Optional error code
            details: Optional error details
        """
        super().__init__(message)
        self.message = message
        self.code = code
        self.details = details or {}


class CacheError(InfrastructureError):
    """Base class for cache errors."""

    pass


class CacheConnectionError(CacheError):
    """Error connecting to cache backend."""

    pass


class CacheSerializationError(CacheError):
    """Error serializing/deserializing cached data."""

    pass


class CacheCapacityError(CacheError):
    """Cache capacity exceeded."""

    pass


class RateLimitError(InfrastructureError):
    """Base class for rate limit errors."""

    pass


class RateLimitExceededError(RateLimitError):
    """Rate limit exceeded."""

    def __init__(self, message: str = "Rate limit exceeded", limit: int = 0, remaining: int = 0, reset_time: Optional[int] = None, **kwargs: Any):
        """Initialize error.

        Args:
            message: Error message
            limit: Rate limit
            remaining: Remaining tokens
            reset_time: Time until limit reset
            **kwargs: Additional error details
        """
        super().__init__(message, details={"limit": limit, "remaining": remaining, "reset_time": reset_time, **kwargs})
        self.limit = limit
        self.remaining = remaining
        self.reset_time = reset_time


class CircuitBreakerError(InfrastructureError):
    """Base class for circuit breaker errors."""

    pass


class CircuitOpenError(CircuitBreakerError):
    """Circuit breaker is open."""

    def __init__(self, message: str = "Circuit breaker is open", name: str = "", failures: int = 0, reset_timeout: Optional[int] = None, **kwargs: Any):
        """Initialize error.

        Args:
            message: Error message
            name: Circuit breaker name
            failures: Number of failures
            reset_timeout: Time until reset
            **kwargs: Additional error details
        """
        super().__init__(message, details={"name": name, "failures": failures, "reset_timeout": reset_timeout, **kwargs})
        self.name = name
        self.failures = failures
        self.reset_timeout = reset_timeout


class CircuitHalfOpenError(CircuitBreakerError):
    """Circuit breaker is half-open."""

    pass


class MonitorError(InfrastructureError):
    """Base class for monitoring errors."""

    pass


class MetricError(MonitorError):
    """Error recording/retrieving metrics."""

    pass


class AlertError(MonitorError):
    """Error managing alerts."""

    pass


class StorageError(MonitorError):
    """Error accessing metric storage."""

    pass


class ValidationError(InfrastructureError):
    """Error validating configuration or parameters."""

    def __init__(self, message: str, field: Optional[str] = None, value: Optional[Any] = None, constraint: Optional[str] = None, **kwargs: Any):
        """Initialize error.

        Args:
            message: Error message
            field: Invalid field name
            value: Invalid value
            constraint: Validation constraint
            **kwargs: Additional error details
        """
        super().__init__(message, details={"field": field, "value": value, "constraint": constraint, **kwargs})
        self.field = field
        self.value = value
        self.constraint = constraint
