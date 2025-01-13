"""Base interfaces for infrastructure components.

This module defines the core interfaces for:
1. Cache systems
2. Rate limiters
3. Circuit breakers
4. Monitoring systems

These interfaces provide a consistent API across different implementations.
"""

from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, TypeVar, Union

T = TypeVar("T")  # Generic type for return values


class BaseCache(ABC):
    """Base interface for all cache implementations."""

    @abstractmethod
    async def get(self, key: str) -> Optional[Any]:
        """Get value from cache.

        Args:
            key: Cache key

        Returns:
            Cached value if found, None otherwise
        """
        pass

    @abstractmethod
    async def set(self, key: str, value: Any, ttl: Optional[int] = None) -> None:
        """Set value in cache.

        Args:
            key: Cache key
            value: Value to cache
            ttl: Optional time-to-live in seconds
        """
        pass

    @abstractmethod
    async def delete(self, key: str) -> None:
        """Delete value from cache.

        Args:
            key: Cache key to delete
        """
        pass

    @abstractmethod
    async def clear(self) -> None:
        """Clear all values from cache."""
        pass

    @abstractmethod
    async def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics.

        Returns:
            Dictionary of cache statistics
        """
        pass


class BaseRateLimiter(ABC):
    """Base interface for rate limiters."""

    @abstractmethod
    async def acquire(self, key: str, tokens: int = 1) -> bool:
        """Acquire rate limit tokens.

        Args:
            key: Rate limit key
            tokens: Number of tokens to acquire

        Returns:
            True if tokens were acquired, False otherwise
        """
        pass

    @abstractmethod
    async def release(self, key: str, tokens: int = 1) -> None:
        """Release rate limit tokens.

        Args:
            key: Rate limit key
            tokens: Number of tokens to release
        """
        pass

    @abstractmethod
    async def get_limit(self, key: str) -> Dict[str, Any]:
        """Get rate limit info.

        Args:
            key: Rate limit key

        Returns:
            Dictionary with rate limit information
        """
        pass


class BaseCircuitBreaker(ABC):
    """Base interface for circuit breakers."""

    @abstractmethod
    async def execute(
        self,
        func: Callable[..., T],
        fallback: Optional[Callable[..., T]] = None,
        *args: Any,
        **kwargs: Any,
    ) -> Optional[T]:
        """Execute function with circuit breaking.

        Args:
            func: Function to execute
            fallback: Optional fallback function
            *args: Positional arguments
            **kwargs: Keyword arguments

        Returns:
            Function result or fallback result

        Raises:
            CircuitBreakerError: If circuit is open and no fallback
        """
        pass

    @abstractmethod
    async def can_execute(self) -> bool:
        """Check if request can be executed.

        Returns:
            True if request can be executed, False otherwise
        """
        pass

    @abstractmethod
    async def record_success(self) -> None:
        """Record successful execution."""
        pass

    @abstractmethod
    async def record_failure(self, error: Exception) -> None:
        """Record failed execution.

        Args:
            error: Exception that occurred
        """
        pass

    @abstractmethod
    async def get_metrics(self) -> Dict[str, Any]:
        """Get circuit breaker metrics.

        Returns:
            Dictionary of circuit breaker metrics
        """
        pass


class BaseMonitor(ABC):
    """Base interface for monitoring."""

    @abstractmethod
    async def start(self) -> None:
        """Start monitoring."""
        pass

    @abstractmethod
    async def stop(self) -> None:
        """Stop monitoring."""
        pass

    @abstractmethod
    async def record_metric(self, name: str, value: Any, tags: Optional[Dict[str, str]] = None) -> None:
        """Record metric value.

        Args:
            name: Metric name
            value: Metric value
            tags: Optional metric tags
        """
        pass

    @abstractmethod
    async def get_metrics(
        self,
        names: Optional[List[str]] = None,
        tags: Optional[Dict[str, str]] = None,
        start_time: Optional[datetime] = None,
        end_time: Optional[datetime] = None,
    ) -> Dict[str, Any]:
        """Get monitoring metrics.

        Args:
            names: Optional list of metric names to retrieve
            tags: Optional tags to filter metrics
            start_time: Optional start time for time range
            end_time: Optional end time for time range

        Returns:
            Dictionary of monitoring metrics
        """
        pass

    @abstractmethod
    async def add_alert(
        self,
        name: str,
        condition: str,
        threshold: float,
        window: int,
        handler: Callable[[str, float], None],
    ) -> None:
        """Add metric alert.

        Args:
            name: Metric name to monitor
            condition: Alert condition (>, <, >=, <=, ==)
            threshold: Alert threshold value
            window: Time window in seconds
            handler: Alert handler function
        """
        pass
