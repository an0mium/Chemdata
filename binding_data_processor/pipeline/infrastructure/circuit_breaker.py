"""Circuit breaker pattern implementation.

This module provides circuit breaker functionality to:
1. Detect failing external services
2. Prevent cascading failures
3. Enable graceful degradation
4. Allow service recovery
5. Track detailed metrics

Example:
    ```python
    breaker = CircuitBreaker("my-service")
    
    # Context manager usage
    try:
        with breaker:
            result = make_external_request()
    except CircuitBreakerError:
        # Handle service unavailable
        pass
        
    # Function wrapper usage
    result = breaker.execute(
        make_external_request,
        fallback=handle_failure
    )
    ```
"""

import logging
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from typing import Dict, Optional, Any, Callable, TypeVar, Generic, List

from ...models.core import ErrorCode, McpError
from ...logger import LogManager
from ...utils.config import CIRCUIT_BREAKER


T = TypeVar("T")  # Generic type for circuit breaker return value


class CircuitState(Enum):
    """Circuit breaker states."""
    CLOSED = "closed"     # Normal operation
    OPEN = "open"         # Failing, reject requests
    HALF_OPEN = "half_open"  # Testing if service recovered


@dataclass
class CircuitConfig:
    """Circuit breaker configuration."""
    
    # Failure settings
    failure_threshold: int = CIRCUIT_BREAKER.get("failure_threshold", 5)
    failure_timeout: int = CIRCUIT_BREAKER.get("failure_timeout", 60)
    min_throughput: int = CIRCUIT_BREAKER.get("min_throughput", 10)
    
    # Recovery settings
    reset_timeout: int = CIRCUIT_BREAKER.get("reset_timeout", 60)
    half_open_timeout: int = CIRCUIT_BREAKER.get("half_open_timeout", 30)
    test_requests: int = CIRCUIT_BREAKER.get("test_requests", 1)
    
    # Monitoring settings
    track_metrics: bool = True
    metric_window: int = 60  # seconds
    metric_buckets: int = 10
    
    def __post_init__(self):
        """Validate configuration values."""
        if self.failure_threshold < 1:
            raise ValueError("failure_threshold must be >= 1")
        if self.reset_timeout < 1:
            raise ValueError("reset_timeout must be >= 1")
        if self.half_open_timeout < 1:
            raise ValueError("half_open_timeout must be >= 1")


@dataclass
class CircuitMetrics:
    """Circuit breaker metrics."""
    
    # Request counts
    total_requests: int = 0
    successful_requests: int = 0
    failed_requests: int = 0
    rejected_requests: int = 0
    
    # State changes
    total_trips: int = 0
    total_resets: int = 0
    last_failure: Optional[datetime] = None
    last_success: Optional[datetime] = None
    
    # Time windows
    request_counts: Dict[int, int] = field(default_factory=dict)
    failure_counts: Dict[int, int] = field(default_factory=dict)
    latencies: Dict[int, float] = field(default_factory=dict)
    
    # Error tracking
    error_types: Dict[str, int] = field(default_factory=dict)
    error_history: List[Dict[str, Any]] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert metrics to dictionary format."""
        return {
            "requests": {
                "total": self.total_requests,
                "successful": self.successful_requests,
                "failed": self.failed_requests,
                "rejected": self.rejected_requests,
                "success_rate": self._get_success_rate(),
            },
            "state": {
                "trips": self.total_trips,
                "resets": self.total_resets,
                "last_failure": (
                    self.last_failure.isoformat()
                    if self.last_failure
                    else None
                ),
                "last_success": (
                    self.last_success.isoformat()
                    if self.last_success
                    else None
                ),
            },
            "windows": {
                "requests": self.request_counts,
                "failures": self.failure_counts,
                "latencies": self.latencies,
            },
            "errors": {
                "types": self.error_types,
                "history": self.error_history[-10:],  # Last 10 errors
            },
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get request success rate."""
        if not self.total_requests:
            return None
        return self.successful_requests / self.total_requests


class CircuitBreakerError(Exception):
    """Exception raised when circuit breaker prevents execution."""
    
    def __init__(self, circuit_name: str):
        """Initialize error.
        
        Args:
            circuit_name: Name of the circuit breaker
        """
        self.circuit_name = circuit_name
        super().__init__(f"Circuit breaker {circuit_name} is open")


class CircuitBreaker(Generic[T]):
    """Circuit breaker implementation."""

    def __init__(
        self,
        name: str,
        config: Optional[CircuitConfig] = None,
        logger: Optional[LogManager] = None,
    ):
        """Initialize circuit breaker.
        
        Args:
            name: Circuit breaker name
            config: Optional configuration
            logger: Optional logger instance
        """
        self.name = name
        self.config = config or CircuitConfig()
        self.logger = logger or LogManager().get_logger(f"circuit_breaker.{name}")
        
        # Initialize state
        self.state = CircuitState.CLOSED
        self._lock = threading.Lock()
        self._last_failure_time = 0
        self._test_requests = 0
        
        # Initialize metrics
        self.metrics = CircuitMetrics()
        self._last_window = self._get_window()
        self._init_metrics()

    def __enter__(self):
        """Context manager entry."""
        if not self.can_execute():
            self.metrics.rejected_requests += 1
            raise CircuitBreakerError(self.name)
        self.metrics.total_requests += 1
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if exc_type is None:
            self.record_success()
        else:
            self.record_failure(exc_val)
        return False  # Don't suppress exceptions

    def _init_metrics(self) -> None:
        """Initialize metric tracking."""
        if not self.config.track_metrics:
            return
            
        # Initialize time windows
        window = self._get_window()
        for offset in range(self.config.metric_buckets):
            bucket = window - offset
            self.metrics.request_counts[bucket] = 0
            self.metrics.failure_counts[bucket] = 0
            self.metrics.latencies[bucket] = 0.0

    def execute(
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
            Exception: Any exception from func if no fallback
        """
        try:
            # Check if request allowed
            if not self.can_execute():
                self.metrics.rejected_requests += 1
                self.logger.warning(f"Circuit {self.name} is {self.state.value}")
                return self._handle_open_circuit(fallback, *args, **kwargs)
            
            # Track request
            start_time = datetime.now()
            self._track_request()
            
            try:
                # Execute function
                result = func(*args, **kwargs)
                
                # Track success
                self._track_success(start_time)
                return result
                
            except Exception as e:
                # Track failure
                self._track_failure(start_time, e)
                raise
            
        except Exception as e:
            self.logger.error(
                f"Circuit {self.name} execution failed: {str(e)}"
            )
            if fallback:
                return self._handle_error(e, fallback, *args, **kwargs)
            raise

    def can_execute(self) -> bool:
        """Check if request can be executed."""
        with self._lock:
            now = time.time()
            
            if self.state == CircuitState.CLOSED:
                return True
                
            if self.state == CircuitState.OPEN:
                # Check reset timeout
                if now - self._last_failure_time >= self.config.reset_timeout:
                    self.state = CircuitState.HALF_OPEN
                    self.logger.info(f"Circuit {self.name} entering half-open state")
                    return True
                return False
                
            # HALF_OPEN state
            if self._test_requests >= self.config.test_requests:
                return False
                
            if now - self._last_failure_time >= self.config.half_open_timeout:
                self._test_requests += 1
                return True
                
            return False

    def _handle_open_circuit(
        self,
        fallback: Optional[Callable[..., T]],
        *args: Any,
        **kwargs: Any,
    ) -> Optional[T]:
        """Handle open circuit state."""
        if fallback:
            try:
                return fallback(*args, **kwargs)
            except Exception as e:
                self.logger.error(
                    f"Circuit {self.name} fallback failed: {str(e)}"
                )
                raise CircuitBreakerError(self.name)
        raise CircuitBreakerError(self.name)

    def _handle_error(
        self,
        error: Exception,
        fallback: Callable[..., T],
        *args: Any,
        **kwargs: Any,
    ) -> T:
        """Handle execution error."""
        try:
            return fallback(*args, **kwargs)
        except Exception as e:
            self.logger.error(
                f"Circuit {self.name} fallback failed: {str(e)}"
            )
            raise CircuitBreakerError(self.name)

    def _track_request(self) -> None:
        """Track request metrics."""
        if not self.config.track_metrics:
            return
            
        with self._lock:
            # Update counts
            self.metrics.total_requests += 1
            
            # Update time window
            window = self._get_window()
            self._roll_windows(window)
            self.metrics.request_counts[window] += 1

    def _track_success(
        self,
        start_time: datetime,
    ) -> None:
        """Track successful request."""
        if not self.config.track_metrics:
            return
            
        with self._lock:
            # Update counts
            self.metrics.successful_requests += 1
            self.metrics.last_success = datetime.now()
            
            # Update latency
            window = self._get_window()
            latency = (datetime.now() - start_time).total_seconds() * 1000
            self.metrics.latencies[window] = (
                (self.metrics.latencies[window] * 
                (self.metrics.request_counts[window] - 1) +
                latency) /
                self.metrics.request_counts[window]
            )
            
            # Check half-open state
            if self.state == CircuitState.HALF_OPEN:
                if self._test_requests >= self.config.test_requests:
                    self._close_circuit()

    def _track_failure(
        self,
        start_time: datetime,
        error: Exception,
    ) -> None:
        """Track failed request."""
        if not self.config.track_metrics:
            return
            
        with self._lock:
            # Update counts
            self.metrics.failed_requests += 1
            self.metrics.last_failure = datetime.now()
            
            # Track error
            error_type = error.__class__.__name__
            self.metrics.error_types[error_type] = (
                self.metrics.error_types.get(error_type, 0) + 1
            )
            self.metrics.error_history.append({
                "timestamp": datetime.now().isoformat(),
                "type": error_type,
                "message": str(error),
            })
            
            # Update time window
            window = self._get_window()
            self.metrics.failure_counts[window] += 1
            
            # Check failure threshold
            if self.state == CircuitState.CLOSED:
                recent_failures = sum(
                    count
                    for bucket, count in self.metrics.failure_counts.items()
                    if bucket > window - self.config.failure_timeout
                )
                if recent_failures >= self.config.failure_threshold:
                    self._open_circuit()
            
            # Check half-open state
            elif self.state == CircuitState.HALF_OPEN:
                self._open_circuit()

    def _open_circuit(self) -> None:
        """Open circuit breaker."""
        self.state = CircuitState.OPEN
        self._last_failure_time = time.time()
        self.metrics.total_trips += 1
        self.logger.warning(f"Circuit {self.name} opened")

    def _close_circuit(self) -> None:
        """Close circuit breaker."""
        self.state = CircuitState.CLOSED
        self._last_failure_time = 0
        self._test_requests = 0
        self.metrics.total_resets += 1
        self.logger.info(f"Circuit {self.name} closed")

    def _roll_windows(
        self,
        current_window: int,
    ) -> None:
        """Roll metric time windows."""
        if current_window == self._last_window:
            return
            
        # Remove old windows
        for bucket in list(self.metrics.request_counts.keys()):
            if bucket <= current_window - self.config.metric_buckets:
                del self.metrics.request_counts[bucket]
                del self.metrics.failure_counts[bucket]
                del self.metrics.latencies[bucket]
        
        # Add new windows
        for bucket in range(
            self._last_window + 1,
            current_window + 1
        ):
            self.metrics.request_counts[bucket] = 0
            self.metrics.failure_counts[bucket] = 0
            self.metrics.latencies[bucket] = 0.0
        
        self._last_window = current_window

    def _get_window(self) -> int:
        """Get current time window."""
        return int(
            time.time() /
            (self.config.metric_window / self.config.metric_buckets)
        )

    def get_metrics(self) -> Dict[str, Any]:
        """Get circuit breaker metrics."""
        return {
            "name": self.name,
            "state": self.state.value,
            "metrics": self.metrics.to_dict(),
        }

    def reset(self) -> None:
        """Reset circuit breaker to initial state."""
        with self._lock:
            self.state = CircuitState.CLOSED
            self._last_failure_time = 0
            self._test_requests = 0
            self.metrics = CircuitMetrics()
            self._last_window = self._get_window()
            self._init_metrics()
            self.logger.info(f"Circuit {self.name} reset to initial state")
