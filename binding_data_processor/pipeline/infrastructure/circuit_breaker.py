"""Circuit breaker pattern implementation for handling external service failures."""

import time
from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, Dict, Optional, TypeVar

from ...models.core import ErrorCode, McpError
from ...logger import LogManager

T = TypeVar('T')  # Generic type for circuit breaker return value


class CircuitState(Enum):
    """Circuit breaker states."""
    CLOSED = 'closed'     # Normal operation
    OPEN = 'open'         # Failing, reject requests
    HALF_OPEN = 'half_open'  # Testing if service recovered


@dataclass
class CircuitBreakerConfig:
    """Configuration for circuit breaker behavior."""
    
    failure_threshold: int = 5      # Number of failures before opening circuit
    reset_timeout: int = 60         # Seconds to wait before attempting reset
    half_open_timeout: int = 30     # Seconds between half-open retry attempts
    
    def __post_init__(self):
        """Validate configuration values."""
        if self.failure_threshold < 1:
            raise ValueError("failure_threshold must be >= 1")
        if self.reset_timeout < 1:
            raise ValueError("reset_timeout must be >= 1")
        if self.half_open_timeout < 1:
            raise ValueError("half_open_timeout must be >= 1")


class CircuitBreaker:
    """
    Implements circuit breaker pattern for protecting external service calls.
    
    The circuit breaker pattern prevents cascading failures by stopping requests
    to a failing service until it recovers. It has three states:
    
    - CLOSED: Normal operation, requests pass through
    - OPEN: Service is failing, requests are rejected
    - HALF_OPEN: Testing if service has recovered
    
    Example:
        ```python
        breaker = CircuitBreaker("my-service")
        
        try:
            with breaker:
                result = make_external_request()
        except CircuitBreakerError:
            # Handle service unavailable
            pass
        ```
    """
    
    def __init__(
        self,
        name: str,
        config: Optional[CircuitBreakerConfig] = None
    ):
        """
        Initialize circuit breaker.
        
        Args:
            name: Name for this circuit breaker instance
            config: Optional configuration, uses defaults if not provided
        """
        self.name = name
        self.config = config or CircuitBreakerConfig()
        self.state = CircuitState.CLOSED
        self.failures = 0
        self.last_failure_time = 0
        self.last_success_time = time.time()
        self.logger = LogManager().get_logger(f"circuit_breaker.{name}")
        self.stats: Dict[str, int] = {
            "requests": 0,
            "failures": 0,
            "rejects": 0,
            "successes": 0,
        }

    def __enter__(self):
        """Context manager entry."""
        if not self.can_execute():
            self.stats["rejects"] += 1
            raise McpError(
                ErrorCode.CircuitBreakerOpen,
                f"Circuit breaker {self.name} is {self.state.value}"
            )
        self.stats["requests"] += 1
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if exc_type is None:
            self.record_success()
        else:
            self.record_failure()
        return False  # Don't suppress exceptions

    def execute(self, func: Callable[..., T], *args, **kwargs) -> T:
        """
        Execute a function with circuit breaker protection.
        
        Args:
            func: Function to execute
            *args: Positional arguments for func
            **kwargs: Keyword arguments for func
            
        Returns:
            Result from func
            
        Raises:
            McpError: If circuit is open
            Exception: Any exception from func
        """
        with self:
            return func(*args, **kwargs)

    def can_execute(self) -> bool:
        """
        Check if request can be executed based on circuit state.
        
        Returns:
            True if request can proceed, False if it should be rejected
        """
        if self.state == CircuitState.CLOSED:
            return True
            
        if self.state == CircuitState.OPEN:
            if time.time() - self.last_failure_time >= self.config.reset_timeout:
                self.state = CircuitState.HALF_OPEN
                self.logger.info(f"Circuit {self.name} entering half-open state")
                return True
            return False
            
        # HALF_OPEN state
        if time.time() - self.last_failure_time >= self.config.half_open_timeout:
            return True
        return False

    def record_failure(self):
        """Record a failure and potentially open the circuit."""
        self.failures += 1
        self.last_failure_time = time.time()
        self.stats["failures"] += 1
        
        if self.failures >= self.config.failure_threshold:
            self.state = CircuitState.OPEN
            self.logger.warning(
                f"Circuit {self.name} opened after {self.failures} failures"
            )

    def record_success(self):
        """Record a success and potentially close the circuit."""
        self.failures = 0
        self.last_success_time = time.time()
        self.stats["successes"] += 1
        
        if self.state in (CircuitState.OPEN, CircuitState.HALF_OPEN):
            self.state = CircuitState.CLOSED
            self.logger.info(f"Circuit {self.name} closed after success")

    def get_stats(self) -> Dict[str, Any]:
        """
        Get circuit breaker statistics.
        
        Returns:
            Dictionary of statistics including:
            - current state
            - failure count
            - request counts
            - timing information
        """
        return {
            "state": self.state.value,
            "failures": self.failures,
            "last_failure": self.last_failure_time,
            "last_success": self.last_success_time,
            "stats": self.stats.copy()
        }

    def reset(self):
        """Reset circuit breaker to initial state."""
        self.state = CircuitState.CLOSED
        self.failures = 0
        self.last_failure_time = 0
        self.last_success_time = time.time()
        self.stats = {k: 0 for k in self.stats}
        self.logger.info(f"Circuit {self.name} reset to initial state")
