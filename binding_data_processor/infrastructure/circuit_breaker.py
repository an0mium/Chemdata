"""Circuit breaker pattern implementation for handling failures and retries."""

from dataclasses import dataclass
from typing import Optional


@dataclass
class CircuitConfig:
    """Configuration for circuit breaker pattern."""

    failure_threshold: int = 5  # Number of failures before opening circuit
    recovery_timeout: int = 60  # Seconds to wait before attempting recovery
    reset_timeout: int = 300  # Seconds to wait before resetting failure count
    concurrency_limit: int = 10  # Maximum number of concurrent requests
    max_retries: int = 3  # Maximum number of retry attempts
    retry_delay: float = 1.0  # Base delay between retries in seconds
    timeout: float = 30.0  # Request timeout in seconds
    jitter: float = 0.1  # Random jitter factor for retry delays

    def __post_init__(self):
        """Validate configuration values."""
        if self.failure_threshold < 1:
            raise ValueError("failure_threshold must be at least 1")
        if self.recovery_timeout < 0:
            raise ValueError("recovery_timeout must be non-negative")
        if self.reset_timeout < 0:
            raise ValueError("reset_timeout must be non-negative")
        if self.concurrency_limit < 1:
            raise ValueError("concurrency_limit must be at least 1")
        if self.max_retries < 0:
            raise ValueError("max_retries must be non-negative")
        if self.retry_delay < 0:
            raise ValueError("retry_delay must be non-negative")
        if self.timeout <= 0:
            raise ValueError("timeout must be positive")
        if self.jitter < 0 or self.jitter > 1:
            raise ValueError("jitter must be between 0 and 1")


class CircuitBreaker:
    """Circuit breaker implementation for handling failures and retries."""

    def __init__(self, config: Optional[CircuitConfig] = None):
        """Initialize circuit breaker with configuration.

        Args:
            config: Optional circuit breaker configuration. If not provided,
                   default configuration will be used.
        """
        self.config = config or CircuitConfig()
        self._failure_count = 0
        self._last_failure_time = 0
        self._is_open = False
        self._concurrent_requests = 0

    @property
    def is_open(self) -> bool:
        """Check if circuit breaker is open (failing)."""
        return self._is_open

    def increment_failure(self) -> None:
        """Increment failure count and potentially open circuit."""
        self._failure_count += 1
        if self._failure_count >= self.config.failure_threshold:
            self._is_open = True

    def reset(self) -> None:
        """Reset circuit breaker state."""
        self._failure_count = 0
        self._is_open = False
        self._concurrent_requests = 0

    def acquire(self) -> bool:
        """Attempt to acquire permission to make a request.

        Returns:
            bool: True if request is allowed, False if circuit is open
                 or concurrency limit reached
        """
        if self._is_open:
            return False
        if self._concurrent_requests >= self.config.concurrency_limit:
            return False
        self._concurrent_requests += 1
        return True

    def release(self) -> None:
        """Release a request slot."""
        if self._concurrent_requests > 0:
            self._concurrent_requests -= 1

    def __enter__(self):
        """Context manager entry."""
        if not self.acquire():
            raise RuntimeError("Circuit breaker is open")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.release()
        if exc_type is not None:
            self.increment_failure()
