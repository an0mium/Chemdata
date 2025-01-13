"""Configuration settings for utility modules."""

from dataclasses import dataclass
from typing import Optional


@dataclass
class CircuitBreakerConfig:
    """Configuration for circuit breaker pattern."""

    failure_threshold: int = 5  # Number of failures before opening circuit
    reset_timeout: float = 60.0  # Seconds to wait before attempting reset
    half_open_timeout: float = 30.0  # Seconds to wait in half-open state
    excluded_exceptions: Optional[tuple] = None  # Exceptions that don't count as failures


# Default circuit breaker configuration
CIRCUIT_BREAKER = CircuitBreakerConfig(
    failure_threshold=5,
    reset_timeout=60.0,
    half_open_timeout=30.0,
    excluded_exceptions=(ValueError, KeyError),  # Common exceptions to ignore
)
