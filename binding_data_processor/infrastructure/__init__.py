"""Infrastructure components for binding data processor.

This package provides core infrastructure components:
1. Cache systems
2. Rate limiters
3. Circuit breakers 
4. Monitoring systems

These components provide essential infrastructure capabilities like:
- Caching with different backends (memory, file, redis)
- Rate limiting with token bucket and window algorithms
- Circuit breaking for fault tolerance
- Metrics collection and monitoring
"""

from .base import (
    BaseCache,
    BaseRateLimiter,
    BaseCircuitBreaker,
    BaseMonitor,
)

from .config import (
    CacheConfig,
    RateLimitConfig,
    CircuitBreakerConfig,
    MonitorConfig,
)

from .errors import (
    InfrastructureError,
    CacheError,
    CacheConnectionError,
    CacheSerializationError,
    CacheCapacityError,
    RateLimitError,
    RateLimitExceededError,
    CircuitBreakerError,
    CircuitOpenError,
    CircuitHalfOpenError,
    MonitorError,
    MetricError,
    AlertError,
    StorageError,
    ValidationError,
)

from .stats import (
    MetricType,
    TimeWindow,
    MetricValue,
    RollingMetric,
    StatsTracker,
    StatsManager,
)

__all__ = [
    # Base interfaces
    "BaseCache",
    "BaseRateLimiter",
    "BaseCircuitBreaker",
    "BaseMonitor",
    # Configuration
    "CacheConfig",
    "RateLimitConfig",
    "CircuitBreakerConfig",
    "MonitorConfig",
    # Errors
    "InfrastructureError",
    "CacheError",
    "CacheConnectionError",
    "CacheSerializationError",
    "CacheCapacityError",
    "RateLimitError",
    "RateLimitExceededError",
    "CircuitBreakerError",
    "CircuitOpenError",
    "CircuitHalfOpenError",
    "MonitorError",
    "MetricError",
    "AlertError",
    "StorageError",
    "ValidationError",
    # Stats
    "MetricType",
    "TimeWindow",
    "MetricValue",
    "RollingMetric",
    "StatsTracker",
    "StatsManager",
]
