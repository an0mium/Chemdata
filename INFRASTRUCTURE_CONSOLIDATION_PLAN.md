# Infrastructure Consolidation Plan
# Integration Plan

## Overview
This document outlines planned integration improvements that have been identified but deferred for future implementation.

## Pending Integration Tasks

### 1. Configuration Integration
- Consolidate duplicate enrichment configurations between AppConfig and EnrichmentConfig
- Review and optimize validation logic in MLSystemConfig.validate_integration()
- Consider moving shared validation logic to a common utility

### 2. Component Integration
- Evaluate shared GPU usage between ML components
- Review batch processing coordination
- Optimize CPU thread allocation across components

### 3. Resource Management
- Implement coordinated cache management
- Add centralized rate limiting
- Develop unified monitoring strategy

### 4. Validation Strategy
- Create consistent validation patterns
- Implement cross-component validation
- Add integration tests

## Implementation Notes

### Configuration Consolidation
Current AppConfig has duplicate enrichment fields that should be consolidated:
```python
@dataclass
class AppConfig:
    # Data sources and enrichment
    data_sources: DataSourceConfig
    web_enrichment: WebEnrichmentConfig  
    enrichment: EnrichmentConfig  # <-- Duplicate with web_enrichment
```

### MLSystem Integration
Current validation has redundant checks that should be optimized:
```python
def validate_integration(self):
    # Resource validation appears in multiple places
    if not 0 < self.gpu_memory_fraction <= 1:
        raise ValueError("gpu_memory_fraction must be in (0, 1]")
    if self.cpu_threads < -1:
        raise ValueError("cpu_threads must be >= -1")
```

## Next Steps
1. Create detailed design docs for each integration area
2. Add integration tests to verify behavior
3. Implement changes incrementally with proper testing
4. Update documentation to reflect new integration patterns

## Timeline
- Phase 1: Design & Planning (Q1 2025)
- Phase 2: Implementation (Q2 2025) 
- Phase 3: Testing & Validation (Q3 2025)
- Phase 4: Documentation & Release (Q4 2025)

## Current State Analysis

We have multiple infrastructure implementations spread across different parts of the system:

### 1. Cache Systems
- Pipeline Cache (CacheManager)
  - File-based storage
  - JSON serialization
  - Global expiry
  - Basic stats
- Web Cache (Cache)
  - In-memory storage
  - LRU eviction
  - Async support
  - Per-item TTL

### 2. Rate Limiters
- Pipeline Rate Limiter
  - Token bucket algorithm
  - Comprehensive stats
  - Configurable limits
  - Stats tracking
- Web Rate Limiter
  - Singleton pattern
  - Async support
  - Simple stats
  - Token-based limiting

### 3. Monitoring
- Pipeline Monitor
  - Comprehensive metrics
  - Log management
  - Report generation
  - Alert system
- Document Monitor
  - File system events
  - Directory watching
  - Async processing
  - Callback system

## Consolidation Strategy

### 1. Create Base Infrastructure Package

```
binding_data_processor/infrastructure/
  ├── __init__.py
  ├── base.py          # Base interfaces
  ├── config.py        # Configuration
  ├── errors.py        # Error handling
  ├── stats.py         # Statistics
  ├── cache/          # Cache implementations
  │   ├── __init__.py
  │   ├── base.py
  │   ├── memory.py
  │   ├── file.py
  │   └── redis.py
  ├── rate_limit/     # Rate limiting
  │   ├── __init__.py
  │   ├── base.py
  │   ├── token.py
  │   └── window.py
  ├── circuit/        # Circuit breaker
  │   ├── __init__.py
  │   ├── base.py
  │   ├── metrics.py
  │   └── breaker.py
  ├── monitoring/     # Monitoring
  │   ├── __init__.py
  │   ├── base.py
  │   ├── metrics.py
  │   ├── alerts.py
  │   └── logging.py
  └── utils/          # Shared utilities
      ├── __init__.py
      ├── async_utils.py
      └── serialization.py
```

### 2. Define Core Interfaces

1. Base Cache Interface:
```python
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

class BaseCache(ABC):
    """Base interface for all cache implementations."""
    
    @abstractmethod
    async def get(self, key: str) -> Optional[Any]:
        """Get value from cache."""
        pass
        
    @abstractmethod
    async def set(
        self, 
        key: str, 
        value: Any, 
        ttl: Optional[int] = None
    ) -> None:
        """Set value in cache."""
        pass
        
    @abstractmethod
    async def delete(self, key: str) -> None:
        """Delete value from cache."""
        pass
        
    @abstractmethod
    async def clear(self) -> None:
        """Clear all values from cache."""
        pass
        
    @abstractmethod
    async def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        pass
```

2. Base Rate Limiter Interface:
```python
class BaseRateLimiter(ABC):
    """Base interface for rate limiters."""
    
    @abstractmethod
    async def acquire(self, key: str, tokens: int = 1) -> bool:
        """Acquire rate limit tokens."""
        pass
        
    @abstractmethod
    async def release(self, key: str, tokens: int = 1) -> None:
        """Release rate limit tokens."""
        pass
        
    @abstractmethod
    async def get_limit(self, key: str) -> Dict[str, Any]:
        """Get rate limit info."""
        pass
```

3. Base Circuit Breaker Interface:
```python
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
        """Execute function with circuit breaking."""
        pass
        
    @abstractmethod
    async def can_execute(self) -> bool:
        """Check if request can be executed."""
        pass
        
    @abstractmethod
    async def record_success(self) -> None:
        """Record successful execution."""
        pass
        
    @abstractmethod
    async def record_failure(self, error: Exception) -> None:
        """Record failed execution."""
        pass
        
    @abstractmethod
    async def get_metrics(self) -> Dict[str, Any]:
        """Get circuit breaker metrics."""
        pass
```

4. Base Monitor Interface:
```python
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
    async def record_metric(
        self, 
        name: str, 
        value: Any
    ) -> None:
        """Record metric value."""
        pass
        
    @abstractmethod
    async def get_metrics(self) -> Dict[str, Any]:
        """Get monitoring metrics."""
        pass
```

### 3. Implement Core Components

1. Cache Components:
```python
class MemoryCache(BaseCache):
    """In-memory cache with LRU eviction."""
    
    def __init__(
        self,
        max_size: int = 1000,
        ttl: int = 300
    ):
        self._cache = OrderedDict()
        self._max_size = max_size
        self._ttl = ttl
        self._lock = asyncio.Lock()

class FileCache(BaseCache):
    """File-based persistent cache."""
    
    def __init__(
        self,
        cache_dir: Path,
        ttl: int = 300
    ):
        self.cache_dir = cache_dir
        self._ttl = ttl
```

2. Rate Limiter Components:
```python
class TokenBucketLimiter(BaseRateLimiter):
    """Token bucket rate limiter."""
    
    def __init__(
        self,
        limit: int,
        window: int,
        burst: Optional[int] = None
    ):
        self._limit = limit
        self._window = window
        self._burst = burst or limit
        self._buckets = {}
        self._lock = asyncio.Lock()

class WindowLimiter(BaseRateLimiter):
    """Sliding window rate limiter."""
    
    def __init__(
        self,
        limit: int,
        window: int
    ):
        self._limit = limit
        self._window = window
        self._windows = {}
        self._lock = asyncio.Lock()
```

3. Circuit Breaker Components:
```python
@dataclass
class CircuitConfig:
    """Circuit breaker configuration."""
    failure_threshold: int = 5
    failure_timeout: int = 60
    reset_timeout: int = 60
    half_open_timeout: int = 30
    test_requests: int = 1
    track_metrics: bool = True

@dataclass
class CircuitMetrics:
    """Circuit breaker metrics."""
    total_requests: int = 0
    successful_requests: int = 0
    failed_requests: int = 0
    rejected_requests: int = 0
    total_trips: int = 0
    total_resets: int = 0
    last_failure: Optional[datetime] = None
    last_success: Optional[datetime] = None

class StandardCircuitBreaker(BaseCircuitBreaker):
    """Standard circuit breaker implementation."""
    
    def __init__(
        self,
        name: str,
        config: Optional[CircuitConfig] = None,
        logger: Optional[logging.Logger] = None
    ):
        self.name = name
        self.config = config or CircuitConfig()
        self.logger = logger or logging.getLogger(__name__)
        self.state = CircuitState.CLOSED
        self.metrics = CircuitMetrics()
        self._lock = asyncio.Lock()
```

4. Monitor Components:
```python
class MetricsMonitor(BaseMonitor):
    """Metrics collection and monitoring."""
    
    def __init__(
        self,
        config: MonitorConfig
    ):
        self.config = config
        self.metrics = MonitoringMetrics()
        self.alerts = AlertComponent()
        self.logger = LoggingComponent()

class FileSystemMonitor(BaseMonitor):
    """File system monitoring."""
    
    def __init__(
        self,
        watch_dirs: Dict[Path, Set[str]],
        processor: Optional[Callable] = None
    ):
        self.watch_dirs = watch_dirs
        self.processor = processor
        self.observer = Observer()
```

### 4. Circuit Breaker States

```python
class CircuitState(Enum):
    """Circuit breaker states."""
    CLOSED = "closed"      # Normal operation
    OPEN = "open"         # Failing, reject requests
    HALF_OPEN = "half_open"  # Testing if service recovered
```

### 5. Configuration System

```python
@dataclass
class CacheConfig:
    """Cache configuration."""
    backend: str = "memory"
    ttl: int = 300
    max_size: int = 1000
    cache_dir: Optional[Path] = None

@dataclass
class RateLimitConfig:
    """Rate limit configuration."""
    algorithm: str = "token"
    limit: int = 100
    window: int = 60
    burst: Optional[int] = None

@dataclass
class CircuitConfig:
    """Circuit breaker configuration."""
    failure_threshold: int = 5
    failure_timeout: int = 60
    reset_timeout: int = 60
    half_open_timeout: int = 30
    test_requests: int = 1
    track_metrics: bool = True

@dataclass
class MonitorConfig:
    """Monitor configuration."""
    metrics_enabled: bool = True
    alerts_enabled: bool = True
    log_dir: Optional[Path] = None
    report_format: str = "json"
```

### 6. Factory System

```python
class InfrastructureFactory:
    """Factory for infrastructure components."""
    
    @staticmethod
    def create_cache(
        config: CacheConfig
    ) -> BaseCache:
        """Create cache instance."""
        if config.backend == "memory":
            return MemoryCache(
                max_size=config.max_size,
                ttl=config.ttl
            )
        elif config.backend == "file":
            return FileCache(
                cache_dir=config.cache_dir,
                ttl=config.ttl
            )
        raise ValueError(f"Unknown cache backend: {config.backend}")
    
    @staticmethod
    def create_rate_limiter(
        config: RateLimitConfig
    ) -> BaseRateLimiter:
        """Create rate limiter instance."""
        if config.algorithm == "token":
            return TokenBucketLimiter(
                limit=config.limit,
                window=config.window,
                burst=config.burst
            )
        elif config.algorithm == "window":
            return WindowLimiter(
                limit=config.limit,
                window=config.window
            )
        raise ValueError(f"Unknown rate limit algorithm: {config.algorithm}")
    
    @staticmethod
    def create_monitor(
        config: MonitorConfig
    ) -> BaseMonitor:
        """Create monitor instance."""
        return MetricsMonitor(config)
        
    @staticmethod
    def create_circuit_breaker(
        name: str,
        config: Optional[CircuitConfig] = None
    ) -> BaseCircuitBreaker:
        """Create circuit breaker instance."""
        return StandardCircuitBreaker(name, config)
```

### 7. Migration Steps

1. Create New Structure:
- Set up infrastructure package
- Implement base interfaces
- Add core components
- Create configuration system
- Implement factory system

2. Update Dependencies:
- Update pipeline code to use new infrastructure
- Update web code to use new infrastructure
- Add tests for new implementation
- Update documentation

3. Deprecate Old Implementations:
- Mark old classes as deprecated
- Provide migration guide
- Remove after transition period

### 8. Testing Strategy

1. Unit Tests:
```python
class TestCache:
    """Test cache implementations."""
    
    async def test_basic_operations(self, cache: BaseCache):
        """Test basic cache operations."""
        await cache.set("key", "value")
        assert await cache.get("key") == "value"
        
    async def test_ttl(self, cache: BaseCache):
        """Test TTL expiration."""
        await cache.set("key", "value", ttl=1)
        await asyncio.sleep(2)
        assert await cache.get("key") is None
```

2. Circuit Breaker Tests:
```python
class TestCircuitBreaker:
    """Test circuit breaker implementations."""
    
    async def test_execution(self, breaker: BaseCircuitBreaker):
        """Test normal execution."""
        result = await breaker.execute(lambda: "success")
        assert result == "success"
        
    async def test_failure_threshold(self, breaker: BaseCircuitBreaker):
        """Test failure threshold."""
        for _ in range(5):
            with pytest.raises(Exception):
                await breaker.execute(lambda: 1/0)
        
        assert not await breaker.can_execute()
        
    async def test_recovery(self, breaker: BaseCircuitBreaker):
        """Test circuit recovery."""
        # Force open state
        for _ in range(5):
            with pytest.raises(Exception):
                await breaker.execute(lambda: 1/0)
                
        # Wait for reset
        await asyncio.sleep(60)
        
        # Should allow test request
        result = await breaker.execute(lambda: "recovered")
        assert result == "recovered"
```

3. Integration Tests:
```python
class TestInfrastructure:
    """Test infrastructure integration."""
    
    async def test_cache_rate_limit(
        self,
        cache: BaseCache,
        limiter: BaseRateLimiter
    ):
        """Test cache with rate limiting."""
        if await limiter.acquire():
            await cache.set("key", "value")
            assert await cache.get("key") == "value"
```

### 9. Benefits

1. Consistency:
- Single interface for each component
- Unified configuration
- Standard error handling
- Consistent stats collection

2. Flexibility:
- Multiple implementations
- Pluggable architecture
- Configurable features
- Easy to extend

3. Reliability:
- Better error handling
- Comprehensive testing
- Clear documentation
- Migration support

4. Performance:
- Optimized implementations
- Proper concurrency
- Efficient resource usage
- Better monitoring

### 10. Next Steps

1. Implementation:
- Create new package structure
- Implement base interfaces
- Add core components
- Write tests

2. Migration:
- Update pipeline code
- Update web code
- Run integration tests
- Update docs

3. Deployment:
- Stage rollout
- Monitor performance
- Collect feedback
- Make adjustments

4. Cleanup:
- Remove old code
- Update dependencies
- Final testing
- Release notes
