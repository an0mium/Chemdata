# Cache Consolidation Plan

## Current State

We have two cache implementations:

### 1. Pipeline Cache (CacheManager)
- Location: binding_data_processor/pipeline/infrastructure/cache.py
- Features:
  - File-based persistent storage
  - JSON serialization
  - Global expiry time
  - Basic size stats
  - Pipeline data focus

### 2. Web Cache (Cache)
- Location: binding_data_processor/web/cache.py  
- Features:
  - In-memory storage
  - LRU eviction
  - Async support
  - Per-item TTL
  - Detailed stats
  - Web request focus

## Consolidation Strategy

### 1. Create Unified Interface

```python
class CacheBackend(ABC):
    """Abstract base class for cache backends."""
    
    @abstractmethod
    async def get(self, key: str) -> Optional[Any]:
        """Get value from cache."""
        pass
        
    @abstractmethod
    async def set(self, key: str, value: Any, ttl: Optional[int] = None) -> None:
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

### 2. Implement Storage Backends

1. Memory Backend:
```python
class MemoryCache(CacheBackend):
    """In-memory cache with LRU eviction."""
    def __init__(self, max_size: int = 1000, ttl: int = 300):
        self._cache = OrderedDict()
        self._max_size = max_size
        self._ttl = ttl
        self._lock = asyncio.Lock()
```

2. File Backend:
```python
class FileCache(CacheBackend):
    """File-based persistent cache."""
    def __init__(self, cache_dir: Path, ttl: int = 300):
        self.cache_dir = cache_dir
        self._ttl = ttl
```

3. Redis Backend (Future):
```python
class RedisCache(CacheBackend):
    """Redis-based distributed cache."""
    def __init__(self, redis_url: str, ttl: int = 300):
        self._redis = Redis.from_url(redis_url)
        self._ttl = ttl
```

### 3. Create Factory

```python
class CacheFactory:
    """Factory for creating cache instances."""
    
    @staticmethod
    def create(
        backend: str = "memory",
        config: Optional[Dict[str, Any]] = None
    ) -> CacheBackend:
        """Create cache instance."""
        if backend == "memory":
            return MemoryCache(**config or {})
        elif backend == "file":
            return FileCache(**config or {})
        elif backend == "redis":
            return RedisCache(**config or {})
        raise ValueError(f"Unknown backend: {backend}")
```

### 4. Enhance Features

1. Stats Collection:
```python
@dataclass
class CacheStats:
    """Cache statistics."""
    hits: int = 0
    misses: int = 0
    size: int = 0
    oldest_entry: float = 0
    newest_entry: float = 0
```

2. Error Handling:
```python
class CacheError(Exception):
    """Base class for cache errors."""
    pass

class CacheKeyError(CacheError):
    """Error for invalid cache keys."""
    pass

class CacheValueError(CacheError):
    """Error for invalid cache values."""
    pass
```

3. Configuration:
```python
@dataclass
class CacheConfig:
    """Cache configuration."""
    backend: str = "memory"
    ttl: int = 300
    max_size: int = 1000
    cache_dir: Optional[Path] = None
    redis_url: Optional[str] = None
```

### 5. Migration Steps

1. Create New Structure:
```
binding_data_processor/infrastructure/cache/
  ├── __init__.py
  ├── base.py        # Base classes and interfaces
  ├── memory.py      # Memory backend
  ├── file.py        # File backend
  ├── redis.py       # Redis backend (future)
  ├── factory.py     # Cache factory
  ├── errors.py      # Error classes
  ├── config.py      # Configuration
  └── stats.py       # Statistics
```

2. Update Dependencies:
- Update pipeline code to use new cache
- Update web code to use new cache
- Add tests for new implementation
- Update documentation

3. Deprecate Old Implementations:
- Mark old classes as deprecated
- Provide migration guide
- Remove after transition period

### 6. Testing Strategy

1. Unit Tests:
```python
class TestCacheBackend:
    """Test cache backend implementation."""
    
    async def test_basic_operations(self, cache: CacheBackend):
        """Test basic cache operations."""
        await cache.set("key", "value")
        assert await cache.get("key") == "value"
        
    async def test_ttl(self, cache: CacheBackend):
        """Test TTL expiration."""
        await cache.set("key", "value", ttl=1)
        await asyncio.sleep(2)
        assert await cache.get("key") is None
```

2. Integration Tests:
```python
class TestCacheIntegration:
    """Test cache integration."""
    
    async def test_pipeline_integration(self, cache: CacheBackend):
        """Test pipeline cache integration."""
        # Test pipeline specific use cases
        
    async def test_web_integration(self, cache: CacheBackend):
        """Test web cache integration."""
        # Test web specific use cases
```

### 7. Documentation

1. API Documentation:
```python
class Cache:
    """Unified cache interface.
    
    This class provides a consistent interface for caching data across
    the application. It supports multiple backend storage options and
    provides configurable TTL, stats collection, and error handling.
    
    Examples:
        >>> cache = Cache(backend="memory")
        >>> await cache.set("key", "value")
        >>> value = await cache.get("key")
    """
```

2. Migration Guide:
- Document changes from old implementations
- Provide examples of updating code
- List breaking changes
- Include troubleshooting tips

### 8. Benefits

1. Consistency:
- Single interface for all caching
- Consistent error handling
- Unified configuration
- Standard stats collection

2. Flexibility:
- Multiple storage backends
- Configurable TTL
- Pluggable architecture
- Easy to extend

3. Reliability:
- Better error handling
- Comprehensive testing
- Clear documentation
- Migration support

4. Performance:
- Optimized implementations
- Proper concurrency
- Efficient storage
- Better monitoring

### 9. Next Steps

1. Implementation:
- Create new package structure
- Implement base classes
- Add storage backends
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
