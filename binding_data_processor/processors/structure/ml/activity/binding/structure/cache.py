"""Cache module for structure analysis system.

This module provides:
1. Result caching
2. Cache management
3. Cache invalidation
4. Caching utilities
"""

import logging
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Any, TypeVar, Generic, Callable
from pathlib import Path
import pickle
import hashlib
import functools

from .logging import get_logger
from .errors import ValidationError

logger = get_logger(__name__)

T = TypeVar("T")  # Generic type for cached values


@dataclass
class CacheEntry(Generic[T]):
    """Cache entry containing value and metadata."""

    value: T
    timestamp: float
    ttl: Optional[float] = None
    metadata: Optional[Dict[str, Any]] = None

    @property
    def is_expired(self) -> bool:
        """Check if entry has expired.

        Returns:
            True if entry has expired
        """
        if self.ttl is None:
            return False
        return time.time() - self.timestamp > self.ttl


class Cache(Generic[T]):
    """Generic cache implementation."""

    def __init__(
        self,
        ttl: Optional[float] = None,
        max_size: Optional[int] = None,
        persist: bool = False,
        cache_dir: Optional[Path] = None,
    ):
        """Initialize cache.

        Args:
            ttl: Time-to-live in seconds for cache entries
            max_size: Maximum number of entries to store
            persist: Whether to persist cache to disk
            cache_dir: Directory for persistent cache
        """
        self.ttl = ttl
        self.max_size = max_size
        self.persist = persist
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self._cache: Dict[str, CacheEntry[T]] = {}

        if self.persist:
            if not self.cache_dir:
                raise ValueError("Cache directory required for persistence")
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self._load_persistent()

    def get(self, key: str) -> Optional[T]:
        """Get value from cache.

        Args:
            key: Cache key

        Returns:
            Cached value or None if not found
        """
        try:
            entry = self._cache.get(key)
            if entry is None:
                return None

            if entry.is_expired:
                self.invalidate(key)
                return None

            return entry.value

        except Exception as e:
            logger.error(f"Error getting cache entry: {str(e)}")
            return None

    def set(
        self,
        key: str,
        value: T,
        ttl: Optional[float] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> bool:
        """Set cache entry.

        Args:
            key: Cache key
            value: Value to cache
            ttl: Optional TTL override
            metadata: Optional metadata

        Returns:
            True if successful
        """
        try:
            # Check cache size limit
            if self.max_size and len(self._cache) >= self.max_size:
                self._evict_oldest()

            # Create entry
            entry = CacheEntry(
                value=value,
                timestamp=time.time(),
                ttl=ttl or self.ttl,
                metadata=metadata,
            )
            self._cache[key] = entry

            # Persist if enabled
            if self.persist:
                self._persist_entry(key, entry)

            return True

        except Exception as e:
            logger.error(f"Error setting cache entry: {str(e)}")
            return False

    def invalidate(self, key: str) -> bool:
        """Invalidate cache entry.

        Args:
            key: Cache key

        Returns:
            True if entry was removed
        """
        try:
            if key in self._cache:
                del self._cache[key]
                if self.persist:
                    self._remove_persistent(key)
                return True
            return False

        except Exception as e:
            logger.error(f"Error invalidating cache entry: {str(e)}")
            return False

    def clear(self):
        """Clear all cache entries."""
        try:
            self._cache.clear()
            if self.persist:
                for path in self.cache_dir.glob("*.cache"):
                    path.unlink()

        except Exception as e:
            logger.error(f"Error clearing cache: {str(e)}")

    def _evict_oldest(self):
        """Evict oldest cache entry."""
        if not self._cache:
            return

        oldest_key = min(
            self._cache.keys(),
            key=lambda k: self._cache[k].timestamp,
        )
        self.invalidate(oldest_key)

    def _persist_entry(self, key: str, entry: CacheEntry[T]):
        """Persist cache entry to disk.

        Args:
            key: Cache key
            entry: Cache entry
        """
        try:
            path = self._get_cache_path(key)
            with open(path, "wb") as f:
                pickle.dump(entry, f)

        except Exception as e:
            logger.error(f"Error persisting cache entry: {str(e)}")

    def _remove_persistent(self, key: str):
        """Remove persistent cache entry.

        Args:
            key: Cache key
        """
        try:
            path = self._get_cache_path(key)
            if path.exists():
                path.unlink()

        except Exception as e:
            logger.error(f"Error removing persistent cache entry: {str(e)}")

    def _load_persistent(self):
        """Load persistent cache entries."""
        try:
            for path in self.cache_dir.glob("*.cache"):
                try:
                    with open(path, "rb") as f:
                        entry: CacheEntry[T] = pickle.load(f)
                        key = path.stem
                        if not entry.is_expired:
                            self._cache[key] = entry
                        else:
                            path.unlink()
                except Exception as e:
                    logger.error(f"Error loading cache entry {path}: {str(e)}")

        except Exception as e:
            logger.error(f"Error loading persistent cache: {str(e)}")

    def _get_cache_path(self, key: str) -> Path:
        """Get path for persistent cache entry.

        Args:
            key: Cache key

        Returns:
            Cache file path
        """
        return self.cache_dir / f"{key}.cache"


def cache_key(*args, **kwargs) -> str:
    """Generate cache key from arguments.

    Args:
        *args: Positional arguments
        **kwargs: Keyword arguments

    Returns:
        Cache key string
    """
    # Convert args/kwargs to strings
    key_parts = [str(arg) for arg in args]
    key_parts.extend(f"{k}={v}" for k, v in sorted(kwargs.items()))

    # Generate hash
    key = "|".join(key_parts)
    return hashlib.sha256(key.encode()).hexdigest()


def cached(
    cache: Cache,
    key_fn: Optional[Callable] = None,
    ttl: Optional[float] = None,
):
    """Decorator for caching function results.

    Args:
        cache: Cache instance
        key_fn: Optional function to generate cache key
        ttl: Optional TTL override

    Returns:
        Decorated function
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Generate cache key
            if key_fn:
                key = key_fn(*args, **kwargs)
            else:
                key = cache_key(*args, **kwargs)

            # Check cache
            result = cache.get(key)
            if result is not None:
                return result

            # Call function
            result = func(*args, **kwargs)

            # Cache result
            cache.set(key, result, ttl=ttl)
            return result

        return wrapper

    return decorator


class StructureCache(Cache[Any]):
    """Cache for structure analysis results."""

    def __init__(
        self,
        ttl: Optional[float] = None,
        max_size: Optional[int] = 1000,
        persist: bool = True,
        cache_dir: Optional[Path] = None,
    ):
        """Initialize structure cache.

        Args:
            ttl: Time-to-live in seconds for cache entries
            max_size: Maximum number of entries to store
            persist: Whether to persist cache to disk
            cache_dir: Directory for persistent cache
        """
        super().__init__(
            ttl=ttl,
            max_size=max_size,
            persist=persist,
            cache_dir=cache_dir,
        )

    def get_structure_key(self, structure_id: str, **kwargs) -> str:
        """Generate cache key for structure.

        Args:
            structure_id: Structure identifier
            **kwargs: Additional key components

        Returns:
            Cache key string
        """
        return cache_key(structure_id, **kwargs)


def create_structure_cache(
    cache_dir: Optional[Path] = None,
    ttl: Optional[float] = None,
    max_size: Optional[int] = 1000,
) -> StructureCache:
    """Create structure cache.

    Args:
        cache_dir: Directory for persistent cache
        ttl: Time-to-live in seconds for cache entries
        max_size: Maximum number of entries to store

    Returns:
        StructureCache instance
    """
    return StructureCache(
        ttl=ttl,
        max_size=max_size,
        persist=True,
        cache_dir=cache_dir,
    )
