"""Cache implementation for web components."""

from typing import Any, Dict, Optional
import time
import asyncio
from collections import OrderedDict


class Cache:
    """Simple in-memory cache with TTL and LRU eviction."""

    def __init__(self, max_size: int = 1000, ttl: int = 300):
        """Initialize cache.

        Args:
            max_size: Maximum number of items to store
            ttl: Time to live in seconds
        """
        self._cache: OrderedDict[str, Dict[str, Any]] = OrderedDict()
        self._max_size = max_size
        self._ttl = ttl
        self._hits = 0
        self._misses = 0
        self._lock = asyncio.Lock()

    async def get(self, key: str) -> Optional[Any]:
        """Get value from cache.

        Args:
            key: Cache key

        Returns:
            Cached value if found and not expired, None otherwise
        """
        async with self._lock:
            if key not in self._cache:
                self._misses += 1
                return None

            entry = self._cache[key]
            if time.time() > entry["expires"]:
                self._cache.pop(key)
                self._misses += 1
                return None

            # Move to end (most recently used)
            self._cache.move_to_end(key)
            self._hits += 1
            return entry["value"]

    async def set(self, key: str, value: Any, ttl: Optional[int] = None) -> None:
        """Set value in cache.

        Args:
            key: Cache key
            value: Value to cache
            ttl: Optional TTL override
        """
        async with self._lock:
            # Evict oldest if at max size
            while len(self._cache) >= self._max_size:
                self._cache.popitem(last=False)

            expires = time.time() + (ttl or self._ttl)
            self._cache[key] = {
                "value": value,
                "expires": expires,
            }

    async def delete(self, key: str) -> None:
        """Delete value from cache.

        Args:
            key: Cache key
        """
        async with self._lock:
            self._cache.pop(key, None)

    async def clear(self) -> None:
        """Clear all values from cache."""
        async with self._lock:
            self._cache.clear()
            self._hits = 0
            self._misses = 0

    async def size(self) -> int:
        """Get current cache size.

        Returns:
            Number of items in cache
        """
        async with self._lock:
            return len(self._cache)

    async def hits(self) -> int:
        """Get cache hit count.

        Returns:
            Number of cache hits
        """
        async with self._lock:
            return self._hits

    async def misses(self) -> int:
        """Get cache miss count.

        Returns:
            Number of cache misses
        """
        async with self._lock:
            return self._misses

    @property
    def ttl(self) -> int:
        """Get cache TTL.

        Returns:
            Cache TTL in seconds
        """
        return self._ttl

    @property
    def max_size(self) -> int:
        """Get cache max size.

        Returns:
            Maximum cache size
        """
        return self._max_size
