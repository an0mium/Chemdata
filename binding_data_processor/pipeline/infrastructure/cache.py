"""Cache management system for API responses and processed data.

This module provides a centralized caching system for the pipeline infrastructure,
handling temporary storage of API responses, processed data, and intermediate results.
"""

import json
import time
from pathlib import Path
from typing import Any, Dict, Optional

from binding_data_processor.config import CACHE_DIR, CACHE_EXPIRY


class CacheManager:
    """Manages caching of API responses and processed data in the pipeline infrastructure.

    This class provides a file-based caching system with automatic expiration and
    cleanup capabilities. It is used throughout the pipeline to cache expensive
    operations like API calls and data processing results.

    Attributes:
        cache_dir: Path to the cache directory
    """

    def __init__(self) -> None:
        """Initialize cache manager and ensure cache directory exists."""
        self.cache_dir = CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_path(self, key: str) -> Path:
        """Get the file path for a cache key.

        Args:
            key: Unique identifier for the cached data

        Returns:
            Path object pointing to the cache file location
        """
        # Use hash of key to avoid filesystem issues with long/invalid characters
        safe_key = str(hash(key))
        return self.cache_dir / f"{safe_key}.json"

    def get(self, key: str) -> Optional[Dict[str, Any]]:
        """Retrieve data from cache if it exists and hasn't expired.

        Args:
            key: Unique identifier for the cached data

        Returns:
            Cached data if valid and not expired, None otherwise

        Example:
            >>> cache = CacheManager()
            >>> data = cache.get("my_expensive_operation")
            >>> if data is None:
            ...     data = perform_expensive_operation()
            ...     cache.set("my_expensive_operation", data)
        """
        cache_path = self._get_cache_path(key)

        if not cache_path.exists():
            return None

        try:
            with cache_path.open("r") as f:
                cached = json.load(f)

            # Check if cache has expired
            if time.time() - cached["timestamp"] > CACHE_EXPIRY:
                cache_path.unlink()  # Remove expired cache
                return None

            return cached["data"]

        except (json.JSONDecodeError, KeyError, OSError) as e:
            print(f"Cache read error for {key}: {str(e)}")
            return None

    def set(self, key: str, data: Dict[str, Any]) -> bool:
        """Store data in cache with timestamp.

        Args:
            key: Unique identifier for the data
            data: Data to cache

        Returns:
            True if successful, False otherwise

        Example:
            >>> cache = CacheManager()
            >>> result = expensive_operation()
            >>> success = cache.set("expensive_operation_result", result)
        """
        cache_path = self._get_cache_path(key)

        try:
            cache_data = {"timestamp": time.time(), "data": data}

            with cache_path.open("w") as f:
                json.dump(cache_data, f)

            return True

        except (OSError, TypeError) as e:
            print(f"Cache write error for {key}: {str(e)}")
            return False

    def invalidate(self, key: str) -> bool:
        """Remove item from cache.

        Args:
            key: Cache key to invalidate

        Returns:
            True if successful or file didn't exist, False on error

        Example:
            >>> cache = CacheManager()
            >>> cache.invalidate("stale_data")
        """
        cache_path = self._get_cache_path(key)

        try:
            if cache_path.exists():
                cache_path.unlink()
            return True

        except OSError as e:
            print(f"Cache invalidation error for {key}: {str(e)}")
            return False

    def clear(self) -> bool:
        """Clear all cached data.

        Returns:
            True if successful, False on error

        Example:
            >>> cache = CacheManager()
            >>> cache.clear()  # Clear all cached data
        """
        try:
            for cache_file in self.cache_dir.glob("*.json"):
                cache_file.unlink()
            return True

        except OSError as e:
            print(f"Cache clear error: {str(e)}")
            return False

    def get_cache_size(self) -> int:
        """Get total size of cached data in bytes.

        Returns:
            Total size of cache in bytes

        Example:
            >>> cache = CacheManager()
            >>> size = cache.get_cache_size()
            >>> print(f"Cache size: {size} bytes")
        """
        return sum(f.stat().st_size for f in self.cache_dir.glob("*.json"))

    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics.

        Returns:
            Dictionary containing:
                - total_entries: Number of cache entries
                - total_size_bytes: Total size of cache in bytes
                - oldest_entry: Timestamp of oldest entry
                - newest_entry: Timestamp of newest entry

        Example:
            >>> cache = CacheManager()
            >>> stats = cache.get_cache_stats()
            >>> print(f"Cache entries: {stats['total_entries']}")
        """
        cache_files = list(self.cache_dir.glob("*.json"))

        return {
            "total_entries": len(cache_files),
            "total_size_bytes": self.get_cache_size(),
            "oldest_entry": min((f.stat().st_mtime for f in cache_files), default=0),
            "newest_entry": max((f.stat().st_mtime for f in cache_files), default=0),
        }
