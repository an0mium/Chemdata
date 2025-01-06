"""Tests for the cache management system."""

import json
import time
from pathlib import Path
from unittest import mock

import pytest

from binding_data_processor.pipeline.infrastructure.cache import CacheManager

# Create tests directory if it doesn't exist
Path(__file__).parent.mkdir(parents=True, exist_ok=True)


@pytest.fixture
def cache_dir(tmp_path: Path) -> Path:
    """Create a temporary directory for cache files."""
    cache_path = tmp_path / "cache"
    cache_path.mkdir()
    return cache_path


@pytest.fixture
def cache_manager(cache_dir: Path) -> CacheManager:
    """Create a CacheManager instance with a temporary cache directory."""
    with mock.patch("binding_data_processor.config.CACHE_DIR", cache_dir):
        return CacheManager()


def test_cache_initialization(cache_dir: Path) -> None:
    """Test that cache manager creates cache directory if it doesn't exist."""
    cache_path = cache_dir / "new_cache"
    with mock.patch("binding_data_processor.config.CACHE_DIR", cache_path):
        CacheManager()
        assert cache_path.exists()
        assert cache_path.is_dir()


def test_cache_set_and_get(cache_manager: CacheManager) -> None:
    """Test setting and retrieving data from cache."""
    test_data = {"key": "value"}
    assert cache_manager.set("test_key", test_data)

    retrieved = cache_manager.get("test_key")
    assert retrieved == test_data


def test_cache_expiry(cache_manager: CacheManager) -> None:
    """Test that expired cache entries are removed."""
    test_data = {"key": "value"}

    # Set cache with very short expiry
    with mock.patch("binding_data_processor.config.CACHE_EXPIRY", 0.1):
        assert cache_manager.set("test_key", test_data)

        # Should be available immediately
        assert cache_manager.get("test_key") == test_data

        # Wait for expiry
        time.sleep(0.2)

        # Should be None after expiry
        assert cache_manager.get("test_key") is None


def test_cache_invalidate(cache_manager: CacheManager) -> None:
    """Test manually invalidating cache entries."""
    test_data = {"key": "value"}
    assert cache_manager.set("test_key", test_data)

    # Verify data was cached
    assert cache_manager.get("test_key") == test_data

    # Invalidate the entry
    assert cache_manager.invalidate("test_key")

    # Verify data is gone
    assert cache_manager.get("test_key") is None


def test_cache_clear(cache_manager: CacheManager) -> None:
    """Test clearing all cache entries."""
    # Add multiple entries
    assert cache_manager.set("key1", {"data": 1})
    assert cache_manager.set("key2", {"data": 2})

    # Verify entries exist
    assert cache_manager.get("key1") == {"data": 1}
    assert cache_manager.get("key2") == {"data": 2}

    # Clear cache
    assert cache_manager.clear()

    # Verify all entries are gone
    assert cache_manager.get("key1") is None
    assert cache_manager.get("key2") is None


def test_cache_size(cache_manager: CacheManager) -> None:
    """Test getting cache size."""
    # Add entries
    assert cache_manager.set("key1", {"data": "small"})
    assert cache_manager.set("key2", {"data": "larger data"})

    # Get size
    size = cache_manager.get_cache_size()
    assert size > 0

    # Clear and verify size is 0
    assert cache_manager.clear()
    assert cache_manager.get_cache_size() == 0


def test_cache_stats(cache_manager: CacheManager) -> None:
    """Test getting cache statistics."""
    # Add entries
    assert cache_manager.set("key1", {"data": 1})
    assert cache_manager.set("key2", {"data": 2})

    # Get stats
    stats = cache_manager.get_cache_stats()

    # Verify stats structure
    assert isinstance(stats, dict)
    assert "total_entries" in stats
    assert "total_size_bytes" in stats
    assert "oldest_entry" in stats
    assert "newest_entry" in stats

    # Verify values
    assert stats["total_entries"] == 2
    assert stats["total_size_bytes"] > 0
    assert stats["oldest_entry"] > 0
    assert stats["newest_entry"] >= stats["oldest_entry"]


def test_cache_error_handling(cache_manager: CacheManager) -> None:
    """Test error handling in cache operations."""
    # Test invalid JSON data
    invalid_data = {"data": object()}
    # Try to serialize invalid data to verify json usage
    with pytest.raises(TypeError):
        json.dumps(invalid_data)

    # Test cache handling of invalid data
    with mock.patch("json.dump", side_effect=TypeError("Invalid data")):
        assert not cache_manager.set("key", invalid_data)

    # Test file read error
    with mock.patch("pathlib.Path.open", side_effect=OSError("Read error")):
        assert cache_manager.get("key") is None

    # Test file write error
    with mock.patch("pathlib.Path.open", side_effect=OSError("Write error")):
        assert not cache_manager.set("key", {"data": "value"})

    # Test clear error
    with mock.patch("pathlib.Path.unlink", side_effect=OSError("Delete error")):
        assert not cache_manager.clear()


def test_cache_path_generation(cache_manager: CacheManager) -> None:
    """Test cache path generation for different keys."""
    # Test different key types
    key1 = "simple_key"
    key2 = "complex/key/with/slashes"
    key3 = "key with spaces and $pecial ch@rs"

    # Get paths
    path1 = cache_manager._get_cache_path(key1)
    path2 = cache_manager._get_cache_path(key2)
    path3 = cache_manager._get_cache_path(key3)

    # Verify paths are valid and unique
    assert path1.suffix == ".json"
    assert path2.suffix == ".json"
    assert path3.suffix == ".json"
    assert path1 != path2
    assert path2 != path3
    assert path1 != path3
