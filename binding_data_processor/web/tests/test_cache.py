"""Tests for cache implementation."""

import pytest
import asyncio
import time
from ..cache import Cache


@pytest.fixture
async def cache():
    """Create test cache instance."""
    cache = Cache(max_size=3, ttl=1)
    yield cache
    await cache.clear()


@pytest.mark.asyncio
async def test_set_get():
    """Test setting and getting cache values."""
    cache = Cache()
    await cache.set("key1", "value1")
    value = await cache.get("key1")
    assert value == "value1"


@pytest.mark.asyncio
async def test_ttl_expiration(cache):
    """Test cache TTL expiration."""
    await cache.set("key1", "value1")
    value = await cache.get("key1")
    assert value == "value1"

    # Wait for TTL to expire
    await asyncio.sleep(1.1)
    value = await cache.get("key1")
    assert value is None


@pytest.mark.asyncio
async def test_max_size(cache):
    """Test cache max size enforcement."""
    # Add items up to max size
    await cache.set("key1", "value1")
    await cache.set("key2", "value2")
    await cache.set("key3", "value3")

    # Verify all items present
    assert await cache.get("key1") == "value1"
    assert await cache.get("key2") == "value2"
    assert await cache.get("key3") == "value3"

    # Add one more item (should evict oldest)
    await cache.set("key4", "value4")

    # Verify oldest item evicted
    assert await cache.get("key1") is None
    assert await cache.get("key2") == "value2"
    assert await cache.get("key3") == "value3"
    assert await cache.get("key4") == "value4"


@pytest.mark.asyncio
async def test_lru_eviction(cache):
    """Test LRU eviction policy."""
    # Add items
    await cache.set("key1", "value1")
    await cache.set("key2", "value2")
    await cache.set("key3", "value3")

    # Access key1 (should move to end)
    await cache.get("key1")

    # Add new item (should evict key2)
    await cache.set("key4", "value4")

    # Verify key2 was evicted
    assert await cache.get("key1") == "value1"
    assert await cache.get("key2") is None
    assert await cache.get("key3") == "value3"
    assert await cache.get("key4") == "value4"


@pytest.mark.asyncio
async def test_clear():
    """Test clearing cache."""
    cache = Cache()
    await cache.set("key1", "value1")
    await cache.set("key2", "value2")

    await cache.clear()
    assert await cache.get("key1") is None
    assert await cache.get("key2") is None
    assert await cache.size() == 0


@pytest.mark.asyncio
async def test_delete():
    """Test deleting cache entry."""
    cache = Cache()
    await cache.set("key1", "value1")
    await cache.set("key2", "value2")

    await cache.delete("key1")
    assert await cache.get("key1") is None
    assert await cache.get("key2") == "value2"


@pytest.mark.asyncio
async def test_hit_miss_counts():
    """Test cache hit/miss counting."""
    cache = Cache()
    await cache.set("key1", "value1")

    # Test hits
    await cache.get("key1")
    await cache.get("key1")
    assert await cache.hits() == 2

    # Test misses
    await cache.get("key2")
    await cache.get("key3")
    assert await cache.misses() == 2


@pytest.mark.asyncio
async def test_custom_ttl():
    """Test custom TTL override."""
    cache = Cache(ttl=10)
    await cache.set("key1", "value1", ttl=1)

    value = await cache.get("key1")
    assert value == "value1"

    # Wait for custom TTL to expire
    await asyncio.sleep(1.1)
    value = await cache.get("key1")
    assert value is None


@pytest.mark.asyncio
async def test_size():
    """Test cache size tracking."""
    cache = Cache()
    assert await cache.size() == 0

    await cache.set("key1", "value1")
    assert await cache.size() == 1

    await cache.set("key2", "value2")
    assert await cache.size() == 2

    await cache.delete("key1")
    assert await cache.size() == 1

    await cache.clear()
    assert await cache.size() == 0


@pytest.mark.asyncio
async def test_concurrent_access():
    """Test concurrent cache access."""
    cache = Cache()

    async def writer():
        for i in range(100):
            await cache.set(f"key{i}", f"value{i}")
            await asyncio.sleep(0.001)

    async def reader():
        for i in range(100):
            await cache.get(f"key{i}")
            await asyncio.sleep(0.001)

    # Run concurrent writers and readers
    writers = [writer() for _ in range(3)]
    readers = [reader() for _ in range(3)]
    await asyncio.gather(*writers, *readers)

    # Verify cache integrity
    size = await cache.size()
    assert 0 <= size <= cache.max_size


@pytest.mark.asyncio
async def test_properties():
    """Test cache properties."""
    cache = Cache(max_size=100, ttl=60)
    assert cache.max_size == 100
    assert cache.ttl == 60
