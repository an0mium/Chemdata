"""Tests for rate limiter implementation."""

import pytest
import asyncio
import time
from ..rate_limiter import RateLimiter


@pytest.fixture
async def rate_limiter():
    """Create test rate limiter instance."""
    limiter = RateLimiter(key="test", limit=10, window=1, burst=20)
    yield limiter
    # Clear singleton instance
    RateLimiter._instances.pop("test", None)


@pytest.mark.asyncio
async def test_singleton():
    """Test singleton pattern."""
    limiter1 = RateLimiter(key="singleton_test", limit=10)
    limiter2 = RateLimiter(key="singleton_test", limit=20)

    # Should be same instance
    assert limiter1 is limiter2
    assert limiter1.limit == 10  # First initialization wins

    # Clean up
    RateLimiter._instances.pop("singleton_test", None)


@pytest.mark.asyncio
async def test_acquire_release(rate_limiter):
    """Test basic token acquisition and release."""
    # Should be able to acquire tokens
    assert await rate_limiter.acquire(5)
    assert rate_limiter.remaining == 15

    # Release tokens
    rate_limiter.release(5)
    assert rate_limiter.remaining == 20


@pytest.mark.asyncio
async def test_burst_limit(rate_limiter):
    """Test burst limit handling."""
    # Try to acquire more than burst
    assert not await rate_limiter.acquire(21)
    assert rate_limiter.remaining == 20

    # Should be able to acquire up to burst
    assert await rate_limiter.acquire(20)
    assert rate_limiter.remaining == 0


@pytest.mark.asyncio
async def test_token_refill():
    """Test token refill behavior."""
    limiter = RateLimiter(key="refill_test", limit=10, window=1, burst=10)

    # Use all tokens
    assert await limiter.acquire(10)
    assert limiter.remaining == 0

    # Wait for refill
    await asyncio.sleep(1.1)
    await limiter._add_tokens()
    assert limiter.remaining == 10

    # Clean up
    RateLimiter._instances.pop("refill_test", None)


@pytest.mark.asyncio
async def test_concurrent_access():
    """Test concurrent access handling."""
    limiter = RateLimiter(key="concurrent_test", limit=100, window=1, burst=100)
    results = []

    async def worker():
        for _ in range(20):
            if await limiter.acquire(1):
                results.append(1)
            await asyncio.sleep(0.001)

    # Run concurrent workers
    workers = [worker() for _ in range(5)]
    await asyncio.gather(*workers)

    # Should have acquired exactly 100 tokens
    assert sum(results) <= 100

    # Clean up
    RateLimiter._instances.pop("concurrent_test", None)


@pytest.mark.asyncio
async def test_properties(rate_limiter):
    """Test rate limiter properties."""
    assert rate_limiter.limit == 10
    assert rate_limiter.window == 1
    assert rate_limiter.burst == 20
    assert rate_limiter.remaining == 20


@pytest.mark.asyncio
async def test_reset_time(rate_limiter):
    """Test reset time calculation."""
    # Use half the tokens
    await rate_limiter.acquire(10)

    # Should take 1 second to refill 10 tokens at rate of 10/second
    reset_time = rate_limiter.reset_time
    assert 0.9 <= reset_time <= 1.1

    # Use all tokens
    await rate_limiter.acquire(10)
    reset_time = rate_limiter.reset_time
    assert 1.9 <= reset_time <= 2.1


@pytest.mark.asyncio
async def test_multiple_instances():
    """Test multiple rate limiter instances."""
    limiter1 = RateLimiter(key="test1", limit=10)
    limiter2 = RateLimiter(key="test2", limit=20)

    assert await limiter1.acquire(5)
    assert await limiter2.acquire(10)

    assert limiter1.remaining == 5
    assert limiter2.remaining == 10

    # Clean up
    RateLimiter._instances.pop("test1", None)
    RateLimiter._instances.pop("test2", None)


@pytest.mark.asyncio
async def test_zero_window():
    """Test handling of zero window."""
    with pytest.raises(ValueError):
        RateLimiter(key="zero_test", limit=10, window=0)


@pytest.mark.asyncio
async def test_negative_values():
    """Test handling of negative values."""
    with pytest.raises(ValueError):
        RateLimiter(key="negative_test", limit=-10)

    with pytest.raises(ValueError):
        RateLimiter(key="negative_test", limit=10, window=-1)

    with pytest.raises(ValueError):
        RateLimiter(key="negative_test", limit=10, burst=-20)


@pytest.mark.asyncio
async def test_burst_defaults_to_limit():
    """Test burst defaulting to limit."""
    limiter = RateLimiter(key="burst_test", limit=10)
    assert limiter.burst == 10
    RateLimiter._instances.pop("burst_test", None)
