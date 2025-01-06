"""Rate limiting implementation.

This module provides the RateLimiter class that:
1. Enforces API rate limits
2. Implements token bucket algorithm
3. Handles burst limits
4. Tracks request counts
5. Provides backoff when limits are hit
"""

import time
import threading
from typing import Dict, Optional, Any
from dataclasses import dataclass
from datetime import datetime


@dataclass
class RateLimit:
    """Rate limit configuration."""

    requests: int  # Number of requests allowed
    period: int  # Time period in seconds
    burst: Optional[int] = None  # Optional burst limit


@dataclass
class RateLimitStats:
    """Rate limit statistics."""

    # Request counts
    total_requests: int = 0
    allowed_requests: int = 0
    throttled_requests: int = 0

    # Time tracking
    last_request: Optional[datetime] = None
    total_wait_time: float = 0.0

    # Rate tracking
    current_rate: float = 0.0
    peak_rate: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "requests": {
                "total": self.total_requests,
                "allowed": self.allowed_requests,
                "throttled": self.throttled_requests,
            },
            "timing": {
                "last_request": (self.last_request.isoformat() if self.last_request else None),
                "total_wait": self.total_wait_time,
            },
            "rates": {
                "current": self.current_rate,
                "peak": self.peak_rate,
            },
        }


class TokenBucket:
    """Token bucket rate limiter implementation."""

    def __init__(
        self,
        rate_limit: RateLimit,
    ):
        """Initialize token bucket.

        Args:
            rate_limit: Rate limit configuration
        """
        self.rate_limit = rate_limit
        self.tokens = rate_limit.requests
        self.last_update = time.time()
        self._lock = threading.Lock()

    def _add_tokens(self) -> None:
        """Add new tokens based on elapsed time."""
        now = time.time()
        elapsed = now - self.last_update

        new_tokens = elapsed * (self.rate_limit.requests / self.rate_limit.period)

        self.tokens = min(
            self.tokens + new_tokens,
            self.rate_limit.burst or self.rate_limit.requests,
        )
        self.last_update = now

    def get_token(self) -> float:
        """Get token from bucket.

        Returns:
            Wait time in seconds (0 if token available)
        """
        with self._lock:
            self._add_tokens()

            if self.tokens >= 1:
                self.tokens -= 1
                return 0.0

            # Calculate wait time
            tokens_needed = 1 - self.tokens
            wait_time = tokens_needed / (self.rate_limit.requests / self.rate_limit.period)

            return wait_time


class RateLimiter:
    """Rate limiter implementation."""

    def __init__(
        self,
        default_limit: Optional[RateLimit] = None,
    ):
        """Initialize rate limiter.

        Args:
            default_limit: Optional default rate limit
        """
        self.default_limit = default_limit or RateLimit(
            requests=30,
            period=60,
            burst=60,
        )

        # Initialize state
        self._buckets: Dict[str, TokenBucket] = {}
        self._limits: Dict[str, RateLimit] = {}
        self._stats: Dict[str, RateLimitStats] = {}
        self._lock = threading.Lock()

    def add_limit(
        self,
        key: str,
        limit: RateLimit,
    ) -> None:
        """Add rate limit for key.

        Args:
            key: Rate limit key (e.g. API endpoint)
            limit: Rate limit configuration
        """
        with self._lock:
            self._limits[key] = limit
            self._buckets[key] = TokenBucket(limit)
            self._stats[key] = RateLimitStats()

    def wait(
        self,
        key: Optional[str] = None,
    ) -> None:
        """Wait for rate limit if needed.

        Args:
            key: Optional rate limit key
        """
        stats = self._get_stats(key)
        stats.total_requests += 1

        bucket = self._get_bucket(key)
        wait_time = bucket.get_token()

        if wait_time > 0:
            stats.throttled_requests += 1
            stats.total_wait_time += wait_time
            time.sleep(wait_time)

        stats.allowed_requests += 1
        stats.last_request = datetime.now()

        # Update rate stats
        if stats.last_request:
            elapsed = (datetime.now() - stats.last_request).total_seconds()
            if elapsed > 0:
                stats.current_rate = 1 / elapsed
                stats.peak_rate = max(
                    stats.peak_rate,
                    stats.current_rate,
                )

    def _get_bucket(
        self,
        key: Optional[str] = None,
    ) -> TokenBucket:
        """Get token bucket for key.

        Args:
            key: Optional rate limit key

        Returns:
            Token bucket instance
        """
        with self._lock:
            if not key:
                # Use default bucket
                if "default" not in self._buckets:
                    self._buckets["default"] = TokenBucket(self.default_limit)
                return self._buckets["default"]

            if key not in self._buckets:
                # Create new bucket with default limit
                self._buckets[key] = TokenBucket(self._limits.get(key, self.default_limit))

            return self._buckets[key]

    def _get_stats(
        self,
        key: Optional[str] = None,
    ) -> RateLimitStats:
        """Get stats for key.

        Args:
            key: Optional rate limit key

        Returns:
            Rate limit statistics
        """
        with self._lock:
            key = key or "default"
            if key not in self._stats:
                self._stats[key] = RateLimitStats()
            return self._stats[key]

    def get_stats(
        self,
        key: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Get rate limit statistics.

        Args:
            key: Optional rate limit key

        Returns:
            Statistics dictionary
        """
        stats = self._get_stats(key)
        return {
            "key": key or "default",
            "limit": self._limits.get(
                key,
                self.default_limit,
            ).__dict__,
            "stats": stats.to_dict(),
        }

    def reset_stats(
        self,
        key: Optional[str] = None,
    ) -> None:
        """Reset statistics for key.

        Args:
            key: Optional rate limit key
        """
        with self._lock:
            if key:
                if key in self._stats:
                    self._stats[key] = RateLimitStats()
            else:
                # Reset all stats
                self._stats.clear()
