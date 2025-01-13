"""Rate limiter implementation for web components."""

import time
import asyncio
from typing import Dict, Optional


class RateLimiter:
    """Token bucket rate limiter."""

    _instances: Dict[str, "RateLimiter"] = {}
    _lock = asyncio.Lock()

    def __new__(cls, key: str = "default") -> "RateLimiter":
        """Get or create rate limiter instance.

        Args:
            key: Rate limiter instance key

        Returns:
            Rate limiter instance
        """
        if key not in cls._instances:
            cls._instances[key] = super().__new__(cls)
        return cls._instances[key]

    def __init__(
        self,
        key: str = "default",
        limit: int = 100,
        window: int = 60,
        burst: Optional[int] = None,
    ):
        """Initialize rate limiter.

        Args:
            key: Rate limiter instance key
            limit: Number of tokens per window
            window: Window size in seconds
            burst: Optional burst size (defaults to limit)

        Raises:
            ValueError: If limit, window or burst are invalid
        """
        if limit <= 0:
            raise ValueError("Limit must be positive")
        if window <= 0:
            raise ValueError("Window must be positive")
        if burst is not None and burst <= 0:
            raise ValueError("Burst must be positive")

        # Skip initialization if already initialized
        if hasattr(self, "_initialized"):
            return

        self._initialized = True
        self._key = key
        self._limit = limit
        self._window = window
        self._burst = burst or limit
        self._tokens = self._burst
        self._last_update = time.time()
        self._lock = asyncio.Lock()

    async def acquire(self, tokens: int = 1) -> bool:
        """Acquire tokens from bucket.

        Args:
            tokens: Number of tokens to acquire

        Returns:
            True if tokens were acquired, False otherwise
        """
        async with self._lock:
            await self._add_tokens()

            if tokens > self._burst:
                return False

            if self._tokens >= tokens:
                self._tokens -= tokens
                return True

            return False

    async def _add_tokens(self) -> None:
        """Add tokens based on time elapsed."""
        now = time.time()
        elapsed = now - self._last_update
        tokens_to_add = int((elapsed * self._limit) / self._window)

        if tokens_to_add > 0:
            self._tokens = min(self._burst, self._tokens + tokens_to_add)
            self._last_update = now

    def release(self, tokens: int = 1) -> None:
        """Release tokens back to bucket.

        Args:
            tokens: Number of tokens to release
        """
        self._tokens = min(self._burst, self._tokens + tokens)

    @property
    def limit(self) -> int:
        """Get rate limit.

        Returns:
            Rate limit per window
        """
        return self._limit

    @property
    def window(self) -> int:
        """Get rate limit window.

        Returns:
            Window size in seconds
        """
        return self._window

    @property
    def burst(self) -> int:
        """Get burst size.

        Returns:
            Maximum burst size
        """
        return self._burst

    @property
    def remaining(self) -> int:
        """Get remaining tokens.

        Returns:
            Number of tokens remaining
        """
        return self._tokens

    @property
    def reset_time(self) -> float:
        """Get time until next token.

        Returns:
            Seconds until next token
        """
        if self._tokens >= self._burst:
            return 0.0

        tokens_needed = self._burst - self._tokens
        return (tokens_needed * self._window) / self._limit
