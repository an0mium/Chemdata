"""FastAPI middleware implementations."""

from fastapi import Request, Response
from starlette.middleware.base import BaseHTTPMiddleware
from .rate_limiter import RateLimiter
from .cache import Cache


class RateLimitMiddleware(BaseHTTPMiddleware):
    """Rate limit middleware."""

    def __init__(self, app, limit: int = 100, window: int = 60):
        """Initialize middleware.

        Args:
            app: FastAPI application
            limit: Number of requests per window
            window: Window size in seconds
        """
        super().__init__(app)
        self.limiter = RateLimiter(
            key="api",
            limit=limit,
            window=window,
        )

    async def dispatch(self, request: Request, call_next):
        """Process request through rate limiter.

        Args:
            request: FastAPI request
            call_next: Next middleware/handler

        Returns:
            Response
        """
        if not await self.limiter.acquire():
            return Response(
                status_code=429,
                content="Rate limit exceeded",
                headers={
                    "retry-after": str(int(self.limiter.reset_time)),
                },
            )

        try:
            response = await call_next(request)
            response.headers["x-ratelimit-limit"] = str(self.limiter.limit)
            response.headers["x-ratelimit-remaining"] = str(self.limiter.remaining)
            response.headers["x-ratelimit-reset"] = str(int(self.limiter.reset_time))
            return response
        except:
            self.limiter.release()
            raise


class CacheMiddleware(BaseHTTPMiddleware):
    """Cache middleware."""

    def __init__(self, app, max_size: int = 1000, ttl: int = 300):
        """Initialize middleware.

        Args:
            app: FastAPI application
            max_size: Maximum cache size
            ttl: Cache TTL in seconds
        """
        super().__init__(app)
        self.cache = Cache(max_size=max_size, ttl=ttl)

    async def dispatch(self, request: Request, call_next):
        """Process request through cache.

        Args:
            request: FastAPI request
            call_next: Next middleware/handler

        Returns:
            Response
        """
        if request.method != "GET":
            return await call_next(request)

        cache_key = f"{request.method}:{request.url.path}:{request.query_params}"
        cached = await self.cache.get(cache_key)

        if cached is not None:
            return Response(
                content=cached["content"],
                status_code=cached["status_code"],
                headers={
                    **cached["headers"],
                    "x-cache": "HIT",
                },
            )

        response = await call_next(request)
        content = await response.body()

        await self.cache.set(
            cache_key,
            {
                "content": content,
                "status_code": response.status_code,
                "headers": dict(response.headers),
            },
        )

        return Response(
            content=content,
            status_code=response.status_code,
            headers={
                **response.headers,
                "x-cache": "MISS",
            },
        )
