"""FastAPI application for compound data API."""

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError

from .compounds import CompoundAPI
from .search import SearchAPI
from .export import ExportAPI
from ..cache import Cache
from ..rate_limiter import RateLimiter
from ..middleware import CacheMiddleware, RateLimitMiddleware


def create_app() -> FastAPI:
    """Create FastAPI application.

    Returns:
        FastAPI application
    """
    # Create FastAPI app
    app = FastAPI(
        title="Compound Data API",
        description="API for accessing and searching compound data",
        version="1.0.0",
        docs_url="/api/docs",
        redoc_url="/api/redoc",
        openapi_url="/api/openapi.json",
    )

    # Add middleware (order matters - cache should be first to avoid caching rate limit headers)
    app.add_middleware(CacheMiddleware, max_size=1000, ttl=300)
    app.add_middleware(RateLimitMiddleware, limit=100, window=60)
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],  # TODO: Configure allowed origins
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Set up error handlers
    @app.exception_handler(RequestValidationError)
    async def validation_exception_handler(request: Request, exc: RequestValidationError):
        """Handle request validation errors.

        Args:
            request: FastAPI request
            exc: Validation exception

        Returns:
            JSON response with error details
        """
        return JSONResponse(
            status_code=422,
            content={
                "detail": exc.errors(),
                "body": exc.body,
            },
        )

    # Set up API routes
    compound_api = CompoundAPI()
    search_api = SearchAPI()
    export_api = ExportAPI()

    app.include_router(
        compound_api.router,
        prefix="/api",
        tags=["compounds"],
    )
    app.include_router(
        search_api.router,
        prefix="/api",
        tags=["search"],
    )
    app.include_router(
        export_api.router,
        prefix="/api",
        tags=["export"],
    )

    # Add health check endpoint
    @app.get("/api/health")
    async def health_check():
        """Health check endpoint.

        Returns:
            Health status
        """
        return {"status": "healthy"}

    # Add rate limit info endpoint
    @app.get("/api/rate-limit")
    async def rate_limit_info():
        """Get rate limit information.

        Returns:
            Rate limit status
        """
        rate_limiter = RateLimiter()
        return {
            "limit": rate_limiter.limit,
            "remaining": rate_limiter.remaining,
            "reset": rate_limiter.reset_time,
        }

    # Add cache info endpoint
    @app.get("/api/cache")
    async def cache_info():
        """Get cache information.

        Returns:
            Cache status
        """
        cache = Cache()
        return {
            "size": await cache.size(),
            "hits": await cache.hits(),
            "misses": await cache.misses(),
            "ttl": cache.ttl,
        }

    return app


# Create FastAPI application instance
app = create_app()

if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
