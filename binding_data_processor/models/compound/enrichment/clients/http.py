"""HTTP client for web data enrichment.

This module provides a robust HTTP client with:
- Circuit breaker pattern
- Rate limiting
- Retries with exponential backoff
- Caching
- User agent rotation
- Proxy support
"""

import logging
import time
import json
from pathlib import Path
from typing import Optional, Dict, Any, TypeVar, Callable
from datetime import datetime, timedelta

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from cachetools import TTLCache
from fake_useragent import UserAgent

from ....pipeline.infrastructure.circuit_breaker import (
    CircuitBreaker,
    CircuitConfig,
)

T = TypeVar("T")  # Generic type for circuit breaker return value


class HTTPClient:
    """HTTP client with circuit breaker and caching."""

    def __init__(
        self,
        name: str,
        cache_dir: Optional[Path] = None,
        rate_limit: float = 1.0,  # Requests per second
        max_retries: int = 3,
        timeout: float = 30.0,
        cache_ttl: int = 86400,  # 24 hours
        cache_size: int = 1000,
        proxies: Optional[Dict[str, str]] = None,
        circuit_config: Optional[CircuitConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize HTTP client.
        
        Args:
            name: Client name for circuit breaker
            cache_dir: Optional directory for persistent cache
            rate_limit: Maximum requests per second
            max_retries: Maximum number of retries
            timeout: Request timeout in seconds
            cache_ttl: Cache TTL in seconds
            cache_size: Maximum cache size
            proxies: Optional proxy configuration
            circuit_config: Optional circuit breaker configuration
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        self.cache_dir = cache_dir
        self.rate_limit = rate_limit
        self.timeout = timeout
        self.last_request_time = 0.0

        # Create cache directory
        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)

        # Initialize session
        self.session = requests.Session()

        # Configure retries
        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=0.5,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

        # Configure proxies
        if proxies:
            self.session.proxies.update(proxies)

        # Initialize caches
        self.memory_cache = TTLCache(maxsize=cache_size, ttl=cache_ttl)
        self.user_agents = UserAgent()

        # Initialize circuit breaker
        self.circuit = CircuitBreaker(
            name=f"http_client.{name}",
            config=circuit_config,
            logger=self.logger,
        )

    def get(
        self,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        use_cache: bool = True,
        cache_key: Optional[str] = None,
        fallback: Optional[Callable[..., requests.Response]] = None,
    ) -> requests.Response:
        """Make GET request with circuit breaker and caching.
        
        Args:
            url: Request URL
            params: Optional query parameters
            headers: Optional request headers
            use_cache: Whether to use cache
            cache_key: Optional cache key override
            fallback: Optional fallback function if circuit is open
            
        Returns:
            Response object
            
        Raises:
            requests.RequestException: For request errors
        """
        # Check cache
        if use_cache:
            cached = self._check_cache(cache_key or self._make_cache_key(url, params))
            if cached:
                return cached

        # Prepare request
        headers = self._prepare_headers(headers)
        self._wait_for_rate_limit()

        # Make request
        try:
            response = self._execute_request(
                lambda: self.session.get(
                    url,
                    params=params,
                    headers=headers,
                    timeout=self.timeout,
                ),
                fallback,
            )

            # Cache successful response
            if use_cache:
                self._cache_response(
                    cache_key or self._make_cache_key(url, params),
                    response,
                )

            return response

        except Exception as e:
            self.logger.error(f"Error fetching {url}: {str(e)}")
            raise

    def post(
        self,
        url: str,
        data: Optional[Dict[str, Any]] = None,
        json_data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        fallback: Optional[Callable[..., requests.Response]] = None,
    ) -> requests.Response:
        """Make POST request with circuit breaker.
        
        Args:
            url: Request URL
            data: Optional form data
            json_data: Optional JSON data
            headers: Optional request headers
            fallback: Optional fallback function if circuit is open
            
        Returns:
            Response object
            
        Raises:
            requests.RequestException: For request errors
        """
        # Prepare request
        headers = self._prepare_headers(headers)
        self._wait_for_rate_limit()

        # Make request
        try:
            response = self._execute_request(
                lambda: self.session.post(
                    url,
                    data=data,
                    json=json_data,
                    headers=headers,
                    timeout=self.timeout,
                ),
                fallback,
            )
            return response

        except Exception as e:
            self.logger.error(f"Error posting to {url}: {str(e)}")
            raise

    def _check_cache(self, cache_key: str) -> Optional[requests.Response]:
        """Check memory and disk cache for response."""
        cached = self.memory_cache.get(cache_key)
        if cached:
            return cached

        if self.cache_dir:
            disk_cached = self._load_from_disk(cache_key)
            if disk_cached:
                self.memory_cache[cache_key] = disk_cached
                return disk_cached

        return None

    def _prepare_headers(self, headers: Optional[Dict[str, str]]) -> Dict[str, str]:
        """Prepare request headers with user agent."""
        headers = headers or {}
        headers["User-Agent"] = self.user_agents.random
        return headers

    def _execute_request(
        self,
        request_fn: Callable[[], requests.Response],
        fallback: Optional[Callable[..., requests.Response]] = None,
    ) -> requests.Response:
        """Execute request with circuit breaker."""
        response = self.circuit.execute(request_fn, fallback=fallback)
        response.raise_for_status()
        return response

    def _cache_response(self, cache_key: str, response: requests.Response) -> None:
        """Cache response in memory and disk."""
        self.memory_cache[cache_key] = response
        if self.cache_dir:
            self._save_to_disk(cache_key, response)

    def _wait_for_rate_limit(self) -> None:
        """Wait to respect rate limit."""
        if self.rate_limit > 0:
            now = time.time()
            time_since_last = now - self.last_request_time
            if time_since_last < (1.0 / self.rate_limit):
                time.sleep((1.0 / self.rate_limit) - time_since_last)
            self.last_request_time = time.time()

    def _make_cache_key(self, url: str, params: Optional[Dict[str, Any]]) -> str:
        """Generate cache key from URL and params."""
        key = url
        if params:
            key += "_" + json.dumps(params, sort_keys=True)
        return key

    def _load_from_disk(self, cache_key: str) -> Optional[requests.Response]:
        """Load cached response from disk."""
        try:
            cache_file = self.cache_dir / f"{cache_key}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                data = json.load(f)

            # Check TTL
            cached_time = datetime.fromisoformat(data["timestamp"])
            if datetime.now() - cached_time > timedelta(seconds=self.memory_cache.ttl):
                return None

            # Reconstruct response
            response = requests.Response()
            response.status_code = data["status_code"]
            response._content = json.dumps(data["content"]).encode()
            response.headers.update(data["headers"])
            return response

        except Exception as e:
            self.logger.error(f"Error loading cache: {str(e)}")
            return None

    def _save_to_disk(self, cache_key: str, response: requests.Response) -> None:
        """Save response to disk cache."""
        try:
            cache_file = self.cache_dir / f"{cache_key}.json"
            data = {
                "timestamp": datetime.now().isoformat(),
                "status_code": response.status_code,
                "content": response.json(),
                "headers": dict(response.headers),
            }
            with cache_file.open("w") as f:
                json.dump(data, f)

        except Exception as e:
            self.logger.error(f"Error saving cache: {str(e)}")

    def clear_cache(self) -> None:
        """Clear all caches."""
        self.memory_cache.clear()
        if self.cache_dir:
            for cache_file in self.cache_dir.glob("*.json"):
                try:
                    cache_file.unlink()
                except Exception as e:
                    self.logger.error(f"Error deleting cache file: {str(e)}")

    def get_metrics(self) -> Dict[str, Any]:
        """Get client metrics including circuit breaker state."""
        return {
            "circuit_breaker": self.circuit.get_metrics(),
            "cache": {
                "memory_size": len(self.memory_cache),
                "disk_size": len(list(self.cache_dir.glob("*.json")))
                if self.cache_dir else 0,
            },
        }

    def close(self) -> None:
        """Close session and cleanup."""
        self.session.close()
        self.clear_cache()
