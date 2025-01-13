"""HTTP client for web data enrichment with async support.

This module provides a robust HTTP client with:
- Async operations with aiohttp
- Rate limiting and exponential backoff retries
- Memory and disk caching
- User agent rotation
- Proxy support
- Progress tracking
- Special handling for problematic domains
- Circuit breaker pattern
"""

import logging
import time
import json
import asyncio
from pathlib import Path
from typing import Optional, Dict, Any, Union
from datetime import datetime, timedelta
from urllib.parse import urljoin

import aiohttp
from cachetools import TTLCache
from fake_useragent import UserAgent
from tqdm.asyncio import tqdm

from ..infrastructure.circuit_breaker import CircuitConfig
from binding_data_processor.logger import LogManager

logger = LogManager().get_logger("web_enrichment.http_client")


class HTTPClient:
    """HTTP client with advanced caching, rate limiting and retry capabilities."""

    def __init__(
        self,
        base_url: str = "",
        cache_dir: Optional[Path] = None,
        rate_limit: float = 1.0,  # Requests per second
        cache_ttl: int = 86400,  # 24 hours
        cache_size: int = 1000,
        proxies: Optional[Dict[str, str]] = None,
        circuit_config: Optional[CircuitConfig] = None,
    ):
        """Initialize HTTP client.

        Args:
            base_url: Base URL for all requests
            cache_dir: Optional directory for persistent cache
            rate_limit: Maximum requests per second
            cache_ttl: Cache TTL in seconds
            cache_size: Maximum cache size
            proxies: Optional proxy configuration
            circuit_config: Circuit breaker configuration
        """
        self.base_url = base_url
        self.cache_dir = cache_dir
        self.rate_limit = rate_limit
        self.circuit_config = circuit_config or CircuitConfig()
        self.last_request_time = 0.0
        self._session: Optional[aiohttp.ClientSession] = None
        self._special_sessions: Dict[str, aiohttp.ClientSession] = {}

        # Create cache directory
        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)

        # Initialize caches
        self.memory_cache = TTLCache(maxsize=cache_size, ttl=cache_ttl)
        self.user_agents = UserAgent()

        # Configure proxies
        self.proxies = proxies

    async def _get_session(self, domain: Optional[str] = None) -> aiohttp.ClientSession:
        """Get or create aiohttp session for domain."""
        if domain and domain in self._special_sessions:
            session = self._special_sessions[domain]
            if not session.closed:
                return session

        if domain == "erowid.org":
            # Special session for Erowid (SSL issues)
            session = aiohttp.ClientSession(headers={"User-Agent": self.user_agents.random}, connector=aiohttp.TCPConnector(verify_ssl=False))
            self._special_sessions[domain] = session
            return session

        if self._session is None or self._session.closed:
            # Configure timeout and proxy
            timeout = aiohttp.ClientTimeout(total=self.circuit_config.timeout)
            connector = aiohttp.TCPConnector(ssl=None)

            self._session = aiohttp.ClientSession(headers={"User-Agent": self.user_agents.random}, timeout=timeout, connector=connector)

            if self.proxies:
                self._session._proxy = self.proxies.get("http") or self.proxies.get("https")

        return self._session

    async def close(self):
        """Close all sessions."""
        if self._session and not self._session.closed:
            await self._session.close()
        for session in self._special_sessions.values():
            if not session.closed:
                await session.close()

    def _update_headers(self, session: aiohttp.ClientSession) -> None:
        """Update session headers with random user agent."""
        session._default_headers.update({"User-Agent": self.user_agents.random})

    async def get(
        self,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        use_cache: bool = True,
        cache_key: Optional[str] = None,
        show_progress: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Make GET request with caching and rate limiting.

        Args:
            url: Request URL (appended to base_url if relative)
            params: Query parameters
            headers: Request headers
            use_cache: Whether to use cache
            cache_key: Optional cache key override
            show_progress: Whether to show progress bar

        Returns:
            Response JSON data or None if failed
        """
        if not url.startswith(("http://", "https://")):
            url = urljoin(self.base_url, url)

        # Generate cache key
        if cache_key is None:
            cache_key = self._make_cache_key(url, params)

        # Check memory cache
        if use_cache:
            cached = self.memory_cache.get(cache_key)
            if cached:
                return cached

            # Check disk cache
            if self.cache_dir:
                disk_cached = await self._load_from_disk(cache_key)
                if disk_cached:
                    self.memory_cache[cache_key] = disk_cached
                    return disk_cached

        # Get appropriate session
        domain = url.split("/")[2] if "://" in url else None
        session = await self._get_session(domain)

        if headers:
            session._default_headers.update(headers)
        else:
            self._update_headers(session)

        # Set up progress tracking
        if show_progress:
            progress = tqdm(total=self.circuit_config.max_retries, desc=f"Requesting {url[:50]}...", unit="tries")

        # Apply rate limiting
        await self._wait_for_rate_limit()

        # Make request with retries
        for attempt in range(self.circuit_config.max_retries):
            if show_progress:
                progress.update(1)

            try:
                async with session.get(url, params=params) as response:
                    response.raise_for_status()
                    data = await response.json()

                    # Cache successful response
                    if use_cache:
                        self.memory_cache[cache_key] = data
                        if self.cache_dir:
                            await self._save_to_disk(cache_key, data)

                    if show_progress:
                        progress.close()
                    return data

            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                if attempt == self.circuit_config.max_retries - 1:
                    logger.error(f"Error making request to {url}: {str(e)}")
                    if show_progress:
                        progress.close()
                    return None
                await asyncio.sleep(self.circuit_config.retry_delay * (attempt + 1))

        return None

    async def post(
        self,
        url: str,
        json_data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        show_progress: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Make POST request with retry logic.

        Args:
            url: Request URL (appended to base_url if relative)
            json_data: JSON request body
            headers: Request headers
            show_progress: Whether to show progress bar

        Returns:
            Response JSON data or None if failed
        """
        if not url.startswith(("http://", "https://")):
            url = urljoin(self.base_url, url)

        # Get appropriate session
        domain = url.split("/")[2] if "://" in url else None
        session = await self._get_session(domain)

        if headers:
            session._default_headers.update(headers)
        else:
            self._update_headers(session)

        # Set up progress tracking
        if show_progress:
            progress = tqdm(total=self.circuit_config.max_retries, desc=f"Posting to {url[:50]}...", unit="tries")

        # Apply rate limiting
        await self._wait_for_rate_limit()

        # Make request with retries
        for attempt in range(self.circuit_config.max_retries):
            if show_progress:
                progress.update(1)

            try:
                async with session.post(url, json=json_data) as response:
                    response.raise_for_status()
                    data = await response.json()

                    if show_progress:
                        progress.close()
                    return data

            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                if attempt == self.circuit_config.max_retries - 1:
                    logger.error(f"Error posting to {url}: {str(e)}")
                    if show_progress:
                        progress.close()
                    return None
                await asyncio.sleep(self.circuit_config.retry_delay * (attempt + 1))

        return None

    async def _wait_for_rate_limit(self) -> None:
        """Wait to respect rate limit."""
        if self.rate_limit > 0:
            now = time.time()
            time_since_last = now - self.last_request_time
            if time_since_last < (1.0 / self.rate_limit):
                await asyncio.sleep((1.0 / self.rate_limit) - time_since_last)
            self.last_request_time = time.time()

    def _make_cache_key(self, url: str, params: Optional[Dict[str, Any]]) -> str:
        """Generate cache key from URL and params."""
        key = url
        if params:
            key += "_" + json.dumps(params, sort_keys=True)
        return key

    async def _load_from_disk(self, cache_key: str) -> Optional[Dict[str, Any]]:
        """Load cached response from disk."""
        try:
            cache_file = self.cache_dir / f"{cache_key}.json"
            if not cache_file.exists():
                return None

            async with aiohttp.ClientSession() as session:
                async with session.get(f"file://{cache_file}") as response:
                    data = await response.json()

            # Check TTL
            cached_time = datetime.fromisoformat(data["timestamp"])
            if datetime.now() - cached_time > timedelta(seconds=self.memory_cache.ttl):
                return None

            return data["content"]

        except Exception as e:
            logger.error(f"Error loading cache: {str(e)}")
            return None

    async def _save_to_disk(self, cache_key: str, content: Dict[str, Any]) -> None:
        """Save response to disk cache."""
        try:
            cache_file = self.cache_dir / f"{cache_key}.json"
            data = {"timestamp": datetime.now().isoformat(), "content": content}

            # Write file asynchronously
            async with aiofiles.open(cache_file, "w") as f:
                await f.write(json.dumps(data))

        except Exception as e:
            logger.error(f"Error saving cache: {str(e)}")

    async def clear_cache(self) -> None:
        """Clear all caches."""
        self.memory_cache.clear()
        if self.cache_dir:
            for cache_file in self.cache_dir.glob("*.json"):
                try:
                    cache_file.unlink()
                except Exception as e:
                    logger.error(f"Error deleting cache file: {str(e)}")
