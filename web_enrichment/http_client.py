"""HTTP client functionality for web enrichment.

This module handles:
1. Rate-limited HTTP requests
2. Request caching
3. Session management
4. SSL verification handling
"""

import json
from typing import Dict, Optional
import requests
from ratelimit import limits, sleep_and_retry

from logger import LogManager

logger = LogManager().get_logger("web_enrichment.http_client")


class HttpClient:
    """Handles HTTP requests with rate limiting and caching."""
    
    def __init__(self):
        """Initialize HTTP client."""
        # Main session for most requests
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'ChemDataCollector/0.1 (Research Project)'
        })
        
        # Separate session for Erowid (SSL issues)
        self.erowid_session = requests.Session()
        self.erowid_session.verify = False
        self.erowid_session.headers.update({
            'User-Agent': 'ChemDataCollector/0.1 (Research Project)'
        })
        
        # Cache for web requests
        self._cache = {}

    @sleep_and_retry
    @limits(calls=3, period=1)  # Rate limit: 3 requests per second
    def make_request(
        self,
        url: str,
        params: Optional[Dict] = None,
        verify: bool = True
    ) -> Optional[requests.Response]:
        """
        Make a rate-limited HTTP request with caching.
        
        Args:
            url: URL to request
            params: Optional query parameters
            verify: Whether to verify SSL certificates
            
        Returns:
            Response object or None if failed
        """
        # Generate cache key
        cache_key = f"{url}?{json.dumps(params or {})}"
        
        # Check cache first
        if cache_key in self._cache:
            return self._cache[cache_key]
            
        try:
            # Choose appropriate session
            session = self.erowid_session if 'erowid.org' in url else self.session
            
            # Make request
            response = session.get(
                url,
                params=params,
                timeout=10,
                verify=verify
            )
            response.raise_for_status()
            
            # Cache successful responses
            self._cache[cache_key] = response
            return response
            
        except Exception as e:
            logger.error(f"Error making request to {url}: {str(e)}")
            return None
