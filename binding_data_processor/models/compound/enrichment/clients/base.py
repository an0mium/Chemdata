"""Base client for web data sources.

This module provides a base client class that specific data source clients
can inherit from to get common functionality:
- HTTP client management
- Caching
- Rate limiting
- Error handling
- Data validation
- Circuit breaker pattern
"""

from ....web_enrichment.clients.base import WebClient as BaseWebClient
from ....web_enrichment.clients.base import ValidationError, WebClientError, RateLimitError

__all__ = [
    "BaseWebClient",
    "WebClientError",
    "ValidationError",
    "RateLimitError",
]
