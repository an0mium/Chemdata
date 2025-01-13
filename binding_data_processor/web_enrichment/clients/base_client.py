"""Base client for web data sources.

This module provides a base client class that specific data source clients
can inherit from to get common functionality:
- HTTP client management
- Caching
- Rate limiting
- Error handling
- Data validation
- Compound processing
- Text extraction and cleaning
"""

import logging
import re
from pathlib import Path
from typing import Optional, Dict, Any, List
from abc import ABC, abstractmethod

from ..http_client_enhanced import HTTPClientEnhanced as HTTPClient
from ...models.compound import Compound

logger = logging.getLogger(__name__)


class WebClientError(Exception):
    """Base exception for web client errors."""

    pass


class RateLimitError(WebClientError):
    """Raised when rate limit is exceeded."""

    pass


class ValidationError(WebClientError):
    """Raised when data validation fails."""

    pass


class BaseWebClient(ABC):
    """Base client for web data sources."""

    def __init__(
        self,
        http_client: Optional[HTTPClient] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
        base_url: Optional[str] = None,
        api_key: Optional[str] = None,
    ):
        """Initialize base web client.

        Args:
            http_client: Optional HTTP client to use
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            logger: Optional logger instance
            base_url: Optional base URL for API requests
            api_key: Optional API key for authentication
        """
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        self.model_dir = model_dir
        self.cache_dir = cache_dir
        self.base_url = base_url.rstrip("/") if base_url else None
        self.api_key = api_key

        # Create or use HTTP client
        if http_client:
            self.http = http_client
        else:
            self.http = HTTPClient(name=self.__class__.__name__, cache_dir=cache_dir / "http" if cache_dir else None, logger=self.logger, base_url=base_url, api_key=api_key)

    @abstractmethod
    def process_compounds(
        self,
        compounds: List[Compound],
        skip_predictions: bool = False,
        use_cache: bool = True,
    ) -> None:
        """Process list of compounds.

        Args:
            compounds: List of compounds to process
            skip_predictions: Whether to skip ML predictions
            use_cache: Whether to use cached results
        """
        pass

    @abstractmethod
    def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get data for a single compound.

        Args:
            name: Compound name
            cas_number: Optional CAS number
            use_cache: Whether to use cached results

        Returns:
            Dictionary of compound data or None if not found
        """
        pass

    def _validate_response(
        self,
        response: Dict[str, Any],
        required_fields: List[str],
    ) -> None:
        """Validate response data.

        Args:
            response: Response data to validate
            required_fields: List of required field names

        Raises:
            ValidationError if validation fails
        """
        missing = []
        for field in required_fields:
            if field not in response:
                missing.append(field)
            elif response[field] is None:
                missing.append(field)

        if missing:
            raise ValidationError(f"Missing required fields: {', '.join(missing)}")

    def _clean_text(self, text: str) -> str:
        """Clean text data.

        Args:
            text: Text to clean

        Returns:
            Cleaned text
        """
        if not text:
            return ""

        # Remove extra whitespace
        text = " ".join(text.split())

        # Remove non-printable characters
        text = "".join(char for char in text if char.isprintable())

        return text

    def _extract_cas(self, text: str) -> Optional[str]:
        """Extract CAS number from text.

        Args:
            text: Text to extract from

        Returns:
            CAS number if found, None otherwise
        """
        # CAS number pattern
        pattern = r"\b\d{1,7}-\d{2}-\d\b"

        match = re.search(pattern, text)
        if match:
            return match.group(0)

        return None

    def _extract_doi(self, text: str) -> Optional[str]:
        """Extract DOI from text.

        Args:
            text: Text to extract from

        Returns:
            DOI if found, None otherwise
        """
        # DOI pattern
        pattern = r"\b10\.\d{4,}/[-._;()/:\w]+\b"

        match = re.search(pattern, text)
        if match:
            return match.group(0)

        return None

    def get(self, endpoint: str, **kwargs) -> Dict[str, Any]:
        """Make a GET request.

        Args:
            endpoint: API endpoint path
            **kwargs: Additional arguments to pass to HTTP client

        Returns:
            Parsed JSON response
        """
        return self.http.get(endpoint, **kwargs)

    def post(self, endpoint: str, **kwargs) -> Dict[str, Any]:
        """Make a POST request.

        Args:
            endpoint: API endpoint path
            **kwargs: Additional arguments to pass to HTTP client

        Returns:
            Parsed JSON response
        """
        return self.http.post(endpoint, **kwargs)

    def put(self, endpoint: str, **kwargs) -> Dict[str, Any]:
        """Make a PUT request.

        Args:
            endpoint: API endpoint path
            **kwargs: Additional arguments to pass to HTTP client

        Returns:
            Parsed JSON response
        """
        return self.http.put(endpoint, **kwargs)

    def delete(self, endpoint: str, **kwargs) -> Dict[str, Any]:
        """Make a DELETE request.

        Args:
            endpoint: API endpoint path
            **kwargs: Additional arguments to pass to HTTP client

        Returns:
            Parsed JSON response
        """
        return self.http.delete(endpoint, **kwargs)

    def close(self) -> None:
        """Close client and cleanup."""
        self.http.close()
