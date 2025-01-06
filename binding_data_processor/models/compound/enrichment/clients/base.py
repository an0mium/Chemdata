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

import logging
import re
from pathlib import Path
from typing import Optional, Dict, Any, List
from abc import ABC, abstractmethod

from ....pipeline.infrastructure.circuit_breaker import CircuitConfig
from ..validation import ValidationError
from .http import HTTPClient
from ...base.core import Compound


class WebClientError(Exception):
    """Base exception for web client errors."""
    pass


class RateLimitError(WebClientError):
    """Raised when rate limit is exceeded."""
    pass


class BaseWebClient(ABC):
    """Base client for web data sources."""

    def __init__(
        self,
        name: str,
        http_client: Optional[HTTPClient] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        circuit_config: Optional[CircuitConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize base web client.
        
        Args:
            name: Client name for circuit breaker
            http_client: Optional HTTP client to use
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            circuit_config: Optional circuit breaker configuration
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        self.model_dir = model_dir
        self.cache_dir = cache_dir

        # Create or use HTTP client
        if http_client:
            self.http = http_client
        else:
            self.http = HTTPClient(
                name=name,
                cache_dir=cache_dir / "http" if cache_dir else None,
                circuit_config=circuit_config,
                logger=self.logger,
            )

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
            raise ValidationError(
                f"Missing required fields: {', '.join(missing)}"
            )

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

    def get_metrics(self) -> Dict[str, Any]:
        """Get client metrics including circuit breaker state."""
        return {
            "http": self.http.get_metrics(),
        }

    def close(self) -> None:
        """Close client and cleanup."""
        self.http.close()
