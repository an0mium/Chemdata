"""Web enrichment clients for compound data.

This module provides clients for enriching compound data from various web sources:
- Base client with common functionality
- HTTP client for web requests
- Swiss tools (SwissTargetPrediction, SwissADME)
- Community sources (PsychonautWiki, Erowid, TripSit)
- Social media (Reddit, Twitter)
"""

from .base import (
    WebClientError,
    RateLimitError,
    ValidationError,
    BaseWebClient,
)
from .http import HTTPClient
from .swiss import SwissClient
from .community import CommunityClient
from .social import SocialClient

__all__ = [
    # Base client
    "BaseWebClient",
    
    # HTTP client
    "HTTPClient",
    
    # Exceptions
    "WebClientError",
    "RateLimitError", 
    "ValidationError",
    
    # Specialized clients
    "SwissClient",
    "CommunityClient",
    "SocialClient",
]
