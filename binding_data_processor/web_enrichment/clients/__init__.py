"""Web enrichment client implementations."""

from .base import WebClient
from .community import CommunityClient
from .swiss import SwissClient
from .reddit import RedditClient
from .bluelight import BluelightClient
from .pubmed import PubmedClient
from .sciencedirect import EnhancedScienceDirectClient as ScienceDirectClient
from .scholar import EnhancedScholarClient as ScholarClient
from .patents import EnhancedPatentClient as PatentsClient

__all__ = [
    "WebClient",
    "CommunityClient",
    "SwissClient",
    "RedditClient",
    "BluelightClient",
    "PubmedClient",
    "ScienceDirectClient",
    "ScholarClient",
    "PatentsClient",
]
