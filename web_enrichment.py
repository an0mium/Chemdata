"""Web scraping and data enrichment from various sources.

This module provides functionality for:
1. Web scraping from multiple chemical databases
2. Data enrichment from regulatory sources
3. Chemical name and identifier extraction
4. Rate-limited and cached web requests
5. Content validation and parsing
"""

from web_enrichment import WebEnrichment

# Re-export main class
__all__ = ['WebEnrichment']
