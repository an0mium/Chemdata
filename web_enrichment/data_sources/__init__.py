"""Data source modules for web enrichment.

This module provides clients for interacting with various chemical data sources:
1. PubChem - Chemical structure and property data
2. ChEMBL - Bioactivity data and drug targets
3. Regulatory - DEA, EMCDDA, WHO scheduling information
4. Community - PsychonautWiki, Erowid
5. Web Search - Generic web search and content analysis
"""

from .pubchem import PubChemClient
from .chembl import ChEMBLClient
from .regulatory import RegulatoryClient
from .community import CommunityClient
from .web_search import WebSearchClient

__all__ = [
    'PubChemClient',
    'ChEMBLClient',
    'RegulatoryClient',
    'CommunityClient',
    'WebSearchClient'
]
