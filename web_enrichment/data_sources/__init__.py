"""Data source modules for web enrichment.

This module provides clients for interacting with various chemical data sources:
1. PubChem - Chemical structure and property data
2. ChEMBL - Bioactivity data and drug targets
3. Regulatory - DEA, EMCDDA, WHO scheduling information
4. Community - PsychonautWiki, Erowid
5. Web Search - Generic web search and content analysis
6. Swiss* Services:
   - SwissTargetPrediction - Target prediction
   - SwissADME - ADME properties and pharmacokinetics
   - SwissSimilarity - Chemical similarity search
"""

from .pubchem import PubChemClient
from .chembl import ChEMBLClient
from .regulatory import RegulatoryClient
from .community import CommunityClient
from .web_search import WebSearchClient
from .swiss import SwissClient

__all__ = [
    'PubChemClient',
    'ChEMBLClient',
    'RegulatoryClient',
    'CommunityClient',
    'WebSearchClient',
    'SwissClient'
]
