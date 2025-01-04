"""Web data enrichment manager for pipeline.

This module provides the WebManager class that handles:
1. Web data collection
2. Social media monitoring
3. Community data integration
4. Literature mining
5. Patent analysis

The manager supports multiple data sources:
- ChEMBL database
- PubChem database
- Swiss* services
- Community sites (PsychonautWiki, Erowid, TripSit)
- Social media (Reddit, Twitter, Discord, Bluesky)
- Literature databases
- Patent databases
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
import json
from datetime import datetime

from ..models import CompoundData
from ..web_enrichment.data_sources import (
    ChEMBLClient,
    PubChemClient,
    SwissClient,
    CommunityDataClient,
    SocialDataClient,
    LiteratureClient,
)
from ..web_enrichment.http_client import HttpClient
from ..web_enrichment.llm_utils import analyze_content_with_llm


@dataclass
class WebConfig:
    """Web enrichment configuration."""
    
    # API credentials
    reddit_client_id: Optional[str] = None
    reddit_client_secret: Optional[str] = None
    twitter_api_key: Optional[str] = None
    twitter_api_secret: Optional[str] = None
    discord_token: Optional[str] = None
    bluesky_handle: Optional[str] = None
    bluesky_password: Optional[str] = None
    llm_api_key: Optional[str] = None
    
    # Enrichment settings
    max_retries: int = 3
    timeout: int = 30
    batch_size: int = 100
    use_cache: bool = True
    
    # Search settings
    search_terms: List[str] = field(default_factory=lambda: [
        "psychoactive",
        "hallucinogen",
        "antidepressant",
        "nootropic",
        "dissociative",
        "stimulant",
        "research chemical",
        "novel compound",
    ])
    
    # Target patterns
    target_patterns: Dict[str, str] = field(default_factory=lambda: {
        "5-HT2": r"5-HT2[A-C]?",
        "NMDA": r"NMDA",
        "nootropic": r"nootropic|cognitive enhancer",
        "psychedelic": r"psychedelic|hallucinogen",
        "dissociative": r"dissociative|NMDA",
        "stimulant": r"stimulant|dopamine",
    })


@dataclass
class EnrichmentStats:
    """Web enrichment statistics."""
    
    # Enrichment counts
    total_enrichments: int = 0
    successful_enrichments: int = 0
    failed_enrichments: int = 0
    cached_enrichments: int = 0
    
    # Source stats
    source_stats: Dict[str, Dict[str, int]] = field(default_factory=dict)
    
    # Error tracking
    errors: List[Dict[str, Any]] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "enrichments": {
                "total": self.total_enrichments,
                "successful": self.successful_enrichments,
                "failed": self.failed_enrichments,
                "cached": self.cached_enrichments,
                "success_rate": self._get_success_rate(),
            },
            "sources": self.source_stats,
            "errors": self.errors,
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get enrichment success rate."""
        if not self.total_enrichments:
            return None
        return self.successful_enrichments / self.total_enrichments


class WebManager:
    """Manager for web data enrichment."""

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[WebConfig] = None,
    ):
        """Initialize web manager.
        
        Args:
            cache_dir: Optional directory for caching
            config: Optional web configuration
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.config = config or WebConfig()
        
        # Initialize clients
        self._init_clients()
        
        # Initialize cache
        self._enrichment_cache: Dict[str, Dict] = {}
        
        # Initialize stats
        self.stats = EnrichmentStats()

    def _init_clients(self) -> None:
        """Initialize web clients."""
        try:
            # Initialize HTTP client
            self.http = HttpClient(
                max_retries=self.config.max_retries,
                timeout=self.config.timeout,
            )
            
            # Database clients
            self.chembl = ChEMBLClient(
                http_client=self.http,
                cache_dir=self.cache_dir,
            )
            self.pubchem = PubChemClient(
                http_client=self.http,
                cache_dir=self.cache_dir,
            )
            self.swiss = SwissClient(
                http_client=self.http,
                cache_dir=self.cache_dir,
            )
            
            # Community client
            self.community = CommunityDataClient(
                http_client=self.http,
                cache_dir=self.cache_dir,
            )
            
            # Social client
            self.social = SocialDataClient(
                http_client=self.http,
                cache_dir=self.cache_dir,
                reddit_client_id=self.config.reddit_client_id,
                reddit_client_secret=self.config.reddit_client_secret,
                twitter_api_key=self.config.twitter_api_key,
                twitter_api_secret=self.config.twitter_api_secret,
                discord_token=self.config.discord_token,
                bluesky_handle=self.config.bluesky_handle,
                bluesky_password=self.config.bluesky_password,
            )
            
            # Literature client
            self.literature = LiteratureClient(
                http_client=self.http,
                cache_dir=self.cache_dir,
            )
            
            self.logger.info("Successfully initialized web clients")
            
        except Exception as e:
            self.logger.error(f"Failed to initialize web clients: {str(e)}")
            raise

    def enrich_compound(
        self,
        compound: CompoundData,
        use_cache: Optional[bool] = None,
    ) -> CompoundData:
        """Enrich compound with web data.
        
        Args:
            compound: CompoundData instance to enrich
            use_cache: Whether to use enrichment cache
            
        Returns:
            CompoundData with web data
        """
        try:
            self.stats.total_enrichments += 1
            
            # Check cache
            if (use_cache if use_cache is not None else self.config.use_cache):
                cached = self._get_cached_enrichment(compound)
                if cached:
                    self.stats.cached_enrichments += 1
                    return self._apply_cached_enrichment(compound, cached)
            
            # Database enrichment
            self._enrich_from_databases(compound)
            
            # Community enrichment
            self._enrich_from_community(compound)
            
            # Social enrichment
            self._enrich_from_social(compound)
            
            # Literature enrichment
            self._enrich_from_literature(compound)
            
            # Cache enrichment
            if self.config.use_cache:
                self._cache_enrichment(compound)
            
            self.stats.successful_enrichments += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to enrich compound {compound.name}: {str(e)}"
            )
            self.stats.failed_enrichments += 1
            self.stats.errors.append({
                "type": "enrichment_error",
                "compound": compound.name,
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            })
            raise

    def _enrich_from_databases(self, compound: CompoundData) -> None:
        """Enrich compound from chemical databases."""
        try:
            # ChEMBL data
            chembl_data = self.chembl.get_compound_data(compound.smiles)
            if chembl_data:
                compound.merge_web_data(chembl_data)
                self._update_source_stats("chembl", True)
            
            # PubChem data
            pubchem_data = self.pubchem.get_compound_data(compound.smiles)
            if pubchem_data:
                compound.merge_web_data(pubchem_data)
                self._update_source_stats("pubchem", True)
            
            # Swiss data
            swiss_data = self.swiss.get_compound_data(compound.smiles)
            if swiss_data:
                compound.merge_web_data(swiss_data)
                self._update_source_stats("swiss", True)
                
        except Exception as e:
            self.logger.error(f"Database enrichment error: {str(e)}")
            self._update_source_stats("databases", False)
            raise

    def _enrich_from_community(self, compound: CompoundData) -> None:
        """Enrich compound from community sources."""
        try:
            community_data = self.community.get_compound_data(
                compound.name,
                search_terms=self.config.search_terms,
            )
            if community_data:
                # Analyze with LLM if available
                if self.config.llm_api_key:
                    community_data = analyze_content_with_llm(
                        community_data,
                        self.config.llm_api_key,
                    )
                compound.merge_web_data(community_data)
                self._update_source_stats("community", True)
                
        except Exception as e:
            self.logger.error(f"Community enrichment error: {str(e)}")
            self._update_source_stats("community", False)
            raise

    def _enrich_from_social(self, compound: CompoundData) -> None:
        """Enrich compound from social media."""
        try:
            social_data = self.social.get_compound_data(
                compound.name,
                search_terms=self.config.search_terms,
            )
            if social_data:
                # Analyze with LLM if available
                if self.config.llm_api_key:
                    social_data = analyze_content_with_llm(
                        social_data,
                        self.config.llm_api_key,
                    )
                compound.merge_web_data(social_data)
                self._update_source_stats("social", True)
                
        except Exception as e:
            self.logger.error(f"Social enrichment error: {str(e)}")
            self._update_source_stats("social", False)
            raise

    def _enrich_from_literature(self, compound: CompoundData) -> None:
        """Enrich compound from literature."""
        try:
            literature_data = self.literature.get_compound_data(
                compound.name,
                search_terms=self.config.search_terms,
            )
            if literature_data:
                # Analyze with LLM if available
                if self.config.llm_api_key:
                    literature_data = analyze_content_with_llm(
                        literature_data,
                        self.config.llm_api_key,
                    )
                compound.merge_web_data(literature_data)
                self._update_source_stats("literature", True)
                
        except Exception as e:
            self.logger.error(f"Literature enrichment error: {str(e)}")
            self._update_source_stats("literature", False)
            raise

    def _get_cached_enrichment(
        self,
        compound: CompoundData,
    ) -> Optional[Dict]:
        """Get cached enrichment for compound."""
        return self._enrichment_cache.get(compound.smiles)

    def _apply_cached_enrichment(
        self,
        compound: CompoundData,
        cached: Dict,
    ) -> CompoundData:
        """Apply cached enrichment to compound."""
        compound.merge_web_data(cached)
        return compound

    def _cache_enrichment(
        self,
        compound: CompoundData,
    ) -> None:
        """Cache enrichment for compound."""
        self._enrichment_cache[compound.smiles] = compound.get_web_data()

    def _update_source_stats(
        self,
        source: str,
        success: bool,
    ) -> None:
        """Update source statistics."""
        if source not in self.stats.source_stats:
            self.stats.source_stats[source] = {
                "total": 0,
                "successful": 0,
                "failed": 0,
            }
        
        stats = self.stats.source_stats[source]
        stats["total"] += 1
        if success:
            stats["successful"] += 1
        else:
            stats["failed"] += 1

    def get_client_info(self) -> Dict[str, Any]:
        """Get web client information."""
        return {
            "chembl": self.chembl.get_info(),
            "pubchem": self.pubchem.get_info(),
            "swiss": self.swiss.get_info(),
            "community": self.community.get_info(),
            "social": self.social.get_info(),
            "literature": self.literature.get_info(),
        }

    def clear_cache(self) -> None:
        """Clear enrichment cache."""
        self._enrichment_cache.clear()
        self.stats.cached_enrichments = 0
