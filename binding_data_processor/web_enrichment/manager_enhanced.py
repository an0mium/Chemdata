"""Enhanced web enrichment manager with Crawl4AI integration.

This module provides an enhanced manager class that coordinates web enrichment clients
with improved error handling, resilience through circuit breakers, and Crawl4AI integration:
- Swiss tools (SwissTargetPrediction, SwissADME)
- Community sources (PsychonautWiki, Erowid, TripSit)
- Social media (Reddit, Twitter)
- Research papers and patents (via Crawl4AI)

The manager handles:
1. Client initialization and configuration
2. Coordinated data enrichment
3. Error handling and recovery with circuit breakers
4. Progress tracking and metrics
5. Result aggregation
"""

import logging
import asyncio
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime

from tqdm import tqdm

from .http_client_enhanced import HTTPClientEnhanced
from .swiss_client import SwissClient
from .clients.community import CommunityClient
from .social_client import SocialClient
from .crawl4ai_client import Crawl4AIClient, ResearchData, PatentData, CommunityData
from ..pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig


@dataclass
class EnrichmentConfig:
    """Configuration for web enrichment."""

    # API credentials
    reddit_client_id: Optional[str] = None
    reddit_client_secret: Optional[str] = None
    twitter_bearer_token: Optional[str] = None

    # Processing options
    skip_predictions: bool = False
    skip_web_data: bool = False
    use_cache: bool = True
    n_workers: int = 4
    batch_size: int = 100

    # Directories
    model_dir: Optional[Path] = None
    cache_dir: Optional[Path] = None

    # Circuit breaker settings
    circuit_config: Optional[CircuitBreakerConfig] = None


@dataclass
class EnrichmentResult:
    """Web enrichment result."""

    research_papers: List[ResearchData]
    patents: List[PatentData]
    community_posts: List[CommunityData]
    swiss_data: Dict[str, Any]
    social_data: Dict[str, Any]
    metadata: Dict[str, Any]


class WebEnrichmentManagerEnhanced:
    """Enhanced web enrichment manager."""

    def __init__(
        self,
        config: EnrichmentConfig,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize web enrichment manager.

        Args:
            config: Enrichment configuration
            logger: Optional logger instance
        """
        self.config = config
        self.logger = logger or logging.getLogger(self.__class__.__name__)

        # Initialize shared HTTP client
        self.http = HTTPClientEnhanced(
            name="web_enrichment",
            cache_dir=config.cache_dir / "http" if config.cache_dir else None,
            circuit_config=config.circuit_config,
            logger=self.logger,
        )

        # Initialize clients
        self.swiss_client = SwissClient(
            http_client=self.http,
            model_dir=config.model_dir,
            cache_dir=config.cache_dir / "swiss" if config.cache_dir else None,
            logger=self.logger,
        )

        self.community_client = CommunityClient(
            http_client=self.http,
            model_dir=config.model_dir,
            cache_dir=config.cache_dir / "community" if config.cache_dir else None,
            logger=self.logger,
        )

        self.social_client = SocialClient(
            reddit_client_id=config.reddit_client_id,
            reddit_client_secret=config.reddit_client_secret,
            twitter_bearer_token=config.twitter_bearer_token,
            http_client=self.http,
            model_dir=config.model_dir,
            cache_dir=config.cache_dir / "social" if config.cache_dir else None,
            logger=self.logger,
        )

        self.crawl4ai_client = Crawl4AIClient()

        # Initialize thread pool
        self.executor = ThreadPoolExecutor(
            max_workers=config.n_workers,
            thread_name_prefix="enrichment",
        )

    def _get_config_defaults(
        self,
        skip_predictions: Optional[bool],
        skip_web_data: Optional[bool],
        use_cache: Optional[bool],
    ) -> Tuple[bool, bool, bool]:
        """Get config defaults for optional parameters."""
        if skip_predictions is None:
            skip_pred = self.config.skip_predictions
        else:
            skip_pred = skip_predictions

        if skip_web_data is None:
            skip_web = self.config.skip_web_data
        else:
            skip_web = skip_web_data

        if use_cache is None:
            use_cache_val = self.config.use_cache
        else:
            use_cache_val = use_cache

        return skip_pred, skip_web, use_cache_val

    def _build_tasks(
        self,
        query: str,
        name: str,
        skip_web_data: bool,
        skip_predictions: bool,
    ) -> List[asyncio.Task]:
        """Build list of enrichment tasks."""
        tasks = []
        if not skip_web_data:
            tasks.extend(
                [
                    self.crawl4ai_client.scrape_research(query),
                    self.crawl4ai_client.scrape_patents(query),
                    self.crawl4ai_client.scrape_community(query),
                    self.social_client.get_compound_mentions(name),
                ]
            )
        if not skip_predictions:
            tasks.append(self.swiss_client.get_compound_data(name))
        return tasks

    def _parse_results(
        self,
        results: List[Any],
        skip_web_data: bool,
        skip_predictions: bool,
    ) -> Tuple[List[Any], List[Any], List[Any], Dict[str, Any], Dict[str, Any]]:
        """Parse enrichment results."""
        research_papers = []
        patents = []
        community_posts = []
        swiss_data = {}
        social_data = {}

        for i, result in enumerate(results):
            if isinstance(result, Exception):
                self.logger.error(f"Error in task {i}: {str(result)}")
                continue

            if not skip_web_data:
                if i == 0:
                    research_papers = result
                elif i == 1:
                    patents = result
                elif i == 2:
                    community_posts = result
                elif i == 3:
                    social_data = result
            elif not skip_predictions and i == 0:
                swiss_data = result

        return research_papers, patents, community_posts, swiss_data, social_data

    def _build_metadata(
        self,
        query: str,
        use_cache: bool,
        research_papers: List[Any],
        patents: List[Any],
        community_posts: List[Any],
        swiss_data: Dict[str, Any],
        social_data: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Build metadata for enrichment result."""
        metadata = {
            "sources": [],
            "query": query,
            "timestamp": datetime.now().isoformat(),
            "cache_used": use_cache,
        }
        if research_papers:
            metadata["sources"].append("research_papers")
        if patents:
            metadata["sources"].append("patents")
        if community_posts:
            metadata["sources"].append("community_posts")
        if swiss_data:
            metadata["sources"].append("swiss_data")
        if social_data:
            metadata["sources"].append("social_data")
        return metadata

    async def enrich_compound(
        self,
        name: str,
        smiles: Optional[str] = None,
        skip_predictions: Optional[bool] = None,
        skip_web_data: Optional[bool] = None,
        use_cache: Optional[bool] = None,
    ) -> EnrichmentResult:
        """Enrich compound with web data.

        Args:
            name: Compound name
            smiles: Optional SMILES structure
            skip_predictions: Whether to skip predictions (overrides config)
            skip_web_data: Whether to skip web data (overrides config)
            use_cache: Whether to use cached results (overrides config)

        Returns:
            Enrichment result with data from all sources
        """
        try:
            # Get config defaults
            config_args = (skip_predictions, skip_web_data, use_cache)
            defaults = self._get_config_defaults(*config_args)
            skip_pred, skip_web, use_cache_val = defaults

            # Build search query
            query = name if not smiles else f"{name} {smiles}"

            # Run enrichment tasks
            tasks = self._build_tasks(query, name, skip_web, skip_pred)
            results = await asyncio.gather(*tasks, return_exceptions=True)

            # Parse results
            result_data = self._parse_results(results, skip_web, skip_pred)
            (
                research_papers,
                patents,
                community_posts,
                swiss_data,
                social_data,
            ) = result_data

            # Build metadata
            metadata = self._build_metadata(
                query,
                use_cache_val,
                research_papers,
                patents,
                community_posts,
                swiss_data,
                social_data,
            )

            return EnrichmentResult(
                research_papers=research_papers,
                patents=patents,
                community_posts=community_posts,
                swiss_data=swiss_data,
                social_data=social_data,
                metadata=metadata,
            )

        except Exception as e:
            self.logger.error(f"Error enriching {name}: {str(e)}")
            raise

    async def enrich_compounds(
        self,
        compounds: List[Dict[str, str]],
        skip_predictions: Optional[bool] = None,
        skip_web_data: Optional[bool] = None,
        use_cache: Optional[bool] = None,
    ) -> List[EnrichmentResult]:
        """Enrich multiple compounds with web data.

        Args:
            compounds: List of compounds with name and optional SMILES
            skip_predictions: Whether to skip predictions (overrides config)
            skip_web_data: Whether to skip web data (overrides config)
            use_cache: Whether to use cached results (overrides config)

        Returns:
            List of enrichment results
        """
        # Process compounds in batches
        results = []
        for i in range(0, len(compounds), self.config.batch_size):
            batch = compounds[i : i + self.config.batch_size]

            # Create tasks for batch
            tasks = []
            for compound in batch:
                name = compound["name"]
                smiles = compound.get("smiles")
                task = self.enrich_compound(
                    name,
                    smiles=smiles,
                    skip_predictions=skip_predictions,
                    skip_web_data=skip_web_data,
                    use_cache=use_cache,
                )
                tasks.append(task)

            # Process batch with progress bar
            batch_results = []
            for result in tqdm(
                asyncio.as_completed(tasks),
                total=len(tasks),
                desc=f"Processing batch {i//self.config.batch_size + 1}",
            ):
                try:
                    batch_results.append(await result)
                except Exception as e:
                    self.logger.error(f"Error in batch: {str(e)}")

            results.extend(batch_results)

        return results

    def get_metrics(self) -> Dict[str, Any]:
        """Get enrichment metrics including circuit breaker states."""
        return {
            "http_client": self.http.get_metrics(),
            "swiss_client": {
                "processed": len(self.swiss_client.processed_compounds),
                "failed": len(self.swiss_client.failed_compounds),
            },
            "community_client": {
                "processed": len(self.community_client.processed_compounds),
                "failed": len(self.community_client.failed_compounds),
            },
            "social_client": {
                "processed": len(self.social_client.processed_compounds),
                "failed": len(self.social_client.failed_compounds),
            },
        }

    def close(self) -> None:
        """Close manager and cleanup."""
        self.executor.shutdown()
        self.http.close()
        self.swiss_client.close()
        self.community_client.close()
        self.social_client.close()
