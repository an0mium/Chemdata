"""Enhanced web enrichment manager.

This module provides an enhanced manager class that coordinates web enrichment clients
with improved error handling and resilience through circuit breakers:
- Swiss tools (SwissTargetPrediction, SwissADME)
- Community sources (PsychonautWiki, Erowid, TripSit)
- Social media (Reddit, Twitter)

The manager handles:
1. Client initialization and configuration
2. Coordinated data enrichment
3. Error handling and recovery with circuit breakers
4. Progress tracking and metrics
5. Result aggregation
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime

from tqdm import tqdm

from .http_client_enhanced import HTTPClientEnhanced
from .swiss_client import SwissClient
from .community_client import CommunityClient
from .social_client import SocialClient
from ..models.compound import Compound
from ..pipeline.infrastructure.circuit_breaker import CircuitConfig


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
    circuit_config: Optional[CircuitConfig] = None


class WebEnrichmentManagerEnhanced:
    """Enhanced manager for web enrichment clients."""

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

        # Initialize thread pool
        self.executor = ThreadPoolExecutor(
            max_workers=config.n_workers,
            thread_name_prefix="enrichment",
        )

    def enrich_compounds(
        self,
        compounds: List[Compound],
        skip_predictions: Optional[bool] = None,
        skip_web_data: Optional[bool] = None,
        use_cache: Optional[bool] = None,
    ) -> None:
        """Enrich compounds with web data.
        
        Args:
            compounds: List of compounds to enrich
            skip_predictions: Whether to skip predictions (overrides config)
            skip_web_data: Whether to skip web data (overrides config)
            use_cache: Whether to use cached results (overrides config)
        """
        # Use config defaults if not specified
        skip_predictions = (
            skip_predictions
            if skip_predictions is not None
            else self.config.skip_predictions
        )
        skip_web_data = (
            skip_web_data
            if skip_web_data is not None
            else self.config.skip_web_data
        )
        use_cache = (
            use_cache
            if use_cache is not None
            else self.config.use_cache
        )

        # Process compounds in batches
        for i in range(0, len(compounds), self.config.batch_size):
            batch = compounds[i:i + self.config.batch_size]
            self._process_batch(
                batch,
                skip_predictions=skip_predictions,
                skip_web_data=skip_web_data,
                use_cache=use_cache,
            )

    def _process_batch(
        self,
        compounds: List[Compound],
        skip_predictions: bool,
        skip_web_data: bool,
        use_cache: bool,
    ) -> None:
        """Process a batch of compounds.
        
        Args:
            compounds: List of compounds to process
            skip_predictions: Whether to skip predictions
            skip_web_data: Whether to skip web data
            use_cache: Whether to use cached results
        """
        # Submit enrichment tasks
        futures = []
        for compound in compounds:
            futures.append(
                self.executor.submit(
                    self._enrich_compound,
                    compound,
                    skip_predictions=skip_predictions,
                    skip_web_data=skip_web_data,
                    use_cache=use_cache,
                )
            )

        # Process results with progress bar
        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Enriching compounds",
        ):
            try:
                future.result()
            except Exception as e:
                self.logger.error(f"Error enriching compound: {str(e)}")

    def _enrich_compound(
        self,
        compound: Compound,
        skip_predictions: bool,
        skip_web_data: bool,
        use_cache: bool,
    ) -> None:
        """Enrich a single compound.
        
        Args:
            compound: Compound to enrich
            skip_predictions: Whether to skip predictions
            skip_web_data: Whether to skip web data
            use_cache: Whether to use cached results
        """
        try:
            # Initialize enrichment data
            compound.swiss_data = {}
            compound.community_data = {}
            compound.social_data = {}

            # Get Swiss data
            if not skip_predictions:
                self.swiss_client.process_compounds(
                    [compound],
                    skip_predictions=skip_predictions,
                    use_cache=use_cache,
                )

            # Get community data
            if not skip_web_data:
                self.community_client.process_compounds(
                    [compound],
                    skip_predictions=skip_predictions,
                    use_cache=use_cache,
                )

            # Get social data
            if not skip_web_data:
                self.social_client.process_compounds(
                    [compound],
                    skip_predictions=skip_predictions,
                    use_cache=use_cache,
                )

            # Add metadata
            compound.enrichment_metadata = {
                "timestamp": datetime.now().isoformat(),
                "sources": [],
            }
            if compound.swiss_data:
                compound.enrichment_metadata["sources"].append("swiss")
            if compound.community_data:
                compound.enrichment_metadata["sources"].append("community")
            if compound.social_data:
                compound.enrichment_metadata["sources"].append("social")

        except Exception as e:
            self.logger.error(
                f"Error enriching {compound.name}: {str(e)}"
            )
            raise

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
