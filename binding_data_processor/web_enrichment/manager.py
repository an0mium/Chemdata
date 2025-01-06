"""Web enrichment manager.

This module provides a manager for web enrichment clients:
- Swiss tools (SwissTargetPrediction, SwissADME)
- Community sources (PsychonautWiki, Erowid, TripSit)
- Social media (Reddit, Twitter)
- Patent data (Google Patents, Lens.org)

The manager handles:
1. Client initialization and configuration
2. Coordinated data enrichment
3. Error handling and recovery
4. Progress tracking
5. Result aggregation
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Type

from tqdm import tqdm

from .base import WebClient, WebClientError
from .clients.swiss import SwissClient
from .clients.community import CommunityClient
from .clients.social import SocialClient
from .clients.patent import PatentClient
from ..models.compound import Compound
from ..pipeline.infrastructure.circuit_breaker import CircuitConfig
from ..pipeline.infrastructure.monitoring import MetricsCollector


@dataclass
class EnrichmentConfig:
    """Configuration for web enrichment."""

    # API credentials
    reddit_client_id: Optional[str] = None
    reddit_client_secret: Optional[str] = None
    twitter_bearer_token: Optional[str] = None
    google_patents_key: Optional[str] = None
    lens_api_key: Optional[str] = None

    # Processing options
    skip_predictions: bool = False
    skip_web_data: bool = False
    use_cache: bool = True
    n_workers: int = 4
    batch_size: int = 100

    # Directories
    model_dir: Optional[Path] = None
    cache_dir: Optional[Path] = None

    # Circuit breaker
    circuit_config: Optional[CircuitConfig] = None

    # Client configs
    client_configs: Dict[str, Dict[str, Any]] = field(default_factory=dict)


class WebEnrichmentError(Exception):
    """Web enrichment error."""

    def __init__(
        self,
        message: str,
        client: str,
        details: Optional[Dict[str, Any]] = None,
    ):
        """Initialize error.

        Args:
            message: Error message
            client: Client name
            details: Optional error details
        """
        super().__init__(message)
        self.client = client
        self.details = details or {}
        self.timestamp = datetime.now().isoformat()


class WebEnrichmentManager:
    """Manager for web enrichment clients."""

    def __init__(
        self,
        config: EnrichmentConfig,
        logger: Optional[logging.Logger] = None,
        metrics_collector: Optional[MetricsCollector] = None,
    ):
        """Initialize manager.

        Args:
            config: Enrichment configuration
            logger: Optional logger instance
            metrics_collector: Optional metrics collector
        """
        self.config = config
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        self.metrics = metrics_collector or MetricsCollector(
            namespace="web_enrichment",
            logger=self.logger,
        )

        # Initialize clients
        self.clients: Dict[str, WebClient] = {}
        self.processed_compounds: Set[str] = set()
        self.failed_compounds: Set[str] = set()

        # Initialize thread pool
        self.executor = ThreadPoolExecutor(
            max_workers=config.n_workers,
            thread_name_prefix="enrichment",
        )

        # Register default clients
        self._register_default_clients()

    def _register_default_clients(self) -> None:
        """Register default web enrichment clients."""
        # Swiss client
        self.register_client(
            SwissClient,
            model_dir=self.config.model_dir,
            cache_dir=(self.config.cache_dir / "swiss" if self.config.cache_dir else None),
            **self.config.client_configs.get("swiss", {}),
        )

        # Community client
        self.register_client(
            CommunityClient,
            model_dir=self.config.model_dir,
            cache_dir=(self.config.cache_dir / "community" if self.config.cache_dir else None),
            **self.config.client_configs.get("community", {}),
        )

        # Social client
        self.register_client(
            SocialClient,
            reddit_client_id=self.config.reddit_client_id,
            reddit_client_secret=self.config.reddit_client_secret,
            twitter_bearer_token=self.config.twitter_bearer_token,
            model_dir=self.config.model_dir,
            cache_dir=(self.config.cache_dir / "social" if self.config.cache_dir else None),
            **self.config.client_configs.get("social", {}),
        )

        # Patent client
        self.register_client(
            PatentClient,
            google_patents_key=self.config.google_patents_key,
            lens_api_key=self.config.lens_api_key,
            model_dir=self.config.model_dir,
            cache_dir=(self.config.cache_dir / "patent" if self.config.cache_dir else None),
            **self.config.client_configs.get("patent", {}),
        )

    def register_client(
        self,
        client_class: Type[WebClient],
        **kwargs: Any,
    ) -> None:
        """Register web enrichment client.

        Args:
            client_class: Client class to register
            **kwargs: Additional client arguments
        """
        client = client_class(
            circuit_config=self.config.circuit_config,
            logger=self.logger,
            **kwargs,
        )
        self.clients[client.name] = client
        self.logger.info(f"Registered client: {client.name}")

    def get_client(self, name: str) -> WebClient:
        """Get client by name.

        Args:
            name: Client name

        Returns:
            Web enrichment client

        Raises:
            WebEnrichmentError: If client not found
        """
        if name not in self.clients:
            raise WebEnrichmentError(
                message=f"Client not found: {name}",
                client=name,
            )
        return self.clients[name]

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
        skip_predictions = skip_predictions if skip_predictions is not None else self.config.skip_predictions
        skip_web_data = skip_web_data if skip_web_data is not None else self.config.skip_web_data
        use_cache = use_cache if use_cache is not None else self.config.use_cache

        # Process compounds in batches
        for i in range(0, len(compounds), self.config.batch_size):
            batch = compounds[i : i + self.config.batch_size]
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
            compound.patent_data = {}

            # Get Swiss data
            if not skip_predictions:
                swiss_client = self.get_client("swiss")
                compound.swiss_data = swiss_client.predict_targets(
                    smiles=compound.smiles,
                    use_cache=use_cache,
                )
                self.processed_compounds.add(compound.smiles)

            # Get community data
            if not skip_web_data:
                community_client = self.get_client("community")
                compound.community_data = community_client.get_data(
                    smiles=compound.smiles,
                    use_cache=use_cache,
                )
                self.processed_compounds.add(compound.smiles)

            # Get social data
            if not skip_web_data:
                social_client = self.get_client("social")
                compound.social_data = social_client.get_data(
                    smiles=compound.smiles,
                    use_cache=use_cache,
                )
                self.processed_compounds.add(compound.smiles)

            # Get patent data
            if not skip_web_data:
                patent_client = self.get_client("patent")
                compound.patent_data = patent_client.search_patents(
                    query=compound.name,
                    chemical_structure=compound.smiles,
                    use_cache=use_cache,
                )
                self.processed_compounds.add(compound.smiles)

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
            if compound.patent_data:
                compound.enrichment_metadata["sources"].append("patent")

        except WebClientError as e:
            self.failed_compounds.add(compound.smiles)
            self.logger.error(f"Error enriching {compound.name}: {str(e)}")
            raise WebEnrichmentError(
                message=str(e),
                client=e.source,
                details=e.details,
            )

        except Exception as e:
            self.failed_compounds.add(compound.smiles)
            self.logger.error(f"Error enriching {compound.name}: {str(e)}")
            raise

    def get_metrics(self) -> Dict[str, Any]:
        """Get manager metrics.

        Returns:
            Manager metrics
        """
        metrics = {
            "processed_compounds": len(self.processed_compounds),
            "failed_compounds": len(self.failed_compounds),
            "clients": {},
        }

        for name, client in self.clients.items():
            metrics["clients"][name] = client.get_metrics()

        return metrics

    def close(self) -> None:
        """Close manager and cleanup."""
        self.executor.shutdown()
        for client in self.clients.values():
            client.close()
