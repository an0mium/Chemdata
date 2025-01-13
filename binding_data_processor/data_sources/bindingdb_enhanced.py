"""Enhanced BindingDB data source integration.

This module enhances the base BindingDB functionality with:
1. Enhanced web enrichment clients for better data harvesting
2. Patent search integration for finding additional compounds
3. Circuit breaker pattern for resilience
4. Better error handling and recovery
5. Improved metrics collection
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

from ..web_enrichment.clients.community import CommunityClient as CommunityClientEnhanced
from ..web_enrichment.swiss_client_enhanced import SwissClientEnhanced
from ..web_enrichment.social_client_enhanced import SocialClientEnhanced
from ..web_enrichment.clients.patents import PatentClient
from ..pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig
from .bindingdb import BindingDBSource


class BindingDBSourceEnhanced(BindingDBSource):
    """Enhanced BindingDB data source with additional data harvesting."""

    def __init__(
        self,
        data_dir: Optional[Path] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        n_workers: int = 4,
        batch_size: int = 1000,
        logger: Optional[logging.Logger] = None,
        circuit_config: Optional[CircuitBreakerConfig] = None,
    ):
        """Initialize enhanced BindingDB source.

        Args:
            data_dir: Optional directory for data files
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            n_workers: Number of worker threads
            batch_size: Batch size for processing
            logger: Optional logger instance
            circuit_config: Optional circuit breaker configuration
        """
        super().__init__(
            data_dir=data_dir,
            model_dir=model_dir,
            cache_dir=cache_dir,
            n_workers=n_workers,
            batch_size=batch_size,
            logger=logger,
        )

        # Initialize enhanced clients
        self.community_client = CommunityClientEnhanced(
            http_client=self.http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            logger=logger,
            circuit_config=circuit_config,
        )
        self.swiss_client = SwissClientEnhanced(
            http_client=self.http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            logger=logger,
            circuit_config=circuit_config,
        )
        self.social_client = SocialClientEnhanced(
            http_client=self.http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            logger=logger,
            circuit_config=circuit_config,
        )
        self.patent_client = PatentClient(
            http_client=self.http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            logger=logger,
            circuit_config=circuit_config,
        )

        # Track additional metrics
        self.stats.update(
            {
                "patent_compounds": 0,
                "community_compounds": 0,
                "swiss_compounds": 0,
                "social_compounds": 0,
            }
        )

    def process_all_sources(
        self,
        bindingdb_file: Optional[str] = None,
        output_dir: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Process compounds from all sources including patents.

        Args:
            bindingdb_file: Optional path to BindingDB file
            output_dir: Optional directory for output files

        Returns:
            Processing statistics
        """
        try:
            # Get base compounds from BindingDB
            stats = super().process_all_sources(
                bindingdb_file=bindingdb_file,
                output_dir=output_dir,
            )

            # Search patents for additional compounds
            self.logger.info("Searching patents for additional compounds...")
            patent_compounds = []
            for pattern_name, pattern in self.TARGET_PATTERNS.items():
                try:
                    compounds = self.patent_client.search_compounds(
                        query=pattern,
                        use_cache=True,
                    )
                    if compounds:
                        self.logger.info(f"Found {len(compounds)} compounds from patents for {pattern_name}")
                        patent_compounds.extend(compounds)
                        self.stats["patent_compounds"] += len(compounds)
                except Exception as e:
                    self.logger.error(f"Error searching patents for {pattern_name}: {str(e)}")

            # Process patent compounds
            if patent_compounds:
                self.logger.info(f"Processing {len(patent_compounds)} patent compounds...")
                processed = self._process_compounds(patent_compounds)
                if processed:
                    self.logger.info(f"Successfully processed {len(processed)} patent compounds")

                    # Save patent compounds
                    if output_dir:
                        output_dir = Path(output_dir)
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                        self._save_results(
                            processed,
                            output_dir / f"patent_compounds_{timestamp}",
                        )

            return self.stats

        except Exception as e:
            self.logger.error(f"Error in enhanced processing: {str(e)}")
            raise

    def _process_compounds(self, compounds: List[Any]) -> List[Any]:
        """Process compounds with enhanced enrichment.

        Args:
            compounds: List of compounds to process

        Returns:
            List of processed compounds
        """
        processed = super()._process_compounds(compounds)

        if not processed:
            return []

        try:
            # Additional enrichment with enhanced clients
            self.logger.info("Performing enhanced enrichment...")

            # Community data
            community_compounds = self.community_client.process_compounds(
                processed,
                skip_predictions=True,
                use_cache=True,
            )
            if community_compounds:
                self.stats["community_compounds"] += len(community_compounds)

            # Swiss data
            swiss_compounds = self.swiss_client.process_compounds(
                processed,
                skip_predictions=True,
                use_cache=True,
            )
            if swiss_compounds:
                self.stats["swiss_compounds"] += len(swiss_compounds)

            # Social data
            social_compounds = self.social_client.process_compounds(
                processed,
                skip_predictions=True,
                use_cache=True,
            )
            if social_compounds:
                self.stats["social_compounds"] += len(social_compounds)

        except Exception as e:
            self.logger.error(f"Error in enhanced enrichment: {str(e)}")

        return processed

    def get_metrics(self) -> Dict[str, Any]:
        """Get enhanced metrics."""
        metrics = super().get_metrics()
        metrics.update(
            {
                "patent_compounds": self.stats["patent_compounds"],
                "community_compounds": self.stats["community_compounds"],
                "swiss_compounds": self.stats["swiss_compounds"],
                "social_compounds": self.stats["social_compounds"],
            }
        )
        return metrics
