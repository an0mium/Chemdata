"""Processing pipeline for compound data enrichment.

This module provides the ProcessingPipeline class that:
1. Coordinates data source loading (BindingDB, ChEMBL, PubChem)
2. Manages web data enrichment (community sources, social media)
3. Handles ML predictions (activity, toxicity, abuse potential)
4. Provides batch processing and checkpointing
5. Exposes results through web interface

The pipeline combines:
- Target pattern matching for receptor ligands
- Activity type classification
- Structure validation and standardization
- Property calculation and enrichment
- ML-based predictions
- Web data harvesting
- Social media monitoring
"""

import json
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Set, Any
from dataclasses import dataclass, field
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

from ...models.compound.enhanced import EnhancedCompound
from ...models.compound.types import CompoundType
from ...data_sources.bindingdb import BindingDBSource
from ...web_enrichment.manager import WebEnrichmentManager
from ...web_enrichment.data_sources.community import CommunityClient
from ...web_enrichment.data_sources.swiss import SwissClient
from ...web_enrichment.llm_utils import LLMProcessor
from ...processors.structure.ml.predictors import (
    ActivityPredictor,
    ToxicityPredictor,
    AbusePotentialPredictor,
)
from ...processors.psychopharm.predictors.bbb.base import BBBPredictor
from ...web.dashboard import Dashboard
from ..infrastructure import InfrastructureManager


@dataclass
class ProcessingConfig:
    """Processing pipeline configuration."""

    # Data directories
    data_dir: Optional[Path] = None
    model_dir: Optional[Path] = None
    cache_dir: Optional[Path] = None

    # Processing settings
    batch_size: int = 100
    n_workers: int = 4

    # Feature flags
    use_ml_predictions: bool = True
    use_web_enrichment: bool = True
    use_social_monitoring: bool = True

    # Web sources
    web_sources: List[str] = field(
        default_factory=lambda: [
            "pubchem",
            "chembl",
            "drugbank",
            "wikipedia",
            "erowid",
            "psychonautwiki",
        ]
    )

    # Social media sources
    subreddits: List[str] = field(
        default_factory=lambda: [
            "researchchemicals",
            "DrugNerds",
            "Nootropics",
            "psychopharmacology",
        ]
    )

    def __post_init__(self):
        """Initialize directories."""
        for path_attr in ["data_dir", "model_dir", "cache_dir"]:
            path = getattr(self, path_attr)
            if path:
                path.parent.mkdir(parents=True, exist_ok=True)


class ProcessingPipeline:
    """Pipeline for processing compound data from multiple sources."""

    def __init__(
        self,
        config: Optional[ProcessingConfig] = None,
        infrastructure: Optional[InfrastructureManager] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize processing pipeline.

        Args:
            config: Optional processing configuration
            infrastructure: Optional infrastructure manager
            logger: Optional logger instance
        """
        self.config = config or ProcessingConfig()
        self.infrastructure = infrastructure or InfrastructureManager()
        self.logger = logger or logging.getLogger(self.__class__.__name__)

        # Initialize components
        self._init_components()

        # Initialize stats
        self.stats = {
            "total_compounds": 0,
            "bindingdb_compounds": 0,
            "web_compounds": 0,
            "social_compounds": 0,
            "with_predictions": 0,
            "with_web_data": 0,
            "errors": [],
        }

    def _init_components(self) -> None:
        """Initialize pipeline components."""
        try:
            # Initialize data sources
            self.bindingdb = BindingDBSource(
                data_dir=self.config.data_dir,
                model_dir=self.config.model_dir,
                cache_dir=self.config.cache_dir,
                n_workers=self.config.n_workers,
                batch_size=self.config.batch_size,
            )

            # Initialize web enrichment
            self.web_enrichment = WebEnrichmentManager()
            self.community_client = CommunityClient()
            self.swiss_client = SwissClient()

            # Initialize ML predictors
            if self.config.use_ml_predictions:
                self.activity_predictor = ActivityPredictor()
                self.toxicity_predictor = ToxicityPredictor()
                self.abuse_predictor = AbusePotentialPredictor()
                self.bbb_predictor = BBBPredictor()

            # Initialize LLM processor
            self.llm_processor = LLMProcessor()

            # Initialize dashboard
            self.dashboard = Dashboard()

        except Exception as e:
            self.logger.error(f"Failed to initialize components: {str(e)}")
            raise

    def process_compounds(
        self,
        input_file: Optional[Path] = None,
        output_dir: Optional[Path] = None,
        checkpoint_key: Optional[str] = None,
    ) -> List[EnhancedCompound]:
        """Process compounds from all sources.

        Args:
            input_file: Optional path to BindingDB file
            output_dir: Optional directory for output files
            checkpoint_key: Optional key for checkpointing

        Returns:
            List of processed compounds
        """
        try:
            # Start infrastructure
            self.infrastructure.start()

            # Load checkpoint if available
            if checkpoint_key:
                checkpoint = self.infrastructure.load_checkpoint(checkpoint_key)
                if checkpoint:
                    self.stats = checkpoint["stats"]
                    return checkpoint["compounds"]

            # Process compounds from all sources
            compounds = []

            # Process BindingDB compounds
            bindingdb_compounds = self._process_bindingdb(input_file)
            compounds.extend(bindingdb_compounds)
            self.stats["bindingdb_compounds"] = len(bindingdb_compounds)

            # Process web sources
            if self.config.use_web_enrichment:
                web_compounds = self._process_web_sources()
                compounds.extend(web_compounds)
                self.stats["web_compounds"] = len(web_compounds)

            # Process social media
            if self.config.use_social_monitoring:
                social_compounds = self._process_social_media()
                compounds.extend(social_compounds)
                self.stats["social_compounds"] = len(social_compounds)

            # Update total count
            self.stats["total_compounds"] = len(compounds)

            # Process all compounds
            processed = self._process_compound_batch(compounds)

            # Save checkpoint
            if checkpoint_key:
                self.infrastructure.save_checkpoint(
                    checkpoint_key,
                    {
                        "compounds": processed,
                        "stats": self.stats,
                        "timestamp": datetime.now().isoformat(),
                    },
                )

            # Save results
            if output_dir:
                self._save_results(processed, output_dir)

            # Update dashboard
            self.dashboard.update_data(processed)

            return processed

        except Exception as e:
            self.logger.error(f"Error processing compounds: {str(e)}")
            raise

        finally:
            # Stop infrastructure
            self.infrastructure.stop()

    def _process_bindingdb(
        self,
        input_file: Optional[Path] = None,
    ) -> List[EnhancedCompound]:
        """Process BindingDB compounds.

        Args:
            input_file: Optional path to BindingDB file

        Returns:
            List of processed compounds
        """
        try:
            # Load compounds from BindingDB
            compounds = self.bindingdb.load_compounds(
                input_file=input_file,
                skip_predictions=True,
                skip_web_data=True,
            )

            self.logger.info(f"Loaded {len(compounds)} compounds from BindingDB")
            return compounds

        except Exception as e:
            self.logger.error(f"Error processing BindingDB: {str(e)}")
            return []

    def _process_web_sources(self) -> List[EnhancedCompound]:
        """Process compounds from web sources.

        Returns:
            List of compounds from web sources
        """
        compounds = []

        try:
            for source in self.config.web_sources:
                try:
                    # Get compounds from source
                    source_compounds = self.community_client.get_compounds(source=source)
                    compounds.extend(source_compounds)

                except Exception as e:
                    self.logger.error(f"Error processing {source}: {str(e)}")
                    continue

            self.logger.info(f"Found {len(compounds)} compounds from web sources")
            return compounds

        except Exception as e:
            self.logger.error(f"Error processing web sources: {str(e)}")
            return []

    def _process_social_media(self) -> List[EnhancedCompound]:
        """Process compounds from social media.

        Returns:
            List of compounds from social media
        """
        compounds = []

        try:
            # Process Reddit
            for subreddit in self.config.subreddits:
                try:
                    reddit_compounds = self.web_enrichment.get_reddit_compounds(subreddit=subreddit)
                    compounds.extend(reddit_compounds)

                except Exception as e:
                    self.logger.error(f"Error processing r/{subreddit}: {str(e)}")
                    continue

            # Process Twitter
            try:
                twitter_compounds = self.web_enrichment.get_twitter_compounds()
                compounds.extend(twitter_compounds)

            except Exception as e:
                self.logger.error(f"Error processing Twitter: {str(e)}")

            self.logger.info(f"Found {len(compounds)} compounds from social media")
            return compounds

        except Exception as e:
            self.logger.error(f"Error processing social media: {str(e)}")
            return []

    def _process_compound_batch(
        self,
        compounds: List[EnhancedCompound],
    ) -> List[EnhancedCompound]:
        """Process a batch of compounds.

        Args:
            compounds: List of compounds to process

        Returns:
            List of processed compounds
        """
        processed = []

        # Process in batches
        for i in range(0, len(compounds), self.config.batch_size):
            batch = compounds[i : i + self.config.batch_size]

            # Process batch in parallel
            with ThreadPoolExecutor(max_workers=self.config.n_workers) as executor:
                futures = []
                for compound in batch:
                    future = executor.submit(
                        self._process_single_compound,
                        compound,
                    )
                    futures.append(future)

                # Collect results
                for future in executor.as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            processed.append(result)

                    except Exception as e:
                        self.logger.error(f"Error processing compound: {str(e)}")
                        continue

            # Save checkpoint after each batch
            self.infrastructure.update_progress(len(batch))

        return processed

    def _process_single_compound(
        self,
        compound: EnhancedCompound,
    ) -> Optional[EnhancedCompound]:
        """Process a single compound.

        Args:
            compound: Compound to process

        Returns:
            Processed compound or None if failed
        """
        try:
            # Add ML predictions
            if self.config.use_ml_predictions and compound.smiles:
                # Activity predictions
                activity = self.activity_predictor.predict(compound.smiles)
                if activity:
                    compound.activity_predictions = activity

                # Toxicity predictions
                toxicity = self.toxicity_predictor.predict(compound.smiles)
                if toxicity:
                    compound.toxicity_predictions = toxicity

                # Abuse potential predictions
                abuse = self.abuse_predictor.predict(compound.smiles)
                if abuse:
                    compound.abuse_predictions = abuse

                # BBB predictions
                bbb = self.bbb_predictor.predict(compound.smiles)
                if bbb:
                    compound.bbb_predictions = bbb

                self.stats["with_predictions"] += 1

            # Add web enrichment
            if self.config.use_web_enrichment:
                # Get community data
                community_data = self.community_client.get_compound_data(
                    compound.name,
                    compound.smiles,
                )
                if community_data:
                    compound.community_data = community_data

                # Get Swiss data
                swiss_data = self.swiss_client.get_compound_data(
                    compound.name,
                    compound.smiles,
                )
                if swiss_data:
                    compound.swiss_data = swiss_data

                self.stats["with_web_data"] += 1

            return compound

        except Exception as e:
            self.logger.error(f"Error processing compound {compound.name}: {str(e)}")
            return None

    def _save_results(
        self,
        compounds: List[EnhancedCompound],
        output_dir: Path,
    ) -> None:
        """Save processed compounds.

        Args:
            compounds: List of compounds to save
            output_dir: Output directory
        """
        try:
            # Create output directory
            output_dir.mkdir(parents=True, exist_ok=True)

            # Save JSON data
            json_file = output_dir / "compounds.json"
            with open(json_file, "w") as f:
                json.dump(
                    [c.to_dict() for c in compounds],
                    f,
                    indent=2,
                )

            # Save TSV data
            tsv_file = output_dir / "compounds.tsv"
            df = pd.DataFrame(
                [
                    {
                        "name": c.name,
                        "smiles": c.smiles,
                        "cas_number": c.cas_number,
                        "source": c.source,
                        "compound_type": c.compound_type.name,
                        "activity_predictions": str(c.activity_predictions),
                        "toxicity_predictions": str(c.toxicity_predictions),
                        "abuse_predictions": str(c.abuse_predictions),
                        "bbb_predictions": str(c.bbb_predictions),
                    }
                    for c in compounds
                ]
            )
            df.to_csv(tsv_file, sep="\t", index=False)

            self.logger.info(f"Saved {len(compounds)} compounds to {output_dir}")

        except Exception as e:
            self.logger.error(f"Error saving results: {str(e)}")
            raise
