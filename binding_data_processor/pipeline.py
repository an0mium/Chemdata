"""Pipeline manager for compound data processing.

This module provides:
1. Automated BindingDB data processing
2. Web data enrichment and social media harvesting
3. ML predictions and analysis
4. Data validation and deduplication
5. Progress tracking and error handling
6. Checkpointing and caching
7. Parallel processing and error recovery
8. Model ensembling and confidence scoring
"""

import logging
import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Set, Tuple
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict

import pandas as pd
import numpy as np
from tqdm import tqdm
import torch
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from transformers import pipeline
from sentence_transformers import SentenceTransformer, util

from binding_data_processor.models.compound import CompoundData
from binding_data_processor.data_sources.bindingdb import BindingDBProcessor
from binding_data_processor.processors.enrichment import EnrichmentProcessor
from binding_data_processor.processors.structure.similarity import (
    calculate_similarity,
    get_fingerprint,
)
from binding_data_processor.processors.structure.ml.predictors.psychoactive import (
    PsychoactivePredictor,
)
from binding_data_processor.processors.psychopharm.predictors.nootropic.base import (
    NootropicPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.activity import (
    ActivityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.toxicity import (
    ToxicityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.abuse import (
    AbusePotentialPredictor,
)
from binding_data_processor.processors.structure.ml.models.ensemble import (
    EnsembleModel,
)
from web_enrichment.data_sources.social import SocialDataHarvester
from web_enrichment.data_sources.community import CommunityClient
from web_enrichment.data_sources.swiss import SwissClient
from web_enrichment.data_sources.chembl import ChEMBLClient
from web_enrichment.data_sources.pubchem import PubChemClient
from web_enrichment.data_sources.regulatory import RegulatoryClient
from web_enrichment.data_sources.web_search import WebSearchClient
from web_enrichment.http_client import HttpClient
from web_enrichment.llm_utils import analyze_content_with_llm


class PipelineManager:
    """Manager for compound data processing pipeline."""

    def __init__(
        self,
        data_dir: str = "../data",
        model_dir: str = "../models",
        n_workers: int = 4,
        batch_size: int = 100,
        checkpoint_interval: int = 1000,
        max_retries: int = 3,
        cache_dir: Optional[str] = None,
        reddit_client_id: Optional[str] = None,
        reddit_client_secret: Optional[str] = None,
        twitter_api_key: Optional[str] = None,
        twitter_api_secret: Optional[str] = None,
        discord_token: Optional[str] = None,
        bluesky_handle: Optional[str] = None,
        bluesky_password: Optional[str] = None,
        llm_api_key: Optional[str] = None,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
    ):
        """Initialize pipeline manager.

        Args:
            data_dir: Directory for data files
            model_dir: Directory containing ML models
            n_workers: Number of worker threads
            batch_size: Batch size for processing
            checkpoint_interval: Save checkpoint every N compounds
            max_retries: Maximum number of retries for failed operations
            cache_dir: Optional custom cache directory
            reddit_client_id: Optional Reddit API client ID
            reddit_client_secret: Optional Reddit API client secret
            twitter_api_key: Optional Twitter API key
            twitter_api_secret: Optional Twitter API secret
            discord_token: Optional Discord bot token
            bluesky_handle: Optional Bluesky handle
            bluesky_password: Optional Bluesky password
            llm_api_key: Optional API key for LLM analysis
            device: Device to run ML models on
        """
        # Initialize base settings
        self.data_dir = Path(data_dir)
        self.model_dir = Path(model_dir)
        self.checkpoint_interval = checkpoint_interval
        self.max_retries = max_retries
        self.n_workers = n_workers
        self.llm_api_key = llm_api_key
        self.device = device

        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        self.logger.addHandler(handler)

        # Initialize ML models
        self.logger.info("Initializing ML models...")

        # Base predictors
        self.activity_predictor = ActivityPredictor(model_dir=f"{model_dir}/activity")
        self.toxicity_predictor = ToxicityPredictor(model_dir=f"{model_dir}/toxicity")
        self.abuse_predictor = AbusePotentialPredictor(model_dir=f"{model_dir}/abuse")
        self.psychoactive_predictor = PsychoactivePredictor(model_dir=f"{model_dir}/psychoactive")
        self.nootropic_predictor = NootropicPredictor(model_dir=f"{model_dir}/nootropic")

        # Classification models
        self.text_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/pipeline/compound_classifier",
            device=device,
        )
        self.activity_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/pipeline/activity_classifier",
            device=device,
        )
        self.safety_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/pipeline/safety_classifier",
            device=device,
        )
        self.mechanism_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/pipeline/mechanism_classifier",
            device=device,
        )

        # Ensemble models
        self.activity_ensemble = EnsembleModel(
            [
                self.activity_predictor,
                self.psychoactive_predictor,
                self.nootropic_predictor,
            ],
            weights=[0.4, 0.3, 0.3],
        )
        self.safety_ensemble = EnsembleModel(
            [self.toxicity_predictor, self.abuse_predictor],
            weights=[0.5, 0.5],
        )

        # Similarity model
        self.similarity_model = SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2").to(device)

        # Initialize processors
        self.logger.info("Initializing processors...")
        self.bindingdb = BindingDBProcessor(
            data_dir=data_dir,
            model_dir=model_dir,
            n_workers=n_workers,
            batch_size=batch_size,
            checkpoint_interval=checkpoint_interval,
            device=device,
        )
        self.enrichment = EnrichmentProcessor(
            model_dir=model_dir,
            n_workers=n_workers,
            batch_size=batch_size,
            checkpoint_interval=checkpoint_interval,
            max_retries=max_retries,
            cache_dir=cache_dir,
            reddit_client_id=reddit_client_id,
            reddit_client_secret=reddit_client_secret,
            twitter_api_key=twitter_api_key,
            twitter_api_secret=twitter_api_secret,
            discord_token=discord_token,
            bluesky_handle=bluesky_handle,
            bluesky_password=bluesky_password,
            device=device,
        )

        # Initialize web clients
        self.logger.info("Initializing web clients...")
        self.http_client = HttpClient()

        # Community data
        self.community_client = CommunityClient(
            http_client=self.http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            n_workers=n_workers,
            max_retries=max_retries,
            device=device,
        )

        # Swiss data
        self.swiss_client = SwissClient(
            http_client=self.http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            n_workers=n_workers,
            max_retries=max_retries,
            device=device,
        )

        # Social data
        self.social_harvester = SocialDataHarvester(
            reddit_client_id=reddit_client_id,
            reddit_client_secret=reddit_client_secret,
            twitter_api_key=twitter_api_key,
            twitter_api_secret=twitter_api_secret,
            discord_token=discord_token,
            bluesky_handle=bluesky_handle,
            bluesky_password=bluesky_password,
            http_client=self.http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            n_workers=n_workers,
            max_retries=max_retries,
            device=device,
        )

        # Database clients
        self.chembl_client = ChEMBLClient(
            http_client=self.http_client,
            cache_dir=cache_dir,
            n_workers=n_workers,
            max_retries=max_retries,
        )
        self.pubchem_client = PubChemClient(
            http_client=self.http_client,
            cache_dir=cache_dir,
            n_workers=n_workers,
            max_retries=max_retries,
        )

        # Additional data sources
        self.regulatory_client = RegulatoryClient(
            http_client=self.http_client,
            cache_dir=cache_dir,
            n_workers=n_workers,
            max_retries=max_retries,
        )
        self.web_search = WebSearchClient(
            http_client=self.http_client,
            cache_dir=cache_dir,
            n_workers=n_workers,
            max_retries=max_retries,
        )

        # Initialize thread pools
        self.prediction_executor = ThreadPoolExecutor(max_workers=n_workers, thread_name_prefix="prediction")
        self.web_executor = ThreadPoolExecutor(max_workers=n_workers * 2, thread_name_prefix="web")  # More threads for web IO

        # Initialize storage
        self.compounds: Dict[str, CompoundData] = {}
        self.processed_ids: Set[str] = set()
        self.stats: Dict[str, Any] = {
            "total_compounds": 0,
            "bindingdb_compounds": 0,
            "web_compounds": 0,
            "social_compounds": 0,
            "with_predictions": 0,
            "with_web_data": 0,
            "from_cache": 0,
            "errors": [],
            "start_time": datetime.now().isoformat(),
        }

    def run_pipeline(
        self,
        bindingdb_file: Optional[str] = None,
        target_patterns: Optional[Dict[str, str]] = None,
        skip_predictions: bool = False,
        skip_web_data: bool = False,
        output_dir: Optional[str] = None,
        checkpoint_file: Optional[str] = None,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Run the complete data processing pipeline.

        Args:
            bindingdb_file: Optional path to BindingDB file
            target_patterns: Optional dict of target regex patterns
            skip_predictions: Whether to skip ML predictions
            skip_web_data: Whether to skip web data enrichment
            output_dir: Optional directory for output files
            checkpoint_file: Optional checkpoint file to resume from
            use_cache: Whether to use cached results

        Returns:
            Processing statistics
        """
        try:
            output_dir = Path(output_dir or self.data_dir)
            os.makedirs(output_dir, exist_ok=True)

            # Load checkpoint if exists
            if checkpoint_file and os.path.exists(checkpoint_file):
                self._load_checkpoint(checkpoint_file)

            # Process BindingDB compounds
            self.logger.info("Processing BindingDB data...")
            bindingdb_compounds = self.bindingdb.process_data(
                input_file=bindingdb_file,
                target_patterns=target_patterns,
                checkpoint_file=checkpoint_file,
                skip_predictions=skip_predictions,
                skip_web_data=skip_web_data,
                use_cache=use_cache,
            )
            self.stats["bindingdb_compounds"] = len(bindingdb_compounds)
            self._add_compounds(bindingdb_compounds)

            # Search web sources
            if not skip_web_data:
                # Search Swiss data
                self.logger.info("Searching Swiss data...")
                swiss_compounds = self._search_swiss_data(use_cache=use_cache)
                self.stats["swiss_compounds"] = len(swiss_compounds)
                self._add_compounds(swiss_compounds)

                # Search community data
                self.logger.info("Searching community data...")
                community_compounds = self._search_community_data(use_cache=use_cache)
                self.stats["community_compounds"] = len(community_compounds)
                self._add_compounds(community_compounds)

                # Search social media
                self.logger.info("Searching social media...")
                social_compounds = self._search_social_media(use_cache=use_cache)
                self.stats["social_compounds"] = len(social_compounds)
                self._add_compounds(social_compounds)

                # Search databases
                self.logger.info("Searching databases...")
                database_compounds = self._search_databases(use_cache=use_cache)
                self.stats["database_compounds"] = len(database_compounds)
                self._add_compounds(database_compounds)

                # Search web
                self.logger.info("Searching web...")
                web_compounds = self._search_web(use_cache=use_cache)
                self.stats["web_compounds"] = len(web_compounds)
                self._add_compounds(web_compounds)

            # Enrich compounds
            self.logger.info("Enriching compounds...")
            enrichment_stats = self.enrichment.enrich_compounds(
                list(self.compounds.values()),
                skip_predictions=skip_predictions,
                skip_web_data=skip_web_data,
                checkpoint_file=checkpoint_file,
                use_cache=use_cache,
            )
            self.stats.update(enrichment_stats)

            # Save results
            self.logger.info("Saving results...")
            self._save_results(output_dir)

            # Update final stats
            self.stats["end_time"] = datetime.now().isoformat()
            return self.stats

        except KeyboardInterrupt:
            self.logger.info("Pipeline interrupted by user")
            if checkpoint_file:
                self._save_checkpoint(checkpoint_file)
            raise
        except Exception as e:
            self.logger.error(f"Error in pipeline: {str(e)}")
            if checkpoint_file:
                self._save_checkpoint(checkpoint_file)
            raise

    def _search_swiss_data(self, use_cache: bool = True) -> List[CompoundData]:
        """Search Swiss data sources.

        Args:
            use_cache: Whether to use cached results

        Returns:
            List of compounds found from Swiss data
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                # Search Swiss data
                swiss_compounds = self.swiss_client.search_compounds(
                    target_types=[
                        "5-HT2",
                        "NMDA",
                        "nootropic",
                        "psychedelic",
                        "dissociative",
                        "stimulant",
                    ],
                    max_compounds=1000,
                    use_cache=use_cache,
                )

                # Filter and validate compounds
                valid_compounds = []
                for compound in tqdm(
                    swiss_compounds,
                    desc="Validating Swiss compounds",
                    unit="compounds",
                ):
                    if self._validate_compound(compound):
                        valid_compounds.append(compound)
                compounds.extend(valid_compounds)
                self.logger.info(f"Found {len(valid_compounds)} valid compounds from Swiss data")

                break  # Success, exit retry loop

            except Exception as e:
                self.logger.error(f"Error searching Swiss data (attempt {retries + 1}/{self.max_retries}): {e}")
                retries += 1
                if retries == self.max_retries:
                    self.logger.error(f"Failed to search Swiss data after {self.max_retries} attempts")
                    self.stats["errors"].append(f"Swiss data search error: {str(e)}")

        return self._deduplicate_compounds(compounds)

    def _search_community_data(self, use_cache: bool = True) -> List[CompoundData]:
        """Search community data sources.

        Args:
            use_cache: Whether to use cached results

        Returns:
            List of compounds found from community data
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                # Search community data
                community_compounds = self.community_client.search_compounds(
                    keywords=[
                        "psychoactive",
                        "hallucinogen",
                        "antidepressant",
                        "nootropic",
                        "dissociative",
                        "stimulant",
                        "research chemical",
                        "novel compound",
                    ],
                    max_compounds=1000,
                    use_cache=use_cache,
                )

                # Filter and validate compounds
                valid_compounds = []
                for compound in tqdm(
                    community_compounds,
                    desc="Validating community compounds",
                    unit="compounds",
                ):
                    if self._validate_compound(compound):
                        valid_compounds.append(compound)
                compounds.extend(valid_compounds)
                self.logger.info(f"Found {len(valid_compounds)} valid compounds from community data")

                break  # Success, exit retry loop

            except Exception as e:
                self.logger.error(f"Error searching community data (attempt {retries + 1}/{self.max_retries}): {e}")
                retries += 1
                if retries == self.max_retries:
                    self.logger.error(f"Failed to search community data after {self.max_retries} attempts")
                    self.stats["errors"].append(f"Community data search error: {str(e)}")

        return self._deduplicate_compounds(compounds)

    def _search_social_media(self, use_cache: bool = True) -> List[CompoundData]:
        """Search social media sources.

        Args:
            use_cache: Whether to use cached results

        Returns:
            List of compounds found from social media
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                # Search all social media sources
                mentions = self.social_harvester.search_all_sources(
                    days=365,
                    limit=10000,
                    use_cache=use_cache,
                )

                # Validate and standardize compounds
                valid_compounds = []
                for mention in tqdm(
                    mentions,
                    desc="Validating social media compounds",
                    unit="mentions",
                ):
                    try:
                        # Validate mention with ML
                        if self._validate_mention(mention):
                            compound = self.social_harvester.validate_compound(
                                mention["identifier"],
                                mention["type"],
                                use_cache=use_cache,
                            )
                            if compound and self._validate_compound(compound):
                                valid_compounds.append(compound)
                    except Exception as e:
                        self.logger.error(f"Error validating compound {mention['identifier']}: {e}")
                        self.stats["errors"].append(f"Compound validation error: {str(e)}")

                compounds.extend(valid_compounds)
                self.logger.info(f"Found {len(valid_compounds)} valid compounds from social media")

                break  # Success, exit retry loop

            except Exception as e:
                self.logger.error(f"Error searching social media (attempt {retries + 1}/{self.max_retries}): {e}")
                retries += 1
                if retries == self.max_retries:
                    self.logger.error(f"Failed to search social media after {self.max_retries} attempts")
                    self.stats["errors"].append(f"Social media search error: {str(e)}")

        return self._deduplicate_compounds(compounds)

    def _search_databases(self, use_cache: bool = True) -> List[CompoundData]:
        """Search chemical databases.

        Args:
            use_cache: Whether to use cached results

        Returns:
            List of compounds found from databases
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                # Search ChEMBL
                self.logger.info("Searching ChEMBL...")
                chembl_compounds = self.chembl_client.search_compounds(
                    target_types=[
                        "5-HT2",
                        "NMDA",
                        "nootropic",
                        "psychedelic",
                        "dissociative",
                        "stimulant",
                    ],
                    max_compounds=1000,
                    use_cache=use_cache,
                )
                valid_chembl = [c for c in chembl_compounds if self._validate_compound(c)]
                compounds.extend(valid_chembl)
                self.logger.info(f"Found {len(valid_chembl)} valid compounds from ChEMBL")

                # Search PubChem
                self.logger.info("Searching PubChem...")
                pubchem_compounds = self.pubchem_client.search_compounds(
                    keywords=[
                        "psychoactive",
                        "hallucinogen",
                        "antidepressant",
                        "nootropic",
                        "dissociative",
                        "stimulant",
                    ],
                    max_compounds=1000,
                    use_cache=use_cache,
                )
                valid_pubchem = [c for c in pubchem_compounds if self._validate_compound(c)]
                compounds.extend(valid_pubchem)
                self.logger.info(f"Found {len(valid_pubchem)} valid compounds from PubChem")

                break  # Success, exit retry loop

            except Exception as e:
                self.logger.error(f"Error searching databases (attempt {retries + 1}/{self.max_retries}): {e}")
                retries += 1
                if retries == self.max_retries:
                    self.logger.error(f"Failed to search databases after {self.max_retries} attempts")
                    self.stats["errors"].append(f"Database search error: {str(e)}")

        return self._deduplicate_compounds(compounds)

    def _search_web(self, use_cache: bool = True) -> List[CompoundData]:
        """Search web sources.

        Args:
            use_cache: Whether to use cached results

        Returns:
            List of compounds found from web sources
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                # Search web
                web_results = self.web_search.search(
                    keywords=[
                        "novel psychoactive substance",
                        "research chemical",
                        "new psychedelic compound",
                        "emerging nootropic",
                        "novel dissociative",
                        "new stimulant compound",
                    ],
                    max_results=1000,
                    use_cache=use_cache,
                )

                # Extract and validate compounds
                valid_compounds = []
                for result in tqdm(
                    web_results,
                    desc="Processing web results",
                    unit="results",
                ):
                    try:
                        compounds = self.web_search.extract_compounds(
                            result["url"],
                            result["text"],
                            use_cache=use_cache,
                        )
                        for compound in compounds:
                            if self._validate_compound(compound):
                                valid_compounds.append(compound)
                    except Exception as e:
                        self.logger.error(f"Error processing web result: {e}")
                        self.stats["errors"].append(f"Web result error: {str(e)}")

                compounds.extend(valid_compounds)
                self.logger.info(f"Found {len(valid_compounds)} valid compounds from web")

                break  # Success, exit retry loop

            except Exception as e:
                self.logger.error(f"Error searching web (attempt {retries + 1}/{self.max_retries}): {e}")
                retries += 1
                if retries == self.max_retries:
                    self.logger.error(f"Failed to search web after {self.max_retries} attempts")
                    self.stats["errors"].append(f"Web search error: {str(e)}")

        return self._deduplicate_compounds(compounds)

    def _validate_compound(self, compound: CompoundData) -> bool:
        """Validate a compound.

        Args:
            compound: Compound to validate

        Returns:
            Whether the compound is valid
        """
        try:
            # Check required fields
            if not compound.smiles:
                return False

            # Validate structure
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                return False

            # Check molecular weight
            if not compound.molecular_weight:
                try:
                    compound.molecular_weight = Chem.Descriptors.ExactMolWt(mol)
                except (ValueError, RuntimeError) as e:
                    self.logger.debug(f"Error calculating molecular weight: {e}")
                    return False

            # Validate with ML
            if compound.name:
                try:
                    # Check compound relevance
                    classification = self.text_classifier(compound.name)[0]
                    if classification["label"] != "relevant" or classification["score"] < 0.8:
                        return False

                    # Check activity type
                    activity = self.activity_classifier(compound.name)[0]
                    if activity["score"] < 0.6:
                        return False

                    # Check safety
                    safety = self.safety_classifier(compound.name)[0]
                    if safety["label"] == "unsafe" and safety["score"] > 0.8:
                        return False

                except (KeyError, IndexError, RuntimeError) as e:
                    self.logger.debug(f"Error running ML validation: {e}")

            return True

        except Exception as e:
            self.logger.error(f"Error validating compound: {e}")
            return False

    def _validate_mention(self, mention: Dict[str, Any]) -> bool:
        """Validate a social media mention.

        Args:
            mention: Mention to validate

        Returns:
            Whether the mention is valid
        """
        try:
            # Check required fields
            if not mention.get("identifier") or not mention.get("context"):
                return False

            # Validate with ML
            try:
                # Check relevance
                classification = self.text_classifier(mention["context"])[0]
                if classification["label"] != "relevant" or classification["score"] < 0.8:
                    return False

                # Check activity type
                activity = self.activity_classifier(mention["context"])[0]
                if activity["score"] < 0.6:
                    return False

                # Check safety
                safety = self.safety_classifier(mention["context"])[0]
                if safety["label"] == "unsafe" and safety["score"] > 0.8:
                    return False

                # Check mechanism
                mechanism = self.mechanism_classifier(mention["context"])[0]
                if mechanism["score"] < 0.6:
                    return False

                return True

            except (KeyError, IndexError, RuntimeError) as e:
                self.logger.debug(f"Error running ML validation: {e}")
                return False

        except Exception as e:
            self.logger.error(f"Error validating mention: {e}")
            return False

    def _deduplicate_compounds(self, compounds: List[CompoundData]) -> List[CompoundData]:
        """Deduplicate compounds based on structure similarity.

        Args:
            compounds: List of compounds to deduplicate

        Returns:
            Deduplicated list of compounds
        """
        try:
            # Sort compounds by molecular weight
            compounds.sort(key=lambda x: x.molecular_weight or 0)

            # Initialize result list
            unique = []
            seen_smiles = set()

            # Process compounds
            for compound in compounds:
                if not compound.smiles:
                    continue

                # Skip exact duplicates
                if compound.smiles in seen_smiles:
                    continue

                # Check similarity to existing compounds
                is_unique = True
                mol = Chem.MolFromSmiles(compound.smiles)
                if not mol:
                    continue

                fp = get_fingerprint(mol)
                for existing in unique:
                    try:
                        existing_mol = Chem.MolFromSmiles(existing.smiles)
                        if not existing_mol:
                            continue
                        existing_fp = get_fingerprint(existing_mol)
                        similarity = calculate_similarity(fp, existing_fp)
                        if similarity > 0.9:  # High similarity threshold
                            is_unique = False
                            break
                    except (ValueError, RuntimeError) as e:
                        self.logger.debug(f"Error calculating similarity: {e}")

                if is_unique:
                    unique.append(compound)
                    seen_smiles.add(compound.smiles)

            return unique

        except Exception as e:
            self.logger.error(f"Error deduplicating compounds: {e}")
            return compounds

    def _add_compounds(self, new_compounds: List[CompoundData]) -> None:
        """Add new compounds to storage.

        Args:
            new_compounds: List of compounds to add
        """
        for compound in new_compounds:
            if compound.smiles not in self.compounds:
                self.compounds[compound.smiles] = compound
                if compound.source_id:
                    self.processed_ids.add(compound.source_id)

    def _load_checkpoint(self, checkpoint_file: str) -> None:
        """Load processing checkpoint.

        Args:
            checkpoint_file: Path to checkpoint file
        """
        try:
            self.logger.info(f"Loading checkpoint from {checkpoint_file}")
            checkpoint = pd.read_pickle(checkpoint_file)
            self.compounds = {c.smiles: c for c in checkpoint["compounds"]}
            self.processed_ids = set(checkpoint["processed_ids"])
            self.stats.update(checkpoint["stats"])
            self.logger.info(f"Loaded {len(self.compounds)} compounds from checkpoint")
        except Exception as e:
            self.logger.error(f"Error loading checkpoint: {str(e)}")
            self.stats["errors"].append(f"Checkpoint load error: {str(e)}")

    def _save_checkpoint(self, checkpoint_file: str) -> None:
        """Save processing checkpoint.

        Args:
            checkpoint_file: Path to checkpoint file
        """
        try:
            checkpoint = {
                "compounds": list(self.compounds.values()),
                "processed_ids": list(self.processed_ids),
                "stats": self.stats,
                "timestamp": datetime.now().isoformat(),
            }
            pd.to_pickle(checkpoint, checkpoint_file)
            self.logger.info(f"Saved checkpoint with {len(self.compounds)} compounds")
        except Exception as e:
            self.logger.error(f"Error saving checkpoint: {str(e)}")
            self.stats["errors"].append(f"Checkpoint save error: {str(e)}")

    def _save_results(self, output_dir: Path) -> None:
        """Save processed compounds.

        Args:
            output_dir: Directory to save results in
        """
        try:
            # Convert compounds to records
            records = []
            for compound in tqdm(
                self.compounds.values(),
                desc="Preparing records",
                unit="compounds",
            ):
                try:
                    record = compound.to_dict()
                    # Convert nested data to JSON strings
                    for key, value in record.items():
                        if isinstance(value, (dict, list)):
                            record[key] = json.dumps(value)
                    records.append(record)
                except Exception as e:
                    self.logger.error(f"Error converting compound {compound.name}: {e}")
                    self.stats["errors"].append(f"Record conversion error: {str(e)}")

            # Create DataFrame
            df = pd.DataFrame.from_records(records)

            # Save to TSV
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = output_dir / f"compounds_{timestamp}.tsv"
            df.to_csv(output_file, sep="\t", index=False)

            self.logger.info(f"Saved {len(records)} compounds to {output_file}")

            # Save pipeline info
            info = {
                "timestamp": timestamp,
                "stats": self.stats,
                "predictor_versions": self.enrichment.get_predictor_versions(),
                "web_clients": self.enrichment.get_web_client_info(),
            }
            info_file = output_dir / f"pipeline_info_{timestamp}.json"
            with open(info_file, "w") as f:
                json.dump(info, f, indent=2)

            self.logger.info(f"Saved pipeline info to {info_file}")

        except Exception as e:
            self.logger.error(f"Error saving results: {str(e)}")
            self.stats["errors"].append(f"Results save error: {str(e)}")
            raise
