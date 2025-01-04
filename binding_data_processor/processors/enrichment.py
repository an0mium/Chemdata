"""Data enrichment processor for compounds.

This module provides:
1. ML predictions for compounds
2. Web data enrichment from multiple sources
3. Social media data harvesting
4. Caching and checkpointing
5. Parallel processing and error handling
6. Model ensembling and confidence scoring
7. Data validation and deduplication
8. Progress tracking and reporting
"""

import logging
import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Set, Tuple
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict

import pandas as pd
import numpy as np
from tqdm import tqdm
import torch
from rdkit import Chem
from transformers import pipeline
from sentence_transformers import SentenceTransformer, util

from binding_data_processor.models.compound import CompoundData
from binding_data_processor.processors.structure.ml.predictors.abuse import (
    AbusePotentialPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.toxicity import (
    ToxicityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.activity import (
    ActivityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.affinity import (
    AffinityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.psychoactive import (
    PsychoactivePredictor,
)
from binding_data_processor.processors.structure.ml.predictors.nootropic import (
    NootropicPredictor,
)
from binding_data_processor.processors.structure.ml.models.ensemble import (
    EnsembleModel,
)
from web_enrichment.data_sources.community import CommunityClient
from web_enrichment.data_sources.swiss import SwissClient
from web_enrichment.data_sources.social import SocialDataHarvester
from web_enrichment.data_sources.chembl import ChEMBLClient
from web_enrichment.data_sources.pubchem import PubChemClient
from web_enrichment.data_sources.regulatory import RegulatoryClient
from web_enrichment.data_sources.web_search import WebSearchClient
from web_enrichment.http_client import HttpClient
from web_enrichment.llm_utils import analyze_content_with_llm


class EnrichmentProcessor:
    """Processor for enriching compound data with predictions and web data."""

    # Cache settings
    CACHE_DIR = Path(".cache/enrichment")
    CACHE_DURATIONS = {
        "predictions": timedelta(days=90),  # Cache predictions longer
        "web_data": timedelta(days=7),  # Refresh web data more frequently
        "social_data": timedelta(days=1),  # Refresh social data daily
        "database_data": timedelta(days=30),  # Refresh database data monthly
        "regulatory_data": timedelta(days=14),  # Refresh regulatory data biweekly
    }

    def __init__(
        self,
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
        """Initialize processor.

        Args:
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
        self.model_dir = model_dir
        self.n_workers = n_workers
        self.batch_size = batch_size
        self.checkpoint_interval = checkpoint_interval
        self.max_retries = max_retries
        self.llm_api_key = llm_api_key
        self.device = device

        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        self.logger.addHandler(handler)

        # Initialize cache
        if cache_dir:
            self.CACHE_DIR = Path(cache_dir)
        for cache_type in self.CACHE_DURATIONS.keys():
            (self.CACHE_DIR / cache_type).mkdir(parents=True, exist_ok=True)

        # Initialize ML predictors
        self.logger.info("Initializing ML predictors...")

        # Base predictors
        self.predictors = {
            "toxicity": ToxicityPredictor(model_dir=f"{model_dir}/toxicity"),
            "abuse": AbusePotentialPredictor(model_dir=f"{model_dir}/abuse"),
            "activity": ActivityPredictor(model_dir=f"{model_dir}/activity"),
            "affinity": AffinityPredictor(model_dir=f"{model_dir}/affinity"),
            "psychoactive": PsychoactivePredictor(
                model_dir=f"{model_dir}/psychoactive"
            ),
            "nootropic": NootropicPredictor(model_dir=f"{model_dir}/nootropic"),
        }

        # Ensemble models
        self.activity_ensemble = EnsembleModel(
            [
                self.predictors["activity"],
                self.predictors["psychoactive"],
                self.predictors["nootropic"],
            ],
            weights=[0.4, 0.3, 0.3],
        )
        self.safety_ensemble = EnsembleModel(
            [self.predictors["toxicity"], self.predictors["abuse"]],
            weights=[0.5, 0.5],
        )

        # Classification models
        self.text_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/enrichment/compound_classifier",
            device=device,
        )
        self.safety_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/enrichment/safety_classifier",
            device=device,
        )
        self.mechanism_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/enrichment/mechanism_classifier",
            device=device,
        )
        self.effect_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/enrichment/effect_classifier",
            device=device,
        )

        # Similarity model
        self.similarity_model = SentenceTransformer(
            "sentence-transformers/all-MiniLM-L6-v2"
        ).to(device)

        # Initialize web clients
        self.logger.info("Initializing web clients...")
        self.http_client = HttpClient()

        # Community data
        self.web_clients = {
            "community": CommunityClient(
                http_client=self.http_client,
                model_dir=model_dir,
                cache_dir=cache_dir,
                n_workers=n_workers,
                max_retries=max_retries,
                device=device,
            ),
            "swiss": SwissClient(
                http_client=self.http_client,
                model_dir=model_dir,
                cache_dir=cache_dir,
                n_workers=n_workers,
                max_retries=max_retries,
                device=device,
            ),
            "social": SocialDataHarvester(
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
            ),
            "chembl": ChEMBLClient(
                http_client=self.http_client,
                cache_dir=cache_dir,
                n_workers=n_workers,
                max_retries=max_retries,
            ),
            "pubchem": PubChemClient(
                http_client=self.http_client,
                cache_dir=cache_dir,
                n_workers=n_workers,
                max_retries=max_retries,
            ),
            "regulatory": RegulatoryClient(
                http_client=self.http_client,
                cache_dir=cache_dir,
                n_workers=n_workers,
                max_retries=max_retries,
            ),
            "web_search": WebSearchClient(
                http_client=self.http_client,
                cache_dir=cache_dir,
                n_workers=n_workers,
                max_retries=max_retries,
            ),
        }

        # Initialize thread pools
        self.prediction_executor = ThreadPoolExecutor(
            max_workers=n_workers, thread_name_prefix="prediction"
        )
        self.web_executor = ThreadPoolExecutor(
            max_workers=n_workers * 2, thread_name_prefix="web"
        )  # More threads for web IO

    def enrich_compounds(
        self,
        compounds: List[CompoundData],
        skip_predictions: bool = False,
        skip_web_data: bool = False,
        checkpoint_file: Optional[str] = None,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Enrich compounds with predictions and web data.

        Args:
            compounds: List of compounds to enrich
            skip_predictions: Whether to skip ML predictions
            skip_web_data: Whether to skip web data enrichment
            checkpoint_file: Optional checkpoint file to resume from
            use_cache: Whether to use cached results

        Returns:
            Dictionary with enrichment statistics
        """
        stats = {
            "total": len(compounds),
            "with_predictions": 0,
            "with_web_data": 0,
            "from_cache": 0,
            "errors": [],
            "start_time": datetime.now().isoformat(),
        }

        # Load checkpoint if exists
        processed_ids = set()
        if checkpoint_file and os.path.exists(checkpoint_file):
            self.logger.info(f"Loading checkpoint from {checkpoint_file}")
            checkpoint = pd.read_pickle(checkpoint_file)
            processed_ids = set(checkpoint["processed_ids"])
            stats.update(checkpoint["stats"])
            self.logger.info(f"Loaded {len(processed_ids)} processed compounds")

        # Process in batches
        for i in range(0, len(compounds), self.batch_size):
            batch = compounds[i : i + self.batch_size]

            # Skip processed compounds
            if processed_ids:
                batch = [c for c in batch if c.source_id not in processed_ids]

            if not batch:
                continue

            # Process batch in parallel
            with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
                futures = []

                # Submit enrichment tasks
                for compound in batch:
                    if not skip_predictions:
                        futures.append(
                            executor.submit(
                                self._add_predictions_with_retry,
                                compound,
                                use_cache=use_cache,
                            )
                        )
                    if not skip_web_data:
                        futures.append(
                            executor.submit(
                                self._add_web_data_with_retry,
                                compound,
                                use_cache=use_cache,
                            )
                        )

                # Collect results with progress bar
                with tqdm(
                    total=len(futures),
                    desc=f"Enriching batch {i // self.batch_size + 1}",
                ) as pbar:
                    for future in as_completed(futures):
                        try:
                            result = future.result()
                            if result:
                                result_type, success, from_cache = result
                                if success:
                                    stats[f"with_{result_type}"] += 1
                                    if from_cache:
                                        stats["from_cache"] += 1
                                    if compound.source_id:
                                        processed_ids.add(compound.source_id)
                        except Exception as e:
                            stats["errors"].append(str(e))
                            self.logger.error(f"Error in enrichment: {e}")
                        pbar.update(1)

            # Save checkpoint
            if checkpoint_file and len(processed_ids) % self.checkpoint_interval == 0:
                self._save_checkpoint(processed_ids, stats, checkpoint_file)

        # Save final checkpoint
        if checkpoint_file:
            self._save_checkpoint(processed_ids, stats, checkpoint_file)

        # Update final stats
        stats["end_time"] = datetime.now().isoformat()
        return stats

    def _add_predictions_with_retry(
        self, compound: CompoundData, use_cache: bool = True
    ) -> Optional[Tuple[str, bool, bool]]:
        """Add ML predictions with retry logic and caching.

        Args:
            compound: Compound to add predictions to
            use_cache: Whether to use cached results

        Returns:
            Tuple of (prediction_type, success, from_cache) or None if failed
        """
        # Check cache first
        if use_cache:
            cache_key = f"predictions_{compound.smiles}"
            cached = self._get_cached_result(cache_key, "predictions")
            if cached:
                for pred_type, prediction in cached.items():
                    setattr(compound, f"{pred_type}_predictions", prediction)
                return "predictions", True, True

        # Try with retries
        for attempt in range(self.max_retries):
            try:
                result = self._add_predictions(compound)
                if result:
                    pred_type, success = result
                    if success and use_cache:
                        # Cache successful predictions
                        predictions = {
                            pred_type: getattr(compound, f"{pred_type}_predictions")
                            for pred_type in self.predictors.keys()
                            if hasattr(compound, f"{pred_type}_predictions")
                        }
                        self._cache_result(cache_key, predictions, "predictions")
                    return pred_type, success, False
            except Exception as e:
                if attempt == self.max_retries - 1:
                    raise
                self.logger.warning(
                    f"Prediction attempt {attempt + 1} failed for {compound.name}: {e}"
                )

        return None

    def _add_web_data_with_retry(
        self, compound: CompoundData, use_cache: bool = True
    ) -> Optional[Tuple[str, bool, bool]]:
        """Add web data with retry logic and caching.

        Args:
            compound: Compound to add web data to
            use_cache: Whether to use cached results

        Returns:
            Tuple of (web_data_type, success, from_cache) or None if failed
        """
        # Check cache first
        if use_cache:
            cache_key = f"web_data_{compound.name}_{compound.smiles}"
            cached = self._get_cached_result(cache_key, "web_data")
            if cached:
                for data_type, data in cached.items():
                    setattr(compound, f"{data_type}_data", data)
                return "web_data", True, True

        # Try with retries
        for attempt in range(self.max_retries):
            try:
                result = self._add_web_data(compound)
                if result:
                    data_type, success = result
                    if success and use_cache:
                        # Cache successful web data
                        web_data = {
                            client_type: getattr(compound, f"{client_type}_data")
                            for client_type in self.web_clients.keys()
                            if hasattr(compound, f"{client_type}_data")
                        }
                        self._cache_result(cache_key, web_data, "web_data")
                    return data_type, success, False
            except Exception as e:
                if attempt == self.max_retries - 1:
                    raise
                self.logger.warning(
                    f"Web data attempt {attempt + 1} failed for {compound.name}: {e}"
                )

        return None

    def _add_predictions(self, compound: CompoundData) -> Optional[Tuple[str, bool]]:
        """Add ML predictions to compound.

        Args:
            compound: Compound to add predictions to

        Returns:
            Tuple of (prediction_type, success) or None if failed
        """
        try:
            # Validate SMILES
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                return None

            # Add predictions from each predictor
            success = False
            for pred_type, predictor in tqdm(
                self.predictors.items(),
                desc=f"Running predictions for {compound.name}",
                leave=False,
            ):
                try:
                    prediction = predictor.predict(compound.smiles)
                    if prediction:
                        setattr(compound, f"{pred_type}_predictions", prediction)
                        success = True
                except Exception as e:
                    self.logger.error(
                        f"Error in {pred_type} prediction for {compound.name}: {e}"
                    )

            # Add ensemble predictions
            try:
                # Activity ensemble
                activity_pred = self.activity_ensemble.predict(compound.smiles)
                if activity_pred:
                    compound.activity_ensemble_predictions = activity_pred
                    success = True

                # Safety ensemble
                safety_pred = self.safety_ensemble.predict(compound.smiles)
                if safety_pred:
                    compound.safety_ensemble_predictions = safety_pred
                    success = True

            except Exception as e:
                self.logger.error(
                    f"Error in ensemble predictions for {compound.name}: {e}"
                )

            return "predictions", success

        except Exception as e:
            self.logger.error(f"Error adding predictions for {compound.name}: {e}")
            return None

    def _add_web_data(self, compound: CompoundData) -> Optional[Tuple[str, bool]]:
        """Add web data to compound.

        Args:
            compound: Compound to add web data to

        Returns:
            Tuple of (web_data_type, success) or None if failed
        """
        try:
            success = False

            # Add data from each web client
            for client_type, client in tqdm(
                self.web_clients.items(),
                desc=f"Getting web data for {compound.name}",
                leave=False,
            ):
                try:
                    data = client.get_compound_data(compound)
                    if data and self._validate_web_data(data, client_type):
                        setattr(compound, f"{client_type}_data", data)
                        success = True
                except Exception as e:
                    self.logger.error(
                        f"Error getting {client_type} data for {compound.name}: {e}"
                    )

            # Add social media enrichment if available
            if "social" in self.web_clients:
                try:
                    social_data = self.web_clients["social"].enrich_compound_data(
                        compound
                    )
                    if social_data and self._validate_web_data(social_data, "social"):
                        compound.social_data = social_data
                        success = True
                except Exception as e:
                    self.logger.error(
                        f"Error getting social data for {compound.name}: {e}"
                    )

            # Use LLM to analyze web data if available
            if success and self.llm_api_key:
                try:
                    llm_analysis = analyze_content_with_llm(
                        {
                            "name": compound.name,
                            "smiles": compound.smiles,
                            "data": {
                                client_type: getattr(compound, f"{client_type}_data")
                                for client_type in self.web_clients.keys()
                                if hasattr(compound, f"{client_type}_data")
                            },
                        },
                        self.llm_api_key,
                    )
                    if llm_analysis:
                        compound.llm_analysis = llm_analysis
                except Exception as e:
                    self.logger.error(f"Error in LLM analysis for {compound.name}: {e}")

            return "web_data", success

        except Exception as e:
            self.logger.error(f"Error adding web data for {compound.name}: {e}")
            return None

    def _validate_web_data(self, data: Dict[str, Any], source: str) -> bool:
        """Validate web data from a source.

        Args:
            data: Web data to validate
            source: Source of the web data

        Returns:
            Whether the data is valid
        """
        try:
            if not isinstance(data, dict):
                return False

            # Check required fields based on source
            if source == "community":
                return bool(
                    data.get("total_mentions")
                    and isinstance(data.get("reported_effects"), list)
                )
            elif source == "swiss":
                return bool(data.get("swiss_id") and data.get("activity_data"))
            elif source == "social":
                return bool(
                    data.get("total_mentions")
                    and data.get("sentiment_summary")
                    and data.get("reported_effects")
                )
            elif source == "chembl":
                return bool(data.get("chembl_id") and data.get("assays"))
            elif source == "pubchem":
                return bool(data.get("cid") and data.get("properties"))
            elif source == "regulatory":
                return bool(data.get("status") and data.get("references"))
            elif source == "web_search":
                return bool(data.get("sources") and data.get("mentions"))

            return True

        except Exception as e:
            self.logger.error(f"Error validating {source} data: {e}")
            return False

    def _get_cached_result(
        self, cache_key: str, cache_type: str
    ) -> Optional[Dict[str, Any]]:
        """Get cached result if available and not expired.

        Args:
            cache_key: Cache key to look up
            cache_type: Type of cached data

        Returns:
            Cached result or None if not found/expired
        """
        cache_file = self.CACHE_DIR / cache_type / f"{cache_key}.json"
        if cache_file.exists():
            try:
                data = json.loads(cache_file.read_text())
                cached_time = datetime.fromisoformat(data["cached_at"])
                if datetime.now() - cached_time < self.CACHE_DURATIONS.get(
                    cache_type, timedelta(days=7)
                ):
                    return data["result"]
            except Exception as e:
                self.logger.error(f"Error reading cache file {cache_file}: {e}")
        return None

    def _cache_result(
        self, cache_key: str, result: Dict[str, Any], cache_type: str
    ) -> None:
        """Cache result with timestamp.

        Args:
            cache_key: Cache key to store under
            result: Result data to cache
            cache_type: Type of data to cache
        """
        try:
            cache_file = self.CACHE_DIR / cache_type / f"{cache_key}.json"
            data = {
                "cached_at": datetime.now().isoformat(),
                "result": result,
            }
            cache_file.write_text(json.dumps(data, indent=2))
        except Exception as e:
            self.logger.error(f"Error writing cache file {cache_file}: {e}")

    def _save_checkpoint(
        self,
        processed_ids: Set[str],
        stats: Dict[str, Any],
        checkpoint_file: str,
    ) -> None:
        """Save processing checkpoint.

        Args:
            processed_ids: Set of processed compound IDs
            stats: Current processing statistics
            checkpoint_file: Path to checkpoint file
        """
        try:
            checkpoint = {
                "processed_ids": list(processed_ids),
                "stats": stats,
                "timestamp": datetime.now().isoformat(),
            }
            pd.to_pickle(checkpoint, checkpoint_file)
            self.logger.info(
                f"Saved checkpoint with {len(processed_ids)} compounds to {checkpoint_file}"
            )
        except Exception as e:
            self.logger.error(f"Error saving checkpoint: {str(e)}")

    def get_predictor_versions(self) -> Dict[str, str]:
        """Get versions of all ML predictors.

        Returns:
            Dictionary mapping predictor names to versions
        """
        versions = {}
        for name, predictor in self.predictors.items():
            try:
                versions[name] = predictor.version
            except Exception:
                versions[name] = "Unknown"
        return versions

    def get_web_client_info(self) -> Dict[str, Dict[str, Any]]:
        """Get information about web clients.

        Returns:
            Dictionary with web client information
        """
        info = {}
        for name, client in self.web_clients.items():
            try:
                info[name] = {
                    "type": client.__class__.__name__,
                    "available": bool(client),
                }
            except Exception:
                info[name] = {
                    "type": "Unknown",
                    "available": False,
                }
        return info
