#!/usr/bin/env python3
"""Process BindingDB data and harvest additional compounds.

This script:
1. Processes BindingDB data to extract receptor ligands
2. Searches web sources for additional compounds
3. Monitors social media for emerging compounds
4. Generates ML predictions for all compounds
5. Enriches data from multiple sources
6. Exposes data through web interface
"""

import concurrent.futures
import json
import logging
import os
import re
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm
import tweepy
import praw
from bs4 import BeautifulSoup

from binding_data_processor.models.compound import CompoundData
from binding_data_processor.processors.activity_analysis import ActivityAnalysis
from binding_data_processor.processors.activity_types import ActivityTypes
from binding_data_processor.processors.structure import StructureProcessor
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
from binding_data_processor.web.dashboard import Dashboard
from web_enrichment.data_sources.community import CommunityClient
from web_enrichment.data_sources.swiss import SwissClient
from web_enrichment.http_client import HttpClient
from web_enrichment.llm_utils import LLMProcessor
from checkpoint_manager import CheckpointManager
from cache_manager import CacheManager
from config import (
    BATCH_SIZE,
    MAX_WORKERS,
    TWITTER_API_KEY,
    TWITTER_API_SECRET,
    REDDIT_CLIENT_ID,
    REDDIT_CLIENT_SECRET,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("compound_processing.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


class CompoundProcessor:
    """Process compounds from multiple sources with ML predictions."""

    # Target patterns for filtering compounds
    TARGET_PATTERNS = {
        "serotonin": r"5-HT\d*[A-Z]?|serotonin|SLC6A4",
        "dopamine": r"D\d+|dopamine|DAT|SLC6A3",
        "norepinephrine": r"norepinephrine|NET|SLC6A2",
        "gaba": r"GABA[A-Z]?\d*|SLC6A1",
        "glutamate": r"glutamate|NMDA|AMPA|mGluR|GluR",
        "opioid": r"[μμκδ]?-?opioid|MOR|KOR|DOR",
        "cannabinoid": r"cannabinoid|CB[12]",
        "psychedelic": r"psychedelic|hallucinogen|entheogen",
    }

    # Web sources to search
    WEB_SOURCES = [
        "pubchem",
        "chembl",
        "drugbank",
        "wikipedia",
        "erowid",
        "psychonautwiki",
    ]

    # Subreddits to monitor
    SUBREDDITS = [
        "researchchemicals",
        "DrugNerds",
        "Nootropics",
        "psychopharmacology",
    ]

    def __init__(
        self,
        data_dir: str = "../data",
        model_dir: str = "../models",
        n_workers: int = 4,
        batch_size: int = 100,
    ):
        """Initialize processor with ML models and clients."""
        self.data_dir = Path(data_dir)
        self.model_dir = Path(model_dir)
        self.n_workers = n_workers
        self.batch_size = batch_size

        # Initialize managers
        self.cache = CacheManager()
        self.checkpoint = CheckpointManager()

        # Initialize processors
        self.activity_types = ActivityTypes()
        self.activity_analyzer = ActivityAnalysis()
        self.structure_processor = StructureProcessor()
        self.llm_processor = LLMProcessor()

        # Initialize ML predictors
        self.toxicity_predictor = ToxicityPredictor()
        self.abuse_predictor = AbusePotentialPredictor()
        self.activity_predictor = ActivityPredictor()
        self.affinity_predictor = AffinityPredictor()

        # Initialize web clients
        self.http_client = HttpClient()
        self.community_client = CommunityClient(self.http_client)
        self.swiss_client = SwissClient(self.http_client)

        # Initialize social media clients
        self._init_social_media_clients()

        # Initialize dashboard
        self.dashboard = Dashboard()

        # Track statistics
        self.stats = {
            "total_compounds": 0,
            "bindingdb_compounds": 0,
            "web_compounds": 0,
            "social_compounds": 0,
            "with_predictions": 0,
            "with_web_data": 0,
            "errors": [],
        }

    def _init_social_media_clients(self):
        """Initialize Twitter and Reddit API clients."""
        try:
            # Twitter client
            auth = tweepy.OAuthHandler(TWITTER_API_KEY, TWITTER_API_SECRET)
            self.twitter = tweepy.API(auth)

            # Reddit client
            self.reddit = praw.Reddit(
                client_id=REDDIT_CLIENT_ID,
                client_secret=REDDIT_CLIENT_SECRET,
                user_agent="ChemDataProcessor/1.0",
            )
        except Exception as e:
            logger.error(f"Error initializing social media clients: {str(e)}")
            self.twitter = None
            self.reddit = None

    def process_all_sources(
        self,
        bindingdb_file: Optional[str] = None,
        output_dir: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Process compounds from all sources.

        Args:
            bindingdb_file: Optional path to BindingDB file
            output_dir: Optional directory for output files

        Returns:
            Processing statistics
        """
        try:
            output_dir = Path(output_dir or self.data_dir)
            os.makedirs(output_dir, exist_ok=True)

            # Process BindingDB compounds
            bindingdb_compounds = self._process_bindingdb(bindingdb_file)
            self.stats["bindingdb_compounds"] = len(bindingdb_compounds)

            # Search web sources
            web_compounds = self._search_web_sources()
            self.stats["web_compounds"] = len(web_compounds)

            # Monitor social media
            social_compounds = self._monitor_social_media()
            self.stats["social_compounds"] = len(social_compounds)

            # Combine all compounds
            all_compounds = bindingdb_compounds + web_compounds + social_compounds
            self.stats["total_compounds"] = len(all_compounds)

            # Process combined compounds
            processed = self._process_compounds(all_compounds)

            # Save results
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self._save_results(processed, output_dir / f"compounds_{timestamp}")

            # Update dashboard
            self.dashboard.update_data(processed)

            return self.stats

        except Exception as e:
            logger.error(f"Error processing compounds: {str(e)}")
            raise

    def _process_bindingdb(
        self, input_file: Optional[str] = None
    ) -> List[CompoundData]:
        """Process BindingDB compounds.

        Args:
            input_file: Optional path to BindingDB file

        Returns:
            List of processed compounds
        """
        try:
            # Download/load BindingDB file
            if not input_file:
                input_file = self._download_bindingdb()

            # Read compounds
            df = pd.read_csv(input_file, sep="\t")
            logger.info(f"Found {len(df)} compounds in BindingDB")

            # Filter by targets
            target_mask = df["Target Name"].str.contains(
                "|".join(self.TARGET_PATTERNS.values()), case=False, na=False
            )
            df = df[target_mask]
            logger.info(f"Found {len(df)} compounds matching target patterns")

            # Convert to CompoundData objects
            compounds = []
            for _, row in df.iterrows():
                try:
                    compound = CompoundData(
                        name=row["Ligand Name"],
                        smiles=row["Ligand SMILES"],
                        source="BindingDB",
                        binding_data=[
                            {
                                "target": row["Target Name"],
                                "activity_type": row["Measurement Type"],
                                "activity_value": float(row["Measurement Value"]),
                                "activity_unit": row["Measurement Unit"],
                                "reference": row["Article DOI"],
                            }
                        ],
                    )
                    compounds.append(compound)
                except Exception as e:
                    logger.error(f"Error processing row: {str(e)}")
                    continue

            return compounds

        except Exception as e:
            logger.error(f"Error processing BindingDB: {str(e)}")
            return []

    def _search_web_sources(self) -> List[CompoundData]:
        """Search web sources for compounds.

        Returns:
            List of compounds found from web sources
        """
        compounds = []
        for source in self.WEB_SOURCES:
            try:
                # Search each target pattern
                for target, pattern in self.TARGET_PATTERNS.items():
                    logger.info(f"Searching {source} for {target} compounds")

                    # Get compounds from source
                    source_compounds = self.community_client.search_compounds(
                        source=source, query=pattern
                    )

                    # Process found compounds
                    for compound_data in source_compounds:
                        try:
                            compound = CompoundData(
                                name=compound_data["name"],
                                smiles=compound_data.get("smiles"),
                                source=source,
                                web_data=compound_data,
                            )
                            compounds.append(compound)
                        except Exception as e:
                            logger.error(
                                f"Error processing {source} compound: {str(e)}"
                            )
                            continue

            except Exception as e:
                logger.error(f"Error searching {source}: {str(e)}")
                continue

        return compounds

    def _monitor_social_media(self) -> List[CompoundData]:
        """Monitor social media for compounds.

        Returns:
            List of compounds found from social media
        """
        compounds = []

        try:
            # Monitor Reddit
            if self.reddit:
                for subreddit in self.SUBREDDITS:
                    try:
                        logger.info(f"Monitoring r/{subreddit}")

                        # Get recent posts and comments
                        subreddit = self.reddit.subreddit(subreddit)
                        texts = []

                        # Get posts
                        for post in subreddit.new(limit=100):
                            texts.append(post.title)
                            texts.append(post.selftext)

                            # Get comments
                            post.comments.replace_more(limit=0)
                            for comment in post.comments.list():
                                texts.append(comment.body)

                        # Extract compounds using LLM
                        compounds.extend(self.llm_processor.extract_compounds(texts))

                    except Exception as e:
                        logger.error(f"Error monitoring r/{subreddit}: {str(e)}")
                        continue

            # Monitor Twitter
            if self.twitter:
                try:
                    logger.info("Monitoring Twitter")

                    # Search tweets
                    queries = [
                        "research chemical",
                        "novel compound",
                        "new synthesis",
                        "receptor binding",
                    ]

                    texts = []
                    for query in queries:
                        tweets = self.twitter.search_tweets(
                            q=query, lang="en", count=100
                        )
                        texts.extend(tweet.text for tweet in tweets)

                    # Extract compounds using LLM
                    compounds.extend(self.llm_processor.extract_compounds(texts))

                except Exception as e:
                    logger.error(f"Error monitoring Twitter: {str(e)}")

        except Exception as e:
            logger.error(f"Error monitoring social media: {str(e)}")

        return compounds

    def _process_compounds(self, compounds: List[CompoundData]) -> List[CompoundData]:
        """Process compounds with predictions and enrichment.

        Args:
            compounds: List of compounds to process

        Returns:
            List of processed compounds
        """
        processed = []

        # Process in batches
        for i in tqdm(
            range(0, len(compounds), self.batch_size), desc="Processing compounds"
        ):
            batch = compounds[i : i + self.batch_size]

            # Process batch in parallel
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=self.n_workers
            ) as executor:
                futures = []
                for compound in batch:
                    future = executor.submit(self._process_single_compound, compound)
                    futures.append(future)

                # Collect results
                for future in concurrent.futures.as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            processed.append(result)
                    except Exception as e:
                        logger.error(f"Error processing compound: {str(e)}")
                        continue

        return processed

    def _process_single_compound(
        self, compound: CompoundData
    ) -> Optional[CompoundData]:
        """Process a single compound.

        Args:
            compound: Compound to process

        Returns:
            Processed compound or None if failed
        """
        try:
            # Add ML predictions if SMILES available
            if compound.smiles:
                # Toxicity predictions
                toxicity = self.toxicity_predictor.predict(compound.smiles)
                if toxicity:
                    compound.toxicity_predictions = toxicity

                # Abuse potential predictions
                abuse = self.abuse_predictor.predict(compound.smiles)
                if abuse:
                    compound.abuse_predictions = abuse

                # Activity predictions
                activity = self.activity_predictor.predict(compound.smiles)
                if activity:
                    compound.activity_predictions = activity

                # Binding predictions
                binding = self.affinity_predictor.predict(compound.smiles)
                if binding:
                    compound.binding_predictions = binding

                self.stats["with_predictions"] += 1

            # Add web data enrichment
            web_data = self.community_client.get_compound_data(
                compound.name, compound.smiles
            )
            if web_data:
                compound.web_data = web_data
                self.stats["with_web_data"] += 1

            return compound

        except Exception as e:
            logger.error(f"Error processing compound {compound.name}: {str(e)}")
            return None

    def _save_results(self, compounds: List[CompoundData], output_prefix: str) -> None:
        """Save processed compounds.

        Args:
            compounds: List of compounds to save
            output_prefix: Prefix for output files
        """
        try:
            # Save full JSON data
            with open(f"{output_prefix}.json", "w") as f:
                json.dump([c.to_dict() for c in compounds], f, indent=2)

            # Save TSV with main fields
            df = pd.DataFrame(
                [
                    {
                        "name": c.name,
                        "smiles": c.smiles,
                        "source": c.source,
                        "toxicity_score": (
                            c.toxicity_predictions.get("score")
                            if c.toxicity_predictions
                            else None
                        ),
                        "abuse_potential": (
                            c.abuse_predictions.get("score")
                            if c.abuse_predictions
                            else None
                        ),
                        "binding_targets": (
                            ";".join(c.binding_predictions.keys())
                            if c.binding_predictions
                            else None
                        ),
                        "activity_types": (
                            ";".join(c.activity_predictions.keys())
                            if c.activity_predictions
                            else None
                        ),
                    }
                    for c in compounds
                ]
            )
            df.to_csv(f"{output_prefix}.tsv", sep="\t", index=False)

            logger.info(f"Saved {len(compounds)} compounds to {output_prefix}")

        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")

    def _download_bindingdb(self) -> str:
        """Download BindingDB data file.

        Returns:
            Path to downloaded file
        """
        url = "https://bindingdb.org/bind/downloads/BindingDB_All.tsv"
        output_file = self.data_dir / "BindingDB_All.tsv"

        try:
            if not output_file.exists():
                logger.info("Downloading BindingDB data...")

                # Download with progress bar
                response = requests.get(url, stream=True)
                total_size = int(response.headers.get("content-length", 0))

                with (
                    open(output_file, "wb") as f,
                    tqdm(
                        total=total_size, unit="iB", unit_scale=True, desc="Downloading"
                    ) as pbar,
                ):
                    for chunk in response.iter_content(chunk_size=8192):
                        size = f.write(chunk)
                        pbar.update(size)

            return str(output_file)

        except Exception as e:
            logger.error(f"Error downloading BindingDB: {str(e)}")
            raise


def main():
    """Run compound processing pipeline."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Process compounds from multiple sources"
    )
    parser.add_argument("--bindingdb", type=str, help="Optional path to BindingDB file")
    parser.add_argument("--output-dir", type=str, help="Optional output directory")
    parser.add_argument(
        "--workers", type=int, default=4, help="Number of worker threads"
    )
    parser.add_argument(
        "--batch-size", type=int, default=100, help="Batch size for processing"
    )
    args = parser.parse_args()

    processor = CompoundProcessor(n_workers=args.workers, batch_size=args.batch_size)
    stats = processor.process_all_sources(
        bindingdb_file=args.bindingdb, output_dir=args.output_dir
    )

    print("\nProcessing Statistics:")
    print(f"Total compounds: {stats['total_compounds']}")
    print(f"  From BindingDB: {stats['bindingdb_compounds']}")
    print(f"  From web sources: {stats['web_compounds']}")
    print(f"  From social media: {stats['social_compounds']}")
    print(f"With predictions: {stats['with_predictions']}")
    print(f"With web data: {stats['with_web_data']}")
    print(f"Errors: {len(stats['errors'])}")


if __name__ == "__main__":
    main()
