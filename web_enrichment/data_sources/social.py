"""Social media data harvesting for compound discovery.

This module provides:
1. Multi-platform social media monitoring
2. Chemical entity extraction and validation
3. ML-based mention analysis
4. Data enrichment and standardization
5. Caching and error handling
6. Model ensembling and confidence scoring
7. Parallel processing and error recovery
8. Progress tracking and reporting
"""

import re
import json
import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Any, Set, Tuple
from dataclasses import asdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import praw
import tweepy
import discord
from atproto import Client as AtprotoClient
from mastodon import Mastodon
from matrix_client.client import MatrixClient
from bs4 import BeautifulSoup
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from transformers import pipeline, AutoTokenizer, AutoModelForTokenClassification
import numpy as np
from tqdm import tqdm
import torch
from sentence_transformers import SentenceTransformer, util

from logger import LogManager
from web_enrichment.http_client import HttpClient
from web_enrichment.data_sources.community import CommunityClient
from binding_data_processor.models.compound import CompoundData
from binding_data_processor.processors.structure.similarity import (
    calculate_similarity,
    get_fingerprint,
)
from binding_data_processor.processors.structure.ml.models.ensemble import (
    EnsembleModel,
)

logger = LogManager().get_logger("web_enrichment.data_sources.social")


class DiscordClient(discord.Client):
    """Custom Discord client for async message history retrieval."""

    async def get_channel_history(
        self,
        channel_name: str,
        limit: int,
        after: datetime,
    ) -> List[Dict[str, Any]]:
        """Get message history from a channel.

        Args:
            channel_name: Name of channel to search
            limit: Maximum number of messages to retrieve
            after: Only get messages after this time

        Returns:
            List of message data dictionaries
        """
        messages = []
        for guild in self.guilds:
            channel = discord.utils.get(guild.channels, name=channel_name)
            if channel:
                async for message in channel.history(limit=limit, after=after):
                    messages.append(
                        {
                            "content": message.content,
                            "source": f"discord/{channel_name}",
                            "url": message.jump_url,
                            "date": message.created_at,
                        }
                    )
        return messages


class MatrixHandler:
    """Handler for Matrix chat platform."""

    def __init__(self, homeserver_url: str, token: str):
        """Initialize Matrix handler.

        Args:
            homeserver_url: URL of Matrix homeserver
            token: Access token for authentication
        """
        self.client = MatrixClient(homeserver_url)
        self.client.login_with_token(token)

    def get_room_history(
        self,
        room_id: str,
        limit: int = 100,
        since: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Get message history from a Matrix room.

        Args:
            room_id: ID of room to search
            limit: Maximum number of messages to retrieve
            since: Token to start fetching from

        Returns:
            List of message data dictionaries
        """
        room = self.client.join_room(room_id)
        messages = []
        for event in room.get_messages(limit=limit, start=since):
            if event["type"] == "m.room.message":
                messages.append(
                    {
                        "content": event["content"]["body"],
                        "source": f"matrix/{room_id}",
                        "url": f"{room.room_id}/{event['event_id']}",
                        "date": datetime.fromtimestamp(
                            event["origin_server_ts"] / 1000
                        ),
                    }
                )
        return messages


class SocialDataHarvester:
    """Harvester for social media data about compounds."""

    # Cache settings
    CACHE_DIR = Path(".cache/social")
    CACHE_DURATIONS = {
        "mentions": timedelta(days=1),  # Refresh mentions daily
        "analysis": timedelta(days=7),  # Cache analysis for a week
        "validation": timedelta(days=30),  # Cache validation for a month
    }

    # Search settings
    SUBREDDITS = [
        "researchchemicals",
        "Drugs",
        "DrugNerds",
        "Nootropics",
        "Psychonaut",
        "PsychedelicStudies",
        "RationalPsychonaut",
        "AskDrugNerds",
        "Neurochemistry",
        "pharmacology",
        "neuroscience",
        "chemistry",
        "MedicinalChemistry",
        "drugdesign",
        "TheeHive",
        "BabyBees",
        "ObscureDrugs",
        "PsychonautUnion",
        "NootropicsDepot",
        "StackAdvice",
        "Peptides",
        "HPPD",
        "ReagentTesting",
        "DrugAnalysis",
    ]

    TWITTER_QUERIES = [
        "new psychedelic",
        "new antidepressant",
        "research chemical",
        "novel compound",
        "nootropic",
        "new dissociative",
        "new stimulant",
        "5-HT2 ligand",
        "NMDA antagonist",
        "serotonergic",
        "dopaminergic",
        "pharmacology",
        "medicinal chemistry",
        "drug design",
        "psychopharmacology",
        "neuropharmacology",
        "receptor binding",
        "binding affinity",
        "structure activity",
        "SAR study",
        "novel scaffold",
        "drug discovery",
        "lead compound",
        "drug development",
        "clinical candidate",
    ]

    DISCORD_CHANNELS = [
        "research-chemicals",
        "psychedelics",
        "nootropics",
        "harm-reduction",
        "pharmacology",
        "chemistry",
        "drug-design",
        "medicinal-chemistry",
        "synthesis",
        "analysis",
        "neuroscience",
        "pharmacokinetics",
        "drug-development",
        "clinical-trials",
        "safety-discussion",
        "mechanism-of-action",
        "structure-activity",
        "receptor-binding",
    ]

    BLUESKY_QUERIES = [
        "psychedelic",
        "research chemical",
        "nootropic",
        "novel compound",
        "pharmacology",
        "drug design",
        "medicinal chemistry",
        "drug development",
        "neuroscience",
        "receptor binding",
        "clinical trials",
        "drug discovery",
        "lead compound",
        "SAR analysis",
        "binding affinity",
    ]

    MASTODON_INSTANCES = [
        "mastodon.social",
        "scholar.social",
        "fediscience.org",
        "scicomm.xyz",
        "genomic.social",
    ]

    MATRIX_ROOMS = [
        "#psychedelics:matrix.org",
        "#drugnerds:matrix.org",
        "#chemistry:matrix.org",
        "#neuroscience:matrix.org",
        "#pharmacology:matrix.org",
        "#research:matrix.org",
    ]

    def __init__(
        self,
        reddit_client_id: Optional[str] = None,
        reddit_client_secret: Optional[str] = None,
        twitter_api_key: Optional[str] = None,
        twitter_api_secret: Optional[str] = None,
        discord_token: Optional[str] = None,
        bluesky_handle: Optional[str] = None,
        bluesky_password: Optional[str] = None,
        mastodon_tokens: Optional[Dict[str, str]] = None,
        matrix_token: Optional[str] = None,
        http_client: Optional[HttpClient] = None,
        model_dir: str = "../models",
        cache_dir: Optional[str] = None,
        n_workers: int = 4,
        max_retries: int = 3,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
    ):
        """Initialize social data harvester.

        Args:
            reddit_client_id: Optional Reddit API client ID
            reddit_client_secret: Optional Reddit API client secret
            twitter_api_key: Optional Twitter API key
            twitter_api_secret: Optional Twitter API secret
            discord_token: Optional Discord bot token
            bluesky_handle: Optional Bluesky handle
            bluesky_password: Optional Bluesky password
            mastodon_tokens: Optional dict mapping instances to access tokens
            matrix_token: Optional Matrix access token
            http_client: Optional HTTP client for web requests
            model_dir: Directory containing ML models
            cache_dir: Optional custom cache directory
            n_workers: Number of worker threads
            max_retries: Maximum number of retries for failed operations
            device: Device to run ML models on
        """
        self.http_client = http_client or HttpClient()
        self.n_workers = n_workers
        self.max_retries = max_retries
        self.device = device

        # Initialize cache
        if cache_dir:
            self.CACHE_DIR = Path(cache_dir)
        for cache_type in self.CACHE_DURATIONS.keys():
            (self.CACHE_DIR / cache_type).mkdir(parents=True, exist_ok=True)

        # Initialize social media clients
        self.reddit = self._init_reddit(reddit_client_id, reddit_client_secret)
        self.twitter = self._init_twitter(twitter_api_key, twitter_api_secret)
        self.discord = self._init_discord(discord_token)
        self.bluesky = self._init_bluesky(bluesky_handle, bluesky_password)
        self.mastodon = self._init_mastodon(mastodon_tokens)
        self.matrix = self._init_matrix(matrix_token)

        # Initialize community client
        self.community = CommunityClient(
            http_client=self.http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            n_workers=n_workers,
            max_retries=max_retries,
            device=device,
        )

        # Initialize ML models
        logger.info("Loading ML models...")

        # Chemical NER models
        self.ner_tokenizer = AutoTokenizer.from_pretrained(
            "allenai/scibert_scivocab_uncased"
        )
        self.ner_model = AutoModelForTokenClassification.from_pretrained(
            f"{model_dir}/ner/chemical"
        ).to(device)
        self.ner_bio = AutoModelForTokenClassification.from_pretrained(
            f"{model_dir}/ner/bio"
        ).to(device)
        self.ner_drug = AutoModelForTokenClassification.from_pretrained(
            f"{model_dir}/ner/drug"
        ).to(device)

        # Classification models
        self.text_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/social/compound_classifier",
            device=device,
        )
        self.abuse_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/social/abuse_classifier",
            device=device,
        )
        self.effect_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/social/effect_classifier",
            device=device,
        )
        self.safety_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/social/safety_classifier",
            device=device,
        )
        self.mechanism_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/social/mechanism_classifier",
            device=device,
        )
        self.relationship_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/social/relationship_classifier",
            device=device,
        )

        # Sentiment and similarity models
        self.sentiment_model = pipeline(
            "sentiment-analysis",
            model="finiteautomata/bertweet-base-sentiment-analysis",
            device=device,
        )
        self.similarity_model = SentenceTransformer(
            "sentence-transformers/all-MiniLM-L6-v2"
        ).to(device)

        # Ensemble models
        self.ner_ensemble = EnsembleModel(
            [self.ner_model, self.ner_bio, self.ner_drug],
            weights=[0.4, 0.3, 0.3],
        )
        self.safety_ensemble = EnsembleModel(
            [self.safety_classifier, self.abuse_classifier],
            weights=[0.5, 0.5],
        )

        # Initialize chemical standardizer
        self.standardizer = rdMolStandardize.Standardizer()

        # Compile regex patterns
        self._compile_patterns()

    def _init_reddit(
        self, client_id: Optional[str], client_secret: Optional[str]
    ) -> Optional[praw.Reddit]:
        """Initialize Reddit client."""
        if client_id and client_secret:
            try:
                return praw.Reddit(
                    client_id=client_id,
                    client_secret=client_secret,
                    user_agent="ChemData/1.0",
                )
            except Exception as e:
                logger.error(f"Error initializing Reddit client: {e}")
        return None

    def _init_twitter(
        self, api_key: Optional[str], api_secret: Optional[str]
    ) -> Optional[tweepy.API]:
        """Initialize Twitter client."""
        if api_key and api_secret:
            try:
                auth = tweepy.OAuthHandler(api_key, api_secret)
                return tweepy.API(auth)
            except Exception as e:
                logger.error(f"Error initializing Twitter client: {e}")
        return None

    def _init_discord(self, token: Optional[str]) -> Optional[DiscordClient]:
        """Initialize Discord client."""
        if token:
            try:
                intents = discord.Intents.default()
                intents.message_content = True
                client = DiscordClient(intents=intents)
                client.run(token)
                return client
            except Exception as e:
                logger.error(f"Error initializing Discord client: {e}")
        return None

    def _init_bluesky(
        self, handle: Optional[str], password: Optional[str]
    ) -> Optional[AtprotoClient]:
        """Initialize Bluesky client."""
        if handle and password:
            try:
                client = AtprotoClient()
                client.login(handle, password)
                return client
            except Exception as e:
                logger.error(f"Error initializing Bluesky client: {e}")
        return None

    def _init_mastodon(
        self, tokens: Optional[Dict[str, str]]
    ) -> Dict[str, Optional[Mastodon]]:
        """Initialize Mastodon clients.

        Args:
            tokens: Dict mapping instances to access tokens

        Returns:
            Dict mapping instances to Mastodon clients
        """
        clients = {}
        if tokens:
            for instance, token in tokens.items():
                try:
                    clients[instance] = Mastodon(
                        access_token=token,
                        api_base_url=f"https://{instance}",
                    )
                except Exception as e:
                    logger.error(
                        f"Error initializing Mastodon client for {instance}: {e}"
                    )
                    clients[instance] = None
        return clients

    def _init_matrix(self, token: Optional[str]) -> Optional[MatrixHandler]:
        """Initialize Matrix handler."""
        if token:
            try:
                return MatrixHandler("https://matrix.org", token)
            except Exception as e:
                logger.error(f"Error initializing Matrix handler: {e}")
        return None

    def _compile_patterns(self) -> None:
        """Compile regex patterns for matching chemical names and identifiers."""
        # Common chemical name patterns
        self.name_patterns = [
            # Common chemical suffixes
            r"\b[A-Z][A-Za-z\d\-]+(?:acid|amine|ol|one|ate|ide|ium|yl|ene|ane)\b",
            # Common chemical groups
            r"\b\d*[A-Z][A-Za-z\d\-]*(?:phenyl|methyl|ethyl|propyl|butyl|benzyl)\b",
            # Drug suffixes
            r"\b[A-Z][A-Za-z\d\-]+(?:zine|dine|pine|line|pram|lam|pam|ine)\b",
            # Research chemical patterns
            r"\b\d?-?[A-Z]{1,3}[A-Z\d]{1,4}(?:-[A-Z\d]+)*\b",
            # Common psychoactive prefixes
            r"\b(?:DMT|LSD|DOx|NBOMe|2C|5-MeO|4-AcO|3-MMC|4-MMC|α-PHP)\b",
            # Additional research chemical patterns
            r"\b(?:25[A-Z]-[A-Z]+|1P-|1B-|1V-|1cP-|AL-|ETH-|PRO-|LSZ)\b",
            # Tryptamine patterns
            r"\b(?:4-HO-|4-AcO-|5-MeO-|5-HO-|DiPT|MiPT|DPT|MPT|EPT)\b",
            # Phenethylamine patterns
            r"\b(?:2C-[A-Z]|DOB|DOC|DOI|DOM|DOET|TMA|PMA|PMMA)\b",
            # Amphetamine patterns
            r"\b(?:3-FA|4-FA|2-FMA|3-FMA|4-FMA|2-FEA|3-FEA|4-FEA)\b",
            # Cathinone patterns
            r"\b(?:3-MMC|4-MMC|3-CMC|4-CMC|3-MEC|4-MEC|NEP|NEH)\b",
            # Dissociative patterns
            r"\b(?:2F-DCK|3-MeO-PCE|3-HO-PCE|O-PCE|MXPr|MXiPr)\b",
            # Nootropic patterns
            r"\b(?:Piracetam|Aniracetam|Oxiracetam|Pramiracetam|Phenylpiracetam)\b",
            r"\b(?:Noopept|Semax|Selank|Bromantane|Modafinil|Adrafinil)\b",
            # Receptor patterns
            r"\b(?:5-HT[1-7]|D[1-5]|NMDA|AMPA|GABA[AB]?|mGlu[1-8])\b",
            # Additional patterns
            r"\b(?:[1-9]\d*[A-Z]-[A-Z]+|[A-Z]{2,3}-\d{2,4})\b",
            # Peptide patterns
            r"\b(?:BPC-157|TB-500|GHK-Cu|Thymosin|Epithalon)\b",
            # Additional drug classes
            r"\b(?:benzofuran|cathinone|pyrrolidine|thiophene|oxazoline)\b",
            # Chemical element combinations
            r"\b(?:[A-Z][a-z]?\d*){2,}\b",
        ]

        # SMILES pattern
        self.smiles_pattern = r"(?:\[[A-Za-z\d@+\-]{1,10}\]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"

        # InChI pattern
        self.inchi_pattern = r"InChI=1S?/[A-Za-z0-9/.]+(?:[+-][A-Za-z0-9/.]+)*"

        # CAS number pattern
        self.cas_pattern = r"\b\d{1,7}-\d{2}-\d\b"

        # Compile all patterns
        self.compiled_patterns = {
            "names": [re.compile(p) for p in self.name_patterns],
            "smiles": re.compile(self.smiles_pattern),
            "inchi": re.compile(self.inchi_pattern),
            "cas": re.compile(self.cas_pattern),
        }

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
                logger.error(f"Error reading cache file {cache_file}: {e}")
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
            logger.error(f"Error writing cache file {cache_file}: {e}")

    def search_all_sources(
        self,
        days: int = 30,
        limit: int = 1000,
        use_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """Search all available social media sources.

        Args:
            days: Number of days to search back
            limit: Maximum items to retrieve per source
            use_cache: Whether to use cached results

        Returns:
            List of dictionaries containing found compounds and context
        """
        # Check cache first
        if use_cache:
            cache_key = f"all_sources_{days}_{limit}"
            cached = self._get_cached_result(cache_key, "mentions")
            if cached:
                return cached

        compounds = []

        # Search each source in parallel
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = []

            # Reddit search
            if self.reddit:
                futures.append(
                    executor.submit(
                        self.search_reddit,
                        self.SUBREDDITS,
                        time_filter="month" if days <= 30 else "year",
                        limit=limit,
                        use_cache=use_cache,
                    )
                )

            # Twitter search
            if self.twitter:
                futures.append(
                    executor.submit(
                        self.search_twitter,
                        self.TWITTER_QUERIES,
                        days=days,
                        limit=limit,
                        use_cache=use_cache,
                    )
                )

            # Discord search
            if self.discord:
                futures.append(
                    executor.submit(
                        self.search_discord,
                        self.DISCORD_CHANNELS,
                        days=days,
                        limit=limit,
                        use_cache=use_cache,
                    )
                )

            # Bluesky search
            if self.bluesky:
                futures.append(
                    executor.submit(
                        self.search_bluesky,
                        self.BLUESKY_QUERIES,
                        days=days,
                        limit=limit,
                        use_cache=use_cache,
                    )
                )

            # Mastodon search
            if self.mastodon:
                for instance, client in self.mastodon.items():
                    if client:
                        futures.append(
                            executor.submit(
                                self.search_mastodon,
                                instance,
                                client,
                                days=days,
                                limit=limit,
                                use_cache=use_cache,
                            )
                        )

            # Matrix search
            if self.matrix:
                futures.append(
                    executor.submit(
                        self.search_matrix,
                        self.MATRIX_ROOMS,
                        days=days,
                        limit=limit,
                        use_cache=use_cache,
                    )
                )

            # Collect results with progress bar
            with tqdm(
                total=len(futures),
                desc="Searching social media",
            ) as pbar:
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            compounds.extend(result)
                    except Exception as e:
                        logger.error(f"Error searching source: {e}")
                    pbar.update(1)

        # Cache results
        if use_cache:
            self._cache_result(cache_key, compounds, "mentions")

        return compounds

    def search_reddit(
        self,
        subreddits: List[str],
        time_filter: str = "month",
        limit: int = 1000,
        use_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """Search Reddit for compound mentions.

        Args:
            subreddits: List of subreddits to search
            time_filter: Time filter (hour, day, week, month, year, all)
            limit: Maximum number of posts to retrieve
            use_cache: Whether to use cached results

        Returns:
            List of dictionaries containing found compounds and context
        """
        if not self.reddit:
            return []

        # Check cache first
        if use_cache:
            cache_key = f"reddit_{'-'.join(subreddits)}_{time_filter}_{limit}"
            cached = self._get_cached_result(cache_key, "mentions")
            if cached:
                return cached

        compounds = []

        # Search each subreddit in parallel
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = []
            for subreddit_name in subreddits:
                futures.append(
                    executor.submit(
                        self._search_subreddit,
                        subreddit_name,
                        time_filter,
                        limit,
                    )
                )

            # Collect results with progress bar
            with tqdm(
                total=len(futures),
                desc="Searching subreddits",
            ) as pbar:
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            compounds.extend(result)
                    except Exception as e:
                        logger.error(f"Error searching subreddit: {e}")
                    pbar.update(1)

        # Cache results
        if use_cache:
            self._cache_result(cache_key, compounds, "mentions")

        return compounds

    def _search_subreddit(
        self,
        subreddit_name: str,
        time_filter: str,
        limit: int,
    ) -> List[Dict[str, Any]]:
        """Search a subreddit for compound mentions.

        Args:
            subreddit_name: Name of subreddit to search
            time_filter: Time filter
            limit: Maximum number of posts to retrieve

        Returns:
            List of dictionaries containing found compounds and context
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                subreddit = self.reddit.subreddit(subreddit_name)

                # Search posts
                for post in tqdm(
                    subreddit.top(time_filter=time_filter, limit=limit),
                    desc=f"Processing r/{subreddit_name} posts",
                    leave=False,
                ):
                    compounds.extend(
                        self._extract_compounds_from_text(
                            f"{post.title}\n{post.selftext}",
                            source=f"reddit/r/{subreddit_name}",
                            url=f"https://reddit.com{post.permalink}",
                            date=datetime.fromtimestamp(post.created_utc),
                        )
                    )

                    # Search comments
                    post.comments.replace_more(limit=0)
                    for comment in post.comments.list():
                        compounds.extend(
                            self._extract_compounds_from_text(
                                comment.body,
                                source=f"reddit/r/{subreddit_name}",
                                url=f"https://reddit.com{comment.permalink}",
                                date=datetime.fromtimestamp(comment.created_utc),
                            )
                        )

                break  # Success, exit retry loop

            except Exception as e:
                logger.error(
                    f"Error searching subreddit {subreddit_name} (attempt {retries + 1}/{self.max_retries}): {e}"
                )
                retries += 1
                if retries == self.max_retries:
                    logger.error(
                        f"Failed to search subreddit {subreddit_name} after {self.max_retries} attempts"
                    )

        return compounds

    def search_twitter(
        self,
        queries: List[str],
        days: int = 7,
        limit: int = 1000,
        use_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """Search Twitter for compound mentions.

        Args:
            queries: List of search queries
            days: Number of days to search back
            limit: Maximum number of tweets to retrieve
            use_cache: Whether to use cached results

        Returns:
            List of dictionaries containing found compounds and context
        """
        if not self.twitter:
            return []

        # Check cache first
        if use_cache:
            cache_key = f"twitter_{'-'.join(queries)}_{days}_{limit}"
            cached = self._get_cached_result(cache_key, "mentions")
            if cached:
                return cached

        compounds = []
        since_date = datetime.now() - timedelta(days=days)

        # Search each query in parallel
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = []
            for query in queries:
                futures.append(
                    executor.submit(
                        self._search_twitter_query,
                        query,
                        since_date,
                        limit,
                    )
                )

            # Collect results with progress bar
            with tqdm(
                total=len(futures),
                desc="Searching Twitter",
            ) as pbar:
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            compounds.extend(result)
                    except Exception as e:
                        logger.error(f"Error searching Twitter: {e}")
                    pbar.update(1)

        # Cache results
        if use_cache:
            self._cache_result(cache_key, compounds)

        return compounds

    def _search_twitter_query(
        self,
        query: str,
        since_date: datetime,
        limit: int,
    ) -> List[Dict[str, Any]]:
        """Search Twitter for a query.

        Args:
            query: Search query
            since_date: Only get tweets after this date
            limit: Maximum number of tweets to retrieve

        Returns:
            List of dictionaries containing found compounds and context
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                tweets = tweepy.Cursor(
                    self.twitter.search_tweets,
                    q=query,
                    lang="en",
                    tweet_mode="extended",
                ).items(limit)

                for tweet in tweets:
                    if tweet.created_at > since_date:
                        compounds.extend(
                            self._extract_compounds_from_text(
                                tweet.full_text,
                                source="twitter",
                                url=f"https://twitter.com/i/web/status/{tweet.id}",
                                date=tweet.created_at,
                            )
                        )

                break  # Success, exit retry loop

            except Exception as e:
                logger.error(
                    f"Error searching Twitter query '{query}' (attempt {retries + 1}/{self.max_retries}): {e}"
                )
                retries += 1
                if retries == self.max_retries:
                    logger.error(
                        f"Failed to search Twitter query '{query}' after {self.max_retries} attempts"
                    )

        return compounds

    def search_discord(
        self,
        channels: List[str],
        days: int = 7,
        limit: int = 1000,
        use_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """Search Discord channels for compound mentions.

        Args:
            channels: List of channels to search
            days: Number of days to search back
            limit: Maximum number of messages to retrieve
            use_cache: Whether to use cached results

        Returns:
            List of dictionaries containing found compounds and context
        """
        if not self.discord:
            return []

        # Check cache first
        if use_cache:
            cache_key = f"discord_{'-'.join(channels)}_{days}_{limit}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return cached

        compounds = []
        since_date = datetime.now() - timedelta(days=days)

        # Search each channel in parallel
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = []
            for channel in channels:
                futures.append(
                    executor.submit(
                        self._search_discord_channel,
                        channel,
                        since_date,
                        limit,
                    )
                )

            # Collect results with progress bar
            with tqdm(
                total=len(futures),
                desc="Searching Discord",
            ) as pbar:
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            compounds.extend(result)
                    except Exception as e:
                        logger.error(f"Error searching Discord: {e}")
                    pbar.update(1)

        # Cache results
        if use_cache:
            self._cache_result(cache_key, compounds)

        return compounds

    def _search_discord_channel(
        self,
        channel: str,
        since_date: datetime,
        limit: int,
    ) -> List[Dict[str, Any]]:
        """Search a Discord channel for compound mentions.

        Args:
            channel: Channel to search
            since_date: Only get messages after this date
            limit: Maximum number of messages to retrieve

        Returns:
            List of dictionaries containing found compounds and context
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                messages = self.discord.get_channel_history(
                    channel_name=channel,
                    limit=limit,
                    after=since_date,
                )

                for message in messages:
                    compounds.extend(
                        self._extract_compounds_from_text(
                            message["content"],
                            source=message["source"],
                            url=message["url"],
                            date=message["date"],
                        )
                    )

                break  # Success, exit retry loop

            except Exception as e:
                logger.error(
                    f"Error searching Discord channel {channel} (attempt {retries + 1}/{self.max_retries}): {e}"
                )
                retries += 1
                if retries == self.max_retries:
                    logger.error(
                        f"Failed to search Discord channel {channel} after {self.max_retries} attempts"
                    )

        return compounds

    def search_bluesky(
        self,
        queries: List[str],
        days: int = 7,
        limit: int = 1000,
        use_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """Search Bluesky for compound mentions.

        Args:
            queries: List of search queries
            days: Number of days to search back
            limit: Maximum number of posts to retrieve
            use_cache: Whether to use cached results

        Returns:
            List of dictionaries containing found compounds and context
        """
        if not self.bluesky:
            return []

        # Check cache first
        if use_cache:
            cache_key = f"bluesky_{'-'.join(queries)}_{days}_{limit}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return cached

        compounds = []
        since_date = datetime.now() - timedelta(days=days)

        # Search each query in parallel
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = []
            for query in queries:
                futures.append(
                    executor.submit(
                        self._search_bluesky_query,
                        query,
                        since_date,
                        limit,
                    )
                )

            # Collect results with progress bar
            with tqdm(
                total=len(futures),
                desc="Searching Bluesky",
            ) as pbar:
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            compounds.extend(result)
                    except Exception as e:
                        logger.error(f"Error searching Bluesky: {e}")
                    pbar.update(1)

        # Cache results
        if use_cache:
            self._cache_result(cache_key, compounds)

        return compounds

    def _search_bluesky_query(
        self,
        query: str,
        since_date: datetime,
        limit: int,
    ) -> List[Dict[str, Any]]:
        """Search Bluesky for a query.

        Args:
            query: Search query
            since_date: Only get posts after this date
            limit: Maximum number of posts to retrieve

        Returns:
            List of dictionaries containing found compounds and context
        """
        compounds = []
        retries = 0

        while retries < self.max_retries:
            try:
                posts = self.bluesky.search_posts(query, limit=limit)

                for post in posts:
                    if datetime.fromisoformat(post.created_at) > since_date:
                        compounds.extend(
                            self._extract_compounds_from_text(
                                post.text,
                                source="bluesky",
                                url=f"https://bsky.app/profile/{post.author.handle}/post/{post.rkey}",
                                date=datetime.fromisoformat(post.created_at),
                            )
                        )

                break  # Success, exit retry loop

            except Exception as e:
                logger.error(
                    f"Error searching Bluesky query '{query}' (attempt {retries + 1}/{self.max_retries}): {e}"
                )
                retries += 1
                if retries == self.max_retries:
                    logger.error(
                        f"Failed to search Bluesky query '{query}' after {self.max_retries} attempts"
                    )

        return compounds

    def _extract_compounds_from_text(
        self,
        text: str,
        source: str,
        url: str,
        date: datetime,
    ) -> List[Dict[str, Any]]:
        """Extract compound mentions from text.

        Args:
            text: Text to extract from
            source: Source of the text
            url: URL of the source
            date: Date of the text

        Returns:
            List of dictionaries containing found compounds and context
        """
        compounds = []

        try:
            # Extract using regex patterns
            for pattern_type, patterns in self.compiled_patterns.items():
                if isinstance(patterns, list):
                    for pattern in patterns:
                        for match in pattern.finditer(text):
                            compounds.append(
                                {
                                    "identifier": match.group(),
                                    "type": pattern_type,
                                    "context": self._get_context(text, match.span()),
                                    "source": source,
                                    "url": url,
                                    "date": date,
                                }
                            )
                else:
                    for match in patterns.finditer(text):
                        compounds.append(
                            {
                                "identifier": match.group(),
                                "type": pattern_type,
                                "context": self._get_context(text, match.span()),
                                "source": source,
                                "url": url,
                                "date": date,
                            }
                        )

            # Extract using NER model
            try:
                # Split text into chunks to avoid token limits
                chunks = [text[i : i + 512] for i in range(0, len(text), 512)]
                for chunk in chunks:
                    entities = self.ner_model(chunk)
                    for entity in entities:
                        if entity["entity"].startswith("B-CHEMICAL"):
                            # Validate the entity
                            if self._validate_chemical_entity(entity["word"]):
                                compounds.append(
                                    {
                                        "identifier": entity["word"],
                                        "type": "names",
                                        "context": self._get_context(
                                            chunk,
                                            (entity["start"], entity["end"]),
                                        ),
                                        "source": source,
                                        "url": url,
                                        "date": date,
                                        "confidence": entity.get("score", 0.0),
                                    }
                                )

            except Exception as e:
                logger.error(f"Error running NER model: {e}")

            # Remove duplicates while preserving order
            seen = set()
            unique_compounds = []
            for compound in compounds:
                key = (compound["identifier"], compound["type"])
                if key not in seen:
                    seen.add(key)
                    unique_compounds.append(compound)

            return unique_compounds

        except Exception as e:
            logger.error(f"Error extracting compounds from text: {e}")
            return []

    def _validate_chemical_entity(self, entity: str) -> bool:
        """Validate a chemical entity extracted by NER.

        Args:
            entity: Entity text to validate

        Returns:
            Whether the entity appears to be a valid chemical name
        """
        try:
            # Check length
            if len(entity) < 3:
                return False

            # Check for common chemical patterns
            for pattern in self.name_patterns:
                if re.search(pattern, entity, re.I):
                    return True

            # Try to convert to structure
            try:
                mol = AllChem.MolFromIUPACName(entity)
                if mol:
                    return True
            except (ValueError, RuntimeError) as e:
                logger.debug(f"Error parsing IUPAC name: {e}")

            # Check for numbers and special characters
            if re.search(r"[0-9\-]", entity):
                return True

            # Check for common chemical elements
            if re.search(
                r"\b(?:H|He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al|Si|P|S|Cl|Ar)\b", entity
            ):
                return True

            return False

        except Exception as e:
            logger.error(f"Error validating chemical entity {entity}: {e}")
            return False

    def _get_context(
        self,
        text: str,
        span: Tuple[int, int],
        context_chars: int = 200,
    ) -> str:
        """Get context around a match in text.

        Args:
            text: Full text
            span: Start and end indices of match
            context_chars: Number of characters of context to include

        Returns:
            Text with context around match
        """
        start, end = span
        context_start = max(0, start - context_chars)
        context_end = min(len(text), end + context_chars)
        return text[context_start:context_end].strip()

    def validate_compound(
        self,
        identifier: str,
        type: str,
        use_cache: bool = True,
    ) -> Optional[CompoundData]:
        """Validate and standardize a compound identifier.

        Args:
            identifier: Compound identifier (name, SMILES, CAS)
            type: Type of identifier
            use_cache: Whether to use cached results

        Returns:
            CompoundData object if valid, None otherwise
        """
        # Check cache first
        if use_cache:
            cache_key = f"validate_{identifier}_{type}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return CompoundData(**cached)

        try:
            mol = None
            name = None
            smiles = None
            inchi = None
            cas = None

            # Convert identifier to RDKit mol based on type
            if type == "smiles":
                mol = Chem.MolFromSmiles(identifier)
                if mol:
                    smiles = identifier
            elif type == "inchi":
                mol = Chem.MolFromInchi(identifier)
                if mol:
                    inchi = identifier
            elif type == "cas":
                cas = identifier
                # TODO: Look up structure from CAS
            elif type == "names":
                name = identifier
                # Try to convert name to structure using chemical name parsing
                try:
                    mol = AllChem.MolFromIUPACName(identifier)
                except (ImportError, ValueError, RuntimeError) as e:
                    logger.debug(f"Error parsing chemical name: {e}")

            if not mol and not cas:
                return None

            # Standardize structure if available
            if mol:
                mol = self.standardizer.standardize(mol)
                if not smiles:
                    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                if not inchi:
                    inchi = Chem.MolToInchi(mol)

                # Calculate basic properties
                mw = Descriptors.ExactMolWt(mol)
                logp = Descriptors.MolLogP(mol)
                tpsa = Descriptors.TPSA(mol)
                hbd = Descriptors.NumHDonors(mol)
                hba = Descriptors.NumHAcceptors(mol)
                rotatable = Descriptors.NumRotatableBonds(mol)
            else:
                mw = logp = tpsa = hbd = hba = rotatable = None

            compound = CompoundData(
                name=name,
                smiles=smiles,
                inchi=inchi,
                cas=cas,
                molecular_weight=mw,
                logp=logp,
                tpsa=tpsa,
                hbd=hbd,
                hba=hba,
                rotatable_bonds=rotatable,
            )

            # Cache result
            if use_cache:
                self._cache_result(cache_key, asdict(compound))

            return compound

        except Exception as e:
            print(f"Error validating compound {identifier}: {e}")
            return None

    def _is_relevant_mention(
        self,
        mention: Dict[str, Any],
        compound: CompoundData,
    ) -> bool:
        """Check if a mention is relevant to a compound.

        Args:
            mention: Mention data
            compound: Compound to check against

        Returns:
            Whether mention is relevant
        """
        try:
            # Check structure similarity if available
            if compound.smiles and mention.get("smiles"):
                similarity = calculate_similarity(
                    get_fingerprint(Chem.MolFromSmiles(compound.smiles)),
                    get_fingerprint(Chem.MolFromSmiles(mention["smiles"])),
                )
                if similarity > 0.8:
                    return True

            # Check name/synonym matching
            if compound.name and mention.get("identifier"):
                if compound.name.lower() in mention["identifier"].lower():
                    return True

            # Use ML to check context relevance
            context = mention.get("context", "")
            if context:
                # Classify text relevance
                classification = self.text_classifier(context)[0]
                if (
                    classification["label"] == "relevant"
                    and classification["score"] > 0.8
                ):
                    return True

                # Check for chemical entity mentions
                entities = self.ner_model(context)
                chemical_entities = [
                    e for e in entities if e["entity"].startswith("B-CHEMICAL")
                ]
                if chemical_entities:
                    return True

            return False

        except Exception as e:
            print(f"Error checking mention relevance: {e}")
            return False

    def enrich_compound_data(
        self,
        compound: CompoundData,
        min_mentions: int = 3,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Enrich compound data with social media mentions.

        Args:
            compound: Compound to enrich
            min_mentions: Minimum number of mentions required
            use_cache: Whether to use cached results

        Returns:
            Dictionary containing enriched data
        """
        # Check cache first
        if use_cache:
            cache_key = f"enrich_{compound.smiles or compound.name}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return cached

        # Search for compound mentions
        mentions = self.search_all_sources(days=365, limit=1000)

        # Filter relevant mentions
        relevant_mentions = [
            m for m in mentions if self._is_relevant_mention(m, compound)
        ]

        if len(relevant_mentions) < min_mentions:
            return {}

        # Analyze mentions
        try:
            # Sentiment analysis
            sentiments = []
            for mention in relevant_mentions:
                try:
                    sentiment = self.sentiment_model(mention["context"])[0]
                    sentiments.append(
                        {
                            "label": sentiment["label"],
                            "score": sentiment["score"],
                        }
                    )
                except Exception as e:
                    print(f"Error analyzing sentiment: {e}")
                    continue

            # Effect classification
            effects = []
            for mention in relevant_mentions:
                try:
                    effect = self.effect_classifier(mention["context"])[0]
                    if effect["score"] > 0.8:
                        effects.append(effect["label"])
                except Exception as e:
                    print(f"Error classifying effects: {e}")
                    continue

            # Abuse potential classification
            abuse_scores = []
            for mention in relevant_mentions:
                try:
                    abuse = self.abuse_classifier(mention["context"])[0]
                    abuse_scores.append(abuse["score"])
                except Exception as e:
                    print(f"Error classifying abuse potential: {e}")
                    continue

            # Safety classification
            safety_concerns = []
            for mention in relevant_mentions:
                try:
                    safety = self.safety_classifier(mention["context"])[0]
                    if safety["score"] > 0.8 and safety["label"] == "concern":
                        safety_concerns.append(mention["context"])
                except Exception as e:
                    print(f"Error classifying safety: {e}")
                    continue

            result = {
                "total_mentions": len(relevant_mentions),
                "sources": list(set(m["source"] for m in relevant_mentions)),
                "first_mention": min(m["date"] for m in relevant_mentions),
                "latest_mention": max(m["date"] for m in relevant_mentions),
                "sample_contexts": [m["context"] for m in relevant_mentions[:5]],
                "urls": list(set(m["url"] for m in relevant_mentions)),
                "sentiment_summary": {
                    "positive": len(
                        [s for s in sentiments if s["label"] == "POSITIVE"]
                    ),
                    "negative": len(
                        [s for s in sentiments if s["label"] == "NEGATIVE"]
                    ),
                    "neutral": len([s for s in sentiments if s["label"] == "NEUTRAL"]),
                    "average_score": (
                        np.mean([s["score"] for s in sentiments])
                        if sentiments
                        else None
                    ),
                },
                "reported_effects": list(set(effects)),
                "abuse_potential": (np.mean(abuse_scores) if abuse_scores else None),
                "safety_concerns": safety_concerns,
            }

            # Cache result
            if use_cache:
                self._cache_result(cache_key, result)

            return result

        except Exception as e:
            print(f"Error enriching compound data: {e}")
            return {}
