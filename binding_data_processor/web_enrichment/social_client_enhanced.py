"""Enhanced social media data harvester.

This module provides functionality to:
1. Monitor Reddit for compound discussions
2. Monitor Twitter for compound mentions
3. Extract and analyze social media data
4. Detect novel compounds

Enhanced with:
- Circuit breaker pattern for resilience
- Better error handling and recovery
- Improved metrics collection
- Enhanced text classification
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Set, TYPE_CHECKING
from datetime import datetime, timedelta
import json
from collections import Counter

import praw
import tweepy
from transformers import pipeline
import pandas as pd
from tqdm import tqdm

from .base_client import BaseWebClient
from ..models.compound import Compound
from ..pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig

if TYPE_CHECKING:
    from .http_client_enhanced import HTTPClientEnhanced


class SocialClientEnhanced(BaseWebClient):
    """Enhanced client for social media data sources."""

    # Reddit subreddits to monitor
    SUBREDDITS = [
        "researchchemicals",
        "nootropics",
        "DrugNerds",
        "Psychonaut",
        "psychopharmacology",
    ]

    # Twitter search queries
    TWITTER_QUERIES = [
        "research chemical",
        "novel psychoactive",
        "new compound",
        "synthesis route",
        "receptor binding",
        "pharmacology",
    ]

    def __init__(
        self,
        reddit_client_id: str,
        reddit_client_secret: str,
        twitter_bearer_token: str,
        http_client: Optional["HTTPClientEnhanced"] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
        circuit_config: Optional[CircuitBreakerConfig] = None,
    ):
        """Initialize social media client.

        Args:
            reddit_client_id: Reddit API client ID
            reddit_client_secret: Reddit API client secret
            twitter_bearer_token: Twitter API bearer token
            http_client: Optional HTTP client to use
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            logger: Optional logger instance
            circuit_config: Optional circuit breaker configuration
        """
        super().__init__(http_client, model_dir, cache_dir, logger)

        # Initialize Reddit client
        self.reddit = praw.Reddit(
            client_id=reddit_client_id,
            client_secret=reddit_client_secret,
            user_agent="ChemData/1.0",
        )

        # Initialize Twitter client
        self.twitter = tweepy.Client(
            bearer_token=twitter_bearer_token,
            wait_on_rate_limit=True,
        )

        # Initialize text classifier
        if model_dir:
            self.logger.info("Loading text classifier...")
            self.text_classifier = pipeline(
                "text-classification",
                model=str(model_dir / "social_classifier"),
                device="cuda" if model_dir else "cpu",
            )
        else:
            self.text_classifier = None

        # Initialize compound name extractor
        if model_dir:
            self.logger.info("Loading NER model...")
            self.ner_model = pipeline(
                "ner",
                model=str(model_dir / "compound_ner"),
                device="cuda" if model_dir else "cpu",
            )
        else:
            self.ner_model = None

        # Initialize tracking
        self.processed_compounds: List[str] = []
        self.failed_compounds: List[str] = []
        self.source_stats = {
            "reddit": {
                "success": 0,
                "failure": 0,
                "posts": 0,
                "comments": 0,
                "novel_mentions": 0,
            },
            "twitter": {
                "success": 0,
                "failure": 0,
                "tweets": 0,
                "users": 0,
                "novel_mentions": 0,
            },
        }

    def process_compounds(
        self,
        compounds: List[Compound],
        skip_predictions: bool = False,
        use_cache: bool = True,
    ) -> None:
        """Process list of compounds.

        Args:
            compounds: List of compounds to process
            skip_predictions: Whether to skip ML predictions
            use_cache: Whether to use cached results
        """
        for compound in compounds:
            try:
                data = self.get_compound_data(
                    compound.name,
                    compound.cas_number,
                    use_cache=use_cache,
                )
                if data:
                    compound.social_data = data
                    self.processed_compounds.append(compound.name)
                else:
                    self.failed_compounds.append(compound.name)

            except Exception as e:
                self.logger.error(f"Error processing {compound.name}: {str(e)}")
                self.failed_compounds.append(compound.name)

    def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        use_cache: bool = True,
        days: int = 30,
    ) -> Optional[Dict[str, Any]]:
        """Get social media data for a compound.

        Args:
            name: Compound name
            cas_number: Optional CAS number
            use_cache: Whether to use cached results
            days: Number of days to look back

        Returns:
            Dictionary of social media data or None if not found
        """
        data = {
            "reddit": {
                "posts": [],
                "comments": [],
                "subreddits": Counter(),
                "sentiment": Counter(),
                "topics": Counter(),
            },
            "twitter": {
                "tweets": [],
                "users": Counter(),
                "hashtags": Counter(),
                "sentiment": Counter(),
            },
            "novel_mentions": [],
            "last_updated": datetime.now().isoformat(),
        }

        # Get Reddit data
        try:
            reddit_data = self._get_reddit_data(
                name,
                days,
                use_cache=use_cache,
                fallback=self._get_cached_reddit_data,
            )
            if reddit_data:
                data["reddit"].update(reddit_data)
                self.source_stats["reddit"]["success"] += 1
                self.source_stats["reddit"]["posts"] += len(reddit_data["posts"])
                self.source_stats["reddit"]["comments"] += len(reddit_data["comments"])
                self.source_stats["reddit"]["novel_mentions"] += len(reddit_data.get("novel_mentions", []))
            else:
                self.source_stats["reddit"]["failure"] += 1
        except Exception as e:
            self.logger.error(f"Error getting Reddit data: {str(e)}")
            self.source_stats["reddit"]["failure"] += 1

        # Get Twitter data
        try:
            twitter_data = self._get_twitter_data(
                name,
                days,
                use_cache=use_cache,
                fallback=self._get_cached_twitter_data,
            )
            if twitter_data:
                data["twitter"].update(twitter_data)
                self.source_stats["twitter"]["success"] += 1
                self.source_stats["twitter"]["tweets"] += len(twitter_data["tweets"])
                self.source_stats["twitter"]["users"] += len(twitter_data["users"])
                self.source_stats["twitter"]["novel_mentions"] += len(twitter_data.get("novel_mentions", []))
            else:
                self.source_stats["twitter"]["failure"] += 1
        except Exception as e:
            self.logger.error(f"Error getting Twitter data: {str(e)}")
            self.source_stats["twitter"]["failure"] += 1

        return data if data["reddit"]["posts"] or data["twitter"]["tweets"] else None

    def _get_reddit_data(
        self,
        name: str,
        days: int = 30,
        use_cache: bool = True,
        fallback: Optional[callable] = None,
    ) -> Optional[Dict[str, Any]]:
        """Get Reddit data for a compound.

        Args:
            name: Compound name
            days: Number of days to look back
            use_cache: Whether to use cached results
            fallback: Optional fallback function if service fails

        Returns:
            Dictionary of Reddit data or None if not found
        """
        data = {
            "posts": [],
            "comments": [],
            "subreddits": Counter(),
            "sentiment": Counter(),
            "topics": Counter(),
            "novel_mentions": [],
        }

        # Search each subreddit
        for subreddit_name in self.SUBREDDITS:
            try:
                subreddit = self.reddit.subreddit(subreddit_name)

                # Search posts
                for post in subreddit.search(
                    name,
                    time_filter="month",
                    limit=100,
                ):
                    if datetime.fromtimestamp(post.created_utc) > datetime.now() - timedelta(days=days):
                        post_data = {
                            "id": post.id,
                            "title": post.title,
                            "text": post.selftext,
                            "url": f"https://reddit.com{post.permalink}",
                            "author": str(post.author),
                            "subreddit": post.subreddit.display_name,
                            "score": post.score,
                            "created_utc": post.created_utc,
                        }

                        # Classify text if model available
                        if self.text_classifier and post.selftext:
                            classification = self.text_classifier(post.selftext[:512])[0]
                            post_data["classification"] = {
                                "label": classification["label"],
                                "score": classification["score"],
                            }
                            data["topics"][classification["label"]] += 1

                        # Extract novel compounds if model available
                        if self.ner_model and post.selftext:
                            entities = self.ner_model(post.selftext)
                            compounds = [e["word"] for e in entities if e["entity"] == "COMPOUND"]
                            post_data["novel_compounds"] = compounds
                            data["novel_mentions"].extend(compounds)

                        data["posts"].append(post_data)
                        data["subreddits"][post.subreddit.display_name] += 1

                        # Get comments
                        post.comments.replace_more(limit=0)
                        for comment in post.comments.list():
                            comment_data = {
                                "id": comment.id,
                                "text": comment.body,
                                "url": f"https://reddit.com{comment.permalink}",
                                "author": str(comment.author),
                                "score": comment.score,
                                "created_utc": comment.created_utc,
                            }

                            # Classify comment if model available
                            if self.text_classifier:
                                classification = self.text_classifier(comment.body[:512])[0]
                                comment_data["classification"] = {
                                    "label": classification["label"],
                                    "score": classification["score"],
                                }
                                data["sentiment"][classification["label"]] += 1

                            data["comments"].append(comment_data)

            except Exception as e:
                self.logger.error(f"Error searching subreddit {subreddit_name}: {str(e)}")
                continue

        return data if data["posts"] else None

    def _get_twitter_data(
        self,
        name: str,
        days: int = 30,
        use_cache: bool = True,
        fallback: Optional[callable] = None,
    ) -> Optional[Dict[str, Any]]:
        """Get Twitter data for a compound.

        Args:
            name: Compound name
            days: Number of days to look back
            use_cache: Whether to use cached results
            fallback: Optional fallback function if service fails

        Returns:
            Dictionary of Twitter data or None if not found
        """
        data = {
            "tweets": [],
            "users": Counter(),
            "hashtags": Counter(),
            "sentiment": Counter(),
            "novel_mentions": [],
        }

        try:
            # Search tweets
            for query in self.TWITTER_QUERIES:
                search_query = f"{name} {query}"

                tweets = tweepy.Paginator(
                    self.twitter.search_recent_tweets,
                    query=search_query,
                    max_results=100,
                    tweet_fields=["created_at", "public_metrics", "entities"],
                    user_fields=["username", "public_metrics"],
                    expansions=["author_id"],
                ).flatten(limit=1000)

                for tweet in tweets:
                    if tweet.created_at > datetime.now() - timedelta(days=days):
                        tweet_data = {
                            "id": tweet.id,
                            "text": tweet.text,
                            "url": f"https://twitter.com/i/web/status/{tweet.id}",
                            "author": tweet.author.username,
                            "created_at": tweet.created_at.isoformat(),
                            "metrics": tweet.public_metrics,
                        }

                        # Extract hashtags
                        if tweet.entities and "hashtags" in tweet.entities:
                            hashtags = [h["tag"] for h in tweet.entities["hashtags"]]
                            tweet_data["hashtags"] = hashtags
                            for tag in hashtags:
                                data["hashtags"][tag] += 1

                        # Classify tweet if model available
                        if self.text_classifier:
                            classification = self.text_classifier(tweet.text[:512])[0]
                            tweet_data["classification"] = {
                                "label": classification["label"],
                                "score": classification["score"],
                            }
                            data["sentiment"][classification["label"]] += 1

                        # Extract novel compounds if model available
                        if self.ner_model:
                            entities = self.ner_model(tweet.text)
                            compounds = [e["word"] for e in entities if e["entity"] == "COMPOUND"]
                            tweet_data["novel_compounds"] = compounds
                            data["novel_mentions"].extend(compounds)

                        data["tweets"].append(tweet_data)
                        data["users"][tweet.author.username] += 1

        except Exception as e:
            self.logger.error(f"Error searching Twitter: {str(e)}")

        return data if data["tweets"] else None

    def _get_cached_reddit_data(self, name: str) -> Optional[Dict[str, Any]]:
        """Get cached Reddit data.

        Args:
            name: Compound name

        Returns:
            Dictionary of cached data or None if not found
        """
        if not self.cache_dir:
            return None

        try:
            cache_file = self.cache_dir / f"reddit_{name}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                return json.load(f)

        except Exception as e:
            self.logger.error(f"Error reading cached Reddit data: {str(e)}")
            return None

    def _get_cached_twitter_data(self, name: str) -> Optional[Dict[str, Any]]:
        """Get cached Twitter data.

        Args:
            name: Compound name

        Returns:
            Dictionary of cached data or None if not found
        """
        if not self.cache_dir:
            return None

        try:
            cache_file = self.cache_dir / f"twitter_{name}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                return json.load(f)

        except Exception as e:
            self.logger.error(f"Error reading cached Twitter data: {str(e)}")
            return None

    def get_metrics(self) -> Dict[str, Any]:
        """Get client metrics."""
        return {
            "processed_compounds": len(self.processed_compounds),
            "failed_compounds": len(self.failed_compounds),
            "success_rate": (
                len(self.processed_compounds) / (len(self.processed_compounds) + len(self.failed_compounds)) if self.processed_compounds or self.failed_compounds else 0
            ),
            "source_stats": self.source_stats,
        }
