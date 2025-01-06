"""Social media data harvester.

This module provides functionality to:
1. Monitor Reddit for compound discussions
2. Monitor Twitter for compound mentions
3. Extract and analyze social media data
4. Detect novel compounds

Features:
- Circuit breaker pattern for resilience
- Caching with fallback
- Metrics collection
- Error tracking
- Enhanced text classification
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, TYPE_CHECKING
from datetime import datetime, timedelta
import json
from collections import Counter

import praw
import tweepy
from transformers import pipeline

from ....pipeline.infrastructure.circuit_breaker import CircuitConfig
from .base import BaseWebClient
from ...base.core import Compound

if TYPE_CHECKING:
    from .http import HTTPClient


class SocialClient(BaseWebClient):
    """Client for social media data sources."""

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
        name: str = "social",
        reddit_client_id: str = "",
        reddit_client_secret: str = "",
        twitter_bearer_token: str = "",
        http_client: Optional["HTTPClient"] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        circuit_config: Optional[CircuitConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize social media client.
        
        Args:
            name: Client name for circuit breaker
            reddit_client_id: Reddit API client ID
            reddit_client_secret: Reddit API client secret
            twitter_bearer_token: Twitter API bearer token
            http_client: Optional HTTP client to use
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            circuit_config: Optional circuit breaker configuration
            logger: Optional logger instance
        """
        super().__init__(
            name=name,
            http_client=http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            circuit_config=circuit_config,
            logger=logger,
        )

        # Initialize Reddit client
        if reddit_client_id and reddit_client_secret:
            self.reddit = praw.Reddit(
                client_id=reddit_client_id,
                client_secret=reddit_client_secret,
                user_agent="ChemData/1.0",
            )
        else:
            self.reddit = None

        # Initialize Twitter client
        if twitter_bearer_token:
            self.twitter = tweepy.Client(
                bearer_token=twitter_bearer_token,
                wait_on_rate_limit=True,
            )
        else:
            self.twitter = None

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
                fallback=lambda: self._get_cached_reddit_data(name),
            )
            if reddit_data:
                data["reddit"].update(reddit_data)
                self._update_reddit_stats(reddit_data)
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
                fallback=lambda: self._get_cached_twitter_data(name),
            )
            if twitter_data:
                data["twitter"].update(twitter_data)
                self._update_twitter_stats(twitter_data)
            else:
                self.source_stats["twitter"]["failure"] += 1
        except Exception as e:
            self.logger.error(f"Error getting Twitter data: {str(e)}")
            self.source_stats["twitter"]["failure"] += 1

        return data if data["reddit"]["posts"] or data["twitter"]["tweets"] else None

    def _update_reddit_stats(self, data: Dict[str, Any]) -> None:
        """Update Reddit statistics.
        
        Args:
            data: Reddit data to update stats with
        """
        self.source_stats["reddit"]["success"] += 1
        self.source_stats["reddit"]["posts"] += len(data["posts"])
        self.source_stats["reddit"]["comments"] += len(data["comments"])
        self.source_stats["reddit"]["novel_mentions"] += len(
            data.get("novel_mentions", [])
        )

    def _update_twitter_stats(self, data: Dict[str, Any]) -> None:
        """Update Twitter statistics.
        
        Args:
            data: Twitter data to update stats with
        """
        self.source_stats["twitter"]["success"] += 1
        self.source_stats["twitter"]["tweets"] += len(data["tweets"])
        self.source_stats["twitter"]["users"] += len(data["users"])
        self.source_stats["twitter"]["novel_mentions"] += len(
            data.get("novel_mentions", [])
        )

    def _classify_text(self, text: str) -> Optional[Dict[str, Any]]:
        """Classify text using ML model.
        
        Args:
            text: Text to classify
            
        Returns:
            Classification result or None if no model
        """
        if not self.text_classifier:
            return None

        classification = self.text_classifier(text[:512])[0]
        return {
            "label": classification["label"],
            "score": classification["score"],
        }

    def _extract_compounds(self, text: str) -> List[str]:
        """Extract compound names from text using NER model.
        
        Args:
            text: Text to extract from
            
        Returns:
            List of extracted compound names
        """
        if not self.ner_model:
            return []

        entities = self.ner_model(text)
        return [
            e["word"] for e in entities
            if e["entity"] == "COMPOUND"
        ]

    def _process_reddit_post(
        self,
        post: Any,
        data: Dict[str, Any],
    ) -> None:
        """Process a Reddit post.
        
        Args:
            post: Reddit post to process
            data: Data dictionary to update
        """
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
        if post.selftext:
            classification = self._classify_text(post.selftext)
            if classification:
                post_data["classification"] = classification
                data["topics"][classification["label"]] += 1

            # Extract novel compounds
            compounds = self._extract_compounds(post.selftext)
            if compounds:
                post_data["novel_compounds"] = compounds
                data["novel_mentions"].extend(compounds)

        data["posts"].append(post_data)
        data["subreddits"][post.subreddit.display_name] += 1

    def _process_reddit_comment(
        self,
        comment: Any,
        data: Dict[str, Any],
    ) -> None:
        """Process a Reddit comment.
        
        Args:
            comment: Reddit comment to process
            data: Data dictionary to update
        """
        comment_data = {
            "id": comment.id,
            "text": comment.body,
            "url": f"https://reddit.com{comment.permalink}",
            "author": str(comment.author),
            "score": comment.score,
            "created_utc": comment.created_utc,
        }

        # Classify comment if model available
        classification = self._classify_text(comment.body)
        if classification:
            comment_data["classification"] = classification
            data["sentiment"][classification["label"]] += 1

        data["comments"].append(comment_data)

    def _process_tweet(
        self,
        tweet: Any,
        data: Dict[str, Any],
    ) -> None:
        """Process a tweet.
        
        Args:
            tweet: Tweet to process
            data: Data dictionary to update
        """
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

        # Classify tweet
        classification = self._classify_text(tweet.text)
        if classification:
            tweet_data["classification"] = classification
            data["sentiment"][classification["label"]] += 1

        # Extract novel compounds
        compounds = self._extract_compounds(tweet.text)
        if compounds:
            tweet_data["novel_compounds"] = compounds
            data["novel_mentions"].extend(compounds)

        data["tweets"].append(tweet_data)
        data["users"][tweet.author.username] += 1

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
        if not self.reddit:
            self.logger.warning("Reddit client not initialized")
            return None

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
                    if (
                        datetime.fromtimestamp(post.created_utc)
                        > datetime.now() - timedelta(days=days)
                    ):
                        self._process_reddit_post(post, data)

                        # Get comments
                        post.comments.replace_more(limit=0)
                        for comment in post.comments.list():
                            self._process_reddit_comment(comment, data)

            except Exception as e:
                self.logger.error(
                    f"Error searching subreddit {subreddit_name}: {str(e)}"
                )
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
        if not self.twitter:
            self.logger.warning("Twitter client not initialized")
            return None

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
                    if (
                        tweet.created_at
                        > datetime.now() - timedelta(days=days)
                    ):
                        self._process_tweet(tweet, data)

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
        metrics = super().get_metrics()
        metrics.update({
            "processed_compounds": len(self.processed_compounds),
            "failed_compounds": len(self.failed_compounds),
            "success_rate": (
                len(self.processed_compounds) /
                (len(self.processed_compounds) + len(self.failed_compounds))
                if self.processed_compounds or self.failed_compounds
                else 0
            ),
            "source_stats": self.source_stats,
        })
        return metrics
