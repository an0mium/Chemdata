"""Social media monitoring for compounds.

This module provides functionality for:
1. Reddit data collection
2. Twitter monitoring
3. Bluesky integration
4. Discord monitoring
5. LLM-based compound extraction
"""

import logging
from typing import List, Optional
from dataclasses import dataclass
import tweepy
import praw
from atproto import Client as AtprotoClient

from ...models.compound.enhanced import EnhancedCompound
from ...web_enrichment.llm_utils import LLMProcessor
from .config import SocialMedia


@dataclass
class SocialConfig:
    """Social media configuration."""
    
    # Reddit settings
    reddit_client_id: Optional[str] = None
    reddit_client_secret: Optional[str] = None
    reddit_user_agent: str = "ChemDataProcessor/1.0"
    
    # Twitter settings
    twitter_api_key: Optional[str] = None
    twitter_api_secret: Optional[str] = None
    
    # Bluesky settings
    bluesky_handle: Optional[str] = None
    bluesky_password: Optional[str] = None
    
    # Discord settings
    discord_token: Optional[str] = None
    
    # Monitoring settings
    post_limit: int = 100
    comment_limit: int = 1000
    batch_size: int = 50


class SocialMonitor:
    """Monitor social media for compounds."""

    def __init__(
        self,
        config: Optional[SocialConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize social monitor.
        
        Args:
            config: Optional social media configuration
            logger: Optional logger instance
        """
        self.config = config or SocialConfig()
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize clients
        self._init_clients()
        
        # Initialize LLM processor
        self.llm = LLMProcessor()

    def _init_clients(self) -> None:
        """Initialize social media clients."""
        try:
            # Initialize Reddit client
            if self.config.reddit_client_id:
                self.reddit = praw.Reddit(
                    client_id=self.config.reddit_client_id,
                    client_secret=self.config.reddit_client_secret,
                    user_agent=self.config.reddit_user_agent,
                )
            else:
                self.reddit = None
            
            # Initialize Twitter client
            if self.config.twitter_api_key:
                auth = tweepy.OAuthHandler(
                    self.config.twitter_api_key,
                    self.config.twitter_api_secret,
                )
                self.twitter = tweepy.API(auth)
            else:
                self.twitter = None
            
            # Initialize Bluesky client
            if self.config.bluesky_handle:
                self.bluesky = AtprotoClient()
                self.bluesky.login(
                    self.config.bluesky_handle,
                    self.config.bluesky_password,
                )
            else:
                self.bluesky = None
            
        except Exception as e:
            self.logger.error(f"Error initializing social clients: {str(e)}")
            raise

    def monitor_all(self) -> List[EnhancedCompound]:
        """Monitor all social media sources.
        
        Returns:
            List of compounds found from social media
        """
        compounds = []
        
        try:
            # Monitor Reddit
            if self.reddit:
                reddit_compounds = self._monitor_reddit()
                compounds.extend(reddit_compounds)
            
            # Monitor Twitter
            if self.twitter:
                twitter_compounds = self._monitor_twitter()
                compounds.extend(twitter_compounds)
            
            # Monitor Bluesky
            if self.bluesky:
                bluesky_compounds = self._monitor_bluesky()
                compounds.extend(bluesky_compounds)
            
            return compounds
            
        except Exception as e:
            self.logger.error(f"Error monitoring social media: {str(e)}")
            return []

    def _monitor_reddit(self) -> List[EnhancedCompound]:
        """Monitor Reddit for compounds.
        
        Returns:
            List of compounds found from Reddit
        """
        compounds = []
        
        try:
            for subreddit_name in SocialMedia.SUBREDDITS:
                try:
                    self.logger.info(f"Monitoring r/{subreddit_name}")
                    
                    # Get subreddit
                    subreddit = self.reddit.subreddit(subreddit_name)
                    texts = []
                    
                    # Get posts
                    for post in subreddit.new(limit=self.config.post_limit):
                        texts.append(post.title)
                        texts.append(post.selftext)
                        
                        # Get comments
                        post.comments.replace_more(limit=0)
                        for comment in post.comments.list():
                            texts.append(comment.body)
                            
                            if len(texts) >= self.config.comment_limit:
                                break
                    
                    # Process texts in batches
                    for i in range(0, len(texts), self.config.batch_size):
                        batch = texts[i:i + self.config.batch_size]
                        
                        # Extract compounds using LLM
                        batch_compounds = self.llm.extract_compounds(
                            batch,
                            source=f"reddit/r/{subreddit_name}",
                        )
                        compounds.extend(batch_compounds)
                    
                except Exception as e:
                    self.logger.error(f"Error monitoring r/{subreddit_name}: {str(e)}")
                    continue
            
            return compounds
            
        except Exception as e:
            self.logger.error(f"Error monitoring Reddit: {str(e)}")
            return []

    def _monitor_twitter(self) -> List[EnhancedCompound]:
        """Monitor Twitter for compounds.
        
        Returns:
            List of compounds found from Twitter
        """
        compounds = []
        
        try:
            for query in SocialMedia.TWITTER_QUERIES:
                try:
                    self.logger.info(f"Searching Twitter for: {query}")
                    
                    # Search tweets
                    tweets = self.twitter.search_tweets(
                        q=query,
                        lang="en",
                        count=self.config.post_limit,
                    )
                    
                    # Extract text
                    texts = [tweet.text for tweet in tweets]
                    
                    # Process in batches
                    for i in range(0, len(texts), self.config.batch_size):
                        batch = texts[i:i + self.config.batch_size]
                        
                        # Extract compounds using LLM
                        batch_compounds = self.llm.extract_compounds(
                            batch,
                            source="twitter",
                        )
                        compounds.extend(batch_compounds)
                    
                except Exception as e:
                    self.logger.error(f"Error searching Twitter: {str(e)}")
                    continue
            
            return compounds
            
        except Exception as e:
            self.logger.error(f"Error monitoring Twitter: {str(e)}")
            return []

    def _monitor_bluesky(self) -> List[EnhancedCompound]:
        """Monitor Bluesky for compounds.
        
        Returns:
            List of compounds found from Bluesky
        """
        compounds = []
        
        try:
            for query in SocialMedia.BLUESKY_QUERIES:
                try:
                    self.logger.info(f"Searching Bluesky for: {query}")
                    
                    # Search posts
                    posts = self.bluesky.search_posts(
                        query,
                        limit=self.config.post_limit,
                    )
                    
                    # Extract text
                    texts = [post.text for post in posts]
                    
                    # Process in batches
                    for i in range(0, len(texts), self.config.batch_size):
                        batch = texts[i:i + self.config.batch_size]
                        
                        # Extract compounds using LLM
                        batch_compounds = self.llm.extract_compounds(
                            batch,
                            source="bluesky",
                        )
                        compounds.extend(batch_compounds)
                    
                except Exception as e:
                    self.logger.error(f"Error searching Bluesky: {str(e)}")
                    continue
            
            return compounds
            
        except Exception as e:
            self.logger.error(f"Error monitoring Bluesky: {str(e)}")
            return []
