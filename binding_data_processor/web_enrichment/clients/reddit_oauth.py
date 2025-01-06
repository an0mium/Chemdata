"""Enhanced Reddit client with OAuth and content monitoring.

This module provides a comprehensive Reddit client that:
1. Uses OAuth for API access with proper token management
2. Falls back to web scraping when needed
3. Monitors subreddits for compound-related content
4. Integrates with storage for persistence
5. Performs safety analysis and alerting
"""

import asyncio
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any
from dataclasses import dataclass

import aiohttp
import praw
from crawl4ai import AsyncWebCrawler, Config

from ..base_client import BaseWebClient
from ..validation.schema import BaseSchema
from ..storage.reddit_storage import RedditStorage
from ...pipeline.infrastructure.rate_limiter import RateLimit
from ...pipeline.infrastructure.circuit_breaker import CircuitBreaker, ErrorCode, McpError


# API Rate Limits
API_RATE_LIMIT = RateLimit(requests=60, period=60)  # 60 requests per minute
OAUTH_RATE_LIMIT = RateLimit(requests=1, period=1)  # 1 token request per second
SCRAPE_RATE_LIMIT = RateLimit(requests=10, period=60)  # 10 scrapes per minute

# Circuit breaker config
CIRCUIT_FAILURE_THRESHOLD = 5
CIRCUIT_RESET_TIMEOUT = 300  # 5 minutes
CIRCUIT_HALF_OPEN_TIMEOUT = 60  # 1 minute

# Confidence scoring weights
FIELD_PRESENCE_WEIGHT = 0.4
QUALITY_FACTORS_WEIGHT = 0.3
EXTRACTED_DATA_WEIGHT = 0.3


@dataclass
class RedditAPIConfig(BaseSchema):
    """Reddit API configuration."""

    client_id: str
    client_secret: str
    user_agent: str
    username: Optional[str] = None
    password: Optional[str] = None
    refresh_token: Optional[str] = None


@dataclass
class RedditPostData(BaseSchema):
    """Schema for Reddit post data."""

    title: str
    content: Optional[str]
    author: Optional[str]
    subreddit: str
    score: Optional[int]
    num_comments: Optional[int]
    compounds: List[str]
    effects: List[str]
    mechanisms: List[str]
    safety_notes: List[str]
    confidence: float
    metadata: Dict[str, Any]


class RedditOAuthClient(BaseWebClient):
    """Enhanced Reddit client with OAuth and monitoring."""

    OAUTH_URL = "https://www.reddit.com/api/v1/authorize"
    TOKEN_URL = "https://www.reddit.com/api/v1/access_token"
    API_URL = "https://oauth.reddit.com"
    BASE_URL = "https://old.reddit.com"  # For scraping fallback

    SCOPES = [
        "read",
        "history",
        "search",
    ]

    SUBREDDITS = [
        "researchchemicals",
        "nootropics",
        "DrugNerds",
        "Psychonaut",
        "psychopharmacology",
    ]

    def __init__(
        self,
        config: RedditAPIConfig,
        llm_provider: str = "ollama/llama2",
        api_token: Optional[str] = None,
        storage_dir: str = "data/reddit",
        **kwargs,
    ):
        """Initialize client.

        Args:
            config: Reddit API configuration
            llm_provider: LLM provider for text analysis
            api_token: Optional API token for LLM
            storage_dir: Directory for storing data
            **kwargs: Additional arguments passed to BaseWebClient
        """
        super().__init__(**kwargs)
        self.config = config
        self.llm_provider = llm_provider
        self.llm_token = api_token
        self.logger = logging.getLogger(__name__)

        # Initialize PRAW client
        self.reddit = praw.Reddit(
            client_id=config.client_id,
            client_secret=config.client_secret,
            user_agent=config.user_agent,
            username=config.username,
            password=config.password,
            refresh_token=config.refresh_token,
        )

        # Initialize storage
        self.storage = RedditStorage(storage_dir=storage_dir)

        # Initialize session and tokens
        self.session = None
        self.access_token = None
        self.token_expiry = None

        # Initialize rate limiters
        self.api_limiter = RateLimit(
            requests=API_RATE_LIMIT.requests,
            period=API_RATE_LIMIT.period,
        )
        self.oauth_limiter = RateLimit(
            requests=OAUTH_RATE_LIMIT.requests,
            period=OAUTH_RATE_LIMIT.period,
        )
        self.scrape_limiter = RateLimit(
            requests=SCRAPE_RATE_LIMIT.requests,
            period=SCRAPE_RATE_LIMIT.period,
        )

        # Initialize circuit breakers
        self.api_circuit = CircuitBreaker(
            "reddit_api",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )
        self.oauth_circuit = CircuitBreaker(
            "reddit_oauth",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )

        # Configure Crawl4AI for fallback scraping
        self.scraper_config = Config(
            javascript=Config.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    ".thing",
                    ".title",
                    ".entry",
                ],
                stealth_mode=True,
            ),
            screenshot=Config.Screenshot(enabled=True, full_page=True),
            extraction=Config.Extraction(
                llm=Config.LLM(
                    provider=llm_provider,
                    api_token=api_token,
                    prompts={
                        "compounds": "Extract mentioned chemical compounds:",
                        "effects": "Extract described effects and outcomes:",
                        "mechanisms": "Extract mechanisms of action:",
                        "safety": "Extract safety information and warnings:",
                    },
                ),
                css={
                    "title": ".title a.title",
                    "content": ".entry .usertext-body",
                    "author": ".entry .author",
                    "subreddit": ".entry .subreddit",
                    "score": ".score.unvoted",
                    "num_comments": ".comments",
                },
            ),
            proxy=Config.Proxy(
                enabled=True,
                rotation=True,
                retry_count=3,
            ),
            rate_limit=Config.RateLimit(
                requests_per_minute=SCRAPE_RATE_LIMIT.requests,
                delay_after_failure=60,
            ),
        )

    async def __aenter__(self):
        """Enter async context."""
        self.session = aiohttp.ClientSession()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Exit async context."""
        if self.session:
            await self.session.close()

    async def ensure_token(self):
        """Ensure valid access token exists."""
        if not self.access_token or (self.token_expiry and datetime.now() >= self.token_expiry):
            await self.refresh_token()

    async def refresh_token(self):
        """Refresh access token."""
        # Wait for rate limit
        await self.oauth_limiter.acquire()

        # Check circuit breaker
        if not self.oauth_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "Token refresh circuit breaker open",
            )

        try:
            if not self.session:
                self.session = aiohttp.ClientSession()

            auth = aiohttp.BasicAuth(self.config.client_id, self.config.client_secret)

            data = {
                "grant_type": "refresh_token" if self.config.refresh_token else "password",
                "duration": "permanent",
            }

            if self.config.refresh_token:
                data["refresh_token"] = self.config.refresh_token
            else:
                data["username"] = self.config.username
                data["password"] = self.config.password

            async with self.session.post(
                self.TOKEN_URL,
                auth=auth,
                data=data,
            ) as response:
                if response.status != 200:
                    error_text = await response.text()
                    self.oauth_circuit.record_failure()
                    raise McpError(
                        ErrorCode.AuthenticationError,
                        f"Token refresh failed: {error_text}",
                    )

                data = await response.json()
                self.access_token = data["access_token"]
                self.token_expiry = datetime.now() + timedelta(seconds=data["expires_in"])

                if "refresh_token" in data:
                    self.config.refresh_token = data["refresh_token"]

                # Record success
                self.oauth_circuit.record_success()

        except Exception as e:
            self.oauth_circuit.record_failure()
            raise McpError(
                ErrorCode.AuthenticationError,
                f"Token refresh failed: {str(e)}",
            )

    async def search_posts(
        self,
        query: str,
        subreddits: Optional[List[str]] = None,
        max_results: int = 100,
        use_api: bool = True,
    ) -> List[RedditPostData]:
        """Search Reddit posts using API with scraping fallback.

        Args:
            query: Search query
            subreddits: Optional list of subreddits (defaults to self.SUBREDDITS)
            max_results: Maximum number of results
            use_api: Whether to try API first (falls back to scraping on failure)

        Returns:
            List of post data
        """
        results = []
        subreddits = subreddits or self.SUBREDDITS

        if use_api:
            try:
                # Try API first
                results = await self._search_posts_api(query, subreddits, max_results)
                if results:
                    return results
            except Exception as e:
                self.logger.warning(f"API search failed, falling back to scraping: {e}")

        # Fall back to scraping
        return await self._search_posts_scrape(query, subreddits, max_results)

    async def _search_posts_api(
        self,
        query: str,
        subreddits: List[str],
        max_results: int,
    ) -> List[RedditPostData]:
        """Search posts using Reddit API.

        Args:
            query: Search query
            subreddits: List of subreddits to search
            max_results: Maximum number of results

        Returns:
            List of post data
        """
        # Check circuit breaker
        if not self.api_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "API circuit breaker open",
            )

        try:
            await self.ensure_token()
            results = []

            headers = {
                "Authorization": f"Bearer {self.access_token}",
                "User-Agent": self.config.user_agent,
            }

            for subreddit in subreddits:
                if len(results) >= max_results:
                    break

                posts = await self._search_subreddit(
                    subreddit=subreddit,
                    query=query,
                    max_results=max_results,
                    headers=headers,
                )

                for post in posts:
                    result = await self._process_post(post["data"])
                    if result:
                        results.append(result)
                        if len(results) >= max_results:
                            break

            # Record success if we got any results
            if results:
                self.api_circuit.record_success()

            return results

        except Exception as e:
            self.api_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"API search failed: {str(e)}",
            )

    async def _search_subreddit(
        self,
        subreddit: str,
        query: str,
        max_results: int,
        headers: Dict[str, str],
    ) -> List[Dict[str, Any]]:
        """Search a single subreddit.

        Args:
            subreddit: Subreddit name
            query: Search query
            max_results: Maximum results to return
            headers: Request headers

        Returns:
            List of post data
        """
        # Wait for rate limit
        await self.api_limiter.acquire()

        url = f"{self.API_URL}/r/{subreddit}/search"
        params = {
            "q": query,
            "restrict_sr": "true",
            "limit": min(100, max_results),
            "sort": "relevance",
        }

        try:
            async with self.session.get(
                url,
                headers=headers,
                params=params,
            ) as response:
                if response.status != 200:
                    error_text = await response.text()
                    error_msg = f"API search failed for r/{subreddit}: {error_text}"
                    self.logger.warning(error_msg)
                    self.api_circuit.record_failure()
                    return []

                data = await response.json()
                return data["data"]["children"]

        except Exception as e:
            self.logger.error(f"Error searching r/{subreddit}: {str(e)}")
            return []

    async def _search_posts_scrape(
        self,
        query: str,
        subreddits: List[str],
        max_results: int,
    ) -> List[RedditPostData]:
        """Search posts using web scraping.

        Args:
            query: Search query
            subreddits: List of subreddits to search
            max_results: Maximum number of results

        Returns:
            List of post data
        """
        results = []

        try:
            async with AsyncWebCrawler() as crawler:
                for subreddit in subreddits:
                    if len(results) >= max_results:
                        break

                    # Wait for rate limit
                    await self.scrape_limiter.acquire()

                    # Construct search URL
                    url = f"{self.BASE_URL}/r/{subreddit}/search" "?q={query}&restrict_sr=on"

                    # Scrape with Crawl4AI
                    result = await crawler.arun(urls=[url], config=self.scraper_config)

                    if not result.success:
                        self.logger.error(f"Failed to scrape r/{subreddit}")
                        continue

                    # Extract data using both CSS and LLM
                    content = result.extracted_content[0]
                    posts = self._parse_scrape_results(content, subreddit)

                    if not posts:
                        continue

                    for post in posts:
                        # Create post data object
                        post_data = RedditPostData(
                            title=post.get("title", ""),
                            content=post.get("content"),
                            author=post.get("author"),
                            subreddit=post.get("subreddit", subreddit),
                            score=post.get("score"),
                            num_comments=post.get("num_comments"),
                            compounds=post.get("compounds", []),
                            effects=post.get("effects", []),
                            mechanisms=post.get("mechanisms", []),
                            safety_notes=post.get("safety", []),
                            confidence=self._calculate_confidence(post),
                            metadata={
                                "source": "reddit_scrape",
                                "scrape_date": datetime.now().isoformat(),
                                "url": post.get("url"),
                                "screenshot": str(result.screenshot),
                                "javascript_logs": result.javascript_logs,
                            },
                        )

                        results.append(post_data)
                        if len(results) >= max_results:
                            break

        except Exception as e:
            self.logger.error(f"Error scraping Reddit: {str(e)}")

        return results

    async def monitor_subreddits(
        self,
        subreddits: Optional[List[str]] = None,
        interval: int = 3600,
        max_posts: int = 100,
    ):
        """Monitor subreddits for new content.

        Args:
            subreddits: List of subreddits to monitor
            interval: Monitoring interval in seconds
            max_posts: Maximum posts to process per interval
        """
        subreddits = subreddits or self.SUBREDDITS
        last_check = {sub: datetime.now() for sub in subreddits}

        while True:
            try:
                for subreddit in subreddits:
                    # Get new posts since last check
                    posts = await self._get_new_posts(
                        subreddit,
                        last_check[subreddit],
                        max_posts,
                    )

                    if posts:
                        # Process and store new posts
                        await self._process_new_posts(posts)

                    last_check[subreddit] = datetime.now()

            except Exception as e:
                self.logger.error(f"Error monitoring subreddits: {str(e)}")

            await asyncio.sleep(interval)

    async def _get_new_posts(
        self,
        subreddit: str,
        since: datetime,
        limit: int,
    ) -> List[Dict[str, Any]]:
        """Get new posts from subreddit.

        Args:
            subreddit: Subreddit name
            since: Get posts newer than this time
            limit: Maximum posts to return

        Returns:
            List of new posts
        """
        # Check circuit breaker
        if not self.api_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "API circuit breaker open",
            )

        try:
            # Wait for rate limit
            await self.api_limiter.acquire()

            await self.ensure_token()
            posts = []

            headers = {
                "Authorization": f"Bearer {self.access_token}",
                "User-Agent": self.config.user_agent,
            }

            url = f"{self.API_URL}/r/{subreddit}/new"
            params = {
                "limit": limit,
            }

            try:
                async with self.session.get(
                    url,
                    headers=headers,
                    params=params,
                ) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        error_msg = f"Failed to get new posts from r/{subreddit}: " f"{error_text}"
                        self.logger.warning(error_msg)
                        self.api_circuit.record_failure()
                        return []

                    data = await response.json()
                    for post in data["data"]["children"]:
                        post_data = post["data"]
                        created = datetime.fromtimestamp(post_data["created_utc"])

                        if created > since:
                            posts.append(post_data)
                        else:
                            break

                # Record success if we got any posts
                if posts:
                    self.api_circuit.record_success()

                return posts

            except Exception as e:
                self.logger.error(f"Error getting posts from r/{subreddit}: {str(e)}")
                self.api_circuit.record_failure()
                return []

        except Exception as e:
            self.api_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"Failed to get new posts: {str(e)}",
            )

    async def _process_new_posts(self, posts: List[Dict[str, Any]]):
        """Process new posts.

        Args:
            posts: List of new posts to process
        """
        for post in posts:
            try:
                # Extract text
                text = f"{post['title']}\n\n{post.get('selftext', '')}"

                # Analyze content
                analysis = await self._analyze_content(text)

                # Check for compounds and safety issues
                if analysis["compounds"] or analysis["safety"]:
                    await self._handle_relevant_post(post, analysis)

            except Exception as e:
                self.logger.error(f"Error processing post: {str(e)}")

    async def _handle_relevant_post(
        self,
        post: Dict[str, Any],
        analysis: Dict[str, List[str]],
    ):
        """Handle post containing relevant content.

        Args:
            post: Post data
            analysis: Content analysis results
        """
        # Store post data
        post_id = post["id"]
        self.storage.store_post(
            post_id,
            {
                "title": post["title"],
                "subreddit": post["subreddit"],
                "url": f"https://reddit.com{post['permalink']}",
                "author": post.get("author"),
                "score": post.get("score"),
                "num_comments": post.get("num_comments"),
                "created_utc": post.get("created_utc"),
                "compounds": analysis["compounds"],
                "effects": analysis.get("effects", []),
                "mechanisms": analysis.get("mechanisms", []),
                "safety": analysis.get("safety", []),
            },
        )

        # Create alert if safety concerns exist
        if analysis["safety"]:
            severity = "high" if len(analysis["safety"]) > 2 else "medium"
            alert = self.storage.create_alert(
                post_id=post_id,
                compounds=analysis["compounds"],
                safety_notes=analysis["safety"],
                severity=severity,
                metadata={
                    "effects": analysis.get("effects", []),
                    "mechanisms": analysis.get("mechanisms", []),
                    "score": post.get("score"),
                    "num_comments": post.get("num_comments"),
                },
            )

            # Log the alert
            self.logger.warning(
                f"Created {severity} severity alert for post in r/{post['subreddit']}:"
                f"\nTitle: {post['title']}"
                f"\nURL: https://reddit.com{post['permalink']}"
                f"\nCompounds: {alert.compounds}"
                f"\nSafety Notes: {alert.safety_notes}"
                f"\nBBB Predictions: {alert.bbb_predictions}"
            )
        else:
            # Log informational finding
            self.logger.info(
                f"Found relevant post in r/{post['subreddit']}:"
                f"\nTitle: {post['title']}"
                f"\nURL: https://reddit.com{post['permalink']}"
                f"\nCompounds: {analysis['compounds']}"
            )

    async def _process_post(self, post_data: Dict[str, Any]) -> Optional[RedditPostData]:
        """Process a single Reddit post.

        Args:
            post_data: Raw post data from API

        Returns:
            Processed post data
        """
        try:
            # Extract text for analysis
            text = f"{post_data['title']}\n\n{post_data.get('selftext', '')}"

            # Analyze content with LLM
            analysis = await self._analyze_content(text)

            # Create post data object
            return RedditPostData(
                title=post_data["title"],
                content=post_data.get("selftext"),
                author=post_data.get("author"),
                subreddit=post_data["subreddit"],
                score=post_data.get("score"),
                num_comments=post_data.get("num_comments"),
                compounds=analysis.get("compounds", []),
                effects=analysis.get("effects", []),
                mechanisms=analysis.get("mechanisms", []),
                safety_notes=analysis.get("safety", []),
                confidence=self._calculate_confidence(post_data),
                metadata={
                    "source": "reddit_api",
                    "extraction_date": datetime.now().isoformat(),
                    "url": f"https://reddit.com{post_data['permalink']}",
                    "created_utc": post_data.get("created_utc"),
                    "upvote_ratio": post_data.get("upvote_ratio"),
                    "is_self": post_data.get("is_self"),
                    "domain": post_data.get("domain"),
                },
            )

        except Exception as e:
            self.logger.error(f"Error processing post: {str(e)}")
            return None

    async def _analyze_content(self, text: str) -> Dict[str, List[str]]:
        """Analyze post content using LLM.

        Args:
            text: Text content to analyze

        Returns:
            Dictionary with analysis results
        """
        try:
            from ..llm_utils import analyze_text

            return await analyze_text(
                text,
                provider=self.llm_provider,
                api_token=self.llm_token,
            )
        except Exception as e:
            self.logger.error(f"Error analyzing content: {str(e)}")
            return {
                "compounds": [],
                "effects": [],
                "mechanisms": [],
                "safety": [],
            }

    def _parse_scrape_results(
        self,
        content: Dict[str, Any],
        subreddit: str,
    ) -> List[Dict[str, Any]]:
        """Parse search results from scraped content.

        Args:
            content: Extracted content from Crawl4AI
            subreddit: Current subreddit being scraped

        Returns:
            List of parsed post data
        """
        posts = []

        # Get post elements
        titles = content.get("title", [])
        contents = content.get("content", [])
        authors = content.get("author", [])
        subreddits = content.get("subreddit", [])
        scores = content.get("score", [])
        num_comments = content.get("num_comments", [])

        # Get LLM-extracted data
        compounds = content.get("compounds", [])
        effects = content.get("effects", [])
        mechanisms = content.get("mechanisms", [])
        safety = content.get("safety", [])

        # Combine data for each post
        for i in range(len(titles)):
            post = {
                "title": titles[i] if i < len(titles) else None,
                "content": contents[i] if i < len(contents) else None,
                "author": authors[i] if i < len(authors) else None,
                "subreddit": subreddits[i] if i < len(subreddits) else subreddit,
                "score": self._parse_score(scores[i]) if i < len(scores) else None,
                "num_comments": self._parse_num_comments(num_comments[i] if i < len(num_comments) else None),
                "url": self._extract_url(titles[i]) if i < len(titles) else None,
                "compounds": compounds,
                "effects": effects,
                "mechanisms": mechanisms,
                "safety": safety,
            }
            posts.append(post)

        return posts

    def _parse_score(self, score_text: Optional[str]) -> Optional[int]:
        """Parse score from Reddit format.

        Args:
            score_text: Score text to parse

        Returns:
            Score as integer or None
        """
        if not score_text:
            return None

        try:
            # Extract number from score text
            import re

            match = re.search(r"\b(\d+)\b", score_text)
            if match:
                return int(match.group(1))
        except Exception:
            pass

        return None

    def _parse_num_comments(self, comments_text: Optional[str]) -> Optional[int]:
        """Parse number of comments from Reddit format.

        Args:
            comments_text: Comments text to parse

        Returns:
            Number of comments as integer or None
        """
        if not comments_text:
            return None

        try:
            # Extract number from comments text
            import re

            match = re.search(r"\b(\d+)\b", comments_text)
            if match:
                return int(match.group(1))
        except Exception:
            pass

        return None

    def _extract_url(self, title_element: Optional[str]) -> Optional[str]:
        """Extract post URL from title element.

        Args:
            title_element: Title element HTML

        Returns:
            Post URL or None
        """
        if not title_element:
            return None

        try:
            # Parse URL from title link
            from bs4 import BeautifulSoup

            soup = BeautifulSoup(title_element, "html.parser")
            link = soup.find("a")
            if link:
                url = link.get("href")
                if not url.startswith("http"):
                    url = f"{self.BASE_URL}{url}"
                return url
        except Exception:
            pass

        return None

    def _calculate_confidence(self, post: Dict[str, Any]) -> float:
        """Calculate confidence score for extracted data.

        Args:
            post: Extracted post data

        Returns:
            Confidence score between 0 and 1
        """
        # Calculate field presence score
        field_score = self._calculate_field_presence(post)

        # Calculate quality factors score
        quality_score = self._calculate_quality_factors(post)

        # Calculate extracted data score
        data_score = self._calculate_extracted_data(post)

        # Combine scores with weights
        total_score = (
            field_score * FIELD_PRESENCE_WEIGHT
            + quality_score * QUALITY_FACTORS_WEIGHT
            + data_score * EXTRACTED_DATA_WEIGHT
        )

        return total_score

    def _calculate_field_presence(self, post: Dict[str, Any]) -> float:
        """Calculate field presence score.

        Args:
            post: Post data

        Returns:
            Score between 0 and 1
        """
        fields = ["title", "content", "author", "subreddit", "score", "num_comments"]
        present = sum(1 for field in fields if post.get(field))
        return present / len(fields)

    def _calculate_quality_factors(self, post: Dict[str, Any]) -> float:
        """Calculate quality factors score.

        Args:
            post: Post data

        Returns:
            Score between 0 and 1
        """
        factors = [
            post.get("score", 0) > 10,
            post.get("num_comments", 0) > 5,
            len(post.get("selftext", "")) > 100,
            post.get("upvote_ratio", 0) > 0.75,
        ]
        return sum(1 for factor in factors if factor) / len(factors)

    def _calculate_extracted_data(self, post: Dict[str, Any]) -> float:
        """Calculate extracted data score.

        Args:
            post: Post data

        Returns:
            Score between 0 and 1
        """
        data_fields = ["compounds", "effects", "mechanisms", "safety"]
        present = sum(1 for field in data_fields if post.get(field))
        return present / len(data_fields)
