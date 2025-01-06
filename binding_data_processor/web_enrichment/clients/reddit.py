"""Reddit integration for Crawl4AI.

This module provides a specialized Crawl4AI client for Reddit that:
1. Uses Crawl4AI for web scraping Reddit pages
2. Handles anti-bot detection
3. Extracts community data about compounds
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from datetime import datetime

from crawl4ai import AsyncWebCrawler, Config
from bs4 import BeautifulSoup

from ..base_client import BaseWebClient
from ..validation.schema import BaseSchema


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


class RedditCrawl4AIClient(BaseWebClient):
    """Enhanced Reddit client using Crawl4AI."""

    BASE_URL = "https://old.reddit.com"  # Use old.reddit.com for better scraping
    SUBREDDITS = [
        "researchchemicals",
        "nootropics",
        "DrugNerds",
        "Psychonaut",
        "psychopharmacology",
    ]

    def __init__(
        self,
        llm_provider: str = "ollama/llama2",
        api_token: Optional[str] = None,
        **kwargs,
    ):
        """Initialize client.

        Args:
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            **kwargs: Additional arguments passed to BaseWebClient
        """
        super().__init__(**kwargs)
        self.llm_provider = llm_provider
        self.api_token = api_token

        # Configure Crawl4AI with anti-bot detection avoidance
        self.config = Config(
            javascript=Config.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    ".thing",
                    ".title",
                    ".entry",
                ],
                # Use stealth mode to avoid detection
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
                rotation=True,  # Rotate proxies to avoid blocks
                retry_count=3,
            ),
            rate_limit=Config.RateLimit(
                requests_per_minute=10,  # Conservative rate limiting
                delay_after_failure=60,
            ),
        )

    async def search_posts(
        self,
        query: str,
        subreddits: Optional[List[str]] = None,
        max_results: int = 100,
    ) -> List[RedditPostData]:
        """Search Reddit posts.

        Args:
            query: Search query
            subreddits: Optional list of subreddits to search (defaults to self.SUBREDDITS)
            max_results: Maximum number of results to return

        Returns:
            List of post data
        """
        results = []
        subreddits = subreddits or self.SUBREDDITS

        try:
            async with AsyncWebCrawler() as crawler:
                for subreddit in subreddits:
                    if len(results) >= max_results:
                        break

                    # Construct search URL
                    url = f"{self.BASE_URL}/r/{subreddit}/search?q={query}&restrict_sr=on"

                    # Scrape with Crawl4AI
                    result = await crawler.arun(urls=[url], config=self.config)

                    if not result.success:
                        self.logger.error(f"Failed to scrape r/{subreddit}")
                        continue

                    # Extract data using both CSS and LLM
                    content = result.extracted_content[0]
                    posts = self._parse_search_results(content, subreddit)

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
                                "source": "reddit",
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

    def _parse_search_results(self, content: Dict[str, Any], subreddit: str) -> List[Dict[str, Any]]:
        """Parse search results from extracted content.

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
                "num_comments": self._parse_num_comments(num_comments[i]) if i < len(num_comments) else None,
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
        score = 0.0
        total = 0

        # Check key fields presence
        fields = ["title", "content", "author", "subreddit", "score", "num_comments", "url"]
        for field in fields:
            total += 1
            if post.get(field):
                score += 1.0

        # Check extracted data quality
        if post.get("compounds"):
            score += 1.0
            total += 1
        if post.get("effects"):
            score += 1.0
            total += 1
        if post.get("mechanisms"):
            score += 1.0
            total += 1
        if post.get("safety"):
            score += 1.0
            total += 1

        # Calculate final score
        return score / total if total > 0 else 0.0
