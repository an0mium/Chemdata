"""Bluelight integration for Crawl4AI.

This module provides a specialized Crawl4AI client for Bluelight that:
1. Uses Crawl4AI for web scraping Bluelight pages
2. Handles anti-bot detection
3. Extracts community data about compounds
"""

import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from datetime import datetime

from crawl4ai import AsyncWebCrawler, Config
from bs4 import BeautifulSoup

from ..base_client import BaseWebClient
from ..validation.schema import BaseSchema
from ..storage.bluelight_storage import BluelightStorage

# Configure logging
logger = logging.getLogger(__name__)


@dataclass
class BluelightPostData(BaseSchema):
    """Schema for Bluelight post data."""

    title: str
    content: Optional[str]
    author: Optional[str]
    subforum: str
    date: Optional[str]
    compounds: List[str]
    effects: List[str]
    mechanisms: List[str]
    safety_notes: List[str]
    confidence: float
    metadata: Dict[str, Any]


class BluelightCrawl4AIClient(BaseWebClient):
    """Enhanced Bluelight client using Crawl4AI."""

    BASE_URL = "https://bluelight.org/xf"
    SUBFORUMS = [
        "advanced-drug-discussion",
        "psychedelic-drugs",
        "other-drugs",
        "scientific-research",
        "drug-studies",
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
                    ".message",
                    ".message-content",
                    ".message-attribution",
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
                    "title": ".p-title-value",
                    "content": ".message-content .bbWrapper",
                    "author": ".message-name",
                    "subforum": ".p-breadcrumbs li:last-child",
                    "date": ".message-attribution-main time",
                },
            ),
            proxy=Config.Proxy(
                enabled=True,
                rotation=True,  # Rotate proxies to avoid blocks
                retry_count=3,
            ),
            rate_limit=Config.RateLimit(
                requests_per_minute=5,  # Very conservative rate limiting
                delay_after_failure=120,
            ),
        )

    async def search_posts(
        self,
        query: str,
        subforums: Optional[List[str]] = None,
        max_results: int = 100,
    ) -> List[BluelightPostData]:
        """Search Bluelight posts.

        Args:
            query: Search query
            subforums: Optional list of subforums to search (defaults to self.SUBFORUMS)
            max_results: Maximum number of results to return

        Returns:
            List of post data
        """
        results = []
        subforums = subforums or self.SUBFORUMS

        try:
            async with AsyncWebCrawler() as crawler:
                for subforum in subforums:
                    if len(results) >= max_results:
                        break

                    # Construct search URL
                    url = f"{self.BASE_URL}/search/search?q={query}&c[]={subforum}"

                    # Scrape with Crawl4AI
                    result = await crawler.arun(urls=[url], config=self.config)

                    if not result.success:
                        logger.error(f"Failed to scrape {subforum}")
                        continue

                    # Extract data using both CSS and LLM
                    content = result.extracted_content[0]
                    posts = self._parse_search_results(content, subforum)

                    if not posts:
                        continue

                    for post in posts:
                        # Create post data object
                        post_data = BluelightPostData(
                            title=post.get("title", ""),
                            content=post.get("content"),
                            author=post.get("author"),
                            subforum=post.get("subforum", subforum),
                            date=post.get("date"),
                            compounds=post.get("compounds", []),
                            effects=post.get("effects", []),
                            mechanisms=post.get("mechanisms", []),
                            safety_notes=post.get("safety", []),
                            confidence=self._calculate_confidence(post),
                            metadata={
                                "source": "bluelight",
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
            logger.error(f"Error scraping Bluelight: {str(e)}")

            # Store error for monitoring
            storage = BluelightStorage()
            storage.store_error(
                error=str(e),
                metadata={
                    "query": query,
                    "subforums": subforums,
                    "max_results": max_results,
                },
            )

        return results

    def _parse_search_results(
        self,
        content: Dict[str, Any],
        subforum: str,
    ) -> List[Dict[str, Any]]:
        """Parse search results from extracted content.

        Args:
            content: Extracted content from Crawl4AI
            subforum: Current subforum being scraped

        Returns:
            List of parsed post data
        """
        posts = []

        # Get post elements
        titles = content.get("title", [])
        contents = content.get("content", [])
        authors = content.get("author", [])
        subforums = content.get("subforum", [])
        dates = content.get("date", [])

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
                "subforum": subforums[i] if i < len(subforums) else subforum,
                "date": self._parse_date(dates[i]) if i < len(dates) else None,
                "url": self._extract_url(titles[i]) if i < len(titles) else None,
                "compounds": compounds,
                "effects": effects,
                "mechanisms": mechanisms,
                "safety": safety,
            }
            posts.append(post)

        # Store posts in database
        storage = BluelightStorage()
        for post in posts:
            storage.store_post(
                post_id=post["url"].split("/")[-1],
                data=post,
            )

        return posts

    def _parse_date(self, date_text: Optional[str]) -> Optional[str]:
        """Parse date from Bluelight format.

        Args:
            date_text: Date text to parse

        Returns:
            ISO format date string or None
        """
        if not date_text:
            return None

        try:
            # Parse datetime from text
            from dateutil.parser import parse

            dt = parse(date_text)
            return dt.isoformat()
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
        fields = ["title", "content", "author", "subforum", "date", "url"]
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
