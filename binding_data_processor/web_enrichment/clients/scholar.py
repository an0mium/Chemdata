"""Enhanced Google Scholar client combining web scraping and anti-bot detection.

This module provides a comprehensive Google Scholar client that:
1. Uses Crawl4AI for advanced web scraping
2. Handles anti-bot detection with stealth mode
3. Provides proxy rotation and rate limiting
4. Includes validation and caching
5. Uses LLM for enhanced extraction
6. Provides robust error recovery
7. Includes comprehensive metadata extraction
"""

from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
from datetime import datetime, timedelta
import logging
import re
import asyncio
from collections import deque

from crawl4ai import AsyncWebCrawler, BrowserConfig
from bs4 import BeautifulSoup

from ..clients.base import WebClient  # Updated import path
from ..validation.enhanced import EnhancedValidator
from ...models.compound import CompoundData


@dataclass
class ScholarResearchData:
    """Schema for Google Scholar research data."""

    title: str
    abstract: Optional[str]
    authors: List[str]
    journal: Optional[str]
    year: Optional[int]
    doi: Optional[str] = None
    citations: Optional[int] = None
    compounds: List[str] = None
    effects: List[str] = None
    mechanisms: List[str] = None
    safety_notes: List[str] = None
    binding_data: List[Dict[str, Any]] = None
    confidence: float = 0.0
    metadata: Dict[str, Any] = None


class SessionManager:
    """Enhanced session management with error recovery."""

    def __init__(
        self,
        max_retries: int = 3,
        backoff_factor: float = 1.5,
        max_concurrent: int = 5,
        rate_limit: int = 10,
        rate_window: int = 60,
    ):
        """Initialize session manager.

        Args:
            max_retries: Maximum retry attempts
            backoff_factor: Exponential backoff factor
            max_concurrent: Maximum concurrent sessions
            rate_limit: Requests per window
            rate_window: Window size in seconds
        """
        self.max_retries = max_retries
        self.backoff_factor = backoff_factor
        self.semaphore = asyncio.Semaphore(max_concurrent)
        self.rate_limit = rate_limit
        self.rate_window = rate_window
        self.request_times = deque(maxlen=rate_limit)
        self.logger = logging.getLogger(__name__)

    async def execute(self, func, *args, **kwargs) -> Tuple[bool, Any]:
        """Execute function with session management.

        Args:
            func: Async function to execute
            *args: Function arguments
            **kwargs: Function keyword arguments

        Returns:
            Tuple of (success, result)
        """
        async with self.semaphore:
            for attempt in range(self.max_retries):
                try:
                    # Check rate limit
                    await self._check_rate_limit()

                    # Execute function
                    result = await func(*args, **kwargs)
                    self.request_times.append(datetime.now())
                    return True, result

                except Exception as e:
                    delay = self.backoff_factor**attempt
                    self.logger.warning(f"Attempt {attempt + 1} failed: {str(e)}. " f"Retrying in {delay:.1f}s...")
                    await asyncio.sleep(delay)

            self.logger.error(f"All {self.max_retries} attempts failed")
            return False, None

    async def _check_rate_limit(self):
        """Check and enforce rate limiting."""
        if len(self.request_times) < self.rate_limit:
            return

        # Calculate time until next allowed request
        window_start = datetime.now() - timedelta(seconds=self.rate_window)
        while self.request_times and self.request_times[0] <= window_start:
            self.request_times.popleft()

        if len(self.request_times) >= self.rate_limit:
            delay = (self.request_times[0] + timedelta(seconds=self.rate_window) - datetime.now()).total_seconds()
            self.logger.info(f"Rate limit reached. Waiting {delay:.1f}s...")
            await asyncio.sleep(delay)


class EnhancedScholarClient(WebClient):  # Updated parent class
    """Enhanced Google Scholar client with anti-bot detection."""

    BASE_URL = "https://scholar.google.com"

    def __init__(
        self,
        validator: Optional[EnhancedValidator] = None,
        llm_provider: str = "ollama/llama2",
        api_token: Optional[str] = None,
        cache_dir: Optional[str] = None,
        proxy_enabled: bool = True,
        proxy_rotation: bool = True,
        proxy_retry_count: int = 3,
        rate_limit: int = 10,
        rate_limit_delay: int = 60,
        max_concurrent: int = 5,
    ):
        """Initialize client.

        Args:
            validator: Optional data validator
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            cache_dir: Optional cache directory
            proxy_enabled: Whether to use proxies
            proxy_rotation: Whether to rotate proxies
            proxy_retry_count: Number of proxy retries
            rate_limit: Requests per minute
            rate_limit_delay: Delay after failure in seconds
            max_concurrent: Maximum concurrent sessions
        """
        super().__init__(cache_dir=cache_dir)
        self.validator = validator or EnhancedValidator()
        self.logger = logging.getLogger(__name__)
        self.session_manager = SessionManager(
            max_retries=proxy_retry_count,
            rate_limit=rate_limit,
            rate_window=rate_limit_delay,
            max_concurrent=max_concurrent,
        )

        # Configure Crawl4AI with anti-bot detection avoidance
        self.config = BrowserConfig(
            javascript=BrowserConfig.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    ".gs_ri",  # Article container
                    ".gs_a",  # Author info
                    ".gs_rs",  # Abstract
                ],
                stealth_mode=True,  # Enable stealth mode
            ),
            screenshot=BrowserConfig.Screenshot(enabled=True, full_page=True),
            extraction=BrowserConfig.Extraction(
                llm=BrowserConfig.LLM(
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
                    "title": ".gs_rt a",
                    "abstract": ".gs_rs",
                    "authors": ".gs_a a",
                    "journal": ".gs_a",
                    "year": ".gs_a",
                    "citations": ".gs_fl a:contains('Cited by')",
                },
            ),
            proxy=BrowserConfig.Proxy(
                enabled=proxy_enabled,
                rotation=proxy_rotation,
                retry_count=proxy_retry_count,
            ),
            rate_limit=BrowserConfig.RateLimit(
                requests_per_minute=rate_limit,
                delay_after_failure=rate_limit_delay,
            ),
        )

    async def search_articles(
        self,
        query: str,
        max_results: int = 100,
        from_year: Optional[int] = None,
        to_year: Optional[int] = None,
        sort: str = "relevance",
    ) -> List[ScholarResearchData]:
        """Search Google Scholar articles.

        Args:
            query: Search query
            max_results: Maximum number of results
            from_year: Optional start year
            to_year: Optional end year
            sort: Sort order (relevance, date)

        Returns:
            List of research paper data
        """
        results = []
        start = 0
        results_per_page = 10

        try:
            async with AsyncWebCrawler() as crawler:
                while len(results) < max_results:
                    # Check cache
                    cache_key = f"search:{query}:{start}"
                    if cached := self.cache.get(cache_key):
                        results.extend(cached)
                        start += results_per_page
                        continue

                    # Build search URL
                    url = self._build_search_url(
                        query=query,
                        from_year=from_year,
                        to_year=to_year,
                        sort=sort,
                        start=start,
                    )

                    # Execute search with session management
                    success, result = await self.session_manager.execute(
                        crawler.arun,
                        urls=[url],
                        config=self.config,
                    )
                    if not success:
                        self.logger.error(f"Failed to scrape page {start//10 + 1}")
                        break

                    # Extract and process articles
                    content = result.extracted_content[0]
                    papers = self._parse_search_results(content)
                    if not papers:
                        break

                    # Process each paper
                    page_results = []
                    for paper in papers:
                        try:
                            # Extract data
                            article = await self._extract_article_data(paper, result)
                            if not article:
                                continue

                            # Validate
                            if await self._validate_article(article):
                                page_results.append(article)

                        except Exception as e:
                            self.logger.error(f"Error processing paper: {str(e)}")
                            continue

                    # Cache page results
                    if page_results:
                        self.cache.set(cache_key, page_results)
                        results.extend(page_results)

                    if len(results) >= max_results:
                        break

                    start += results_per_page

            return results[:max_results]

        except Exception as e:
            self.logger.error(f"Error searching articles: {str(e)}")
            return []

    async def get_compound_articles(
        self,
        compound: CompoundData,
        max_results: int = 100,
    ) -> List[ScholarResearchData]:
        """Get articles about a compound.

        Args:
            compound: Compound to search for
            max_results: Maximum number of results

        Returns:
            List of relevant articles
        """
        # Check cache
        cache_key = f"compound:{compound.name}"
        if cached := self.cache.get(cache_key):
            return cached

        # Build search query
        query = self._build_compound_query(compound)

        # Search articles
        articles = await self.search_articles(query, max_results=max_results)

        # Filter by relevance
        relevant = []
        for article in articles:
            if await self._is_relevant(article, compound):
                relevant.append(article)

        # Cache results
        self.cache.set(cache_key, relevant)
        return relevant

    async def track_citations(
        self,
        article: ScholarResearchData,
        history_days: int = 30,
    ) -> Dict[str, Any]:
        """Track citation changes over time.

        Args:
            article: Article to track
            history_days: Number of days of history

        Returns:
            Citation tracking data
        """
        if not article.metadata.get("url"):
            return {}

        try:
            async with AsyncWebCrawler() as crawler:
                # Get current citations
                success, result = await self.session_manager.execute(
                    crawler.arun,
                    urls=[article.metadata["url"]],
                    config=self.config,
                )
                if not success:
                    return {}

                current_citations = self._parse_citations(result.extracted_content[0].get("citations", [""])[0])

                # Calculate citation velocity
                if article.citations is not None and current_citations is not None:
                    citation_change = current_citations - article.citations
                    days_elapsed = (datetime.now() - datetime.fromisoformat(article.metadata["extraction_date"])).days
                    velocity = citation_change / max(days_elapsed, 1)
                else:
                    velocity = 0

                return {
                    "initial_citations": article.citations,
                    "current_citations": current_citations,
                    "citation_change": citation_change if article.citations is not None else None,
                    "citation_velocity": velocity,
                    "tracking_period_days": days_elapsed,
                    "last_updated": datetime.now().isoformat(),
                }

        except Exception as e:
            self.logger.error(f"Error tracking citations: {str(e)}")
            return {}

    def _parse_search_results(self, content: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse search results from extracted content.

        Args:
            content: Extracted content from Crawl4AI

        Returns:
            List of parsed paper data
        """
        papers = []

        # Get paper elements
        titles = content.get("title", [])
        abstracts = content.get("abstract", [])
        author_lists = content.get("authors", [])
        journals = content.get("journal", [])
        years = content.get("year", [])
        citations = content.get("citations", [])

        # Get LLM-extracted data
        compounds = content.get("compounds", [])
        effects = content.get("effects", [])
        mechanisms = content.get("mechanisms", [])
        safety = content.get("safety", [])

        # Combine data for each paper
        for i in range(len(titles)):
            paper = {
                "title": titles[i] if i < len(titles) else None,
                "abstract": abstracts[i] if i < len(abstracts) else None,
                "authors": author_lists[i] if i < len(author_lists) else [],
                "journal": journals[i] if i < len(journals) else None,
                "year": self._parse_year(years[i]) if i < len(years) else None,
                "citations": self._parse_citations(citations[i]) if i < len(citations) else None,
                "url": self._extract_url(titles[i]) if i < len(titles) else None,
                "compounds": compounds,
                "effects": effects,
                "mechanisms": mechanisms,
                "safety": safety,
            }
            papers.append(paper)

        return papers

    async def _extract_article_data(
        self,
        paper: Dict[str, Any],
        result: Any,
    ) -> Optional[ScholarResearchData]:
        """Extract article data.

        Args:
            paper: Paper data
            result: Crawl result

        Returns:
            Article data if successful
        """
        try:
            # Extract DOI from paper URL
            doi = self._extract_doi(paper.get("url", ""))

            # Create article object
            article = ScholarResearchData(
                title=paper.get("title", ""),
                abstract=paper.get("abstract"),
                authors=paper.get("authors", []),
                journal=paper.get("journal"),
                year=paper.get("year"),
                doi=doi,
                citations=paper.get("citations"),
                compounds=paper.get("compounds", []),
                effects=paper.get("effects", []),
                mechanisms=paper.get("mechanisms", []),
                safety_notes=paper.get("safety", []),
                binding_data=[],  # No binding data from Google Scholar
                confidence=self._calculate_confidence(paper),
                metadata={
                    "source": "google_scholar",
                    "extraction_date": datetime.now().isoformat(),
                    "url": paper.get("url"),
                    "screenshot": str(result.screenshot),
                    "javascript_logs": result.javascript_logs,
                },
            )

            return article

        except Exception as e:
            self.logger.error(f"Error extracting article data: {str(e)}")
            return None

    async def _validate_article(self, article: ScholarResearchData) -> bool:
        """Validate article data.

        Args:
            article: Article to validate

        Returns:
            True if valid
        """
        try:
            # Convert to dict for validation
            data = {
                "title": article.title,
                "abstract": article.abstract,
                "authors": article.authors,
                "journal": article.journal,
                "year": article.year,
                "doi": article.doi,
                "citations": article.citations,
                "compounds": article.compounds,
                "effects": article.effects,
                "mechanisms": article.mechanisms,
                "safety_notes": article.safety_notes,
                "binding_data": article.binding_data,
            }

            # Validate
            result = await self.validator.validate(
                data=data,
                source_data={},  # No cross-source validation for now
                temporal_data=[],  # No temporal validation for now
                community_data=[],  # No community validation for now
            )

            return result.is_valid

        except Exception as e:
            self.logger.error(f"Error validating article: {str(e)}")
            return False

    def _build_search_url(
        self,
        query: str,
        from_year: Optional[int] = None,
        to_year: Optional[int] = None,
        sort: str = "relevance",
        start: int = 0,
    ) -> str:
        """Build Google Scholar search URL.

        Args:
            query: Search query
            from_year: Optional start year
            to_year: Optional end year
            sort: Sort order
            start: Start index for pagination

        Returns:
            Search URL
        """
        url = f"{self.BASE_URL}/scholar?q={query}&start={start}"

        if from_year:
            url += f"&as_ylo={from_year}"
        if to_year:
            url += f"&as_yhi={to_year}"

        if sort == "date":
            url += "&scisbd=1"  # Sort by date
        else:
            url += "&scisbd=0"  # Sort by relevance

        return url

    def _build_compound_query(self, compound: CompoundData) -> str:
        """Build search query for compound.

        Args:
            compound: Compound to search for

        Returns:
            Search query
        """
        query_parts = []

        # Add identifiers
        if compound.name:
            query_parts.append(compound.name)
        if compound.cas_number:
            query_parts.append(compound.cas_number)
        if compound.iupac_name:
            query_parts.append(compound.iupac_name)

        # Add common names
        if compound.common_name_1 != "N/A":
            query_parts.append(compound.common_name_1)
        if compound.common_name_2 != "N/A":
            query_parts.append(compound.common_name_2)
        if compound.common_name_3 != "N/A":
            query_parts.append(compound.common_name_3)

        # Build query
        return " OR ".join(f'"{term}"' for term in query_parts)

    async def _is_relevant(
        self,
        article: ScholarResearchData,
        compound: CompoundData,
    ) -> bool:
        """Check if article is relevant to compound.

        Args:
            article: Article to check
            compound: Compound to check against

        Returns:
            True if relevant
        """
        try:
            # Check if compound is mentioned
            if not any(
                name.lower() in article.abstract.lower()
                for name in [compound.name, compound.cas_number, compound.iupac_name] + [compound.common_name_1, compound.common_name_2, compound.common_name_3]
                if name and name != "N/A"
            ):
                return False

            # Check for relevant content
            if not article.compounds or not article.effects:
                return False

            # Check for research findings
            if not article.mechanisms:
                return False

            return True

        except Exception as e:
            self.logger.error(f"Error checking relevance: {str(e)}")
            return False

    def _parse_year(self, year_text: Optional[str]) -> Optional[int]:
        """Parse year from Google Scholar format.

        Args:
            year_text: Year text to parse

        Returns:
            Parsed year or None
        """
        if not year_text:
            return None

        try:
            # Extract first 4-digit number
            match = re.search(r"\b\d{4}\b", year_text)
            if match:
                return int(match.group())
        except Exception:
            pass

        return None

    def _parse_citations(self, citations_text: Optional[str]) -> Optional[int]:
        """Parse citation count from Google Scholar format.

        Args:
            citations_text: Citations text to parse

        Returns:
            Number of citations or None
        """
        if not citations_text:
            return None

        try:
            # Extract number after "Cited by"
            match = re.search(r"Cited by (\d+)", citations_text)
            if match:
                return int(match.group(1))
        except Exception:
            pass

        return None

    def _extract_url(self, title_element: Optional[str]) -> Optional[str]:
        """Extract paper URL from title element.

        Args:
            title_element: Title element HTML

        Returns:
            Paper URL or None
        """
        if not title_element:
            return None

        try:
            # Parse URL from title link
            soup = BeautifulSoup(title_element, "html.parser")
            link = soup.find("a")
            if link:
                return link.get("href")
        except Exception:
            pass

        return None

    def _extract_doi(self, url: Optional[str]) -> Optional[str]:
        """Extract DOI from paper URL.

        Args:
            url: Paper URL

        Returns:
            DOI or None
        """
        if not url:
            return None

        try:
            # Look for DOI patterns in URL
            patterns = [
                r"doi\.org/([^/\s]+/[^/\s]+)",
                r"doi:([^/\s]+/[^/\s]+)",
            ]
            for pattern in patterns:
                match = re.search(pattern, url)
                if match:
                    return match.group(1)
        except Exception:
            pass

        return None

    def _calculate_confidence(self, paper: Dict[str, Any]) -> float:
        """Calculate confidence score.

        Args:
            paper: Paper data

        Returns:
            Confidence score between 0 and 1
        """
        score = 0.0
        total = 0

        # Check required fields
        fields = ["title", "abstract", "authors", "journal", "year", "doi"]
        for field in fields:
            total += 1
            if paper.get(field):
                score += 1.0

        # Check extracted data
        if paper.get("compounds"):
            score += 1.0
            total += 1
        if paper.get("effects"):
            score += 1.0
            total += 1
        if paper.get("mechanisms"):
            score += 1.0
            total += 1
        if paper.get("safety"):
            score += 1.0
            total += 1

        # Calculate final score
        return score / total if total > 0 else 0.0
