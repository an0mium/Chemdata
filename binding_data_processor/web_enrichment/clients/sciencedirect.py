"""Enhanced ScienceDirect client with web scraping and API access.

This module provides a comprehensive ScienceDirect client that:
1. Uses Crawl4AI for advanced web scraping
2. Provides LLM-enhanced extraction
3. Includes validation and caching
4. Tracks citations and impact
5. Handles anti-bot detection
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from datetime import datetime
import logging
import re

from crawl4ai import AsyncWebCrawler, BrowserConfig
from bs4 import BeautifulSoup

from ..clients.base import WebClient  # Updated import path
from ..validation.enhanced import EnhancedValidator
from ...models.compound.base.core import CompoundData  # Updated import path
from ...pipeline.infrastructure.rate_limiter import RateLimiter, RateLimit
from ...pipeline.infrastructure.circuit_breaker import CircuitBreaker, ErrorCode, McpError

# API Rate Limits
API_RATE_LIMIT = RateLimit(requests=6, period=1)  # 6 requests per second
SCRAPE_RATE_LIMIT = RateLimit(requests=1, period=1)  # 1 request per second

# Weekly quota tracking (10,000 requests per week for academic)
WEEKLY_QUOTA = 10000

# Circuit Breaker Config
CIRCUIT_FAILURE_THRESHOLD = 5
CIRCUIT_RESET_TIMEOUT = 300  # 5 minutes
CIRCUIT_HALF_OPEN_TIMEOUT = 60  # 1 minute


@dataclass
class ScienceDirectResearchData:
    """Schema for ScienceDirect research data."""

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


class EnhancedScienceDirectClient(WebClient):  # Updated parent class
    """Enhanced ScienceDirect client with anti-bot detection."""

    BASE_URL = "https://www.sciencedirect.com"

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
        weekly_quota: int = WEEKLY_QUOTA,
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
        """
        super().__init__(cache_dir=cache_dir)
        self.validator = validator or EnhancedValidator()
        self.logger = logging.getLogger(__name__)

        # Initialize rate limiters
        self.rate_limiter = RateLimiter()
        self.rate_limiter.add_limit("api", API_RATE_LIMIT)
        self.rate_limiter.add_limit("scrape", SCRAPE_RATE_LIMIT)

        # Initialize circuit breakers
        self.api_circuit = CircuitBreaker(
            "sciencedirect_api",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )
        self.scrape_circuit = CircuitBreaker(
            "sciencedirect_scrape",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )

        # Weekly quota tracking
        self.weekly_quota = weekly_quota
        self.request_count = 0
        self.quota_reset = datetime.now()

        # Configure Crawl4AI with anti-bot detection avoidance
        self.config = BrowserConfig(
            javascript=BrowserConfig.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    ".ResultItem",  # Article container
                    ".Author",  # Author info
                    ".Abstract",  # Abstract
                    ".DOI",  # DOI
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
                    "title": ".ResultItem h2",
                    "abstract": ".Abstract",
                    "authors": ".Author",
                    "journal": ".Journal",
                    "year": ".Year",
                    "citations": ".CitationCount",
                    "doi": ".DOI",
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

    async def _check_quota(self) -> None:
        """Check and update weekly quota."""
        # Reset quota if a week has passed
        now = datetime.now()
        if (now - self.quota_reset).days >= 7:
            self.request_count = 0
            self.quota_reset = now

        # Check quota
        if self.request_count >= self.weekly_quota:
            raise McpError(
                ErrorCode.QuotaExceeded,
                f"Weekly quota of {self.weekly_quota} requests exceeded",
            )

        self.request_count += 1

    async def _process_search_page(
        self,
        crawler: AsyncWebCrawler,
        url: str,
        page_num: int,
    ) -> Optional[List[Dict[str, Any]]]:
        """Process a single search results page.

        Args:
            crawler: AsyncWebCrawler instance
            url: Page URL
            page_num: Page number

        Returns:
            List of papers if successful, None otherwise
        """
        # Wait for scrape rate limit
        self.rate_limiter.wait("scrape")

        # Crawl search results
        result = await crawler.arun(urls=[url], config=self.config)
        if not result.success:
            error_text = str(result.error) if result.error else "Unknown error"
            self.logger.error(f"Failed to scrape page {page_num}: {error_text}")
            self.api_circuit.record_failure()
            return None

        # Extract papers
        content = result.extracted_content[0]
        papers = self._parse_search_results(content)
        if not papers:
            return None

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

        return page_results

    async def _get_cached_results(
        self,
        query: str,
        start: int,
    ) -> Optional[List[ScienceDirectResearchData]]:
        """Get cached search results.

        Args:
            query: Search query
            start: Start index

        Returns:
            Cached results if found, None otherwise
        """
        cache_key = f"search:{query}:{start}"
        return self.cache.get(cache_key)

    async def _cache_results(
        self,
        query: str,
        start: int,
        results: List[ScienceDirectResearchData],
    ) -> None:
        """Cache search results.

        Args:
            query: Search query
            start: Start index
            results: Results to cache
        """
        cache_key = f"search:{query}:{start}"
        self.cache.set(cache_key, results)

    async def search_articles(
        self,
        query: str,
        max_results: int = 100,
        from_year: Optional[int] = None,
        to_year: Optional[int] = None,
        sort: str = "relevance",
    ) -> List[ScienceDirectResearchData]:
        """Search ScienceDirect articles.

        Args:
            query: Search query
            max_results: Maximum number of results
            from_year: Optional start year
            to_year: Optional end year
            sort: Sort order (relevance, date)

        Returns:
            List of research paper data
        """
        # Check circuit breaker
        if not self.api_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "ScienceDirect API circuit breaker open",
            )

        try:
            # Check quota and rate limit
            await self._check_quota()
            self.rate_limiter.wait("api")

            results = []
            start = 0
            results_per_page = 25  # ScienceDirect shows 25 results per page

            async with AsyncWebCrawler() as crawler:
                while len(results) < max_results:
                    # Check cache
                    if cached := await self._get_cached_results(query, start):
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

                    # Process page
                    page_results = await self._process_search_page(
                        crawler=crawler,
                        url=url,
                        page_num=start // results_per_page + 1,
                    )
                    if not page_results:
                        break

                    # Cache and add results
                    await self._cache_results(query, start, page_results)
                    results.extend(page_results)

                    if len(results) >= max_results:
                        break

                    start += results_per_page

            # Record success if we got any results
            if results:
                self.api_circuit.record_success()

            return results[:max_results]

        except Exception as e:
            self.api_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"Error searching ScienceDirect articles: {str(e)}",
            )

    async def get_compound_articles(
        self,
        compound: CompoundData,
        max_results: int = 100,
    ) -> List[ScienceDirectResearchData]:
        """Get articles about a compound.

        Args:
            compound: CompoundData to search for
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
        dois = content.get("doi", [])

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
                "doi": self._extract_doi(dois[i]) if i < len(dois) else None,
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
    ) -> Optional[ScienceDirectResearchData]:
        """Extract article data.

        Args:
            paper: Paper data
            result: Crawl result

        Returns:
            Article data if successful
        """
        try:
            # Create article object
            article = ScienceDirectResearchData(
                title=paper.get("title", ""),
                abstract=paper.get("abstract"),
                authors=paper.get("authors", []),
                journal=paper.get("journal"),
                year=paper.get("year"),
                doi=paper.get("doi"),
                citations=paper.get("citations"),
                compounds=paper.get("compounds", []),
                effects=paper.get("effects", []),
                mechanisms=paper.get("mechanisms", []),
                safety_notes=paper.get("safety", []),
                binding_data=[],  # No binding data from ScienceDirect
                confidence=self._calculate_confidence(paper),
                metadata={
                    "source": "sciencedirect",
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

    async def _validate_article(self, article: ScienceDirectResearchData) -> bool:
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
        """Build ScienceDirect search URL.

        Args:
            query: Search query
            from_year: Optional start year
            to_year: Optional end year
            sort: Sort order
            start: Start index for pagination

        Returns:
            Search URL
        """
        url = f"{self.BASE_URL}/search"
        params = [
            ("qs", query),
            ("offset", str(start)),
        ]

        if from_year:
            params.append(("date", f"{from_year}"))
        if to_year:
            params.append(("end", f"{to_year}"))

        if sort == "date":
            params.append(("sortBy", "date"))
        else:
            params.append(("sortBy", "relevance"))

        return f"{url}?{'&'.join(f'{k}={v}' for k, v in params)}"

    def _build_compound_query(self, compound: CompoundData) -> str:
        """Build search query for compound.

        Args:
            compound: CompoundData to search for

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
        article: ScienceDirectResearchData,
        compound: CompoundData,
    ) -> bool:
        """Check if article is relevant to compound.

        Args:
            article: Article to check
            compound: CompoundData to check against

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
        """Parse year from ScienceDirect format.

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
        """Parse citation count from ScienceDirect format.

        Args:
            citations_text: Citations text to parse

        Returns:
            Number of citations or None
        """
        if not citations_text:
            return None

        try:
            # Extract number from citation count
            match = re.search(r"\b(\d+)\b", citations_text)
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

    def _extract_doi(self, doi_text: Optional[str]) -> Optional[str]:
        """Parse DOI from text.

        Args:
            doi_text: DOI text

        Returns:
            DOI or None
        """
        if not doi_text:
            return None

        try:
            # Look for DOI patterns
            patterns = [
                r"10\.\d{4,}/[-._;()/:\w]+",  # Full DOI pattern
                r"doi\.org/([^/\s]+/[^/\s]+)",  # DOI.org URL pattern
                r"doi:([^/\s]+/[^/\s]+)",  # DOI prefix pattern
            ]
            for pattern in patterns:
                match = re.search(pattern, doi_text)
                if match:
                    return match.group(1) if "/" in match.group(1) else match.group()
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
