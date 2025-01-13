"""Enhanced web scraping with Crawl4AI integration.

This module provides a comprehensive web scraping client that combines:
1. Research paper scraping with LLM extraction
2. Patent document scraping with structure analysis  
3. Community forum scraping with sentiment analysis
4. Session-based crawling for dynamic content
5. Circuit breaker for reliability
6. Caching for efficiency
7. Media and link extraction
8. Metadata enrichment
"""

import asyncio
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
from datetime import datetime
from tenacity import retry, stop_after_attempt, wait_exponential
from pydantic import BaseModel, Field

from crawl4ai import (
    AsyncWebCrawler,
    BrowserConfig,
    CrawlerRunConfig,
    CacheMode,
    JavaScriptConfig,
    ScreenshotConfig,
)
from crawl4ai.extraction_strategy import (
    JsonCssExtractionStrategy,
    LLMExtractionStrategy,
    CosineStrategy,
)
from crawl4ai.chunking import (
    RegexChunking,
    NlpSentenceChunking,
    TopicSegmentationChunking,
    FixedLengthWordChunking,
    SlidingWindowChunking,
)

from .base_client import WebClient, WebClientError
from .validation.schema import DataSource
from ..models.psychopharm.types import (
    ResearchData,
    PatentData,
    CommunityData,
)
from ..pipeline.infrastructure.circuit_breaker import CircuitBreaker, CircuitBreakerConfig
from .llm_utils import (
    extract_chemical_info,
    analyze_patent_text,
    NootropicLLMExtractor,
)


# Custom exception classes
class ExtractionError(Exception):
    """Base class for extraction errors."""

    pass


class ValidationError(ExtractionError):
    """Raised when extracted data fails validation."""

    pass


class ConnectionError(ExtractionError):
    """Raised when connection or network issues occur."""

    pass


@dataclass
class ExtractionResult:
    """Container for extraction results."""

    success: bool
    content: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class WebDataSchema(BaseModel):
    """Schema for web data extraction."""

    research: List[ResearchData]
    patents: List[PatentData]
    community: List[CommunityData]
    confidence: float
    metadata: Dict[str, Any] = Field(
        default_factory=lambda: {
            "extraction_date": datetime.now().isoformat(),
            "extraction_config": {},
            "validation_results": [],
            "clustering_results": [],
            "chunking_stats": {
                "total_chunks": 0,
                "avg_chunk_size": 0,
                "overlap_rate": 0.0,
                "strategy_used": "",
            },
        }
    )


class Crawl4AIClient(WebClient):
    """Enhanced web scraping client using Crawl4AI."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        base_url: str = "https://api.crawl4ai.com/v1",
        requests_per_second: float = 1.0,
        max_retries: int = 3,
        timeout: int = 30,
        cache_ttl: int = 3600,
        browser_config: Optional[Dict[str, Any]] = None,
        crawler_config: Optional[Dict[str, Any]] = None,
        circuit_config: Optional[CircuitBreakerConfig] = None,
        logger: Optional[logging.Logger] = None,
        chunking_config: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
        llm_provider: str = "ollama/llama2",
        api_token: Optional[str] = None,
        proxy_enabled: bool = True,
        proxy_rotation: bool = True,
        proxy_retry_count: int = 3,
        rate_limit: int = 10,
        rate_limit_delay: int = 60,
    ):
        """Initialize client.

        Args:
            api_key: Optional Crawl4AI API key for REST API
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            base_url: Base API URL
            requests_per_second: Maximum requests per second
            max_retries: Maximum number of retries
            timeout: Request timeout in seconds
            cache_ttl: Cache TTL in seconds
            browser_config: Optional browser configuration
            crawler_config: Optional crawler configuration
            circuit_config: Optional circuit breaker config
            logger: Optional logger instance
            chunking_config: Optional chunking configuration
            session_id: Optional session ID
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            proxy_enabled: Whether to use proxies
            proxy_rotation: Whether to rotate proxies
            proxy_retry_count: Number of proxy retries
            rate_limit: Requests per minute
            rate_limit_delay: Delay after failure in seconds
        """
        super().__init__(
            name="crawl4ai",
            base_url=base_url,
            data_source=DataSource.COMMUNITY,
            requests_per_second=requests_per_second,
            max_retries=max_retries,
            timeout=timeout,
            cache_ttl=cache_ttl,
            logger=logger,
        )

        # API config
        self.api_key = api_key or os.getenv("CRAWL4AI_API_KEY")

        # Browser config
        self.browser_config = BrowserConfig(
            headless=True,
            verbose=True,
            screenshot=ScreenshotConfig(
                enabled=True,
                full_page=True,
                quality=80,
            ),
            javascript=JavaScriptConfig(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=True,
                stealth_mode=True,
            ),
            **(browser_config or {}),
        )

        # Crawler config
        self.crawler_config = CrawlerRunConfig(
            cache_mode=CacheMode.ENABLED,
            hooks={
                "before_request": self._before_request_hook,
                "after_request": self._after_request_hook,
                "on_error": self._on_error_hook,
            },
            proxy=CrawlerRunConfig.Proxy(
                enabled=proxy_enabled,
                rotation=proxy_rotation,
                retry_count=proxy_retry_count,
            ),
            rate_limit=CrawlerRunConfig.RateLimit(
                requests_per_minute=rate_limit,
                delay_after_failure=rate_limit_delay,
            ),
            **(crawler_config or {}),
        )

        # Circuit breakers
        self.circuit_config = circuit_config or CircuitBreakerConfig(
            failure_threshold=5,
            reset_timeout=300,
        )
        self.research_breaker = CircuitBreaker(
            name="research_scraper",
            config=self.circuit_config,
            logger=self.logger,
        )
        self.patent_breaker = CircuitBreaker(
            name="patent_scraper",
            config=self.circuit_config,
            logger=self.logger,
        )
        self.community_breaker = CircuitBreaker(
            name="community_scraper",
            config=self.circuit_config,
            logger=self.logger,
        )

        # LLM config
        self.llm_strategy = LLMExtractionStrategy(
            input_format="html",
            provider=llm_provider,
            api_token=api_token,
            instruction="Extract structured data with high accuracy",
            temperature=0,
            chunking=chunking_config
            or {
                "strategy": "sliding_window",
                "window_size": 1000,
                "step": 500,
            },
        )

        # Cosine similarity config
        self.cosine_strategy = CosineStrategy(
            word_count_threshold=15,
            max_dist=0.3,
            sim_threshold=0.2,
            model_name="sentence-transformers/all-MiniLM-L6-v2",
            top_k=2,
        )

        # Session config
        self.session_id = session_id or f"session_{int(time.time())}"

        # Chunking strategies
        self.chunkers = {
            "regex": RegexChunking(patterns=[r"\n\n"]),
            "sentence": NlpSentenceChunking(),
            "topic": TopicSegmentationChunking(num_keywords=3),
            "fixed": FixedLengthWordChunking(chunk_size=1000),
            "sliding": SlidingWindowChunking(window_size=1000, step=500),
        }

        # Initialize LLM extractors
        self.nootropic_extractor = NootropicLLMExtractor()

        # Session-based actions
        self.actions = {
            "scroll": lambda page: page.evaluate("window.scrollTo(0, document.body.scrollHeight)"),
            "click": lambda page, selector: page.click(selector),
            "wait": lambda page, selector: page.wait_for_selector(selector),
            "type": lambda page, selector, text: page.type(selector, text),
        }

        # Source URLs with session handling
        self.research_sources = {
            "pubmed": {
                "url": "https://pubmed.ncbi.nlm.nih.gov/",
                "selectors": {
                    "load_more": ".load-more-papers",
                    "next_page": ".next-page",
                    "results": ".paper-result",
                },
                "actions": [
                    {"type": "wait", "selector": ".results-section"},
                    {"type": "scroll", "amount": "full"},
                    {"type": "click", "selector": ".load-more-papers"},
                ],
            },
            "scholar": {
                "url": "https://scholar.google.com/",
                "selectors": {
                    "results": ".gs_r",
                    "next": ".gs_btnPR",
                    "cite": ".gs_fl a:first-child",
                },
                "actions": [
                    {"type": "wait", "selector": "#gs_res_ccl"},
                    {"type": "scroll", "amount": "viewport"},
                ],
            },
            "sciencedirect": {
                "url": "https://www.sciencedirect.com/",
                "selectors": {
                    "results": ".ResultItem",
                    "load_more": ".show-more-results",
                    "abstract": ".abstract.author",
                },
                "actions": [
                    {"type": "wait", "selector": ".search-results"},
                    {"type": "scroll", "amount": "full"},
                    {"type": "click", "selector": ".show-more-results"},
                ],
            },
        }
        self.community_sources = [
            "reddit.com/r/researchchemicals",
            "reddit.com/r/nootropics",
            "reddit.com/r/DrugNerds",
            "bluelight.org",
            "psychonautwiki.org",
        ]

    def _get_headers(self) -> Dict[str, str]:
        """Get request headers.

        Returns:
            Request headers
        """
        return {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }

    async def _before_request_hook(self, url: str) -> None:
        """Hook called before each request."""
        self.logger.info(f"Starting request to {url}")

    async def _after_request_hook(self, url: str, result: Any) -> None:
        """Hook called after each successful request."""
        self.logger.info(f"Completed request to {url}")

    async def _on_error_hook(self, url: str, error: Exception) -> None:
        """Hook called on request error."""
        self.logger.error(f"Error scraping {url}: {str(error)}")

    async def _execute_session_actions(
        self,
        crawler: AsyncWebCrawler,
        url: str,
        actions: List[Dict[str, Any]],
    ) -> Any:
        """Execute session actions sequentially.

        Args:
            crawler: AsyncWebCrawler instance
            url: URL to execute actions on
            actions: List of actions to execute

        Returns:
            Result of final action
        """
        result = None
        for action in actions:
            action_type = action["type"]
            if action_type == "scroll":
                result = await crawler.arun(
                    url=url,
                    session_id=self.session_id,
                    js_code="window.scrollTo(0, document.body.scrollHeight);",
                    wait_for="networkidle",
                )
            elif action_type == "click":
                result = await crawler.arun(
                    url=url,
                    session_id=self.session_id,
                    js_code=f"document.querySelector('{action['selector']}').click();",
                    wait_for=action.get("wait_for", "networkidle"),
                )
            elif action_type == "wait":
                result = await crawler.arun(
                    url=url,
                    session_id=self.session_id,
                    wait_for=action["selector"],
                )
            elif action_type == "type":
                result = await crawler.arun(
                    url=url,
                    session_id=self.session_id,
                    js_code=f"""
                    const el = document.querySelector('{action['selector']}');
                    el.value = '{action['text']}';
                    el.dispatchEvent(new Event('input'));
                    """,
                    wait_for="networkidle",
                )

        return result

    async def _process_text_chunks(
        self,
        text: str,
        chunking_strategy: str = "sliding",
        **kwargs,
    ) -> List[str]:
        """Process text using specified chunking strategy.

        Args:
            text: Text to process
            chunking_strategy: Chunking strategy to use
            **kwargs: Additional arguments for chunker

        Returns:
            List of text chunks
        """
        chunker = self.chunkers.get(chunking_strategy)
        if not chunker:
            raise ValueError(f"Unknown chunking strategy: {chunking_strategy}")

        return chunker.chunk(text, **kwargs)

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=4, max=10),
        reraise=True,
    )
    async def _scrape_research_api(
        self,
        query: str,
        max_results: int,
        min_date: Optional[str],
        max_date: Optional[str],
    ) -> List[ResearchData]:
        """Scrape research papers using REST API.

        Args:
            query: Search query
            max_results: Maximum number of results
            min_date: Minimum publication date
            max_date: Maximum publication date

        Returns:
            List of research paper data
        """
        params = {
            "query": query,
            "max_results": max_results,
            "min_date": min_date,
            "max_date": max_date,
        }

        data = self.request(
            method="GET",
            endpoint="/research",
            params=params,
            headers=self._get_headers(),
            schema_name="research",
        )

        return [ResearchData(**item) for item in data["results"]]

    async def _extract_research_content(
        self,
        crawler: AsyncWebCrawler,
        url: str,
        config: CrawlerRunConfig,
    ) -> Optional[Dict[str, Any]]:
        """Extract research content using multiple strategies.

        Args:
            crawler: AsyncWebCrawler instance
            url: URL to extract from
            config: Crawler configuration

        Returns:
            Extracted content if successful
        """
        # Get full page text
        text = await crawler.arun(
            url=url,
            config=config,
            session_id=self.session_id,
            return_html=True,
        )

        # Process text in chunks
        chunks = await self._process_text_chunks(
            text=text.html,
            chunking_strategy="sliding",
            window_size=1000,
            step=500,
        )

        # Extract from each chunk
        chunk_results = []
        for chunk in chunks:
            result = await crawler.arun(
                url=url,
                config=config,
                session_id=self.session_id,
                content=chunk,
            )
            if result.extracted_content:
                chunk_results.extend(result.extracted_content)

        # Analyze research text with LLM
        if chunk_results:
            for item in chunk_results:
                analysis = analyze_patent_text(
                    title=item.get("title", ""),
                    abstract=item.get("abstract", ""),
                    description=item.get("full_text", ""),
                    claims="",  # Research papers don't have claims
                )
                item["llm_analysis"] = analysis

        return chunk_results

    async def _validate_research_data(
        self,
        data: List[Dict[str, Any]],
        source_name: str,
    ) -> bool:
        """Validate research data.

        Args:
            data: Data to validate
            source_name: Source name for logging

        Returns:
            True if valid
        """
        validation = self.validator.validate(
            data=data,
            source=DataSource.COMMUNITY,
            schema_name="research",
        )

        if not validation.is_valid:
            msg = f"Validation failed for {source_name}: {validation.errors}"
            self.logger.warning(msg)
            return False

        return True

    async def _scrape_research_browser(
        self,
        crawler: AsyncWebCrawler,
        query: str,
        source_name: str,
        source_config: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Scrape research using browser.

        Args:
            crawler: AsyncWebCrawler instance
            query: Search query
            source_name: Source name
            source_config: Source configuration

        Returns:
            List of research data
        """
        # Configure research extraction
        research_schema = {
            "baseSelector": "article",
            "fields": [
                {
                    "name": "title",
                    "selector": "h1,h2",
                    "type": "text",
                },
                {
                    "name": "abstract",
                    "selector": ".abstract,.summary",
                    "type": "text",
                },
                {
                    "name": "authors",
                    "selector": ".authors,.author-list",
                    "type": "text",
                },
                {
                    "name": "date",
                    "selector": ".date,.pub-date",
                    "type": "text",
                },
                {
                    "name": "doi",
                    "selector": ".doi",
                    "type": "text",
                },
            ],
        }
        css_strategy = JsonCssExtractionStrategy(schema=research_schema)

        # Execute session-based actions
        url = f"{source_config['url']}?q={query}"
        await self._execute_session_actions(
            crawler=crawler,
            url=url,
            actions=source_config["actions"],
        )

        # Try extraction strategies
        config = self.crawler_config
        for strategy in [css_strategy, self.cosine_strategy, self.llm_strategy]:
            config.extraction_strategy = strategy
            content = await self._extract_research_content(crawler, url, config)
            if content:
                if await self._validate_research_data(content, source_name):
                    return content

        return []

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=4, max=10),
        reraise=True,
    )
    async def scrape_research(
        self,
        query: str,
        max_results: int = 100,
        min_date: Optional[str] = None,
        max_date: Optional[str] = None,
        use_api: bool = True,
    ) -> List[ResearchData]:
        """Scrape research papers.

        Args:
            query: Search query
            max_results: Maximum number of results
            min_date: Minimum publication date (YYYY-MM-DD)
            max_date: Maximum publication date (YYYY-MM-DD)
            use_api: Whether to use REST API (vs browser)

        Returns:
            List of research paper data
        """
        try:
            with self.research_breaker:
                if use_api and self.api_key:
                    return await self._scrape_research_api(
                        query=query,
                        max_results=max_results,
                        min_date=min_date,
                        max_date=max_date,
                    )
                else:
                    async with AsyncWebCrawler(config=self.browser_config) as crawler:
                        results = []
                        for source_name, source_config in self.research_sources.items():
                            content = await self._scrape_research_browser(
                                crawler=crawler,
                                query=query,
                                source_name=source_name,
                                source_config=source_config,
                            )
                            results.extend(content)

                        return [ResearchData(**item) for item in results[:max_results]]

        except Exception as e:
            raise WebClientError(f"Research scraping failed: {str(e)}")

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=4, max=10),
        reraise=True,
    )
    async def scrape_patents(
        self,
        query: str,
        max_results: int = 100,
        min_date: Optional[str] = None,
        max_date: Optional[str] = None,
        use_api: bool = True,
    ) -> List[PatentData]:
        """Scrape patent documents.

        Args:
            query: Search query
            max_results: Maximum number of results
            min_date: Minimum publication date (YYYY-MM-DD)
            max_date: Maximum publication date (YYYY-MM-DD)
            use_api: Whether to use REST API (vs MCP)

        Returns:
            List of patent data
        """
        try:
            with self.patent_breaker:
                if use_api and self.api_key:
                    # Use REST API
                    params = {
                        "query": query,
                        "max_results": max_results,
                        "min_date": min_date,
                        "max_date": max_date,
                    }

                    data = self.request(
                        method="GET",
                        endpoint="/patents",
                        params=params,
                        headers=self._get_headers(),
                        schema_name="patent",
                    )

                    return [PatentData(**item) for item in data["results"]]
                else:
                    # Use patent search MCP server
                    from ..pipeline.infrastructure.mcp import use_mcp_tool

                    response = await use_mcp_tool(
                        server_name="patent-search",
                        tool_name="search_patents",
                        arguments={
                            "query": query,
                            "max_results": max_results,
                            "min_date": min_date,
                            "max_date": max_date,
                        },
                    )

                    # Extract chemical info and analyze patent text
                    for item in response["results"]:
                        chemicals = extract_chemical_info(
                            description=item.get("description", ""),
                            claims=item.get("claims", ""),
                        )
                        analysis = analyze_patent_text(
                            title=item.get("title", ""),
                            abstract=item.get("abstract", ""),
                            description=item.get("description", ""),
                            claims=item.get("claims", ""),
                        )
                        item["chemical_info"] = chemicals
                        item["llm_analysis"] = analysis

                    return [PatentData(**item) for item in response["results"]]

        except Exception as e:
            raise WebClientError(f"Patent scraping failed: {str(e)}")

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=4, max=10),
        reraise=True,
    )
    async def scrape_community(
        self,
        query: str,
        sources: Optional[List[str]] = None,
        max_results: int = 100,
        min_date: Optional[str] = None,
        max_date: Optional[str] = None,
        use_api: bool = True,
    ) -> List[CommunityData]:
        """Scrape community forums.

        Args:
            query: Search query
            sources: List of sources to search
            max_results: Maximum number of results
            min_date: Minimum post date (YYYY-MM-DD)
            max_date: Maximum post date (YYYY-MM-DD)
            use_api: Whether to use REST API (vs browser)

        Returns:
            List of community data
        """
        try:
            with self.community_breaker:
                if use_api and self.api_key:
                    # Use REST API
                    params = {
                        "query": query,
                        "sources": sources,
                        "max_results": max_results,
                        "min_date": min_date,
                        "max_date": max_date,
                    }

                    data = self.request(
                        method="GET",
                        endpoint="/community",
                        params=params,
                        headers=self._get_headers(),
                        schema_name="forum",
                    )

                    return [CommunityData(**item) for item in data["results"]]
                else:
                    # Use browser scraping
                    async with AsyncWebCrawler(config=self.browser_config) as crawler:
                        sources = sources or self.community_sources

                        results = []
                        for source in sources:
                            # Configure community extraction
                            community_schema = {
                                "baseSelector": ".post,.comment,.entry",
                                "fields": [
                                    {
                                        "name": "title",
                                        "selector": "h1,h2,.title",
                                        "type": "text",
                                    },
                                    {
                                        "name": "content",
                                        "selector": ".content,.body,.text",
                                        "type": "text",
                                    },
                                    {
                                        "name": "author",
                                        "selector": ".author,.username",
                                        "type": "text",
                                    },
                                    {
                                        "name": "date",
                                        "selector": ".date,.timestamp",
                                        "type": "text",
                                    },
                                    {
                                        "name": "source",
                                        "selector": ".source,.forum-name",
                                        "type": "text",
                                    },
                                ],
                            }
                            css_strategy = JsonCssExtractionStrategy(schema=community_schema)

                            # Try CSS extraction first
                            config = self.crawler_config
                            config.extraction_strategy = css_strategy
                            result = await crawler.arun(
                                url=f"https://{source}/search?q={query}",
                                config=config,
                            )

                            # Fall back to LLM if needed
                            if not result.extracted_content:
                                config.extraction_strategy = self.llm_strategy
                                result = await crawler.arun(
                                    url=f"https://{source}/search?q={query}",
                                    config=config,
                                )

                            # Analyze community text with LLM
                            if result.extracted_content:
                                for item in result.extracted_content:
                                    title = item.get("title", "")
                                    content = item.get("content", "")
                                    text = f"{title} {content}"
                                    extractor = self.nootropic_extractor
                                    item["effects"] = extractor.extract_effects(text)
                                    item["mechanisms"] = extractor.extract_mechanisms(text)
                                    item["safety"] = extractor.extract_safety(text)

                            # Validate data
                            validation = self.validator.validate(
                                data=result.extracted_content,
                                source=DataSource.COMMUNITY,
                                schema_name="forum",
                            )

                            if not validation.is_valid:
                                msg = f"Validation failed for {source}: {validation.errors}"
                                self.logger.warning(msg)
                                continue

                            results.extend(result.extracted_content)

                        return [CommunityData(**item) for item in results[:max_results]]

        except Exception as e:
            raise WebClientError(f"Community scraping failed: {str(e)}")

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=4, max=10),
        reraise=True,
    )
    async def scrape_all(
        self,
        query: str,
        max_results: int = 100,
        min_date: Optional[str] = None,
        max_date: Optional[str] = None,
        use_api: bool = True,
    ) -> Dict[str, Any]:
        """Scrape all data sources.

        Args:
            query: Search query
            max_results: Maximum number of results
            min_date: Minimum date (YYYY-MM-DD)
            max_date: Maximum date (YYYY-MM-DD)
            use_api: Whether to use REST API

        Returns:
            Combined results from all sources
        """
        try:
            # Scrape all sources concurrently
            research_task = self.scrape_research(
                query=query,
                max_results=max_results,
                min_date=min_date,
                max_date=max_date,
                use_api=use_api,
            )
            patents_task = self.scrape_patents(
                query=query,
                max_results=max_results,
                min_date=min_date,
                max_date=max_date,
                use_api=use_api,
            )
            community_task = self.scrape_community(
                query=query,
                max_results=max_results,
                min_date=min_date,
                max_date=max_date,
                use_api=use_api,
            )

            research, patents, community = await asyncio.gather(
                research_task,
                patents_task,
                community_task,
            )

            results = {
                "research": research,
                "patents": patents,
                "community": community,
                "metadata": {
                    "query": query,
                    "timestamp": self._get_timestamp(),
                    "source": self.__class__.__name__,
                },
            }

            # Validate combined results
            validation = self.validator.validate(
                data=results,
                source=DataSource.COMMUNITY,
                schema_name="web",
            )

            if not validation.is_valid:
                msg = f"Combined validation failed: {validation.errors}"
                self.logger.warning(msg)

            return results

        except Exception as e:
            raise WebClientError(f"Combined scraping failed: {str(e)}")

    def process_compounds(
        self,
        compounds: List[Any],
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
            # Get compound data
            data = self.get_compound_data(
                name=compound.name,
                cas_number=compound.cas_number,
                use_cache=use_cache,
            )

            if data:
                # Update compound with web data
                compound.web_data = data

    def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get data for a single compound.

        Args:
            name: Compound name
            cas_number: Optional CAS number
            use_cache: Whether to use cached results

        Returns:
            Dictionary of compound data or None if not found
        """
        # Build query
        query = name
        if cas_number:
            query = f"{query} {cas_number}"

        # Run scraping
        try:
            loop = asyncio.get_event_loop()
            results = loop.run_until_complete(self.scrape_all(query))
            return results
        except Exception as e:
            self.logger.error(f"Error getting compound data: {str(e)}")
            return None

    def extract_media(
        self,
        url: str,
        media_types: Optional[List[str]] = None,
    ) -> Dict[str, List[str]]:
        """Extract media from URL.

        Args:
            url: URL to extract from
            media_types: List of media types (images, videos, etc)

        Returns:
            Dictionary of media URLs by type
        """
        params = {
            "url": url,
            "media_types": media_types,
        }

        return self.request(
            method="GET",
            endpoint="/media",
            params=params,
            headers=self._get_headers(),
        )

    def extract_links(
        self,
        url: str,
        depth: int = 1,
        follow_external: bool = False,
    ) -> List[Dict[str, Any]]:
        """Extract links from URL.

        Args:
            url: URL to extract from
            depth: Link crawl depth
            follow_external: Whether to follow external links

        Returns:
            List of link data
        """
        params = {
            "url": url,
            "depth": depth,
            "follow_external": follow_external,
        }

        return self.request(
            method="GET",
            endpoint="/links",
            params=params,
            headers=self._get_headers(),
        )

    def extract_metadata(
        self,
        url: str,
        extract_schema: bool = True,
    ) -> Dict[str, Any]:
        """Extract metadata from URL.

        Args:
            url: URL to extract from
            extract_schema: Whether to extract schema.org data

        Returns:
            Page metadata
        """
        params = {
            "url": url,
            "extract_schema": extract_schema,
        }

        return self.request(
            method="GET",
            endpoint="/metadata",
            params=params,
            headers=self._get_headers(),
        )

    def execute_javascript(
        self,
        url: str,
        script: str,
        timeout: int = 30,
    ) -> Dict[str, Any]:
        """Execute JavaScript on page.

        Args:
            url: URL to execute on
            script: JavaScript code to execute
            timeout: Script timeout in seconds

        Returns:
            Script execution results
        """
        data = {
            "url": url,
            "script": script,
            "timeout": timeout,
        }

        return self.request(
            method="POST",
            endpoint="/execute",
            json=data,
            headers=self._get_headers(),
        )

    def capture_screenshot(
        self,
        url: str,
        width: int = 1280,
        height: int = 720,
        full_page: bool = False,
    ) -> bytes:
        """Capture page screenshot.

        Args:
            url: URL to capture
            width: Viewport width
            height: Viewport height
            full_page: Whether to capture full page

        Returns:
            Screenshot image bytes
        """
        params = {
            "url": url,
            "width": width,
            "height": height,
            "full_page": full_page,
        }

        response = self.http.request(
            method="GET",
            url=f"{self.base_url}/screenshot",
            params=params,
            headers=self._get_headers(),
        )

        return response.content

    def _get_timestamp(self) -> str:
        """Get current timestamp.

        Returns:
            ISO format timestamp
        """
        return datetime.now().isoformat()
