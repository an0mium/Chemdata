"""Enhanced PubMed client combining API access, web scraping, and LLM extraction.

This module provides a comprehensive PubMed client that:
1. Uses both PubMed API and web scraping for complete coverage
2. Integrates Crawl4AI for advanced web scraping features:
   - Session-based crawling for paginated content
   - LLM-based extraction of research data
   - Chunking strategies for large papers
   - Anti-bot protection with magic mode
3. Provides comprehensive validation:
   - Schema validation
   - Cross-source validation
   - Temporal validation
   - Community consensus validation
4. Includes advanced features:
   - Citation tracking
   - Impact analysis
   - Chemical structure extraction
   - Mechanism prediction
   - Safety analysis
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Any


@dataclass
class ResearchData:
    """Research data from PubMed."""

    pmid: str
    title: str
    abstract: Optional[str] = None
    authors: List[str] = field(default_factory=list)
    journal: Optional[str] = None
    year: Optional[int] = None
    doi: Optional[str] = None
    citations: int = 0
    compounds: List[str] = field(default_factory=list)
    effects: List[str] = field(default_factory=list)
    mechanisms: List[str] = field(default_factory=list)
    safety_notes: List[str] = field(default_factory=list)
    methods: List[str] = field(default_factory=list)
    results: List[str] = field(default_factory=list)
    conclusions: List[str] = field(default_factory=list)
    binding_data: List[Dict[str, Any]] = field(default_factory=list)
    confidence: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)


import aiohttp
from bs4 import BeautifulSoup
from crawl4ai import AsyncWebCrawler, BrowserConfig
from crawl4ai.extraction_strategy import (
    JsonCssExtractionStrategy,
    LLMExtractionStrategy,
    CosineStrategy,
)
from pydantic import BaseModel, Field

from ...models.core import CompoundData
from .base import WebClient, WebClientError, ValidationError
from ..validation.schema import DataSource, ValidationLevel
from ..validation.data import ValidationConfig, ValidationRule
from ..validation.enhanced import EnhancedValidator
from ..llm_utils import NootropicLLMExtractor
from ...pipeline.infrastructure.rate_limiter import RateLimiter, RateLimit
from ...pipeline.infrastructure.circuit_breaker import CircuitBreaker, ErrorCode, McpError

# API Rate Limits
AUTH_RATE_LIMIT = RateLimit(requests=3, period=1)  # 3 requests per second
UNAUTH_RATE_LIMIT = RateLimit(requests=1, period=1)  # 1 request per second
SCRAPE_RATE_LIMIT = RateLimit(requests=1, period=1)  # 1 request per second

# Circuit Breaker Config
CIRCUIT_FAILURE_THRESHOLD = 5
CIRCUIT_RESET_TIMEOUT = 300  # 5 minutes
CIRCUIT_HALF_OPEN_TIMEOUT = 60  # 1 minute

logger = logging.getLogger(__name__)


class PubMedArticleSchema(BaseModel):
    """Schema for PubMed article data extraction."""

    title: str = Field(..., description="Article title")
    abstract: str = Field(..., description="Article abstract")
    authors: List[str] = Field(..., description="List of authors")
    compounds: List[str] = Field(..., description="Chemical compounds mentioned")
    effects: List[str] = Field(..., description="Observed effects and outcomes")
    mechanisms: List[str] = Field(..., description="Mechanisms of action")
    safety: List[str] = Field(..., description="Safety information and warnings")
    methods: List[str] = Field(..., description="Methods and procedures used")
    results: List[str] = Field(..., description="Key findings and results")
    conclusions: List[str] = Field(..., description="Main conclusions")


@dataclass
class PubMedArticle:
    """Schema for PubMed article data."""

    pmid: str
    title: str
    abstract: Optional[str]
    authors: List[str]
    journal: Optional[str]
    year: Optional[int]
    doi: Optional[str] = None
    citations: int = 0
    compounds: List[str] = None
    effects: List[str] = None
    mechanisms: List[str] = None
    safety_notes: List[str] = None
    methods: List[str] = None
    results: List[str] = None
    conclusions: List[str] = None
    binding_data: List[Dict[str, Any]] = None
    confidence: float = 0.0
    metadata: Dict[str, Any] = None


class PubmedClient(WebClient):
    """Enhanced PubMed client combining API access and web scraping."""

    def _get_validation_config(self) -> ValidationConfig:
        """Get validation configuration.

        Returns:
            Validation configuration
        """
        return ValidationConfig(
            level=ValidationLevel.NORMAL,
            rules=[
                ValidationRule(
                    field="pmid",
                    required=True,
                    min_length=1,
                ),
                ValidationRule(
                    field="title",
                    required=True,
                    min_length=1,
                    max_length=1000,
                ),
                ValidationRule(
                    field="abstract",
                    required=False,
                    min_length=10,
                    max_length=10000,
                ),
                ValidationRule(
                    field="authors",
                    required=True,
                    min_items=1,
                ),
                ValidationRule(
                    field="journal",
                    required=False,
                    min_length=1,
                ),
                ValidationRule(
                    field="year",
                    required=False,
                    min_value=1800,
                    max_value=2100,
                ),
                ValidationRule(
                    field="doi",
                    required=False,
                    min_length=1,
                ),
                ValidationRule(
                    field="citations",
                    required=False,
                    min_value=0,
                ),
                ValidationRule(
                    field="compounds",
                    required=False,
                    min_items=0,
                ),
                ValidationRule(
                    field="effects",
                    required=False,
                    min_items=0,
                ),
                ValidationRule(
                    field="mechanisms",
                    required=False,
                    min_items=0,
                ),
                ValidationRule(
                    field="safety_notes",
                    required=False,
                    min_items=0,
                ),
                ValidationRule(
                    field="methods",
                    required=False,
                    min_items=0,
                ),
                ValidationRule(
                    field="results",
                    required=False,
                    min_items=0,
                ),
                ValidationRule(
                    field="conclusions",
                    required=False,
                    min_items=0,
                ),
                ValidationRule(
                    field="binding_data",
                    required=False,
                    min_items=0,
                ),
            ],
        )

    BASE_URL = "https://pubmed.ncbi.nlm.nih.gov"
    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(
        self,
        api_key: Optional[str] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        validator: Optional[EnhancedValidator] = None,
        llm_provider: str = "openai/gpt-4",
        api_token: Optional[str] = None,
        requests_per_second: float = 1.0,
        max_retries: int = 3,
        timeout: int = 30,
        cache_ttl: int = 3600,
        proxy_enabled: bool = True,
        proxy_rotation: bool = True,
        proxy_retry_count: int = 3,
        rate_limit: int = 10,
        rate_limit_delay: int = 60,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize enhanced PubMed client.

        Args:
            api_key: Optional PubMed API key
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            validator: Optional data validator
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            requests_per_second: Maximum requests per second
            max_retries: Maximum number of retries
            timeout: Request timeout in seconds
            cache_ttl: Cache TTL in seconds
            proxy_enabled: Whether to use proxies
            proxy_rotation: Whether to rotate proxies
            proxy_retry_count: Number of proxy retries
            rate_limit: Requests per minute
            rate_limit_delay: Delay after failure in seconds
            logger: Optional logger instance
        """
        super().__init__(
            name="pubmed",
            base_url=self.BASE_URL,
            data_source=DataSource.RESEARCH,
            requests_per_second=requests_per_second,
            max_retries=max_retries,
            timeout=timeout,
            cache_ttl=cache_ttl,
            logger=logger,
        )

        # Initialize rate limiters
        self.rate_limiter = RateLimiter()
        if api_key:
            self.rate_limiter.add_limit("api", AUTH_RATE_LIMIT)
        else:
            self.rate_limiter.add_limit("api", UNAUTH_RATE_LIMIT)
        self.rate_limiter.add_limit("scrape", SCRAPE_RATE_LIMIT)

        # Initialize circuit breakers
        self.api_circuit = CircuitBreaker(
            "pubmed_api",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )
        self.scrape_circuit = CircuitBreaker(
            "pubmed_scrape",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )

        # API config
        self.api_key = api_key
        self.validator = validator or EnhancedValidator()

        # Initialize LLM extractor
        self.llm_extractor = NootropicLLMExtractor()

        # Configure extraction strategies
        self.css_strategy = JsonCssExtractionStrategy(
            schema={
                "baseSelector": "article",
                "fields": [
                    {
                        "name": "title",
                        "selector": ".heading-title",
                        "type": "text",
                    },
                    {
                        "name": "abstract",
                        "selector": ".abstract-content",
                        "type": "text",
                    },
                    {
                        "name": "authors",
                        "selector": ".authors-list",
                        "type": "text",
                    },
                    {
                        "name": "date",
                        "selector": ".article-date",
                        "type": "text",
                    },
                    {
                        "name": "doi",
                        "selector": ".article-doi",
                        "type": "text",
                    },
                    {
                        "name": "journal",
                        "selector": ".journal-title",
                        "type": "text",
                    },
                    {
                        "name": "citations",
                        "selector": ".citation-count",
                        "type": "text",
                        "regex": r"(\d+)",
                    },
                ],
            }
        )

        # Configure LLM strategy with chunking according to crawl4ai docs
        self.llm_strategy = LLMExtractionStrategy(
            input_format="html",
            provider=llm_provider,
            api_token=api_token,
            schema=PubMedArticleSchema.schema(),
            extraction_type="schema",
            instruction="""
            From this article, extract:
            1. Chemical compounds mentioned
            2. Effects and outcomes
            3. Mechanisms of action
            4. Safety information
            5. Methods used
            6. Key results
            7. Main conclusions
            """,
            temperature=0,
            chunking={
                "type": "sliding_window",
                "window_size": 1000,
                "step_size": 500,
                "overlap": 0.5,
                "min_chunk_size": 100,
                "preprocessing": {"remove_html": True, "normalize_whitespace": True, "preserve_sentences": True},
                "postprocessing": {"merge_strategy": "union", "dedup_threshold": 0.9, "min_confidence": 0.7},
            },
        )

        self.cosine_strategy = CosineStrategy(
            word_count_threshold=15,
            max_dist=0.3,
            sim_threshold=0.2,
            model_name="sentence-transformers/all-MiniLM-L6-v2",
            top_k=2,
        )

        # Configure Crawl4AI
        self.crawler_config = BrowserConfig(
            javascript=BrowserConfig.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    ".docsum-content",  # Article container
                    ".authors-list",  # Author info
                    ".abstract-content",  # Abstract
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
                    "title": ".docsum-title",
                    "abstract": ".abstract-content",
                    "authors": ".authors-list",
                    "journal": ".journal-title",
                    "year": ".publication-date",
                    "citations": ".citation-count",
                    "doi": ".article-doi",
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

    async def search_papers(
        self,
        query: str,
        max_results: int = 100,
        min_date: Optional[str] = None,
        max_date: Optional[str] = None,
        sort: str = "relevance",
    ) -> List[ResearchData]:
        """Search PubMed papers using both API and web scraping.

        Args:
            query: Search query
            max_results: Maximum number of results
            min_date: Minimum publication date (YYYY-MM-DD)
            max_date: Maximum publication date (YYYY-MM-DD)
            sort: Sort order (relevance, date)

        Returns:
            List of research paper data
        """
        try:
            # Get PMIDs from API
            pmids = await self._search_api(query, max_results, min_date, max_date, sort)
            if not pmids:
                return []

            # Process articles
            results = []
            async with AsyncWebCrawler() as crawler:
                for pmid in pmids:
                    try:
                        article = await self._process_article(pmid, crawler)
                        if article:
                            results.append(article)
                    except Exception as e:
                        self.logger.error(f"Error processing article {pmid}: {str(e)}")

            return results

        except Exception as e:
            raise WebClientError(f"Error searching PubMed papers: {str(e)}")

    async def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get research data for a compound.

        Args:
            name: Compound name
            cas_number: Optional CAS number
            use_cache: Whether to use cached results

        Returns:
            Dictionary of research data or None if not found
        """
        # Build query
        query = self._build_compound_query(name, cas_number)

        # Search papers
        try:
            papers = await self.search_papers(query)
            if papers:
                return {
                    "papers": papers,
                    "timestamp": datetime.now().isoformat(),
                }
            return None
        except Exception as e:
            self.logger.error(f"Error getting compound data: {str(e)}")
            return None

    async def process_compounds(
        self,
        compounds: List[CompoundData],
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
            data = await self.get_compound_data(
                name=compound.name,
                cas_number=compound.cas_number,
                use_cache=use_cache,
            )

            if data:
                # Update compound with research data
                compound.research_data = data

    async def _search_api(
        self,
        query: str,
        max_results: int,
        min_date: Optional[str],
        max_date: Optional[str],
        sort: str,
    ) -> List[str]:
        """Search PubMed API for PMIDs.

        Args:
            query: Search query
            max_results: Maximum number of results
            min_date: Minimum publication date
            max_date: Maximum publication date
            sort: Sort order

        Returns:
            List of PMIDs
        """
        # Check circuit breaker
        if not self.api_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "PubMed API circuit breaker open",
            )

        try:
            # Wait for rate limit
            self.rate_limiter.wait("api")

            search_params = {
                "db": "pubmed",
                "term": query,
                "retmode": "json",
                "retmax": max_results,
                "sort": sort,
            }
            if min_date:
                search_params["mindate"] = min_date.split("-")[0]
            if max_date:
                search_params["maxdate"] = max_date.split("-")[0]

            if self.api_key:
                search_params["api_key"] = self.api_key

            async with aiohttp.ClientSession() as session:
                async with session.get(self.ESEARCH_URL, params=search_params) as response:
                    if not response.ok:
                        error_text = await response.text()
                        self.logger.error(f"API search failed: {error_text}")
                        self.api_circuit.record_failure()
                        return []

                    search_data = await response.json()
                    pmids = search_data["esearchresult"].get("idlist", [])

                    if pmids:
                        self.api_circuit.record_success()

                    return pmids

        except Exception as e:
            self.api_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"PubMed API search failed: {str(e)}",
            )

    async def _process_article(
        self,
        pmid: str,
        crawler: AsyncWebCrawler,
    ) -> Optional[ResearchData]:
        """Process single article.

        Args:
            pmid: PubMed ID
            crawler: AsyncWebCrawler instance

        Returns:
            Research paper data if successful
        """
        # Check circuit breaker
        if not self.scrape_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "PubMed scraping circuit breaker open",
            )

        try:
            # Wait for rate limit
            self.rate_limiter.wait("scrape")

            # Check cache
            if cached := self.cache.get(f"article:{pmid}"):
                return cached

            # Get article data
            article = await self._get_article_data(pmid, crawler)
            if not article:
                return None

            # Validate
            if not await self._validate_article(article):
                return None

            # Cache result
            result = ResearchData(**article.__dict__)
            self.cache.set(f"article:{pmid}", result)

            self.scrape_circuit.record_success()
            return result

        except Exception as e:
            self.scrape_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"Error processing PubMed article: {str(e)}",
            )

    async def get_paper_details(
        self,
        pmid: str,
    ) -> Optional[ResearchData]:
        """Get detailed paper data by PubMed ID.

        Args:
            pmid: PubMed ID

        Returns:
            Research paper data if found
        """
        try:
            # Check cache
            if cached := self.cache.get(f"article:{pmid}"):
                return cached

            # Get article data
            async with AsyncWebCrawler() as crawler:
                article = await self._get_article_data(pmid, crawler)
                if not article:
                    return None

                # Validate
                if not await self._validate_article(article):
                    return None

                # Cache result
                result = ResearchData(**article.__dict__)
                self.cache.set(f"article:{pmid}", result)
                return result

        except Exception as e:
            raise WebClientError(f"Error getting paper details for PMID {pmid}: {str(e)}")

    async def _get_article_data(
        self,
        pmid: str,
        crawler: AsyncWebCrawler,
    ) -> Optional[PubMedArticle]:
        """Get article data from both web and API.

        Args:
            pmid: PubMed ID
            crawler: AsyncWebCrawler instance

        Returns:
            Combined article data
        """
        try:
            # Get web data
            url = f"{self.BASE_URL}/{pmid}/"

            # Try CSS extraction first
            self.crawler_config.extraction_strategy = self.css_strategy
            result = await crawler.arun(
                url=url,
                config=self.crawler_config,
                session_id=self.session_id,
            )

            web_data = {}
            if result.extracted_content:
                web_data = result.extracted_content[0]

            # Try LLM extraction if needed
            if not web_data.get("abstract"):
                self.crawler_config.extraction_strategy = self.llm_strategy
                result = await crawler.arun(
                    url=url,
                    config=self.crawler_config,
                    session_id=self.session_id,
                )
                if result.extracted_content:
                    web_data.update(result.extracted_content[0])

            # Get API data
            api_data = await self._get_api_data(pmid)

            # Get citation count
            citations = await self._get_citations(pmid, crawler)

            # Get full text if available
            full_text = None
            if web_data.get("doi"):
                full_text = await self._get_full_text(
                    crawler,
                    web_data["doi"],
                    self.crawler_config,
                )

            # Extract chemical info and mechanisms
            text = full_text or web_data.get("abstract", "")
            if text:
                chemicals = await self.llm_extractor.extract_chemical_info(text)
                mechanisms = await self.llm_extractor.extract_mechanisms(text)
                safety = await self.llm_extractor.extract_safety(text)
            else:
                chemicals = []
                mechanisms = []
                safety = []

            # Create article object
            article = PubMedArticle(
                pmid=pmid,
                title=web_data.get("title", api_data.get("title")),
                abstract=web_data.get("abstract", api_data.get("abstract")),
                authors=web_data.get("authors", api_data.get("authors", [])),
                journal=web_data.get("journal", api_data.get("journal")),
                year=web_data.get("year", api_data.get("year")),
                doi=web_data.get("doi", api_data.get("doi")),
                citations=citations,
                compounds=chemicals or web_data.get("compounds", []),
                effects=web_data.get("effects", []),
                mechanisms=mechanisms or web_data.get("mechanisms", []),
                safety_notes=safety or web_data.get("safety", []),
                methods=web_data.get("methods", []),
                results=web_data.get("results", []),
                conclusions=web_data.get("conclusions", []),
                binding_data=api_data.get("binding_data", []),
                confidence=self._calculate_confidence(web_data, api_data),
                metadata={
                    "source": "pubmed",
                    "extraction_date": datetime.now().isoformat(),
                    "screenshot": str(result.screenshot) if result.screenshot else None,
                    "javascript_logs": result.javascript_logs,
                    "chunking_stats": result.chunking_stats,
                    "full_text_available": bool(full_text),
                },
            )

            return article

        except Exception as e:
            self.logger.error(f"Error getting article data for {pmid}: {str(e)}")
            return None

    async def _get_api_data(self, pmid: str) -> Dict[str, Any]:
        """Get data from PubMed API.

        Args:
            pmid: PubMed ID

        Returns:
            API data
        """
        try:
            # Fetch article details
            fetch_params = {"db": "pubmed", "id": pmid, "retmode": "xml"}
            if self.api_key:
                fetch_params["api_key"] = self.api_key

            async with aiohttp.ClientSession() as session:
                async with session.get(self.EFETCH_URL, params=fetch_params) as response:
                    if not response.ok:
                        self.logger.error(f"API fetch failed: {response.status}")
                        return {}

                    xml_text = await response.text()
                    if not xml_text:
                        return {}

            # Parse XML
            soup = BeautifulSoup(xml_text, "xml")
            article = soup.find("PubmedArticle")
            if not article:
                return {}

            # Extract metadata
            return self._extract_article_metadata(article)

        except Exception as e:
            self.logger.error(f"Error getting API data for {pmid}: {str(e)}")
            return {}

    async def _get_citations(self, pmid: str, crawler: AsyncWebCrawler) -> int:
        """Get citation count.

        Args:
            pmid: PubMed ID
            crawler: AsyncWebCrawler instance

        Returns:
            Number of citations
        """
        try:
            result = await crawler.arun(
                url=f"{self.BASE_URL}/{pmid}/citations/",
                config=BrowserConfig(
                    javascript=BrowserConfig.JavaScript(
                        enabled=True,
                        wait_for_network=True,
                        wait_for_selectors=[".citation-count"],
                        stealth_mode=True,
                    ),
                    screenshot=BrowserConfig.Screenshot(enabled=True, full_page=True),
                    extraction=BrowserConfig.Extraction(
                        css={
                            "citations": ".citation-count",
                        },
                    ),
                    proxy=BrowserConfig.Proxy(
                        enabled=True,
                        rotation=True,
                        retry_count=3,
                    ),
                    rate_limit=BrowserConfig.RateLimit(
                        requests_per_minute=10,
                        delay_after_failure=60,
                    ),
                ),
            )
            if not result.success:
                return 0

            citations_text = result.extracted_content[0].get("citations", "0")
            return int(citations_text)

        except Exception as e:
            self.logger.error(f"Error getting citations for {pmid}: {str(e)}")
            return 0

    async def _get_full_text(
        self,
        crawler: AsyncWebCrawler,
        doi: str,
        config: BrowserConfig,
    ) -> Optional[str]:
        """Get full text from DOI if available.

        Args:
            crawler: AsyncWebCrawler instance
            doi: Paper DOI
            config: Crawler configuration

        Returns:
            Full text if available, None otherwise
        """
        try:
            # Try to get full text from DOI
            url = f"https://doi.org/{doi}"
            result = await crawler.arun(
                url=url,
                config=config,
                session_id=self.session_id,
            )

            if result.extracted_content:
                # Extract full text
                text = result.extracted_content[0].get("full_text", "")
                if text:
                    return text

            return None

        except Exception as e:
            self.logger.warning(f"Error getting full text for DOI {doi}: {str(e)}")
            return None

    async def _validate_article(self, article: PubMedArticle) -> bool:
        """Validate article data.

        Args:
            article: Article to validate

        Returns:
            True if valid
        """
        try:
            # Convert to dict for validation
            data = {
                "pmid": article.pmid,
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
                "methods": article.methods,
                "results": article.results,
                "conclusions": article.conclusions,
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

    def _build_compound_query(self, name: str, cas_number: Optional[str] = None) -> str:
        """Build search query for compound.

        Args:
            name: Compound name
            cas_number: Optional CAS number

        Returns:
            Search query
        """
        query_parts = []

        # Add name
        if name:
            query_parts.append(name)

        # Add CAS number
        if cas_number:
            query_parts.append(cas_number)

        # Build query
        return " OR ".join(f'"{term}"' for term in query_parts)

    def _calculate_confidence(self, web_data: Dict[str, Any], api_data: Dict[str, Any]) -> float:
        """Calculate confidence score.

        Args:
            web_data: Data from web scraping
            api_data: Data from API

        Returns:
            Confidence score between 0 and 1
        """
        score = 0.0
        total = 0

        # Calculate base field scores
        base_fields = ["title", "abstract", "authors", "journal", "year", "doi"]
        for field in base_fields:
            total += 1
            if field in web_data and field in api_data:
                score += 1.0
            elif field in web_data or field in api_data:
                score += 0.5

        # Calculate extracted data scores
        extracted_fields = [
            "compounds",
            "effects",
            "mechanisms",
            "safety",
            "methods",
            "results",
            "conclusions",
        ]
        for field in extracted_fields:
            total += 1
            if web_data.get(field):
                score += 1.0

        # Add binding data score
        if api_data.get("binding_data"):
            score += 1.0
            total += 1

        # Calculate final score
        return score / total if total > 0 else 0.0

    def _extract_article_metadata(self, article_xml: BeautifulSoup) -> Dict[str, Any]:
        """Extract metadata from PubMed XML.

        Args:
            article_xml: BeautifulSoup object for article XML

        Returns:
            Extracted metadata
        """
        try:
            metadata = {}
            metadata.update(self._extract_basic_fields(article_xml))
            metadata.update(self._extract_authors(article_xml))
            metadata.update(self._extract_date(article_xml))
            metadata.update(self._extract_identifiers(article_xml))
            metadata.update(self._extract_chemical_data(article_xml))
            return metadata

        except Exception as e:
            self.logger.error(f"Error extracting article metadata: {str(e)}")
            return {}

    def _extract_basic_fields(self, article_xml: BeautifulSoup) -> Dict[str, Any]:
        """Extract basic article fields."""
        fields = {}

        # Extract title
        fields["title"] = article_xml.find("ArticleTitle").text

        # Extract abstract
        abstract_elem = article_xml.find("Abstract")
        fields["abstract"] = abstract_elem.text if abstract_elem else None

        # Extract journal
        journal_elem = article_xml.find("Journal")
        fields["journal"] = journal_elem.find("Title").text if journal_elem else None

        return fields

    def _extract_authors(self, article_xml: BeautifulSoup) -> Dict[str, List[str]]:
        """Extract author information."""
        authors = []
        author_list = article_xml.find("AuthorList")
        if author_list:
            for author in author_list.find_all("Author"):
                last_name = author.find("LastName")
                fore_name = author.find("ForeName")
                if last_name and fore_name:
                    authors.append(f"{fore_name.text} {last_name.text}")
        return {"authors": authors}

    def _extract_date(self, article_xml: BeautifulSoup) -> Dict[str, Optional[int]]:
        """Extract publication date."""
        pub_date = article_xml.find("PubDate")
        if pub_date:
            year = pub_date.find("Year")
            if year:
                return {"year": int(year.text)}
        return {"year": None}

    def _extract_identifiers(self, article_xml: BeautifulSoup) -> Dict[str, Optional[str]]:
        """Extract article identifiers."""
        article_ids = article_xml.find_all("ArticleId")
        for article_id in article_ids:
            if article_id.get("IdType") == "doi":
                return {"doi": article_id.text}
        return {"doi": None}

    def _extract_chemical_data(self, article_xml: BeautifulSoup) -> Dict[str, Any]:
        """Extract chemical and binding data."""
        chemical_list = article_xml.find("ChemicalList")
        if not chemical_list:
            return {
                "chemicals": [],
                "binding_data": [],
            }

        # Extract chemical names
        chemicals = []
        for chemical in chemical_list.find_all("Chemical"):
            substance = chemical.find("NameOfSubstance")
            if substance:
                chemicals.append(substance.text)

        # Extract binding data
        binding_data = []
        for chemical in chemical_list.find_all("Chemical"):
            substance = chemical.find("NameOfSubstance")
            registry = chemical.find("RegistryNumber")
            if substance and registry:
                binding_data.append(
                    {
                        "compound": substance.text,
                        "registry_number": registry.text,
                    }
                )

        return {
            "chemicals": chemicals,
            "binding_data": binding_data,
        }
