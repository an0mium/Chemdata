"""Enhanced patent integration combining Google Patents, USPTO and Espacenet.

This module provides a comprehensive patent client that:
1. Uses Crawl4AI for web scraping Google Patents pages
2. Integrates with USPTO API for additional data
3. Integrates with Espacenet API for European patents and INPADOC data
4. Handles patent families and citations with relevance scoring
5. Extracts chemical compounds and structures
6. Analyzes patent classifications
7. Tracks legal status and lifecycle
"""

from typing import Dict, List, Optional, Any, Tuple, Set
from dataclasses import dataclass
from datetime import datetime
import logging
import re

from crawl4ai import AsyncWebCrawler, Config

from ..base_client import BaseWebClient
from ..validation.schema import BaseSchema


@dataclass
class PatentCitation:
    """Patent citation data."""

    patent_number: str
    title: str
    filing_date: Optional[str]
    relevance: Optional[float]  # Relevance score 0-1
    citation_type: str  # "forward" or "backward"
    metadata: Dict[str, Any]


@dataclass
class PatentFamily:
    """Patent family data."""

    family_id: str
    members: List[str]  # Patent numbers
    priority_date: Optional[str]
    countries: Set[str]
    metadata: Dict[str, Any]


@dataclass
class PatentClassification:
    """Patent classification data."""

    system: str  # "IPC", "CPC", "USPC"
    code: str
    description: str
    level: str  # "section", "class", "subclass", "group"
    relevance: Optional[float]  # Relevance score 0-1


@dataclass
class PatentLegalStatus:
    """Patent legal status data."""

    status: str  # e.g. "granted", "withdrawn", "expired"
    date: str
    country: str
    description: Optional[str]
    metadata: Dict[str, Any]


@dataclass
class PatentData(BaseSchema):
    """Enhanced schema for patent data."""

    title: str
    abstract: Optional[str]
    inventors: List[str]
    assignee: Optional[str]
    filing_date: Optional[str]
    publication_date: Optional[str]
    patent_number: Optional[str]
    compounds: List[str]
    effects: List[str]
    mechanisms: List[str]
    safety_notes: List[str]
    family: Optional[PatentFamily]
    citations: List[PatentCitation]
    classifications: List[PatentClassification]
    legal_status: Optional[PatentLegalStatus]
    confidence: float
    metadata: Dict[str, Any]


class EnhancedPatentClient(BaseWebClient):
    """Enhanced patent client combining Google Patents, USPTO and Espacenet."""

    GOOGLE_PATENTS_URL = "https://patents.google.com"
    USPTO_API_URL = "https://developer.uspto.gov/ibd-api/v1"
    ESPACENET_API_URL = "https://ops.epo.org/3.2/rest-services"

    def __init__(
        self,
        llm_provider: str = "ollama/llama2",
        api_token: Optional[str] = None,
        uspto_api_key: Optional[str] = None,
        espacenet_api_key: Optional[str] = None,
        **kwargs,
    ):
        """Initialize client.

        Args:
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            uspto_api_key: Optional USPTO API key
            espacenet_api_key: Optional Espacenet API key
            **kwargs: Additional arguments passed to BaseWebClient
        """
        super().__init__(**kwargs)
        self.llm_provider = llm_provider
        self.api_token = api_token
        self.uspto_api_key = uspto_api_key
        self.espacenet_api_key = espacenet_api_key
        self.logger = logging.getLogger(__name__)

        # Configure Crawl4AI with anti-bot detection avoidance
        self.config = Config(
            javascript=Config.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    ".patent-result",
                    ".description",
                    ".claims",
                    ".patent-citations",
                    ".patent-family",
                    ".patent-classifications",
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
                    "title": ".patent-title",
                    "abstract": ".abstract",
                    "inventors": ".inventors-line a",
                    "assignee": ".assignee-line a",
                    "filing_date": ".filing-date",
                    "publication_date": ".publication-date",
                    "patent_number": ".patent-number",
                    "citations": ".patent-citation",
                    "family": ".patent-family",
                    "classifications": ".patent-classification",
                },
            ),
            proxy=Config.Proxy(
                enabled=True,
                rotation=True,
                retry_count=3,
            ),
            rate_limit=Config.RateLimit(
                requests_per_minute=10,
                delay_after_failure=60,
            ),
        )

    def _build_search_query(
        self,
        query: str,
        structure: Optional[str] = None,
        date_range: Optional[Tuple[str, str]] = None,
        classification: Optional[str] = None,
    ) -> str:
        """Build search query string.

        Args:
            query: Base search query
            structure: Optional chemical structure
            date_range: Optional date range
            classification: Optional classification code

        Returns:
            Complete search query
        """
        search_query = query
        if structure:
            search_query += f" structure:({structure})"
        if date_range:
            start_date, end_date = date_range
            search_query += f" after:{start_date} before:{end_date}"
        if classification:
            search_query += f" classification:({classification})"
        return search_query

    async def _process_patent(
        self,
        patent: Dict[str, Any],
        include_family: bool = True,
        include_citations: bool = True,
    ) -> Optional[PatentData]:
        """Process single patent data.

        Args:
            patent: Raw patent data
            include_family: Whether to fetch family data
            include_citations: Whether to fetch citation data

        Returns:
            Processed patent data
        """
        try:
            # Get patent details
            patent_data = await self._get_patent_details(
                patent["patent_number"],
                include_family=include_family,
                include_citations=include_citations,
            )
            return patent_data

        except Exception as e:
            self.logger.error(f"Error processing patent: {str(e)}")
            return None

    async def search_patents(
        self,
        query: str,
        structure: Optional[str] = None,
        date_range: Optional[Tuple[str, str]] = None,
        classification: Optional[str] = None,
        max_results: int = 100,
        include_family: bool = True,
        include_citations: bool = True,
    ) -> List[PatentData]:
        """Search patents with enhanced data.

        Args:
            query: Search query
            structure: Optional chemical structure (SMILES/InChI)
            date_range: Optional (start_date, end_date) tuple in YYYY-MM-DD format
            classification: Optional classification code (IPC/CPC/USPC)
            max_results: Maximum number of results
            include_family: Whether to fetch family data
            include_citations: Whether to fetch citation data

        Returns:
            List of patent data
        """
        results = []
        start = 0
        results_per_page = 10

        # Build search query
        search_query = self._build_search_query(query, structure, date_range, classification)

        try:
            async with AsyncWebCrawler() as crawler:
                while len(results) < max_results:
                    # Construct search URL with pagination
                    url = f"{self.GOOGLE_PATENTS_URL}/search?q={search_query}&page={start//10}"

                    # Scrape with Crawl4AI
                    result = await crawler.arun(urls=[url], config=self.config)

                    if not result.success:
                        self.logger.error(f"Failed to scrape page {start//10 + 1}")
                        break

                    # Extract data using both CSS and LLM
                    content = result.extracted_content[0]
                    patents = self._parse_search_results(content)

                    if not patents:
                        break

                    # Process each patent
                    for patent in patents:
                        patent_data = await self._process_patent(
                            patent,
                            include_family=include_family,
                            include_citations=include_citations,
                        )
                        if patent_data:
                            results.append(patent_data)
                            if len(results) >= max_results:
                                break

                    start += results_per_page

        except Exception as e:
            self.logger.error(f"Error searching patents: {str(e)}")

        return results

    async def get_patent_family(self, patent_number: str) -> Optional[PatentFamily]:
        """Get patent family information.

        Args:
            patent_number: Patent number

        Returns:
            Patent family data if found
        """
        try:
            # Try Espacenet first for most comprehensive family data
            if self.espacenet_api_key:
                family = await self._get_espacenet_family(patent_number)
                if family:
                    return family

            # Try Google Patents next
            family = await self._get_google_patents_family(patent_number)
            if family:
                return family

            # Fall back to USPTO API
            if self.uspto_api_key:
                family = await self._get_uspto_family(patent_number)
                if family:
                    return family

        except Exception as e:
            self.logger.error(f"Error getting patent family: {str(e)}")

        return None

    async def get_patent_citations(
        self,
        patent_number: str,
        max_depth: int = 2,
    ) -> List[PatentCitation]:
        """Get patent citation network.

        Args:
            patent_number: Patent number
            max_depth: Maximum citation depth to traverse

        Returns:
            List of patent citations
        """
        citations = []
        seen = set()
        queue = [(patent_number, 0)]

        try:
            while queue and len(citations) < 1000:  # Limit total citations
                current_patent, depth = queue.pop(0)
                if current_patent in seen or depth > max_depth:
                    continue

                seen.add(current_patent)

                # Get citations from all sources
                current_citations = []

                # Get Espacenet citations
                if self.espacenet_api_key:
                    espacenet_citations = await self._get_espacenet_citations(current_patent)
                    if espacenet_citations:
                        current_citations.extend(espacenet_citations)

                # Get Google Patents citations
                google_citations = await self._get_google_patents_citations(current_patent)
                if google_citations:
                    current_citations.extend(google_citations)

                # Get USPTO citations
                if self.uspto_api_key:
                    uspto_citations = await self._get_uspto_citations(current_patent)
                    if uspto_citations:
                        current_citations.extend(uspto_citations)

                # Merge citations by patent number
                by_number = {c.patent_number: c for c in citations}
                for citation in current_citations:
                    existing = by_number.get(citation.patent_number)
                    if existing:
                        # Prefer citation with higher relevance score
                        if citation.relevance and (not existing.relevance or citation.relevance > existing.relevance):
                            by_number[citation.patent_number] = citation
                    else:
                        by_number[citation.patent_number] = citation

                citations = list(by_number.values())

                # Add connected patents to queue
                if depth < max_depth:
                    for citation in current_citations:
                        queue.append((citation.patent_number, depth + 1))

        except Exception as e:
            self.logger.error(f"Error getting citation network: {str(e)}")

        return citations

    async def _get_patent_details(
        self,
        patent_number: str,
        include_family: bool = True,
        include_citations: bool = True,
    ) -> Optional[PatentData]:
        """Get detailed patent data.

        Args:
            patent_number: Patent number
            include_family: Whether to fetch family data
            include_citations: Whether to fetch citation data

        Returns:
            Patent data if found
        """
        try:
            # Get basic patent data
            url = f"{self.GOOGLE_PATENTS_URL}/patent/{patent_number}"
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.config)
                if not result.success:
                    return None

                content = result.extracted_content[0]
                patent = self._parse_patent_details(content)

                # Add family data if requested
                if include_family:
                    family = await self.get_patent_family(patent_number)
                    if family:
                        patent["family"] = family

                # Add citation data if requested
                if include_citations:
                    citations = await self.get_patent_citations(patent_number)
                    if citations:
                        patent["citations"] = citations

                # Get legal status from Espacenet
                legal_status = None
                if self.espacenet_api_key:
                    legal_status = await self._get_espacenet_legal_status(patent_number)

                # Get USPTO data if available
                uspto_data = None
                if self.uspto_api_key:
                    uspto_data = await self._get_uspto_data(patent_number)
                    if uspto_data:
                        patent.update(uspto_data)

                # Create patent data object
                return PatentData(
                    title=patent.get("title", ""),
                    abstract=patent.get("abstract"),
                    inventors=patent.get("inventors", []),
                    assignee=patent.get("assignee"),
                    filing_date=patent.get("filing_date"),
                    publication_date=patent.get("publication_date"),
                    patent_number=patent_number,
                    compounds=patent.get("compounds", []),
                    effects=patent.get("effects", []),
                    mechanisms=patent.get("mechanisms", []),
                    safety_notes=patent.get("safety", []),
                    family=patent.get("family"),
                    citations=patent.get("citations", []),
                    classifications=patent.get("classifications", []),
                    legal_status=legal_status,
                    confidence=self._calculate_confidence(patent),
                    metadata={
                        "sources": [
                            "google_patents",
                            *(["espacenet"] if self.espacenet_api_key else []),
                            *(["uspto"] if uspto_data else []),
                        ],
                        "scrape_date": datetime.now().isoformat(),
                        "url": url,
                        "screenshot": str(result.screenshot),
                        "javascript_logs": result.javascript_logs,
                    },
                )

        except Exception as e:
            self.logger.error(f"Error getting patent details: {str(e)}")
            return None

    async def _get_espacenet_family(self, patent_number: str) -> Optional[PatentFamily]:
        """Get patent family from Espacenet API.

        Args:
            patent_number: Patent number

        Returns:
            Patent family data if found
        """
        if not self.espacenet_api_key:
            return None

        try:
            # Call Espacenet API for INPADOC family data
            url = f"{self.ESPACENET_API_URL}/family/publication/{patent_number}"
            headers = {
                "Authorization": f"Bearer {self.espacenet_api_key}",
                "Accept": "application/json",
            }
            async with self.session.get(url, headers=headers) as response:
                if response.status != 200:
                    return None

                data = await response.json()
                family_data = data.get("ops:world-patent-data", {}).get("ops:patent-family", {})

                # Extract family members
                members = []
                countries = set()
                priority_date = None

                for member in family_data.get("family-member", []):
                    pub_ref = member.get("publication-reference", {})
                    doc_id = pub_ref.get("document-id", {})

                    # Get publication number
                    pub_number = doc_id.get("doc-number")
                    if pub_number:
                        members.append(pub_number)

                    # Get country
                    country = doc_id.get("country")
                    if country:
                        countries.add(country)

                    # Get earliest priority date
                    if not priority_date:
                        priority_date = doc_id.get("date")

                return PatentFamily(
                    family_id=family_data.get("@family-id", ""),
                    members=members,
                    priority_date=priority_date,
                    countries=countries,
                    metadata={
                        "source": "espacenet",
                        "extraction_date": datetime.now().isoformat(),
                    },
                )

        except Exception as e:
            self.logger.error(f"Error getting Espacenet family: {str(e)}")
            return None

    async def _get_google_patents_family(self, patent_number: str) -> Optional[PatentFamily]:
        """Get patent family from Google Patents.

        Args:
            patent_number: Patent number

        Returns:
            Patent family data if found
        """
        try:
            url = f"{self.GOOGLE_PATENTS_URL}/patent/{patent_number}/family"
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.config)
                if not result.success:
                    return None

                content = result.extracted_content[0]
                family_data = self._parse_family_data(content)

                if family_data:
                    return PatentFamily(
                        family_id=family_data["family_id"],
                        members=family_data["members"],
                        priority_date=family_data.get("priority_date"),
                        countries=set(family_data.get("countries", [])),
                        metadata={
                            "source": "google_patents",
                            "extraction_date": datetime.now().isoformat(),
                        },
                    )

        except Exception as e:
            self.logger.error(f"Error getting Google Patents family: {str(e)}")

        return None

    async def _get_uspto_family(self, patent_number: str) -> Optional[PatentFamily]:
        """Get patent family from USPTO API.

        Args:
            patent_number: Patent number

        Returns:
            Patent family data if found
        """
        if not self.uspto_api_key:
            return None

        try:
            # Call USPTO API
            url = f"{self.USPTO_API_URL}/patents/{patent_number}/family"
            headers = {"X-Api-Key": self.uspto_api_key}
            async with self.session.get(url, headers=headers) as response:
                if response.status != 200:
                    return None

                data = await response.json()
                return PatentFamily(
                    family_id=data["family_id"],
                    members=data["members"],
                    priority_date=data.get("priority_date"),
                    countries=set(data.get("countries", [])),
                    metadata={
                        "source": "uspto",
                        "extraction_date": datetime.now().isoformat(),
                    },
                )

        except Exception as e:
            self.logger.error(f"Error getting USPTO family: {str(e)}")
            return None

    async def _get_espacenet_citations(self, patent_number: str) -> List[PatentCitation]:
        """Get citations from Espacenet API.

        Args:
            patent_number: Patent number

        Returns:
            List of patent citations
        """
        if not self.espacenet_api_key:
            return []

        citations = []

        try:
            # Call Espacenet API for citation data
            url = f"{self.ESPACENET_API_URL}/published-data/publication/{patent_number}/citations"
            headers = {
                "Authorization": f"Bearer {self.espacenet_api_key}",
                "Accept": "application/json",
            }
            async with self.session.get(url, headers=headers) as response:
                if response.status != 200:
                    return citations

                data = await response.json()
                citations_data = data.get("ops:world-patent-data", {}).get("ops:citation-list", {})

                for citation in citations_data.get("citation", []):
                    pat_citation = citation.get("patcit", {})
                    doc_id = pat_citation.get("document-id", {})

                    citations.append(
                        PatentCitation(
                            patent_number=doc_id.get("doc-number", ""),
                            title=citation.get("passage", {}).get("text", ""),
                            filing_date=doc_id.get("date"),
                            relevance=self._parse_citation_relevance(citation),
                            citation_type=citation.get("@cited-phase", ""),
                            metadata={
                                "source": "espacenet",
                                "extraction_date": datetime.now().isoformat(),
                                "category": citation.get("@cited-category"),
                            },
                        )
                    )

        except Exception as e:
            self.logger.error(f"Error getting Espacenet citations: {str(e)}")

        return citations

    async def _get_google_patents_citations(self, patent_number: str) -> List[PatentCitation]:
        """Get citations from Google Patents.

        Args:
            patent_number: Patent number

        Returns:
            List of patent citations
        """
        citations = []

        try:
            url = f"{self.GOOGLE_PATENTS_URL}/patent/{patent_number}/citations"
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.config)
                if not result.success:
                    return citations

                content = result.extracted_content[0]
                citations_data = self._parse_citations_data(content)

                for citation in citations_data:
                    citations.append(
                        PatentCitation(
                            patent_number=citation["patent_number"],
                            title=citation["title"],
                            filing_date=citation.get("filing_date"),
                            relevance=citation.get("relevance"),
                            citation_type=citation["type"],
                            metadata={
                                "source": "google_patents",
                                "extraction_date": datetime.now().isoformat(),
                            },
                        )
                    )

        except Exception as e:
            self.logger.error(f"Error getting Google Patents citations: {str(e)}")

        return citations

    async def _get_uspto_citations(self, patent_number: str) -> List[PatentCitation]:
        """Get citations from USPTO API.

        Args:
            patent_number: Patent number

        Returns:
            List of patent citations
        """
        if not self.uspto_api_key:
            return []

        citations = []

        try:
            # Call USPTO API
            url = f"{self.USPTO_API_URL}/patents/{patent_number}/citations"
            headers = {"X-Api-Key": self.uspto_api_key}
            async with self.session.get(url, headers=headers) as response:
                if response.status != 200:
                    return citations

                data = await response.json()
                for citation in data["citations"]:
                    citations.append(
                        PatentCitation(
                            patent_number=citation["patent_number"],
                            title=citation["title"],
                            filing_date=citation.get("filing_date"),
                            relevance=citation.get("relevance"),
                            citation_type=citation["type"],
                            metadata={
                                "source": "uspto",
                                "extraction_date": datetime.now().isoformat(),
                            },
                        )
                    )

        except Exception as e:
            self.logger.error(f"Error getting USPTO citations: {str(e)}")

        return citations

    async def _get_uspto_data(self, patent_number: str) -> Optional[Dict[str, Any]]:
        """Get additional data from USPTO API.

        Args:
            patent_number: Patent number

        Returns:
            USPTO data if available
        """
        if not self.uspto_api_key:
            return None

        try:
            # Call USPTO API
            url = f"{self.USPTO_API_URL}/patents/{patent_number}"
            headers = {"X-Api-Key": self.uspto_api_key}
            async with self.session.get(url, headers=headers) as response:
                if response.status != 200:
                    return None

                return await response.json()

        except Exception as e:
            self.logger.error(f"Error getting USPTO data: {str(e)}")
            return None

    async def _get_espacenet_legal_status(self, patent_number: str) -> Optional[PatentLegalStatus]:
        """Get legal status from Espacenet API.

        Args:
            patent_number: Patent number

        Returns:
            Legal status data if found
        """
        if not self.espacenet_api_key:
            return None

        try:
            url = f"{self.ESPACENET_API_URL}/legal/publication/{patent_number}/status"
            headers = {
                "Authorization": f"Bearer {self.espacenet_api_key}",
                "Accept": "application/json",
            }
            async with self.session.get(url, headers=headers) as response:
                if response.status != 200:
                    return None

                data = await response.json()
                status_data = data.get("ops:world-patent-data", {}).get("ops:legal-status", {})

                if not status_data:
                    return None

                latest_status = status_data[-1] if isinstance(status_data, list) else status_data

                return PatentLegalStatus(
                    status=latest_status.get("@status-code", ""),
                    date=latest_status.get("date", ""),
                    country=latest_status.get("country", ""),
                    description=latest_status.get("text", ""),
                    metadata={
                        "source": "espacenet",
                        "extraction_date": datetime.now().isoformat(),
                        "raw_status": latest_status,
                    },
                )

        except Exception as e:
            self.logger.error(f"Error getting Espacenet legal status: {str(e)}")
            return None

    def _parse_citation_relevance(self, citation: Dict[str, Any]) -> Optional[float]:
        """Parse citation relevance from Espacenet data.

        Args:
            citation: Citation data from Espacenet

        Returns:
            Relevance score between 0 and 1
        """
        # Espacenet uses categories like "X", "Y", "A" to indicate relevance
        category = citation.get("@cited-category", "")

        # Map categories to scores
        category_scores = {
            "X": 1.0,  # Particularly relevant if taken alone
            "Y": 0.8,  # Particularly relevant if combined with other documents
            "A": 0.5,  # Background art
            "O": 0.3,  # Non-written disclosure
            "P": 0.4,  # Intermediate document
            "T": 0.6,  # Theory/principle underlying the invention
            "E": 0.7,  # Earlier patent document
            "D": 0.4,  # Document cited in application
        }

        return category_scores.get(category)

    def _parse_search_results(self, content: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse search results from extracted content.

        Args:
            content: Extracted content from Crawl4AI

        Returns:
            List of parsed patent data
        """
        patents = []

        # Get patent elements
        titles = content.get("title", [])
        abstracts = content.get("abstract", [])
        inventors_list = content.get("inventors", [])
        assignees = content.get("assignee", [])
        filing_dates = content.get("filing_date", [])
        publication_dates = content.get("publication_date", [])
        patent_numbers = content.get("patent_number", [])

        # Get LLM-extracted data
        compounds = content.get("compounds", [])
        effects = content.get("effects", [])
        mechanisms = content.get("mechanisms", [])
        safety = content.get("safety", [])

        # Get classification data
        classifications = self._parse_classifications(content.get("classifications", []))

        # Combine data for each patent
        for i in range(len(titles)):
            patent = {
                "title": titles[i] if i < len(titles) else None,
                "abstract": abstracts[i] if i < len(abstracts) else None,
                "inventors": inventors_list[i] if i < len(inventors_list) else [],
                "assignee": assignees[i] if i < len(assignees) else None,
                "filing_date": filing_dates[i] if i < len(filing_dates) else None,
                "publication_date": publication_dates[i] if i < len(publication_dates) else None,
                "patent_number": patent_numbers[i] if i < len(patent_numbers) else None,
                "url": self._build_patent_url(patent_numbers[i]) if i < len(patent_numbers) else None,
                "compounds": compounds,
                "effects": effects,
                "mechanisms": mechanisms,
                "safety": safety,
                "classifications": classifications[i] if i < len(classifications) else [],
            }
            patents.append(patent)

        return patents

    def _parse_patent_details(self, content: Dict[str, Any]) -> Dict[str, Any]:
        """Parse detailed patent data.

        Args:
            content: Extracted content

        Returns:
            Parsed patent data
        """
        # Extract basic fields
        patent = self._parse_search_results([content])[0]

        # Add classification data
        patent["classifications"] = self._parse_classifications(content.get("classifications", []))

        return patent

    def _parse_family_data(self, content: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Parse patent family data.

        Args:
            content: Extracted content

        Returns:
            Parsed family data
        """
        try:
            family_id = content.get("family_id", [None])[0]
            if not family_id:
                return None

            return {
                "family_id": family_id,
                "members": content.get("family_members", []),
                "priority_date": content.get("priority_date", [None])[0],
                "countries": list(set(content.get("countries", []))),
            }

        except Exception as e:
            self.logger.error(f"Error parsing family data: {str(e)}")
            return None

    def _parse_citations_data(self, content: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse patent citation data.

        Args:
            content: Extracted content

        Returns:
            List of parsed citations
        """
        citations = []

        try:
            # Extract forward citations
            forward = content.get("forward_citations", [])
            for citation in forward:
                citations.append(
                    {
                        "patent_number": citation.get("patent_number"),
                        "title": citation.get("title"),
                        "filing_date": citation.get("filing_date"),
                        "relevance": citation.get("relevance"),
                        "type": "forward",
                    }
                )

            # Extract backward citations
            backward = content.get("backward_citations", [])
            for citation in backward:
                citations.append(
                    {
                        "patent_number": citation.get("patent_number"),
                        "title": citation.get("title"),
                        "filing_date": citation.get("filing_date"),
                        "relevance": citation.get("relevance"),
                        "type": "backward",
                    }
                )

        except Exception as e:
            self.logger.error(f"Error parsing citations data: {str(e)}")

        return citations

    def _parse_classifications(
        self,
        classifications: List[str],
    ) -> List[PatentClassification]:
        """Parse patent classifications.

        Args:
            classifications: Raw classification strings

        Returns:
            List of parsed classifications
        """
        parsed = []

        # Define regex pattern components
        system_pattern = r"(?P<system>[A-Z]+)"
        code_pattern = r"(?P<code>[A-Z0-9/]+)"
        desc_pattern = r"(?P<description>.*)"
        pattern = rf"{system_pattern}\s*" rf"{code_pattern}\s*" rf"{desc_pattern}"

        for raw in classifications:
            try:
                # Extract classification components
                match = re.match(pattern, raw)
                if not match:
                    continue

                system = match.group("system")
                code = match.group("code")
                description = match.group("description")

                # Determine classification level
                if len(code) == 1:
                    level = "section"
                elif len(code) <= 3:
                    level = "class"
                elif "/" not in code:
                    level = "subclass"
                else:
                    level = "group"

                parsed.append(
                    PatentClassification(
                        system=system,
                        code=code,
                        description=description,
                        level=level,
                        relevance=None,  # Could be calculated based on position
                    )
                )

            except Exception as e:
                self.logger.error(f"Error parsing classification: {str(e)}")
                continue

        return parsed

    def _build_patent_url(self, patent_number: Optional[str]) -> Optional[str]:
        """Build patent URL from patent number.

        Args:
            patent_number: Patent number

        Returns:
            Patent URL or None
        """
        if not patent_number:
            return None

        return f"{self.GOOGLE_PATENTS_URL}/patent/{patent_number}"

    def _calculate_confidence(self, patent: Dict[str, Any]) -> float:
        """Calculate confidence score for extracted data.

        Args:
            patent: Extracted patent data

        Returns:
            Confidence score between 0 and 1
        """
        score = 0.0
        total = 0

        # Check key fields presence
        fields = [
            "title",
            "abstract",
            "inventors",
            "assignee",
            "filing_date",
            "publication_date",
            "patent_number",
            "url",
        ]
        for field in fields:
            total += 1
            if patent.get(field):
                score += 1.0

        # Check extracted data quality
        if patent.get("compounds"):
            score += 1.0
            total += 1
        if patent.get("effects"):
            score += 1.0
            total += 1
        if patent.get("mechanisms"):
            score += 1.0
            total += 1
        if patent.get("safety"):
            score += 1.0
            total += 1

        # Check additional data
        if patent.get("family"):
            score += 1.0
            total += 1
        if patent.get("citations"):
            score += 1.0
            total += 1
        if patent.get("classifications"):
            score += 1.0
            total += 1

        # Calculate final score
        return score / total if total > 0 else 0.0
