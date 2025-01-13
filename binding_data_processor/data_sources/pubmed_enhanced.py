"""PubMed API client for literature data enrichment.

This module provides an enhanced client for accessing PubMed data with:
- Binding relevance scoring
- Reference retrieval
- Literature-based analysis
- Circuit breaker pattern
- Metrics collection
- Enhanced error handling
- Rate limiting
- Caching
- Web scraping fallback when API fails
"""

from typing import Dict, Any, List, Optional, Tuple
import re
from bs4 import BeautifulSoup
from crawl4ai import AsyncWebCrawler, BrowserConfig
from crawl4ai.extraction_strategy import JsonCssExtractionStrategy, LLMExtractionStrategy

from ..web_enrichment.base_client_enhanced import BaseWebClientEnhanced
from ..models.compound import Compound
from ..pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig
from ..web_enrichment.llm_utils import extract_literature_info


class PubMedClientEnhanced(BaseWebClientEnhanced):
    """Enhanced client for PubMed API with improved literature analysis."""

    # NCBI E-utilities base URLs
    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(self, llm_provider: str = "ollama/llama2", api_token: Optional[str] = None, api_key: Optional[str] = None, **kwargs):
        """Initialize PubMed client with circuit breaker and metrics.

        Args:
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            api_key: Optional NCBI API key for higher rate limits
            **kwargs: Additional arguments passed to BaseWebClientEnhanced
        """
        # Configure circuit breaker
        circuit_config = CircuitBreakerConfig(
            failure_threshold=5,
            reset_timeout=300,
        )

        # Initialize with name and circuit breaker config
        super().__init__(name="pubmed", circuit_config=circuit_config, **kwargs)
        self.api_key = api_key

        # Update headers
        self.http.session.headers.update({"User-Agent": "ChemDataCollector/0.1 (Research Project)"})

        # Configure rate limits based on API key
        if api_key:
            # With API key: 10 requests per second
            self.http.configure_rate_limit(calls=10, period=1)
        else:
            # Without API key: 3 requests per second with 0.34s delay
            self.http.configure_rate_limit(calls=3, period=1, delay=0.34)

        # Configure crawl4ai for web scraping fallback
        self.crawler_config = BrowserConfig(
            javascript=BrowserConfig.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    "#article-details",
                    "#abstract",
                    "#full-view-heading",
                ],
                stealth_mode=True,
            ),
            screenshot=BrowserConfig.Screenshot(enabled=True, full_page=True),
            extraction=BrowserConfig.Extraction(
                llm=BrowserConfig.LLM(
                    provider=llm_provider,
                    api_token=api_token,
                    prompts={
                        "article": "Extract article details:",
                        "abstract": "Extract abstract and key findings:",
                        "references": "Extract reference information:",
                    },
                ),
                css={
                    "title": "#full-view-heading .heading-title",
                    "authors": "#full-view-heading .authors-list",
                    "abstract": "#abstract",
                    "journal": ".journal-citation",
                    "publication_details": ".publication-details",
                    "references": ".reference-list",
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
        )

    async def _scrape_article(self, pmid: str) -> Optional[Dict[str, Any]]:
        """Scrape article data from PubMed website.

        Args:
            pmid: PubMed ID

        Returns:
            Article data or None if scraping failed
        """
        try:
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            # Scrape with crawl4ai
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.crawler_config)

                if not result.success:
                    self.logger.error(f"Failed to scrape PubMed page for {pmid}")
                    return None

                # Extract data using both CSS and LLM
                content = result.extracted_content[0]

                # Extract article info using LLM
                article_info = extract_literature_info(
                    title=content.get("title", ""),
                    abstract=content.get("abstract", ""),
                )

                # Process authors
                authors = []
                for author in content.get("authors", "").split(";"):
                    author = author.strip()
                    if author:
                        authors.append(author)

                # Extract journal info
                journal_info = {}
                journal_citation = content.get("journal", "")
                if journal_citation:
                    # Parse journal citation (e.g. "Nature. 2021 Jan;589(7841):123-456.")
                    match = re.match(r"([^.]+)\.\s*(\d{4})\s*([^;]+);([^:]+):(.+)", journal_citation)
                    if match:
                        journal_info = {
                            "title": match.group(1).strip(),
                            "year": match.group(2),
                            "month": match.group(3).strip(),
                            "volume": match.group(4).strip(),
                            "pages": match.group(5).strip(),
                        }

                # Combine CSS and LLM extracted data
                data = {
                    "pmid": pmid,
                    "title": content.get("title"),
                    "abstract": content.get("abstract"),
                    "authors": authors,
                    "journal": journal_info.get("title"),
                    "year": journal_info.get("year"),
                    "volume": journal_info.get("volume"),
                    "issue": None,  # Not always available in citation
                    "doi": None,  # Need to extract from page metadata
                    "url": url,
                    **article_info,
                }

                return data

        except Exception as e:
            self.logger.error(f"Error scraping PubMed data: {str(e)}")
            return None

    def _make_request(self, url: str, params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Make a rate-limited request to NCBI E-utilities.

        Args:
            url: E-utilities endpoint URL
            params: Query parameters

        Returns:
            Response data or None if failed
        """
        try:
            if self.api_key:
                params["api_key"] = self.api_key

            def request_operation():
                response = self.http.get(url, params=params)
                return response.json() if "json" in params.get("retmode", "") else response.text

            # Use circuit breaker for request
            return self.http.circuit.execute(request_operation)

        except Exception as e:
            self.logger.error(f"Error making request to {url}: {str(e)}")
            return None

    async def get_binding_relevance(self, compound_name: str, target_name: str) -> int:
        """Get relevance score for compound-target binding pair.

        Args:
            compound_name: Name of compound
            target_name: Name of target

        Returns:
            Relevance score (number of PubMed results)
        """
        try:
            # Construct search query
            query = f'"{compound_name}"[Title/Abstract] AND "{target_name}"[Title/Abstract] AND ("binding" OR "affinity" OR "Ki" OR "IC50" OR "EC50" OR "Kd")'

            # Try API first
            params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": 1000}
            data = self._make_request(self.ESEARCH_URL, params)
            if data:
                return int(data["esearchresult"].get("count", 0))

            # Fall back to web scraping
            self.logger.info(f"Falling back to web scraping for relevance search")
            url = f"https://pubmed.ncbi.nlm.nih.gov/?term={query}"
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.crawler_config)
                if result.success:
                    content = result.extracted_content[0]
                    # Extract result count from page
                    count_text = content.get("search_count", "0")
                    return int(re.search(r"\d+", count_text).group())

        except Exception as e:
            self.logger.error(f"Error getting binding relevance: {str(e)}")

        return 0

    async def sort_names_by_relevance(self, names: List[str], context: str = "") -> List[Tuple[str, int]]:
        """Sort names by PubMed relevance score.

        Args:
            names: List of names to sort
            context: Optional context terms to include in search

        Returns:
            List of (name, score) tuples sorted by score
        """
        scored_names = []

        for name in names:
            try:
                # Construct search query
                query = f'"{name}"[Title/Abstract]'
                if context:
                    query += f" AND ({context})"

                # Try API first
                params = {"db": "pubmed", "term": query, "retmode": "json"}
                data = self._make_request(self.ESEARCH_URL, params)
                if data:
                    score = int(data["esearchresult"].get("count", 0))
                    scored_names.append((name, score))
                    continue

                # Fall back to web scraping
                self.logger.info(f"Falling back to web scraping for name relevance: {name}")
                url = f"https://pubmed.ncbi.nlm.nih.gov/?term={query}"
                async with AsyncWebCrawler() as crawler:
                    result = await crawler.arun(urls=[url], config=self.crawler_config)
                    if result.success:
                        content = result.extracted_content[0]
                        count_text = content.get("search_count", "0")
                        score = int(re.search(r"\d+", count_text).group())
                        scored_names.append((name, score))
                    else:
                        scored_names.append((name, 0))

            except Exception as e:
                self.logger.error(f"Error scoring name {name}: {str(e)}")
                scored_names.append((name, 0))

        # Sort by score
        return sorted(scored_names, key=lambda x: x[1], reverse=True)

    def _extract_article_metadata(self, article) -> Optional[Dict[str, Any]]:
        """Extract metadata from a PubMed article.

        Args:
            article: BeautifulSoup article element

        Returns:
            Dictionary of metadata or None if extraction failed
        """
        try:
            # Extract basic metadata
            pmid = article.find("PMID").text
            title = article.find("ArticleTitle").text
            abstract = article.find("Abstract")
            abstract_text = abstract.find("AbstractText").text if abstract else None

            # Extract authors
            authors = self._extract_authors(article.find("AuthorList"))

            # Extract journal info
            journal_info = self._extract_journal_info(article.find("Journal"))

            # Extract DOI
            doi = self._extract_doi(article.find("ArticleIdList"))

            # Create reference dictionary
            return {
                "pmid": pmid,
                "title": title,
                "abstract": abstract_text,
                "authors": authors,
                "journal": journal_info.get("title"),
                "year": journal_info.get("year"),
                "volume": journal_info.get("volume"),
                "issue": journal_info.get("issue"),
                "doi": doi,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }

        except Exception as e:
            self.logger.error(f"Error extracting article metadata: {str(e)}")
            return None

    def _extract_authors(self, author_list) -> List[str]:
        """Extract author names from author list.

        Args:
            author_list: BeautifulSoup author list element

        Returns:
            List of author names
        """
        authors = []
        if author_list:
            for author in author_list.find_all("Author"):
                last_name = author.find("LastName")
                fore_name = author.find("ForeName")
                if last_name and fore_name:
                    authors.append(f"{last_name.text}, {fore_name.text}")
        return authors

    def _extract_journal_info(self, journal) -> Dict[str, Optional[str]]:
        """Extract journal information.

        Args:
            journal: BeautifulSoup journal element

        Returns:
            Dictionary of journal information
        """
        info = {"title": None, "year": None, "volume": None, "issue": None}
        if journal:
            info["title"] = journal.find("Title").text if journal.find("Title") else None
            info["year"] = journal.find("Year").text if journal.find("Year") else None
            info["volume"] = journal.find("Volume").text if journal.find("Volume") else None
            info["issue"] = journal.find("Issue").text if journal.find("Issue") else None
        return info

    def _extract_doi(self, article_ids) -> Optional[str]:
        """Extract DOI from article IDs.

        Args:
            article_ids: BeautifulSoup article IDs element

        Returns:
            DOI string or None
        """
        if article_ids:
            for id_elem in article_ids.find_all("ArticleId"):
                if id_elem.get("IdType") == "doi":
                    return id_elem.text
        return None

    async def get_compound_references(self, compound_name: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """Get relevant PubMed references for a compound.

        Args:
            compound_name: Name of compound
            max_results: Maximum number of results to return

        Returns:
            List of reference dictionaries
        """
        references = []

        try:
            # Try API first
            search_params = {
                "db": "pubmed",
                "term": f'"{compound_name}"[Title/Abstract]',
                "retmode": "json",
                "retmax": max_results,
                "sort": "relevance",
            }

            search_data = self._make_request(self.ESEARCH_URL, search_params)
            if search_data:
                pmids = search_data["esearchresult"].get("idlist", [])
                if pmids:
                    # Fetch article details
                    fetch_params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
                    xml_text = self._make_request(self.EFETCH_URL, fetch_params)
                    if xml_text:
                        # Parse XML response
                        soup = BeautifulSoup(xml_text, "xml")
                        for article in soup.find_all("PubmedArticle"):
                            reference = self._extract_article_metadata(article)
                            if reference:
                                references.append(reference)
                        return references

            # Fall back to web scraping
            self.logger.info(f"Falling back to web scraping for compound references: {compound_name}")
            url = f"https://pubmed.ncbi.nlm.nih.gov/?term={compound_name}"
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.crawler_config)
                if result.success:
                    content = result.extracted_content[0]
                    # Extract PMIDs from search results
                    pmids = re.findall(r"/(\d+)/", content.get("search_results", ""))[:max_results]

                    # Fetch each article
                    for pmid in pmids:
                        article_data = await self._scrape_article(pmid)
                        if article_data:
                            references.append(article_data)

        except Exception as e:
            self.logger.error(f"Error getting compound references: {str(e)}")

        return references

    def _analyze_abstract(self, abstract_text: str, analysis: Dict[str, Any]) -> None:
        """Analyze abstract text for binding data.

        Args:
            abstract_text: Abstract text to analyze
            analysis: Analysis dictionary to update
        """
        # Check for binding-related content
        if re.search(r"bind|affinity|potency", abstract_text, re.I):
            analysis["binding_articles"] += 1

            # Count affinity types
            affinity_patterns = {"Ki": r"Ki\s*[=~]", "IC50": r"IC50\s*[=~]", "EC50": r"EC50\s*[=~]", "Kd": r"Kd\s*[=~]"}
            for affinity_type, pattern in affinity_patterns.items():
                if re.search(pattern, abstract_text):
                    analysis["affinity_types"][affinity_type] += 1

            # Extract key findings
            findings = re.findall(r"([^.]*(?:Ki|IC50|EC50|Kd)\s*[=~]\s*\d+(?:\.\d+)?\s*(?:nM|µM|pM|mM)[^.]*\.)", abstract_text)
            if findings:
                analysis["key_findings"].extend(findings)

    async def analyze_binding_data(self, compound_name: str, target_name: str) -> Dict[str, Any]:
        """Analyze PubMed articles for binding data between compound and target.

        Args:
            compound_name: Name of compound
            target_name: Name of target

        Returns:
            Dictionary containing binding data analysis
        """
        analysis = {
            "total_articles": 0,
            "binding_articles": 0,
            "affinity_types": {"Ki": 0, "IC50": 0, "EC50": 0, "Kd": 0},
            "key_findings": [],
        }

        try:
            # Search for binding-related articles
            query = f'"{compound_name}"[Title/Abstract] AND "{target_name}"[Title/Abstract]'

            # Try API first
            search_params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": 100}
            search_data = self._make_request(self.ESEARCH_URL, search_params)

            if search_data:
                pmids = search_data["esearchresult"].get("idlist", [])
                analysis["total_articles"] = int(search_data["esearchresult"].get("count", 0))

                if pmids:
                    # Fetch and analyze articles
                    fetch_params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
                    xml_text = self._make_request(self.EFETCH_URL, fetch_params)
                    if xml_text:
                        # Parse XML and analyze content
                        soup = BeautifulSoup(xml_text, "xml")
                        for article in soup.find_all("PubmedArticle"):
                            try:
                                abstract = article.find("Abstract")
                                if not abstract:
                                    continue

                                abstract_text = abstract.find("AbstractText").text
                                self._analyze_abstract(abstract_text, analysis)

                            except Exception as e:
                                self.logger.error(f"Error analyzing article: {str(e)}")
                                continue
                        return analysis

            # Fall back to web scraping
            self.logger.info(f"Falling back to web scraping for binding analysis")
            url = f"https://pubmed.ncbi.nlm.nih.gov/?term={query}"
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.crawler_config)
                if result.success:
                    content = result.extracted_content[0]

                    # Get total count
                    count_text = content.get("search_count", "0")
                    analysis["total_articles"] = int(re.search(r"\d+", count_text).group())

                    # Extract PMIDs from search results
                    pmids = re.findall(r"/(\d+)/", content.get("search_results", ""))[:100]

                    # Analyze each article
                    for pmid in pmids:
                        article_data = await self._scrape_article(pmid)
                        if article_data and article_data.get("abstract"):
                            self._analyze_abstract(article_data["abstract"], analysis)

        except Exception as e:
            self.logger.error(f"Error analyzing binding data: {str(e)}")

        return analysis

    async def process_compounds(
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
            # Get references
            if compound.name:
                references = await self.get_compound_references(compound.name)
                compound.pubmed_references = references

            # Get binding data if target is known
            if compound.name and compound.target_name:
                binding_data = await self.analyze_binding_data(compound.name, compound.target_name)
                compound.binding_data = binding_data

    async def get_compound_data(
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
        # Get references
        references = await self.get_compound_references(name)
        if not references:
            return None

        # Extract data from references
        data = {
            "name": name,
            "cas_number": cas_number,
            "references": references,
            "binding_data": [],
            "activity_data": [],
            "safety_data": [],
        }

        # Analyze abstracts for relevant data
        for ref in references:
            abstract = ref.get("abstract", "")
            if not abstract:
                continue

            # Look for binding data
            if re.search(r"bind|affinity|Ki|IC50|EC50|Kd", abstract, re.I):
                findings = re.findall(r"([^.]*(?:Ki|IC50|EC50|Kd)\s*[=~]\s*\d+(?:\.\d+)?\s*(?:nM|µM|pM|mM)[^.]*\.)", abstract)
                if findings:
                    data["binding_data"].extend(findings)

            # Look for activity data
            if re.search(r"activity|effect|potency|response", abstract, re.I):
                findings = re.findall(r"([^.]*(?:activity|effect|potency|response)[^.]*\.)", abstract)
                if findings:
                    data["activity_data"].extend(findings)

            # Look for safety data
            if re.search(r"toxic|safety|adverse|side effect", abstract, re.I):
                findings = re.findall(r"([^.]*(?:toxic|safety|adverse|side effect)[^.]*\.)", abstract)
                if findings:
                    data["safety_data"].extend(findings)

        return data
