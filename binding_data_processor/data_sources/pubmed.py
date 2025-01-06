"""PubMed API client for literature data enrichment.

This module provides a client for accessing PubMed data with:
- Binding relevance scoring
- Reference retrieval
- Literature-based analysis
- Rate limiting
- Caching
"""

from typing import Dict, Any, List, Optional, Tuple
import re
import time
from bs4 import BeautifulSoup
from ratelimit import limits, sleep_and_retry

from ..web_enrichment.base_client import BaseWebClient
from ..models.compound import Compound


class PubMedClient(BaseWebClient):
    """Client for PubMed API with enhanced literature analysis."""

    # NCBI E-utilities base URLs
    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(self, api_key: Optional[str] = None, **kwargs):
        """Initialize PubMed client.

        Args:
            api_key: Optional NCBI API key for higher rate limits
            **kwargs: Additional arguments passed to BaseWebClient
        """
        super().__init__(**kwargs)
        self.api_key = api_key

        # Update headers
        self.http.session.headers.update({"User-Agent": "ChemDataCollector/0.1 (Research Project)"})

    @sleep_and_retry
    @limits(calls=3, period=1)  # Rate limit: 3 requests per second with API key
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

            response = self.http.get(url, params=params)

            # Required by NCBI: wait 0.34 seconds between requests without API key
            if not self.api_key:
                time.sleep(0.34)

            return response.json() if "json" in params.get("retmode", "") else response.text

        except Exception as e:
            self.logger.error(f"Error making request to {url}: {str(e)}")
            return None

    def get_binding_relevance(self, compound_name: str, target_name: str) -> int:
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

            # Search PubMed
            params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": 1000}

            data = self._make_request(self.ESEARCH_URL, params)
            if data:
                return int(data["esearchresult"].get("count", 0))

        except Exception as e:
            self.logger.error(f"Error getting binding relevance: {str(e)}")

        return 0

    def sort_names_by_relevance(self, names: List[str], context: str = "") -> List[Tuple[str, int]]:
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

                # Search PubMed
                params = {"db": "pubmed", "term": query, "retmode": "json"}

                data = self._make_request(self.ESEARCH_URL, params)
                if data:
                    score = int(data["esearchresult"].get("count", 0))
                    scored_names.append((name, score))

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

    def get_compound_references(self, compound_name: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """Get relevant PubMed references for a compound.

        Args:
            compound_name: Name of compound
            max_results: Maximum number of results to return

        Returns:
            List of reference dictionaries
        """
        references = []

        try:
            # Search PubMed
            search_params = {
                "db": "pubmed",
                "term": f'"{compound_name}"[Title/Abstract]',
                "retmode": "json",
                "retmax": max_results,
                "sort": "relevance",
            }

            search_data = self._make_request(self.ESEARCH_URL, search_params)
            if not search_data:
                return references

            pmids = search_data["esearchresult"].get("idlist", [])
            if not pmids:
                return references

            # Fetch article details
            fetch_params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
            xml_text = self._make_request(self.EFETCH_URL, fetch_params)
            if not xml_text:
                return references

            # Parse XML response
            soup = BeautifulSoup(xml_text, "xml")
            for article in soup.find_all("PubmedArticle"):
                reference = self._extract_article_metadata(article)
                if reference:
                    references.append(reference)

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
            findings = re.findall(
                r"([^.]*(?:Ki|IC50|EC50|Kd)\s*[=~]\s*\d+(?:\.\d+)?\s*(?:nM|µM|pM|mM)[^.]*\.)", abstract_text
            )
            if findings:
                analysis["key_findings"].extend(findings)

    def analyze_binding_data(self, compound_name: str, target_name: str) -> Dict[str, Any]:
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
            search_params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": 100}

            search_data = self._make_request(self.ESEARCH_URL, search_params)
            if not search_data:
                return analysis

            pmids = search_data["esearchresult"].get("idlist", [])
            analysis["total_articles"] = int(search_data["esearchresult"].get("count", 0))

            if not pmids:
                return analysis

            # Fetch and analyze articles
            fetch_params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
            xml_text = self._make_request(self.EFETCH_URL, fetch_params)
            if not xml_text:
                return analysis

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

        except Exception as e:
            self.logger.error(f"Error analyzing binding data: {str(e)}")

        return analysis

    def process_compounds(
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
                references = self.get_compound_references(compound.name)
                compound.pubmed_references = references

            # Get binding data if target is known
            if compound.name and compound.target_name:
                binding_data = self.analyze_binding_data(compound.name, compound.target_name)
                compound.binding_data = binding_data

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
        # Get references
        references = self.get_compound_references(name)
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
                findings = re.findall(
                    r"([^.]*(?:Ki|IC50|EC50|Kd)\s*[=~]\s*\d+(?:\.\d+)?\s*(?:nM|µM|pM|mM)[^.]*\.)", abstract
                )
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
