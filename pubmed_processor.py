"""PubMed data processing and relevance scoring."""

import re
import time
from typing import Dict, Any, List, Optional, Tuple
from urllib.parse import quote_plus
import requests
from bs4 import BeautifulSoup
from ratelimit import limits, sleep_and_retry

from logger import LogManager


class PubMedProcessor:
    """Handles PubMed data processing and relevance scoring."""
    
    # NCBI E-utilities base URLs
    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize PubMed processor.
        
        Args:
            api_key: Optional NCBI API key for higher rate limits
        """
        self.logger = LogManager().get_logger("pubmed_processor")
        self.api_key = api_key
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'ChemDataCollector/0.1 (Research Project)'
        })

    @sleep_and_retry
    @limits(calls=3, period=1)  # Rate limit: 3 requests per second with API key, 1 without
    def _make_request(self, url: str, params: Dict[str, Any]) -> Optional[requests.Response]:
        """
        Make a rate-limited request to NCBI E-utilities.
        
        Args:
            url: E-utilities endpoint URL
            params: Query parameters
            
        Returns:
            Response object or None if failed
        """
        try:
            if self.api_key:
                params['api_key'] = self.api_key
                
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            # Required by NCBI: wait 0.34 seconds between requests without API key
            if not self.api_key:
                time.sleep(0.34)
                
            return response
        except Exception as e:
            self.logger.error(f"Error making request to {url}: {str(e)}")
            return None

    def get_binding_relevance(self, compound_name: str, target_name: str) -> int:
        """
        Get relevance score for compound-target binding pair based on PubMed results.
        
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
            params = {
                'db': 'pubmed',
                'term': query,
                'retmode': 'json',
                'retmax': 1000
            }
            
            response = self._make_request(self.ESEARCH_URL, params)
            if response:
                data = response.json()
                return int(data['esearchresult'].get('count', 0))
                
        except Exception as e:
            self.logger.error(f"Error getting binding relevance: {str(e)}")
            
        return 0

    def sort_names_by_relevance(self, names: List[str], context: str = "") -> List[Tuple[str, int]]:
        """
        Sort names by PubMed relevance score.
        
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
                    query += f' AND ({context})'
                    
                # Search PubMed
                params = {
                    'db': 'pubmed',
                    'term': query,
                    'retmode': 'json'
                }
                
                response = self._make_request(self.ESEARCH_URL, params)
                if response:
                    data = response.json()
                    score = int(data['esearchresult'].get('count', 0))
                    scored_names.append((name, score))
                    
            except Exception as e:
                self.logger.error(f"Error scoring name {name}: {str(e)}")
                scored_names.append((name, 0))
                
        # Sort by score
        return sorted(scored_names, key=lambda x: x[1], reverse=True)

    def get_compound_references(self, compound_name: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Get relevant PubMed references for a compound.
        
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
                'db': 'pubmed',
                'term': f'"{compound_name}"[Title/Abstract]',
                'retmode': 'json',
                'retmax': max_results,
                'sort': 'relevance'
            }
            
            search_response = self._make_request(self.ESEARCH_URL, search_params)
            if not search_response:
                return references
                
            search_data = search_response.json()
            pmids = search_data['esearchresult'].get('idlist', [])
            
            if not pmids:
                return references
                
            # Fetch article details
            fetch_params = {
                'db': 'pubmed',
                'id': ','.join(pmids),
                'retmode': 'xml'
            }
            
            fetch_response = self._make_request(self.EFETCH_URL, fetch_params)
            if not fetch_response:
                return references
                
            # Parse XML response
            soup = BeautifulSoup(fetch_response.text, 'xml')
            for article in soup.find_all('PubmedArticle'):
                try:
                    # Extract basic metadata
                    pmid = article.find('PMID').text
                    title = article.find('ArticleTitle').text
                    abstract = article.find('Abstract')
                    abstract_text = abstract.find('AbstractText').text if abstract else None
                    
                    # Extract authors
                    authors = []
                    author_list = article.find('AuthorList')
                    if author_list:
                        for author in author_list.find_all('Author'):
                            last_name = author.find('LastName')
                            fore_name = author.find('ForeName')
                            if last_name and fore_name:
                                authors.append(f"{last_name.text}, {fore_name.text}")
                    
                    # Extract journal info
                    journal = article.find('Journal')
                    if journal:
                        journal_title = journal.find('Title').text if journal.find('Title') else None
                        year = journal.find('Year').text if journal.find('Year') else None
                        volume = journal.find('Volume').text if journal.find('Volume') else None
                        issue = journal.find('Issue').text if journal.find('Issue') else None
                    
                    # Extract DOI
                    article_ids = article.find('ArticleIdList')
                    doi = None
                    if article_ids:
                        for id_elem in article_ids.find_all('ArticleId'):
                            if id_elem.get('IdType') == 'doi':
                                doi = id_elem.text
                                break
                    
                    # Create reference dictionary
                    reference = {
                        'pmid': pmid,
                        'title': title,
                        'abstract': abstract_text,
                        'authors': authors,
                        'journal': journal_title,
                        'year': year,
                        'volume': volume,
                        'issue': issue,
                        'doi': doi,
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    }
                    
                    references.append(reference)
                    
                except Exception as e:
                    self.logger.error(f"Error parsing article {pmid}: {str(e)}")
                    continue
                    
        except Exception as e:
            self.logger.error(f"Error getting compound references: {str(e)}")
            
        return references

    def analyze_binding_data(self, compound_name: str, target_name: str) -> Dict[str, Any]:
        """
        Analyze PubMed articles for binding data between compound and target.
        
        Args:
            compound_name: Name of compound
            target_name: Name of target
            
        Returns:
            Dictionary containing binding data analysis
        """
        analysis = {
            'total_articles': 0,
            'binding_articles': 0,
            'affinity_types': {
                'Ki': 0,
                'IC50': 0,
                'EC50': 0,
                'Kd': 0
            },
            'key_findings': []
        }
        
        try:
            # Search for binding-related articles
            query = f'"{compound_name}"[Title/Abstract] AND "{target_name}"[Title/Abstract]'
            
            search_params = {
                'db': 'pubmed',
                'term': query,
                'retmode': 'json',
                'retmax': 100
            }
            
            search_response = self._make_request(self.ESEARCH_URL, search_params)
            if not search_response:
                return analysis
                
            search_data = search_response.json()
            pmids = search_data['esearchresult'].get('idlist', [])
            analysis['total_articles'] = int(search_data['esearchresult'].get('count', 0))
            
            if not pmids:
                return analysis
                
            # Fetch and analyze articles
            fetch_params = {
                'db': 'pubmed',
                'id': ','.join(pmids),
                'retmode': 'xml'
            }
            
            fetch_response = self._make_request(self.EFETCH_URL, fetch_params)
            if not fetch_response:
                return analysis
                
            # Parse XML and analyze content
            soup = BeautifulSoup(fetch_response.text, 'xml')
            for article in soup.find_all('PubmedArticle'):
                try:
                    abstract = article.find('Abstract')
                    if not abstract:
                        continue
                        
                    abstract_text = abstract.find('AbstractText').text
                    
                    # Check for binding-related content
                    if re.search(r'bind|affinity|potency', abstract_text, re.I):
                        analysis['binding_articles'] += 1
                        
                        # Count affinity types
                        if re.search(r'Ki\s*[=~]', abstract_text):
                            analysis['affinity_types']['Ki'] += 1
                        if re.search(r'IC50\s*[=~]', abstract_text):
                            analysis['affinity_types']['IC50'] += 1
                        if re.search(r'EC50\s*[=~]', abstract_text):
                            analysis['affinity_types']['EC50'] += 1
                        if re.search(r'Kd\s*[=~]', abstract_text):
                            analysis['affinity_types']['Kd'] += 1
                            
                        # Extract key findings
                        findings = re.findall(
                            r'([^.]*(?:Ki|IC50|EC50|Kd)\s*[=~]\s*\d+(?:\.\d+)?\s*(?:nM|µM|pM|mM)[^.]*\.)',
                            abstract_text
                        )
                        if findings:
                            analysis['key_findings'].extend(findings)
                            
                except Exception as e:
                    self.logger.error(f"Error analyzing article: {str(e)}")
                    continue
                    
        except Exception as e:
            self.logger.error(f"Error analyzing binding data: {str(e)}")
            
        return analysis
