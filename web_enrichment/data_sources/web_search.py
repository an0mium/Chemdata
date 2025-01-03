"""Web search and content analysis functionality.

This module handles:
1. Web search using SERP API
2. Content extraction from generic web pages
3. LLM-based content analysis
4. Structure-based chemical searching
5. Patent document analysis
6. Scientific literature search
"""

from typing import Dict, Any, List, Optional, Tuple, Set
from urllib.parse import quote, urljoin
import os
import re
import requests
from bs4 import BeautifulSoup
import json

from logger import LogManager
from ..http_client import HttpClient
from ..llm_utils import analyze_content_with_llm
from ..name_utils import standardize_chemical_name

logger = LogManager().get_logger("web_enrichment.data_sources.web_search")


class WebSearchClient:
    """Client for web search and content analysis."""
    
    # Chemical information sources
    CHEMICAL_SOURCES = {
        'databases': [
            'pubchem.ncbi.nlm.nih.gov',
            'bindingdb.org',
            'chembl.ebi.ac.uk',
            'drugbank.ca',
            'zinc.docking.org',
            'commonchemistry.cas.org',
            'chemspider.com'
        ],
        'literature': [
            'pubmed.ncbi.nlm.nih.gov',
            'scholar.google.com',
            'europepmc.org',
            'sciencedirect.com',
            'pubs.acs.org',
            'nature.com',
            'science.org'
        ],
        'patents': [
            'patents.google.com',
            'worldwide.espacenet.com',
            'patentscope.wipo.int',
            'lens.org',
            'freepatentsonline.com'
        ],
        'regulatory': [
            'deadiversion.usdoj.gov',
            'who.int',
            'emcdda.europa.eu',
            'ema.europa.eu',
            'fda.gov',
            'isomerdesign.com/Cdsa'
        ],
        'community': [
            'wikipedia.org',
            'erowid.org',
            'psychonautwiki.org',
            'isomerdesign.com',
            'drugbank.ca/drugs'
        ]
    }
    
    # Chemical identifier patterns with improved regex
    IDENTIFIER_PATTERNS = {
        'cas': r'\b\d{1,7}-\d{2}-\d\b',
        'inchi': r'InChI=1S?/[A-Za-z0-9/.]+(?:[+-][A-Za-z0-9/.]+)*',
        'smiles': r'(?:[A-Za-z0-9@\[\](){}\\/=#\-+]+\.?)+',
        'pubchem_cid': r'\bCID:?\s*(\d+)\b',
        'chembl_id': r'\bCHEMBL\d+\b',
        'drugbank_id': r'\bDB\d{5}\b',
        'unii': r'\b[A-Z0-9]{10}\b',  # FDA UNII codes
        'einecs': r'\b[24][0-9]{2}-[0-9]{3}-[0-9]\b'  # European chemical numbers
    }
    
    # Chemical name patterns
    CHEMICAL_NAME_PATTERNS = [
        r'\b[A-Z][a-z]*(?:-[0-9]+[A-Z][a-z]*)+\b',  # e.g., 2C-B, 5-MeO-DMT
        r'\b[0-9]?-[A-Z]{2,4}-[A-Z0-9]+\b',  # e.g., 2-FA, 4-HO-MET
        r'\b[A-Z]{3,5}-[0-9]+\b',  # e.g., NBOMe-2C-C
        r'\b(?:alpha|beta|delta|gamma)-[A-Z][A-Za-z0-9-]+\b',  # e.g., alpha-PVP
        r'\b[0-9][A-Z]{1,2}[0-9]{0,2}-[A-Z][A-Za-z-]+\b'  # e.g., 25I-NBOMe
    ]

    
    def __init__(self, http_client: HttpClient):
        """Initialize web search client."""
        self.http = http_client

    def search_and_analyze(
        self,
        query: str,
        llm_api_key: str,
        compound_data: Dict[str, Any],
        excluded_domains: Optional[List[str]] = None,
        structure_search: bool = True,
        search_patents: bool = True,
        search_literature: bool = True,
        max_results_per_category: int = 5
    ) -> Dict[str, Any]:
        """
        Search for information and analyze content using LLM.
        
        Args:
            query: Search query
            llm_api_key: API key for LLM service
            compound_data: Known compound data for context
            excluded_domains: Optional list of domains to exclude
            structure_search: Whether to include structure-based searching
            search_patents: Whether to search patent documents
            search_literature: Whether to search scientific literature
            max_results_per_category: Maximum results to process per category
            
        Returns:
            Dictionary containing:
            - urls: List of found URLs
            - extracted_data: Dictionary of extracted data
            - identifiers: Dictionary of found chemical identifiers
            - references: Dictionary of references by type
            - chemical_names: List of found chemical names
            - structures: List of found chemical structures
        """
        result = {
            'urls': [],
            'extracted_data': {},
            'identifiers': {
                id_type: set() for id_type in self.IDENTIFIER_PATTERNS
            },
            'references': {
                'patents': [],
                'papers': [],
                'databases': []
            },
            'chemical_names': set(),
            'structures': []
        }
        
        # Build enhanced search queries
        search_queries = self._build_search_queries(
            query,
            compound_data,
            structure_search,
            search_patents,
            search_literature
        )
        
        # Search with each query
        all_results = []
        for search_query in search_queries:
            results = self._search_google(search_query, excluded_domains)
            all_results.extend(results)
        
        # Deduplicate and categorize results
        seen_urls = set()
        categorized_results = {
            category: [] for category in self.CHEMICAL_SOURCES
        }
        
        for search_result in all_results:
            url = search_result['link']
            if url in seen_urls:
                continue
            seen_urls.add(url)
            
            # Categorize URL
            for category, domains in self.CHEMICAL_SOURCES.items():
                if any(domain in url for domain in domains):
                    categorized_results[category].append(search_result)
                    break
        
        # Process each category with appropriate extraction method
        for category, results in categorized_results.items():
            for result_data in results[:max_results_per_category]:
                url = result_data['link']
                
                # Get page content
                is_valid, content = self._get_page_content(url)
                if not is_valid or not content:
                    continue
                
                # Extract chemical identifiers
                self._extract_identifiers(content, result['identifiers'])
                
                # Extract chemical names
                self._extract_chemical_names(content, result['chemical_names'])
                
                # Extract references
                self._extract_references(content, url, result['references'])
                
                # Category-specific processing
                if category == 'patents':
                    extracted = self._process_patent_content(
                        content,
                        url,
                        compound_data,
                        llm_api_key
                    )
                elif category == 'literature':
                    extracted = self._process_literature_content(
                        content,
                        url,
                        compound_data,
                        llm_api_key
                    )
                else:
                    # General content analysis
                    extracted = analyze_content_with_llm(
                        content,
                        compound_data,
                        llm_api_key
                    )
                
                if extracted:
                    result['urls'].append(url)
                    result['extracted_data'][url] = extracted
                    
                    # Extract structures from extracted data
                    if 'structures' in extracted:
                        result['structures'].extend(extracted['structures'])
        
        # Convert sets to lists for JSON serialization
        result['identifiers'] = {
            k: list(v) for k, v in result['identifiers'].items()
        }
        result['chemical_names'] = list(result['chemical_names'])
        
        return result

    def _build_search_queries(
        self,
        base_query: str,
        compound_data: Dict[str, Any],
        include_structure: bool = True,
        include_patents: bool = True,
        include_literature: bool = True
    ) -> List[str]:
        """
        Build multiple search queries using different identifiers.
        
        Args:
            base_query: Base search query
            compound_data: Known compound data
            include_structure: Whether to include structure-based queries
            include_patents: Whether to include patent-specific queries
            include_literature: Whether to include literature-specific queries
            
        Returns:
            List of search queries
        """
        queries = []
        
        # Add base query
        queries.append(base_query)
        
        # Add standardized name query
        std_name = standardize_chemical_name(base_query)
        if std_name != base_query:
            queries.append(std_name)
        
        if include_structure and compound_data:
            # Add structure-based identifiers
            if 'smiles' in compound_data:
                queries.append(f'"{compound_data["smiles"]}" {base_query}')
            if 'inchi' in compound_data:
                queries.append(f'"{compound_data["inchi"]}" {base_query}')
            if 'cas' in compound_data:
                queries.append(f'"{compound_data["cas"]}" {base_query}')
            
            # Add specific database identifiers
            for id_type in ['pubchem_cid', 'chembl_id', 'drugbank_id']:
                if id_type in compound_data:
                    queries.append(f'"{compound_data[id_type]}" {base_query}')
        
        # Add source-specific queries
        for category, domains in self.CHEMICAL_SOURCES.items():
            if category == 'patents' and not include_patents:
                continue
            if category == 'literature' and not include_literature:
                continue
                
            for domain in domains:
                queries.append(f'{base_query} site:{domain}')
                
                # Add advanced queries for certain sources
                if domain == 'patents.google.com' and include_patents:
                    queries.append(f'"{base_query}" synthesis characterization')
                    queries.append(f'"{base_query}" preparation example')
                elif domain == 'pubmed.ncbi.nlm.nih.gov' and include_literature:
                    queries.append(f'"{base_query}" pharmacology mechanism')
                    queries.append(f'"{base_query}" binding affinity')
        
        return queries

    def _search_google(
        self,
        query: str,
        excluded_domains: Optional[List[str]] = None
    ) -> List[Dict[str, str]]:
        """
        Search Google using SERP API.
        
        Args:
            query: Search query
            excluded_domains: Optional list of domains to exclude
            
        Returns:
            List of search results
        """
        try:
            # Build query with exclusions
            if excluded_domains:
                exclusions = ' '.join(f'-site:{domain}' for domain in excluded_domains)
                query = f'{query} {exclusions}'
            
            url = "https://serpapi.com/search"
            params = {
                'q': query,
                'api_key': os.getenv('SERP_API_KEY'),
                'engine': 'google',
                'num': 10,  # Get more results per query
                'gl': 'us',  # Set region to US for consistent results
                'hl': 'en'   # Set language to English
            }
            
            response = self.http.make_request(url, params=params)
            if response:
                data = response.json()
                return data.get('organic_results', [])
                
        except Exception as e:
            logger.error(f"Error searching Google: {str(e)}")
            
        return []

    def _get_page_content(self, url: str) -> Tuple[bool, Optional[str]]:
        """
        Get and validate page content.
        
        Args:
            url: URL to fetch
            
        Returns:
            Tuple of (is_valid, content)
        """
        try:
            response = self.http.make_request(url)
            if not response:
                return False, None
                
            # Extract text content
            soup = BeautifulSoup(response.text, 'html.parser')
            
            # Remove script and style elements
            for script in soup(["script", "style"]):
                script.decompose()
            
            # Get text and normalize whitespace
            content = ' '.join(soup.stripped_strings)
            
            # Basic validation
            if len(content) < 100:  # Too short to be useful
                return False, None
                
            return True, content
            
        except Exception as e:
            logger.error(f"Error getting page content: {str(e)}")
            return False, None

    def _extract_identifiers(
        self,
        content: str,
        identifiers: Dict[str, Set[str]]
    ) -> None:
        """Extract chemical identifiers from text."""
        for id_type, pattern in self.IDENTIFIER_PATTERNS.items():
            matches = re.finditer(pattern, content, re.I)
            for match in matches:
                identifiers[id_type].add(match.group())

    def _extract_references(
        self,
        content: str,
        url: str,
        references: Dict[str, List[Dict[str, str]]]
    ) -> None:
        """Extract references from content."""
        # Extract patent numbers
        patent_patterns = [
            r'US\d{7,8}[A-Z]\d?',
            r'EP\d{7,8}[A-Z]\d?',
            r'WO\d{4}/\d{6}'
        ]
        for pattern in patent_patterns:
            matches = re.finditer(pattern, content)
            for match in matches:
                references['patents'].append({
                    'id': match.group(),
                    'source_url': url
                })
        
        # Extract DOIs
        doi_matches = re.finditer(
            r'10\.\d{4,}/[-._;()/:\w]+',
            content
        )
        for match in doi_matches:
            references['papers'].append({
                'doi': match.group(),
                'source_url': url
            })
        
        # Extract database references
        db_patterns = {
            'pubchem': r'CID:\s*\d+',
            'chembl': r'CHEMBL\d+',
            'drugbank': r'DB\d{5}'
        }
        for db, pattern in db_patterns.items():
            matches = re.finditer(pattern, content)
            for match in matches:
                references['databases'].append({
                    'database': db,
                    'id': match.group(),
                    'source_url': url
                })

    def _extract_chemical_names(
        self,
        content: str,
        chemical_names: Set[str]
    ) -> None:
        """Extract chemical names from text."""
        # Extract using patterns
        for pattern in self.CHEMICAL_NAME_PATTERNS:
            matches = re.finditer(pattern, content)
            for match in matches:
                name = match.group().strip()
                if name:
                    chemical_names.add(name)
        
        # Extract IUPAC-like names
        iupac_pattern = r'\b[1-9]?-?(?:(?:[A-Z][a-z]*)?(?:yl|ylidene|oxy|oxo|thio|amino|imino|hydro|hydroxy|mercapto|phospho|sulfo|nitro|nitroso|azo|diazo|phenyl|methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl|octyl|nonyl|decyl|fluoro|chloro|bromo|iodo))+[a-z]*\b'
        matches = re.finditer(iupac_pattern, content)
        for match in matches:
            name = match.group().strip()
            if name and len(name) > 5:  # Filter out very short matches
                chemical_names.add(name)

    def _process_patent_content(
        self,
        content: str,
        url: str,
        compound_data: Dict[str, Any],
        llm_api_key: str
    ) -> Optional[Dict[str, Any]]:
        """Process patent content with specialized prompts."""
        # Check if URL is a PDF
        if url.endswith('.pdf'):
            content = self._extract_text_from_pdf(url)
            if not content:
                return None
        
        # Build context-aware prompt
        prompt = f"""Analyze this patent content for chemical compound information:

Known compound data:
Name: {compound_data.get('name', 'Unknown')}
SMILES: {compound_data.get('smiles', 'Unknown')}
CAS: {compound_data.get('cas', 'Unknown')}

Extract the following information if present:
1. Synthesis procedures and conditions
2. Physical and chemical properties
3. Biological activity data
4. Structure-activity relationships
5. Example compounds and their properties
6. Chemical structures (SMILES, InChI)
7. Analytical data
8. Pharmacological data
9. Target receptors and binding data

Format the response as a JSON object with these fields."""

        try:
            # Extract example sections
            examples = self._extract_patent_examples(content)
            
            # Process each example
            extracted_data = {
                'examples': [],
                'structures': [],
                'properties': {},
                'synthesis': [],
                'activity': []
            }
            
            for example in examples:
                example_data = analyze_content_with_llm(
                    example,
                    compound_data,
                    llm_api_key,
                    prompt=prompt
                )
                if example_data:
                    extracted_data['examples'].append(example_data)
                    
                    # Extract structures
                    if 'structures' in example_data:
                        extracted_data['structures'].extend(
                            example_data['structures']
                        )
                    
                    # Extract properties
                    if 'properties' in example_data:
                        extracted_data['properties'].update(
                            example_data['properties']
                        )
                    
                    # Extract synthesis
                    if 'synthesis' in example_data:
                        extracted_data['synthesis'].extend(
                            example_data['synthesis']
                        )
                    
                    # Extract activity
                    if 'activity' in example_data:
                        extracted_data['activity'].extend(
                            example_data['activity']
                        )
            
            return extracted_data
            
        except Exception as e:
            logger.error(f"Error processing patent content: {str(e)}")
            return None

    def _extract_patent_examples(self, content: str) -> List[str]:
        """Extract example sections from patent content."""
        examples = []
        
        # Common example section patterns
        patterns = [
            r'Example\s+\d+[.\s]+(.*?)(?=Example\s+\d+[.\s]+|$)',
            r'Example\s+[A-Z][.\s]+(.*?)(?=Example\s+[A-Z][.\s]+|$)',
            r'Preparation\s+\d+[.\s]+(.*?)(?=Preparation\s+\d+[.\s]+|$)'
        ]
        
        for pattern in patterns:
            matches = re.finditer(pattern, content, re.DOTALL)
            for match in matches:
                example = match.group(1).strip()
                if example:
                    examples.append(example)
        
        return examples

    def _extract_text_from_pdf(self, url: str) -> Optional[str]:
        """Extract text from PDF URL."""
        try:
            # Download PDF
            response = self.http.make_request(url)
            if not response:
                return None
            
            # Use PyPDF2 to extract text
            import io
            from PyPDF2 import PdfReader
            
            pdf_file = io.BytesIO(response.content)
            pdf_reader = PdfReader(pdf_file)
            
            # Extract text from all pages
            text = []
            for page in pdf_reader.pages:
                text.append(page.extract_text())
            
            return '\n'.join(text)
            
        except Exception as e:
            logger.error(f"Error extracting PDF text: {str(e)}")
            return None


    def _process_literature_content(
        self,
        content: str,
        compound_data: Dict[str, Any],
        llm_api_key: str
    ) -> Optional[Dict[str, Any]]:
        """Process scientific literature content with specialized prompts."""
        prompt = f"""Analyze this scientific literature for chemical compound information:

Known compound data:
Name: {compound_data.get('name', 'Unknown')}
SMILES: {compound_data.get('smiles', 'Unknown')}
CAS: {compound_data.get('cas', 'Unknown')}

Extract the following information if present:
1. Experimental methods and conditions
2. Analytical data and characterization
3. Pharmacological activity data
4. Structure-activity relationships
5. Binding affinity data
6. Mechanism of action

Format the response as a JSON object with these fields."""

        try:
            extracted = analyze_content_with_llm(
                content,
                compound_data,
                llm_api_key,
                prompt=prompt
            )
            return extracted
        except Exception as e:
            logger.error(f"Error processing literature content: {str(e)}")
            return None
