"""API client with circuit breaker pattern and retry mechanisms."""

import time
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote, urljoin

import requests
from requests.exceptions import RequestException

from cache_manager import CacheManager
from config import (API_RATE_LIMIT, CIRCUIT_BREAKER, MAX_RETRIES, PUBCHEM_BASE_URL,
                   PUBMED_BASE_URL, RETRY_DELAY, SERP_API_KEY)
from logger import LogManager


class CircuitState(Enum):
    """Circuit breaker states."""
    CLOSED = 'closed'  # Normal operation
    OPEN = 'open'     # Failing, reject requests
    HALF_OPEN = 'half_open'  # Testing if service recovered


@dataclass
class CircuitBreaker:
    """Implements circuit breaker pattern for API calls."""
    
    name: str
    failure_threshold: int = CIRCUIT_BREAKER['failure_threshold']
    reset_timeout: int = CIRCUIT_BREAKER['reset_timeout']
    half_open_timeout: int = CIRCUIT_BREAKER['half_open_timeout']
    
    def __post_init__(self):
        """Initialize circuit breaker state."""
        self.state = CircuitState.CLOSED
        self.failures = 0
        self.last_failure_time = 0
        self.logger = LogManager().get_logger(f"circuit_breaker.{self.name}")

    def record_failure(self):
        """Record a failure and potentially open the circuit."""
        self.failures += 1
        self.last_failure_time = time.time()
        
        if self.failures >= self.failure_threshold:
            self.state = CircuitState.OPEN
            self.logger.warning(
                f"Circuit {self.name} opened after {self.failures} failures"
            )

    def record_success(self):
        """Record a success and potentially close the circuit."""
        self.failures = 0
        if self.state in (CircuitState.OPEN, CircuitState.HALF_OPEN):
            self.state = CircuitState.CLOSED
            self.logger.info(f"Circuit {self.name} closed after success")

    def can_execute(self) -> bool:
        """Check if request can be executed based on circuit state."""
        if self.state == CircuitState.CLOSED:
            return True
            
        if self.state == CircuitState.OPEN:
            if time.time() - self.last_failure_time >= self.reset_timeout:
                self.state = CircuitState.HALF_OPEN
                self.logger.info(f"Circuit {self.name} entering half-open state")
                return True
            return False
            
        # HALF_OPEN state
        if time.time() - self.last_failure_time >= self.half_open_timeout:
            return True
        return False


class APIClient:
    """Client for making API requests with caching and circuit breaker."""
    
    def __init__(self, base_url: str, name: str):
        """
        Initialize API client.
        
        Args:
            base_url: Base URL for API requests
            name: Name for this API client instance
        """
        self.base_url = base_url
        self.name = name
        self.cache = CacheManager()
        self.circuit_breaker = CircuitBreaker(name)
        self.logger = LogManager().get_logger(f"api_client.{name}")
        self.last_request_time = 0

    def _wait_for_rate_limit(self):
        """Wait to respect rate limiting."""
        elapsed = time.time() - self.last_request_time
        if elapsed < API_RATE_LIMIT:
            time.sleep(API_RATE_LIMIT - elapsed)
        self.last_request_time = time.time()

    def _make_request(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Make HTTP request with retries and circuit breaker.
        
        Args:
            method: HTTP method (GET, POST, etc.)
            endpoint: API endpoint
            params: Query parameters
            **kwargs: Additional arguments for requests
            
        Returns:
            API response data
            
        Raises:
            McpError: On API or circuit breaker errors
        """
        if not self.circuit_breaker.can_execute():
            raise McpError(
                ErrorCode.CircuitBreakerOpen,
                f"Circuit breaker {self.name} is open"
            )
            
        url = urljoin(self.base_url, endpoint)
        cache_key = f"{method}:{url}:{str(params)}:{str(kwargs)}"
        
        # Check cache first
        cached = self.cache.get(cache_key)
        if cached is not None:
            self.logger.debug(f"Cache hit for {cache_key}")
            return cached
            
        # Make request with retries
        for attempt in range(MAX_RETRIES):
            try:
                self._wait_for_rate_limit()
                
                response = requests.request(
                    method,
                    url,
                    params=params,
                    **kwargs
                )
                response.raise_for_status()
                
                data = response.json()
                self.cache.set(cache_key, data)
                self.circuit_breaker.record_success()
                return data
                
            except RequestException as e:
                self.logger.warning(
                    f"Request failed (attempt {attempt + 1}/{MAX_RETRIES}): {str(e)}"
                )
                self.circuit_breaker.record_failure()
                
                if attempt < MAX_RETRIES - 1:
                    time.sleep(RETRY_DELAY * (attempt + 1))
                    continue
                    
                raise McpError(
                    ErrorCode.APIError,
                    f"API request failed after {MAX_RETRIES} attempts: {str(e)}"
                )


class PubChemClient(APIClient):
    """Client for PubChem API."""
    
    def __init__(self):
        """Initialize PubChem client."""
        super().__init__(PUBCHEM_BASE_URL, "pubchem")

    def get_compound_by_name(self, name: str) -> Dict[str, Any]:
        """
        Get compound data by name.
        
        Args:
            name: Compound name
            
        Returns:
            Compound data
        """
        # First get CID
        endpoint = f"compound/name/{quote(name)}/cids/JSON"
        response = self._make_request('GET', endpoint)
        if not response or 'IdentifierList' not in response:
            return {}
            
        cid = str(response['IdentifierList'].get('CID', [None])[0])
        if not cid:
            return {}
            
        return self.get_compound_by_cid(cid)

    def search_by_cas(self, cas: str) -> Dict[str, Any]:
        """
        Search compound by CAS number.
        
        Args:
            cas: CAS registry number
            
        Returns:
            Compound data
        """
        endpoint = f"compound/fastidentity/cas/{cas}/cids/JSON"
        response = self._make_request('GET', endpoint)
        if not response or 'IdentifierList' not in response:
            return {}
            
        cid = str(response['IdentifierList'].get('CID', [None])[0])
        if not cid:
            return {}
            
        return self.get_compound_by_cid(cid)

    def search_by_inchikey(self, inchikey: str) -> Dict[str, Any]:
        """
        Search compound by InChIKey.
        
        Args:
            inchikey: InChIKey identifier
            
        Returns:
            Compound data
        """
        endpoint = f"compound/inchikey/{inchikey}/cids/JSON"
        response = self._make_request('GET', endpoint)
        if not response or 'IdentifierList' not in response:
            return {}
            
        cid = str(response['IdentifierList'].get('CID', [None])[0])
        if not cid:
            return {}
            
        return self.get_compound_by_cid(cid)

    def get_compound_by_cid(self, cid: str) -> Dict[str, Any]:
        """
        Get compound data by CID.
        
        Args:
            cid: PubChem CID
            
        Returns:
            Compound data
        """
        # Get basic properties
        props_endpoint = f"compound/cid/{cid}/property/IUPACName,MolecularWeight,InChI,InChIKey/JSON"
        props_response = self._make_request('GET', props_endpoint)
        
        # Get synonyms
        synonyms_endpoint = f"compound/cid/{cid}/synonyms/JSON"
        synonyms_response = self._make_request('GET', synonyms_endpoint)
        
        # Get computed properties
        computed_endpoint = f"compound/cid/{cid}/property/XLogP,TPSA/JSON"
        computed_response = self._make_request('GET', computed_endpoint)
        
        data = {}
        
        # Process properties
        if props_response and 'PropertyTable' in props_response:
            props = props_response['PropertyTable'].get('Properties', [])
            if props:
                data.update(props[0])
                
        # Process synonyms
        if synonyms_response and 'InformationList' in synonyms_response:
            info = synonyms_response['InformationList'].get('Information', [])
            if info and 'Synonym' in info[0]:
                data['synonyms'] = info[0]['Synonym']
                
        # Process computed properties
        if computed_response and 'PropertyTable' in computed_response:
            props = computed_response['PropertyTable'].get('Properties', [])
            if props:
                data.update(props[0])
                
        # Add URLs
        data['pubchem_url'] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        
        return data

    def get_compound_properties(
        self,
        cid: str,
        properties: list[str]
    ) -> Dict[str, Any]:
        """
        Get specific compound properties.
        
        Args:
            cid: PubChem CID
            properties: List of property names
            
        Returns:
            Property data
        """
        endpoint = f"compound/cid/{cid}/property/{','.join(properties)}/JSON"
        return self._make_request('GET', endpoint)


class McpError(Exception):
    """Base exception for MCP errors."""
    
    def __init__(self, code: 'ErrorCode', message: str):
        """
        Initialize error.
        
        Args:
            code: Error code
            message: Error message
        """
        self.code = code
        self.message = message
        super().__init__(f"{code.name}: {message}")


class WebSearchClient(APIClient):
    """Client for web search operations using SerpAPI."""
    
    def __init__(self):
        """Initialize web search client."""
        super().__init__("https://serpapi.com", "serp")
        self.api_key = SERP_API_KEY

    def get_search_results(self, query: str, num_results: int = 50) -> List[Dict[str, str]]:
        """
        Get web search results.
        
        Args:
            query: Search query
            num_results: Maximum number of results to return
            
        Returns:
            List of search results with title, link, and snippet
        """
        params = {
            'api_key': self.api_key,
            'engine': 'google',
            'q': query,
            'num': num_results
        }
        
        response = self._make_request('GET', 'search', params=params)
        if not response or 'organic_results' not in response:
            return []
            
        return [
            {
                'title': result.get('title', ''),
                'link': result.get('link', ''),
                'snippet': result.get('snippet', '')
            }
            for result in response['organic_results']
        ]

    def get_search_results_count(self, query: str) -> int:
        """
        Get number of search results for a query.
        
        Args:
            query: Search query
            
        Returns:
            Number of results
        """
        params = {
            'api_key': self.api_key,
            'engine': 'google',
            'q': query,
            'num': 1  # We only need result count
        }
        
        response = self._make_request('GET', 'search', params=params)
        if not response or 'search_information' not in response:
            return 0
            
        return int(response['search_information'].get('total_results', 0))

    def rank_common_names(self, names: list[str]) -> list[tuple[str, int]]:
        """
        Rank common names by search result count.
        
        Args:
            names: List of names to rank
            
        Returns:
            List of (name, result_count) tuples, sorted by count
        """
        results = []
        for name in names:
            count = self.get_search_results_count(name)
            results.append((name, count))
            time.sleep(API_RATE_LIMIT)  # Respect rate limits
            
        return sorted(results, key=lambda x: x[1], reverse=True)


class PubMedClient(APIClient):
    """Client for PubMed API."""
    
    def __init__(self):
        """Initialize PubMed client."""
        super().__init__(PUBMED_BASE_URL, "pubmed")

    def get_result_count(self, query: str) -> int:
        """
        Get number of PubMed results for a query.
        
        Args:
            query: Search query
            
        Returns:
            Number of results
        """
        params = {
            'term': query,
            'retmax': 0,  # We only need count
            'format': 'json'
        }
        
        response = self._make_request('GET', 'esearch.fcgi', params=params)
        if not response or 'esearchresult' not in response:
            return 0
            
        return int(response['esearchresult'].get('count', 0))

    def get_abstract(self, pmid: str) -> Optional[str]:
        """
        Get abstract text for a PubMed ID.
        
        Args:
            pmid: PubMed ID
            
        Returns:
            Abstract text if available, None otherwise
        """
        params = {
            'id': pmid,
            'rettype': 'abstract',
            'retmode': 'text'
        }
        
        try:
            response = self._make_request('GET', 'efetch.fcgi', params=params)
            if response and isinstance(response, str):
                return response
            return None
        except Exception as e:
            self.logger.warning(f"Error fetching abstract for {pmid}: {str(e)}")
            return None

    def get_relevant_pmids(self, query: str, max_results: int = 5) -> List[str]:
        """
        Get most relevant PubMed IDs for a query.
        
        Args:
            query: Search query
            max_results: Maximum number of results to return
            
        Returns:
            List of PubMed IDs
        """
        params = {
            'term': query,
            'retmax': max_results,
            'format': 'json'
        }
        
        response = self._make_request('GET', 'esearch.fcgi', params=params)
        if not response or 'esearchresult' not in response:
            return []
            
        return response['esearchresult'].get('idlist', [])

    def get_binding_relevance(self, compound: str, target: str) -> int:
        """
        Get relevance score for compound-target binding.
        
        Args:
            compound: Compound name/identifier
            target: Target name/identifier
            
        Returns:
            Number of PubMed results
        """
        query = f'"{compound}"[Title/Abstract] AND "{target}"[Title/Abstract] AND binding[Title/Abstract]'
        return self.get_result_count(query)


class ErrorCode(Enum):
    """Error codes for MCP operations."""
    CircuitBreakerOpen = 'circuit_breaker_open'
    APIError = 'api_error'
    CacheError = 'cache_error'
    ValidationError = 'validation_error'
    ConfigError = 'config_error'
