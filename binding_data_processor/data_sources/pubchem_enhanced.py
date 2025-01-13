"""PubChem API client.

This module provides an enhanced client for accessing the PubChem API with:
- Compound lookup by name, CAS, InChIKey
- Property retrieval
- Circuit breaker pattern
- Metrics collection
- Enhanced error handling
- Caching
- Rate limiting
"""

from typing import Dict, Any, Optional, List
from urllib.parse import quote
from crawl4ai import AsyncWebCrawler, BrowserConfig
from crawl4ai.extraction_strategy import JsonCssExtractionStrategy, LLMExtractionStrategy

from ..web_enrichment.base_client_enhanced import BaseWebClientEnhanced
from ..models.compound import Compound
from ..pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig
from ..web_enrichment.llm_utils import extract_chemical_info


class PubChemClientEnhanced(BaseWebClientEnhanced):
    """Enhanced client for PubChem API."""

    def __init__(self, llm_provider: str = "ollama/llama2", api_token: Optional[str] = None, **kwargs):
        """Initialize PubChem client with circuit breaker and metrics.

        Args:
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            **kwargs: Additional arguments passed to BaseWebClientEnhanced
        """
        # Configure circuit breaker - PubChem has good uptime so we can be lenient
        circuit_config = CircuitBreakerConfig(
            failure_threshold=10,
            reset_timeout=300,
        )

        # Initialize with name and circuit breaker config
        super().__init__(name="pubchem", circuit_config=circuit_config, **kwargs)
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

        # Configure crawl4ai for web scraping fallback
        self.crawler_config = BrowserConfig(
            javascript=BrowserConfig.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    "#Summary",
                    "#Structures",
                    "#Chemical-Properties",
                ],
                stealth_mode=True,
            ),
            screenshot=BrowserConfig.Screenshot(enabled=True, full_page=True),
            extraction=BrowserConfig.Extraction(
                llm=BrowserConfig.LLM(
                    provider=llm_provider,
                    api_token=api_token,
                    prompts={
                        "properties": "Extract chemical properties:",
                        "structure": "Extract structural information:",
                        "names": "Extract chemical names and identifiers:",
                    },
                ),
                css={
                    "iupac_name": "#Summary .IUPACName",
                    "molecular_weight": "#Chemical-Properties .Molecular-Weight",
                    "inchi": "#Structures .InChI",
                    "inchikey": "#Structures .InChIKey",
                    "xlogp": "#Chemical-Properties .XLogP",
                    "tpsa": "#Chemical-Properties .TPSA",
                    "synonyms": "#Names-and-Identifiers .Synonym",
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

    async def _scrape_compound_data(self, identifier: str, id_type: str = "name") -> Optional[Dict[str, Any]]:
        """Scrape compound data from PubChem website.

        Args:
            identifier: Compound identifier (name, CAS, InChIKey)
            id_type: Type of identifier (name, cas, inchikey)

        Returns:
            Compound data or None if scraping failed
        """
        try:
            # Construct URL based on identifier type
            if id_type == "cas":
                url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{quote(identifier)}#section=CAS"
            elif id_type == "inchikey":
                url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{quote(identifier)}"
            else:
                url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{quote(identifier)}#section=Names-and-Identifiers"

            # Scrape with crawl4ai
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.crawler_config)

                if not result.success:
                    self.logger.error(f"Failed to scrape PubChem page for {identifier}")
                    return None

                # Extract data using both CSS and LLM
                content = result.extracted_content[0]

                # Extract chemical info using LLM
                chemical_info = extract_chemical_info(
                    title=content.get("title", ""),
                    description=content.get("description", ""),
                )

                # Combine CSS and LLM extracted data
                data = {
                    "IUPACName": content.get("iupac_name"),
                    "MolecularWeight": content.get("molecular_weight"),
                    "InChI": content.get("inchi"),
                    "InChIKey": content.get("inchikey"),
                    "XLogP": content.get("xlogp"),
                    "TPSA": content.get("tpsa"),
                    "synonyms": content.get("synonyms", []),
                    **chemical_info,
                }

                return data

        except Exception as e:
            self.logger.error(f"Error scraping PubChem data: {str(e)}")
            return None

    async def get_compound_by_name(self, name: str) -> Dict[str, Any]:
        """Get compound data by name.

        Args:
            name: Compound name

        Returns:
            Compound data
        """
        try:

            def get_cid():
                endpoint = f"compound/name/{quote(name)}/cids/JSON"
                response = self.http.get(f"{self.base_url}/{endpoint}")
                data = response.json()

                if not data or "IdentifierList" not in data:
                    return None

                return str(data["IdentifierList"].get("CID", [None])[0])

            # Try API first
            cid = self.http.circuit.execute(get_cid)
            if cid:
                return self.get_compound_by_cid(cid)

            # Fall back to web scraping if API fails
            self.logger.info(f"Falling back to web scraping for {name}")
            scraped_data = await self._scrape_compound_data(name)
            if scraped_data:
                return scraped_data

        except Exception as e:
            self.logger.error(f"Error getting compound by name '{name}': {str(e)}")
            return {}

    async def search_by_cas(self, cas: str) -> Dict[str, Any]:
        """Search compound by CAS number.

        Args:
            cas: CAS registry number

        Returns:
            Compound data
        """
        try:

            def get_cid():
                endpoint = f"compound/fastidentity/cas/{cas}/cids/JSON"
                response = self.http.get(f"{self.base_url}/{endpoint}")
                data = response.json()

                if not data or "IdentifierList" not in data:
                    return None

                return str(data["IdentifierList"].get("CID", [None])[0])

            # Try API first
            cid = self.http.circuit.execute(get_cid)
            if cid:
                return self.get_compound_by_cid(cid)

            # Fall back to web scraping if API fails
            self.logger.info(f"Falling back to web scraping for CAS {cas}")
            scraped_data = await self._scrape_compound_data(cas, id_type="cas")
            if scraped_data:
                return scraped_data

        except Exception as e:
            self.logger.error(f"Error searching by CAS '{cas}': {str(e)}")
            return {}

    async def search_by_inchikey(self, inchikey: str) -> Dict[str, Any]:
        """Search compound by InChIKey.

        Args:
            inchikey: InChIKey identifier

        Returns:
            Compound data
        """
        try:

            def get_cid():
                endpoint = f"compound/inchikey/{inchikey}/cids/JSON"
                response = self.http.get(f"{self.base_url}/{endpoint}")
                data = response.json()

                if not data or "IdentifierList" not in data:
                    return None

                return str(data["IdentifierList"].get("CID", [None])[0])

            # Try API first
            cid = self.http.circuit.execute(get_cid)
            if cid:
                return self.get_compound_by_cid(cid)

            # Fall back to web scraping if API fails
            self.logger.info(f"Falling back to web scraping for InChIKey {inchikey}")
            scraped_data = await self._scrape_compound_data(inchikey, id_type="inchikey")
            if scraped_data:
                return scraped_data

        except Exception as e:
            self.logger.error(f"Error searching by InChIKey '{inchikey}': {str(e)}")
            return {}

    def get_compound_by_cid(self, cid: str) -> Dict[str, Any]:
        """Get compound data by CID.

        Args:
            cid: PubChem CID

        Returns:
            Compound data
        """
        try:

            def get_properties():
                # Get basic properties
                props_endpoint = f"compound/cid/{cid}/property/IUPACName,MolecularWeight,InChI,InChIKey/JSON"
                props_response = self.http.get(f"{self.base_url}/{props_endpoint}")
                props_data = props_response.json()

                # Get synonyms
                synonyms_endpoint = f"compound/cid/{cid}/synonyms/JSON"
                synonyms_response = self.http.get(f"{self.base_url}/{synonyms_endpoint}")
                synonyms_data = synonyms_response.json()

                # Get computed properties
                computed_endpoint = f"compound/cid/{cid}/property/XLogP,TPSA/JSON"
                computed_response = self.http.get(f"{self.base_url}/{computed_endpoint}")
                computed_data = computed_response.json()

                return {"props": props_data, "synonyms": synonyms_data, "computed": computed_data}

            # Use circuit breaker to get all properties
            response = self.http.circuit.execute(get_properties)
            props_data = response["props"]
            synonyms_data = response["synonyms"]
            computed_data = response["computed"]

            data = {}

            # Process properties
            if props_data and "PropertyTable" in props_data:
                props = props_data["PropertyTable"].get("Properties", [])
                if props:
                    data.update(props[0])

            # Process synonyms
            if synonyms_data and "InformationList" in synonyms_data:
                info = synonyms_data["InformationList"].get("Information", [])
                if info and "Synonym" in info[0]:
                    data["synonyms"] = info[0]["Synonym"]

            # Process computed properties
            if computed_data and "PropertyTable" in computed_data:
                props = computed_data["PropertyTable"].get("Properties", [])
                if props:
                    data.update(props[0])

            # Add URLs
            data["pubchem_url"] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

            return data

        except Exception as e:
            self.logger.error(f"Error getting compound by CID '{cid}': {str(e)}")
            return {}

    def get_compound_properties(self, cid: str, properties: List[str]) -> Dict[str, Any]:
        """Get specific compound properties.

        Args:
            cid: PubChem CID
            properties: List of property names

        Returns:
            Property data
        """
        try:

            def get_props():
                endpoint = f"compound/cid/{cid}/property/{','.join(properties)}/JSON"
                response = self.http.get(f"{self.base_url}/{endpoint}")
                return response.json()

            # Use circuit breaker to get properties
            return self.http.circuit.execute(get_props)

        except Exception as e:
            self.logger.error(f"Error getting properties for CID '{cid}': {str(e)}")
            return {}

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
            try:
                # Try different identifiers
                data = {}
                if compound.name:
                    data = await self.get_compound_by_name(compound.name)
                if not data and compound.cas_number:
                    data = await self.search_by_cas(compound.cas_number)
                if not data and compound.inchikey:
                    data = await self.search_by_inchikey(compound.inchikey)

                if data:
                    # Update compound with PubChem data
                    if "IUPACName" in data:
                        compound.iupac_name = data["IUPACName"]
                    if "MolecularWeight" in data:
                        compound.molecular_weight = float(data["MolecularWeight"])
                    if "InChI" in data:
                        compound.inchi = data["InChI"]
                    if "InChIKey" in data:
                        compound.inchikey = data["InChIKey"]
                    if "XLogP" in data:
                        compound.xlogp = float(data["XLogP"])
                    if "TPSA" in data:
                        compound.tpsa = float(data["TPSA"])
                    if "synonyms" in data:
                        compound.synonyms = data["synonyms"]
                    if "pubchem_url" in data:
                        compound.pubchem_url = data["pubchem_url"]

            except Exception as e:
                self.logger.error(f"Error processing compound {compound.name}: {str(e)}")
                continue

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
        # Try name first
        data = await self.get_compound_by_name(name)
        if data:
            return data

        # Try CAS if provided
        if cas_number:
            data = await self.search_by_cas(cas_number)
            if data:
                return data

        return None
