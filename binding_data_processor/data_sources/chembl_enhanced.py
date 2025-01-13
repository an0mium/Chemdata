"""ChEMBL API client for retrieving binding and activity data.

This module provides an enhanced client for the ChEMBL API with:
- Compound searching by multiple identifiers
- Binding data retrieval
- Target data retrieval
- Psychopharmacological data enrichment
- Circuit breaker pattern
- Robust error handling and retries
- Caching support
- Web scraping fallback when API fails
"""

from typing import Any, Dict, List, Optional, Set
from chembl_webresource_client.new_client import new_client
import time
from functools import lru_cache
from crawl4ai import AsyncWebCrawler, BrowserConfig
from crawl4ai.extraction_strategy import JsonCssExtractionStrategy, LLMExtractionStrategy

from ..web_enrichment.base_client_enhanced import BaseWebClientEnhanced
from ..models.compound import Compound, BindingData
from ..pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig
from ..web_enrichment.llm_utils import extract_chemical_info


class ChEMBLClientEnhanced(BaseWebClientEnhanced):
    """Enhanced client for ChEMBL API using official Python client."""

    # Activity types of interest for psychopharmacology
    ACTIVITY_TYPES = {"Ki", "IC50", "Kd", "EC50", "potency", "activity"}

    # Target types of interest
    TARGET_TYPES = {"SINGLE PROTEIN", "PROTEIN COMPLEX", "PROTEIN FAMILY", "SELECTIVITY GROUP"}

    # Receptor families of interest
    RECEPTOR_FAMILIES = {"5-HT", "NMDA", "GABA", "dopamine", "opioid", "cannabinoid", "sigma", "adrenergic"}

    def __init__(self, llm_provider: str = "ollama/llama2", api_token: Optional[str] = None, **kwargs):
        """Initialize ChEMBL client with caching and circuit breaker.

        Args:
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            **kwargs: Additional arguments passed to BaseWebClientEnhanced
        """
        # Configure circuit breaker
        circuit_config = CircuitBreakerConfig(
            failure_threshold=5,
            reset_timeout=300,
        )

        # Initialize with name and circuit breaker config
        super().__init__(name="chembl", circuit_config=circuit_config, **kwargs)
        self._init_clients()

        # Configure crawl4ai for web scraping fallback
        self.crawler_config = BrowserConfig(
            javascript=BrowserConfig.JavaScript(
                enabled=True,
                wait_for_network=True,
                wait_for_selectors=[
                    "#CompoundInformation",
                    "#MolecularProperties",
                    "#Bioactivity",
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
                        "bioactivity": "Extract bioactivity data:",
                        "targets": "Extract target information:",
                    },
                ),
                css={
                    "compound_name": "#CompoundInformation .pref-name",
                    "molecular_weight": "#MolecularProperties .mol-weight",
                    "alogp": "#MolecularProperties .alogp",
                    "psa": "#MolecularProperties .psa",
                    "ro5_violations": "#MolecularProperties .ro5-violations",
                    "rotatable_bonds": "#MolecularProperties .rotatable-bonds",
                    "aromatic_rings": "#MolecularProperties .aromatic-rings",
                    "hba": "#MolecularProperties .hba",
                    "hbd": "#MolecularProperties .hbd",
                    "bioactivities": "#Bioactivity .activity-data",
                    "targets": "#Bioactivity .target-data",
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

    def _init_clients(self):
        """Initialize ChEMBL API clients with retry and circuit breaker."""
        max_retries = 3
        for attempt in range(max_retries):
            try:
                # Use the circuit breaker from the enhanced base client
                response = self.http.circuit.execute(
                    lambda: {
                        "molecule": new_client.molecule,
                        "activity": new_client.activity,
                        "target": new_client.target,
                        "mechanism": new_client.mechanism,
                    }
                )
                self.molecule = response["molecule"]
                self.activity = response["activity"]
                self.target = response["target"]
                self.mechanism = response["mechanism"]
                return
            except Exception as e:
                self.logger.warning(f"Error initializing ChEMBL clients (attempt {attempt + 1}/{max_retries}): {str(e)}")
                if attempt < max_retries - 1:
                    time.sleep(2**attempt)  # Exponential backoff
                else:
                    raise

    async def _scrape_compound_data(self, identifier: str, id_type: str = "name") -> Optional[Dict[str, Any]]:
        """Scrape compound data from ChEMBL website.

        Args:
            identifier: Compound identifier (name, CHEMBL ID)
            id_type: Type of identifier (name, chembl_id)

        Returns:
            Compound data or None if scraping failed
        """
        try:
            # Construct URL based on identifier type
            if id_type == "chembl_id":
                url = f"https://www.ebi.ac.uk/chembl/compound_report_card/{identifier}/"
            else:
                url = f"https://www.ebi.ac.uk/chembl/g/#search_results/compounds/query={identifier}"

            # Scrape with crawl4ai
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(urls=[url], config=self.crawler_config)

                if not result.success:
                    self.logger.error(f"Failed to scrape ChEMBL page for {identifier}")
                    return None

                # Extract data using both CSS and LLM
                content = result.extracted_content[0]

                # Extract chemical info using LLM
                chemical_info = extract_chemical_info(
                    title=content.get("title", ""),
                    description=content.get("description", ""),
                )

                # Process bioactivity data
                bioactivities = []
                for activity in content.get("bioactivities", []):
                    try:
                        bioactivity = BindingData(
                            target_common_name=activity.get("target_name", "N/A"),
                            target_protein_name=activity.get("protein_name", "N/A"),
                            target_gene_name=activity.get("gene_name", "N/A"),
                            target_organism=activity.get("organism", "N/A"),
                            affinity_value=float(activity.get("value", 0.0)),
                            affinity_type=activity.get("type", "N/A"),
                            affinity_unit=activity.get("unit", "N/A"),
                            affinity_relation=activity.get("relation", "="),
                            assay_description=activity.get("assay_description", "N/A"),
                            assay_organism=activity.get("assay_organism", "N/A"),
                            assay_type=activity.get("assay_type", "N/A"),
                            reference=activity.get("reference", "N/A"),
                            confidence_score=float(activity.get("confidence", 0.0)),
                        )
                        bioactivities.append(bioactivity)
                    except Exception as e:
                        self.logger.warning(f"Error processing bioactivity: {str(e)}")
                        continue

                # Combine CSS and LLM extracted data
                data = {
                    "name": content.get("compound_name"),
                    "molecular_weight": content.get("molecular_weight"),
                    "alogp": content.get("alogp"),
                    "psa": content.get("psa"),
                    "ro5_violations": content.get("ro5_violations"),
                    "rotatable_bonds": content.get("rotatable_bonds"),
                    "aromatic_rings": content.get("aromatic_rings"),
                    "hba": content.get("hba"),
                    "hbd": content.get("hbd"),
                    "binding_data": bioactivities,
                    **chemical_info,
                }

                return data

        except Exception as e:
            self.logger.error(f"Error scraping ChEMBL data: {str(e)}")
            return None

    @lru_cache(maxsize=1000)
    async def search_compound(self, compound: Dict[str, str], use_cache: bool = True) -> Dict[str, Any]:
        """Search for a compound using multiple identifiers.

        Args:
            compound: Dictionary with available identifiers (name, cas, smiles, inchi, etc.)
            use_cache: Whether to use cached results

        Returns:
            Compound data with binding information
        """
        cache_key = f"compound_search_{compound.get('name', '')}_{compound.get('cas', '')}"

        if use_cache and self.cache_dir:
            cached = self.http.get_cache(cache_key)
            if cached:
                return cached

        try:
            # Try API first
            def search_operation():
                # Try exact name match first
                if compound.get("name"):
                    # Try exact match
                    results = list(self.molecule.filter(pref_name__iexact=compound["name"]))
                    if not results:
                        # Try synonym match
                        results = list(self.molecule.filter(molecule_synonyms__synonym__iexact=compound["name"]))
                    if not results:
                        # Try fuzzy name match
                        results = list(self.molecule.filter(pref_name__icontains=compound["name"]))
                    if results:
                        return self.get_compound_by_chembl_id(results[0]["molecule_chembl_id"])

                # Try CAS number
                if compound.get("cas"):
                    results = list(self.molecule.filter(molecule_synonyms__synonym__iexact=compound["cas"]))
                    if results:
                        return self.get_compound_by_chembl_id(results[0]["molecule_chembl_id"])

                # Try structure search if SMILES is available
                if compound.get("smiles"):
                    results = list(self.molecule.filter(molecule_structures__canonical_smiles__flexmatch=compound["smiles"]))
                    if results:
                        return self.get_compound_by_chembl_id(results[0]["molecule_chembl_id"])

                # Try InChI search
                if compound.get("inchi"):
                    results = list(self.molecule.filter(molecule_structures__standard_inchi=compound["inchi"]))
                    if results:
                        return self.get_compound_by_chembl_id(results[0]["molecule_chembl_id"])

                return {}

            data = self.http.circuit.execute(search_operation)
            if data:
                if use_cache and self.cache_dir:
                    self.http.set_cache(cache_key, data)
                return data

            # Fall back to web scraping if API fails
            self.logger.info(f"Falling back to web scraping for {compound.get('name', '')}")
            scraped_data = await self._scrape_compound_data(compound.get("name", ""))
            if scraped_data:
                if use_cache and self.cache_dir:
                    self.http.set_cache(cache_key, scraped_data)
                return scraped_data

            return {}

        except Exception as e:
            self.logger.warning(f"Error searching compound: {str(e)}")
            return {}

    async def get_compound_by_chembl_id(self, chembl_id: str, use_cache: bool = True) -> Dict[str, Any]:
        """Get compound data by ChEMBL ID with enhanced binding data.

        Args:
            chembl_id: ChEMBL ID
            use_cache: Whether to use cached results

        Returns:
            Compound data with binding information
        """
        cache_key = f"compound_data_{chembl_id}"

        if use_cache and self.cache_dir:
            cached = self.http.get_cache(cache_key)
            if cached:
                return cached

        max_retries = 3
        for attempt in range(max_retries):
            try:

                def fetch_operation():
                    # Re-initialize clients if needed
                    if attempt > 0:
                        self._init_clients()

                    # Get molecule details
                    molecule_data = self.molecule.get(chembl_id)

                    # Get binding data with focus on psychoactive targets
                    activities = list(
                        self.activity.filter(
                            molecule_chembl_id=chembl_id,
                            type__in=list(self.ACTIVITY_TYPES),
                            target_type__in=list(self.TARGET_TYPES),
                            relation__in=["=", "<", ">", "<=", ">="],
                            standard_units__isnull=False,
                        ).order_by("-confidence_score")[:20]
                    )

                    # Get mechanism of action data
                    mechanisms = list(self.mechanism.filter(molecule_chembl_id=chembl_id))
                    return {"molecule_data": molecule_data, "activities": activities, "mechanisms": mechanisms}

                # Try API first
                response = self.http.circuit.execute(fetch_operation)
                if response:
                    molecule_data = response["molecule_data"]
                    activities = response["activities"]
                    mechanisms = response["mechanisms"]
                    break
                else:
                    # Fall back to web scraping if API fails
                    self.logger.info(f"Falling back to web scraping for ChEMBL ID {chembl_id}")
                    scraped_data = await self._scrape_compound_data(chembl_id, id_type="chembl_id")
                    if scraped_data:
                        if use_cache and self.cache_dir:
                            self.http.set_cache(cache_key, scraped_data)
                        return scraped_data
                    return {}

            except Exception as e:
                self.logger.warning(f"Error getting compound data (attempt {attempt + 1}/{max_retries}): {str(e)}")
                if attempt < max_retries - 1:
                    time.sleep(2**attempt)  # Exponential backoff
                else:
                    return {}

        data = {"chembl_id": chembl_id, "binding_data": [], "mechanisms": []}

        try:
            if molecule_data:
                data.update(
                    {
                        "name": molecule_data.get("pref_name", ""),
                        "synonyms": [s.get("synonym", "") for s in molecule_data.get("molecule_synonyms", []) if s.get("synonym")],
                        "smiles": molecule_data.get("molecule_structures", {}).get("canonical_smiles"),
                        "inchi": molecule_data.get("molecule_structures", {}).get("standard_inchi"),
                        "inchi_key": molecule_data.get("molecule_structures", {}).get("standard_inchi_key"),
                        "molecular_weight": molecule_data.get("molecule_properties", {}).get("full_mwt"),
                        "alogp": molecule_data.get("molecule_properties", {}).get("alogp"),
                        "psa": molecule_data.get("molecule_properties", {}).get("psa"),
                        "ro5_violations": molecule_data.get("molecule_properties", {}).get("num_ro5_violations"),
                        "rotatable_bonds": molecule_data.get("molecule_properties", {}).get("rtb"),
                        "aromatic_rings": molecule_data.get("molecule_properties", {}).get("aromatic_rings"),
                        "hba": molecule_data.get("molecule_properties", {}).get("hba"),
                        "hbd": molecule_data.get("molecule_properties", {}).get("hbd"),
                        "chembl_url": f"https://www.ebi.ac.uk/chembl/compound/{chembl_id}",
                    }
                )

            # Process binding data
            for activity in activities:
                try:
                    # Check if target is a receptor of interest
                    target_name = activity.get("target_pref_name", "").lower()
                    if not any(family.lower() in target_name for family in self.RECEPTOR_FAMILIES):
                        continue

                    binding = BindingData(
                        target_common_name=activity.get("target_pref_name") or "N/A",
                        target_protein_name=(activity.get("target_components", [{}]) or [{}])[0].get("protein_name") or "N/A",
                        target_gene_name=(activity.get("target_components", [{}]) or [{}])[0].get("gene_name") or "N/A",
                        target_organism=activity.get("target_organism", "N/A"),
                        affinity_value=float(activity.get("value") or 0.0),
                        affinity_type=activity.get("type") or "N/A",
                        affinity_unit=activity.get("units") or "N/A",
                        affinity_relation=activity.get("relation") or "=",
                        assay_description=activity.get("assay_description") or "N/A",
                        assay_organism=activity.get("assay_organism") or "N/A",
                        assay_type=activity.get("assay_type") or "N/A",
                        reference=activity.get("document_chembl_id") or "N/A",
                        confidence_score=float(activity.get("confidence_score") or 0.0),
                    )
                    data["binding_data"].append(binding)
                except Exception as e:
                    self.logger.warning(f"Error processing activity data: {str(e)}")
                    continue

            # Process mechanism data
            for mech in mechanisms:
                try:
                    mechanism = {
                        "mechanism_of_action": mech.get("mechanism_of_action"),
                        "target_name": mech.get("target_name"),
                        "action_type": mech.get("action_type"),
                        "binding_site": mech.get("binding_site_name"),
                        "mechanism_refs": mech.get("mechanism_refs", []),
                        "confidence_score": float(mech.get("confidence_score") or 0.0),
                    }
                    data["mechanisms"].append(mechanism)
                except Exception as e:
                    self.logger.warning(f"Error processing mechanism data: {str(e)}")
                    continue

        except Exception as e:
            self.logger.warning(f"Error processing molecule data: {str(e)}")

        if use_cache and self.cache_dir:
            self.http.set_cache(cache_key, data)

        return data

    @lru_cache(maxsize=100)
    async def search_targets(self, query: str, target_type: Optional[str] = None) -> List[Dict[str, Any]]:
        """Search for protein targets with filtering.

        Args:
            query: Search query
            target_type: Optional target type filter

        Returns:
            List of matching targets
        """
        try:

            def search_operation():
                filters = {"pref_name__icontains": query}
                if target_type:
                    filters["target_type"] = target_type
                return list(self.target.filter(**filters))

            # Try API first
            data = self.http.circuit.execute(search_operation)
            if data:
                return data

            # Fall back to web scraping if API fails
            self.logger.info(f"Falling back to web scraping for target search: {query}")
            async with AsyncWebCrawler() as crawler:
                url = f"https://www.ebi.ac.uk/chembl/g/#search_results/targets/query={query}"
                result = await crawler.arun(urls=[url], config=self.crawler_config)
                if result.success:
                    return result.extracted_content
                return []

        except Exception as e:
            self.logger.error(f"Error searching targets: {str(e)}")
            return []

    @lru_cache(maxsize=100)
    async def get_target_by_chembl_id(self, target_id: str) -> Dict[str, Any]:
        """Get target data by ChEMBL target ID.

        Args:
            target_id: ChEMBL target ID

        Returns:
            Target data
        """
        try:
            # Try API first
            data = self.http.circuit.execute(lambda: self.target.get(target_id))
            if data:
                return data

            # Fall back to web scraping if API fails
            self.logger.info(f"Falling back to web scraping for target {target_id}")
            async with AsyncWebCrawler() as crawler:
                url = f"https://www.ebi.ac.uk/chembl/target_report_card/{target_id}/"
                result = await crawler.arun(urls=[url], config=self.crawler_config)
                if result.success:
                    return result.extracted_content[0]
                return {}

        except Exception as e:
            self.logger.error(f"Error getting target: {str(e)}")
            return {}

    async def process_compounds(
        self,
        compounds: List[Compound],
        skip_predictions: bool = False,
        use_cache: bool = True,
    ) -> None:
        """Process list of compounds to enrich with ChEMBL data.

        Args:
            compounds: List of compounds to process
            skip_predictions: Whether to skip ML predictions
            use_cache: Whether to use cached results
        """
        for compound in compounds:
            try:
                # Search ChEMBL
                chembl_data = await self.search_compound(
                    {
                        "name": compound.name,
                        "cas": compound.cas_number,
                        "smiles": compound.smiles,
                        "inchi": compound.inchi,
                    },
                    use_cache=use_cache,
                )

                if chembl_data:
                    # Update compound with ChEMBL data
                    compound.chembl_id = chembl_data.get("chembl_id")
                    compound.synonyms.extend(chembl_data.get("synonyms", []))
                    compound.binding_data.extend(chembl_data.get("binding_data", []))
                    compound.mechanisms.extend(chembl_data.get("mechanisms", []))

                    # Update chemical properties
                    if not compound.molecular_weight:
                        compound.molecular_weight = chembl_data.get("molecular_weight")
                    if not compound.alogp:
                        compound.alogp = chembl_data.get("alogp")
                    if not compound.psa:
                        compound.psa = chembl_data.get("psa")
                    if not compound.ro5_violations:
                        compound.ro5_violations = chembl_data.get("ro5_violations")
                    if not compound.rotatable_bonds:
                        compound.rotatable_bonds = chembl_data.get("rotatable_bonds")
                    if not compound.aromatic_rings:
                        compound.aromatic_rings = chembl_data.get("aromatic_rings")
                    if not compound.hba:
                        compound.hba = chembl_data.get("hba")
                    if not compound.hbd:
                        compound.hbd = chembl_data.get("hbd")

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
        # Search ChEMBL with API first, then fallback to web scraping
        data = await self.search_compound(
            {
                "name": name,
                "cas": cas_number,
            },
            use_cache=use_cache,
        )

        if not data:
            return None

        return data
