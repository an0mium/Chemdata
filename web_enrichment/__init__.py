"""Web scraping and data enrichment from various sources.

This module provides functionality for:
1. Web scraping from multiple chemical databases
2. Data enrichment from regulatory sources
3. Chemical name and identifier extraction
4. Rate-limited and cached web requests
5. Content validation and parsing
"""

from typing import Dict, Any, List, Optional
import os

from logger import LogManager
from .http_client import HttpClient
from .name_utils import extract_identifiers, clean_name
from .data_sources import (
    PubChemClient,
    ChEMBLClient,
    RegulatoryClient,
    CommunityClient,
    SwissClient,
)
from .data_sources.web_search import WebSearchClient
from .llm_utils import analyze_content_with_llm, extract_patent_compound

logger = LogManager().get_logger("web_enrichment")


class WebEnrichment:
    """Handles web scraping and data enrichment from various sources."""

    # Known API domains to exclude from general web search
    EXCLUDED_DOMAINS = [
        "pubchem.ncbi.nlm.nih.gov",
        "ebi.ac.uk",
        "psychonautwiki.org",
        "erowid.org",
        "deadiversion.usdoj.gov",
        "emcdda.europa.eu",
        "who.int",
        "bindingdb.org",
        "drugbank.ca",
        "zinc.docking.org",
        "chemspider.com",
        "isomerdesign.com",
    ]

    # Additional search domains for chemical information
    CHEMICAL_DOMAINS = [
        "site:patents.google.com",
        "site:wikipedia.org",
        "site:chemspider.com",
        "site:commonchemistry.cas.org",
        "site:drugbank.ca",
        "site:isomerdesign.com/Cdsa",
        "site:zinc.docking.org",
    ]

    def __init__(self):
        """Initialize web enrichment processor."""
        # Initialize HTTP client
        self.http = HttpClient()

        # Initialize data source clients
        self.pubchem = PubChemClient(self.http)
        self.chembl = ChEMBLClient(self.http)
        self.regulatory = RegulatoryClient(self.http)
        self.community = CommunityClient(self.http)
        self.web_search = WebSearchClient(self.http)
        self.swiss = SwissClient(self.http)

    async def use_mcp_tool(
        self, server_name: str, tool_name: str, arguments: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Use an MCP tool."""
        from cline_utils import use_mcp_tool

        try:
            result = use_mcp_tool(server_name, tool_name, arguments)
            if isinstance(result, dict):
                return result
            return {"error": "Invalid response format"}
        except Exception as e:
            logger.error(f"Error using MCP tool {tool_name}: {str(e)}")
            return {"error": str(e)}

    def get_pubchem_data(
        self,
        cas: Optional[str] = None,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        """Get compound data from PubChem."""
        try:
            logger.info("Fetching PubChem data...")
            data = self.pubchem.get_compound_data(
                cas=cas, smiles=smiles, inchi=inchi, include_bioassays=True
            )
            if data:
                logger.info("Successfully retrieved PubChem data")
                return data
            else:
                logger.warning("No PubChem data found")
                return None
        except Exception as e:
            logger.error(f"Error getting PubChem data: {str(e)}")
            return None

    def get_reference_urls(
        self,
        compound_name: str,
        llm_api_key: Optional[str] = None,
        search_patents: bool = True,
    ) -> Dict[str, Any]:
        """
        Get reference URLs and extracted data for various sources.

        Args:
            compound_name: Name of compound to search for
            llm_api_key: Optional API key for LLM service
            search_patents: Whether to search patents

        Returns:
            Dictionary containing:
            - urls: Dictionary of source names and URLs
            - extracted_data: Dictionary of data extracted by LLM for each source
        """
        result = {"urls": {}, "extracted_data": {}}

        # Extract identifiers
        clean_name_str, chembl_id, cas_number, patent_id = extract_identifiers(
            compound_name
        )
        compound_data = {
            "name": clean_name_str,
            "cas": cas_number,
            "chembl_id": chembl_id,
            "patent_id": patent_id,
        }

        # Get ChEMBL URL and data
        if chembl_id:
            result["urls"]["chembl_url"] = self.chembl.get_compound_url(chembl_id)
            chembl_data = self.chembl.get_compound_data(chembl_id)
            if chembl_data and llm_api_key:
                result["extracted_data"]["chembl"] = analyze_content_with_llm(
                    str(chembl_data), compound_data, llm_api_key
                )

        # Get community source URLs and data
        community_urls = self.community.get_urls(clean_name_str, cas_number)
        result["urls"].update(community_urls)
        if llm_api_key:
            for source, url in community_urls.items():
                content = self.community.get_content(url)
                if content:
                    result["extracted_data"][source] = analyze_content_with_llm(
                        content, compound_data, llm_api_key
                    )

        # Get PubChem data
        pubchem_data = self.pubchem.get_compound_data(
            cas=cas_number,
            smiles=None,  # Add SMILES if available
            inchi=None,  # Add InChI if available
        )
        if pubchem_data:
            result["urls"]["pubchem_url"] = pubchem_data["url"]
            if llm_api_key:
                result["extracted_data"]["pubchem"] = analyze_content_with_llm(
                    str(pubchem_data["data"]), compound_data, llm_api_key
                )

        # Search patents if requested
        if search_patents and llm_api_key:
            # Build patent search query
            patent_terms = [clean_name_str]
            if cas_number:
                patent_terms.append(cas_number)
            if chembl_id:
                patent_terms.append(chembl_id)
            patent_query = " OR ".join(f'"{term}"' for term in patent_terms)

            # Search patents
            patent_results = self.web_search.search_patents(
                patent_query, llm_api_key, compound_data
            )

            # Add patent results
            if patent_results:
                result["urls"]["patents"] = patent_results["urls"]
                result["extracted_data"]["patents"] = patent_results["extracted_data"]

        # Search other web sources if LLM API key is provided
        if llm_api_key:
            # Build search queries for different domains
            for domain in self.CHEMICAL_DOMAINS:
                query = f'{domain} ("{clean_name_str}"'
                if cas_number:
                    query += f' OR "{cas_number}"'
                if chembl_id:
                    query += f' OR "{chembl_id}"'
                query += ")"

                # Get and analyze web search results
                web_results = self.web_search.search_and_analyze(
                    query,
                    llm_api_key,
                    compound_data,
                    excluded_domains=self.EXCLUDED_DOMAINS,
                )

                # Add results
                source = domain.split(":")[1].split(".")[0]  # Extract domain name
                if web_results["urls"]:
                    result["urls"][f"{source}_url"] = web_results["urls"][0]
                if web_results["extracted_data"]:
                    result["extracted_data"][source] = web_results["extracted_data"]

        return result

    def get_common_names(
        self,
        compound_name: str,
        chembl_id: Optional[str] = None,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
        search_patents: bool = True,
    ) -> List[Dict[str, Any]]:
        """
        Get common names for a compound sorted by web search relevance.

        Args:
            compound_name: Name of compound to search for
            chembl_id: Optional ChEMBL ID
            smiles: Optional SMILES string
            inchi: Optional InChI string
            search_patents: Whether to search patents

        Returns:
            List of dictionaries containing name information
        """
        names = []
        clean_name_str, extracted_chembl_id, cas_number, patent_id = (
            extract_identifiers(compound_name)
        )
        chembl_id = chembl_id or extracted_chembl_id

        # Get ChEMBL names
        if chembl_id:
            names.extend(self.chembl.get_compound_names(chembl_id))

        # Get PubChem names
        pubchem_data = self.pubchem.get_compound_data(
            cas=cas_number, smiles=smiles, inchi=inchi
        )
        if pubchem_data:
            names.extend(self.pubchem.get_compound_names(pubchem_data["cid"]))

        # Get community source names
        community_names = self.community.get_compound_names(clean_name_str, cas_number)
        if community_names:
            names.extend(community_names)

        # Search patents if requested
        if search_patents and patent_id:
            patent_names = self.web_search.get_patent_names(patent_id)
            if patent_names:
                names.extend(patent_names)

        # Add identifiers with high priority
        if cas_number:
            names.append({"name": cas_number, "source": "CAS", "relevance": 100})
        if chembl_id:
            names.append({"name": chembl_id, "source": "ChEMBL", "relevance": 95})

        # Remove duplicates while preserving order
        seen = set()
        unique_names = []
        for name_data in names:
            clean = clean_name(name_data["name"])
            if clean and clean not in seen:  # Skip empty names
                seen.add(clean)
                unique_names.append(name_data)

        # Sort by relevance score
        unique_names.sort(key=lambda x: x.get("relevance", 0), reverse=True)
        return unique_names[:3]  # Return top 3

    def get_legal_status(self, compound_name: str) -> Dict[str, Any]:
        """
        Get legal status information from various sources.

        Args:
            compound_name: Name of compound to search for

        Returns:
            Dictionary containing legal status information
        """
        clean_name_str, _, cas_number, _ = extract_identifiers(compound_name)

        # Get regulatory data
        regulatory_data = self.regulatory.get_legal_status(clean_name_str, cas_number)

        # Get community data
        community_data = self.community.get_legal_status(clean_name_str, cas_number)

        # Merge data
        merged = {"scheduling": [], "sources": set()}

        if regulatory_data:
            merged["scheduling"].extend(regulatory_data.get("scheduling", []))
            merged["sources"].update(regulatory_data.get("sources", []))

        if community_data:
            merged["scheduling"].extend(community_data.get("scheduling", []))
            merged["sources"].update(community_data.get("sources", []))

        # Remove duplicates while preserving order
        seen = set()
        unique_scheduling = []
        for schedule in merged["scheduling"]:
            key = (schedule["jurisdiction"], schedule["schedule"])
            if key not in seen:
                seen.add(key)
                unique_scheduling.append(schedule)

        merged["scheduling"] = unique_scheduling
        merged["sources"] = list(merged["sources"])

        return merged

    def get_pharmacology(
        self, compound_name: str, smiles: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Get pharmacological information from various sources.

        Args:
            compound_name: Name of compound to search for

        Returns:
            Dictionary containing pharmacological information
        """
        info = {
            "mechanism_of_action": [],
            "primary_targets": [],
            "metabolism": [],
            "toxicity": [],
            "pharmacokinetics": [],  # Added pharmacokinetics
            "clinical_data": [],  # Added clinical data
            "drug_interactions": [],  # Added drug interactions
            "contraindications": [],  # Added contraindications
            "sources": [],
        }

        clean_name_str, chembl_id, cas_number, _ = extract_identifiers(compound_name)

        # Get Swiss* data if SMILES is available
        if smiles:
            # Get target predictions
            target_data = self.swiss.get_target_predictions(smiles)
            if target_data:
                info["primary_targets"].extend(
                    [
                        f"{pred['target']} (probability: {pred['probability']:.2f})"
                        for pred in target_data["predictions"][:5]  # Top 5 predictions
                    ]
                )
                info["sources"].append("SwissTargetPrediction")

            # Get ADME properties
            adme_data = self.swiss.get_adme_properties(smiles)
            if adme_data:
                # Add pharmacokinetics data
                pk_info = []
                if adme_data["absorption"]["gi_absorption"]:
                    pk_info.append(
                        f"GI absorption: {adme_data['absorption']['gi_absorption']}"
                    )
                if adme_data["absorption"]["bbb_permeant"]:
                    pk_info.append("BBB permeant")
                if adme_data["absorption"]["pgp_substrate"]:
                    pk_info.append("P-gp substrate")

                # Add metabolism data
                cyp_info = []
                for cyp, is_inhibitor in adme_data["metabolism"][
                    "cyp_inhibition"
                ].items():
                    if is_inhibitor:
                        cyp_info.append(f"{cyp.upper()} inhibitor")

                if pk_info:
                    info["pharmacokinetics"].extend(pk_info)
                if cyp_info:
                    info["metabolism"].extend(cyp_info)

                info["sources"].append("SwissADME")

            # Get similar compounds
            similar_data = self.swiss.search_similar_compounds(
                smiles, similarity_threshold=0.7, max_results=5
            )
            if similar_data:
                info["similar_compounds"] = [
                    f"{cmpd['name']} (similarity: {cmpd['similarity']:.2f})"
                    for cmpd in similar_data["similar_compounds"]
                ]
                info["sources"].append("SwissSimilarity")

        # Try ChEMBL first if available
        if chembl_id:
            chembl_info = self.chembl.get_pharmacology(chembl_id)
            if chembl_info:
                for key in info:
                    if key in chembl_info and isinstance(chembl_info[key], list):
                        info[key].extend(chembl_info[key])
                info["sources"].append("ChEMBL")

        # Get PubChem data
        pubchem_info = self.pubchem.get_pharmacology(clean_name_str, cas_number)
        if pubchem_info:
            for key in info:
                if key in pubchem_info and isinstance(pubchem_info[key], list):
                    info[key].extend(pubchem_info[key])
            info["sources"].append("PubChem")

        # Get community data
        community_info = self.community.get_pharmacology(clean_name_str, cas_number)
        if community_info:
            for key in info:
                if key in community_info and isinstance(community_info[key], list):
                    info[key].extend(community_info[key])
            info["sources"].extend(community_info.get("sources", []))

        # Remove duplicates while preserving order
        for key in info:
            if isinstance(info[key], list):
                info[key] = list(dict.fromkeys(info[key]))

        return info

    def fill_empty_fields(
        self, tsv_path: str, output_path: str, llm_api_key: str
    ) -> None:
        """
        Fill empty fields in TSV file using LLM analysis of reference URLs.

        Args:
            tsv_path: Path to input TSV file
            output_path: Path to write enriched TSV file
            llm_api_key: API key for LLM service
        """
        try:
            # Read TSV file
            with open(tsv_path, "r") as f:
                lines = f.readlines()

            # Parse header and data
            header = lines[0].strip().split("\t")
            data = [line.strip().split("\t") for line in lines[1:]]

            # Process each compound
            enriched_data = []
            for row in data:
                compound_dict = dict(zip(header, row))

                # Get reference URLs and extracted data
                name = compound_dict.get("name", "")
                if name:
                    result = self.get_reference_urls(name, llm_api_key)

                    # Update empty fields with extracted data
                    if result["extracted_data"]:
                        for source, extracted in result["extracted_data"].items():
                            for field, value in extracted.items():
                                # Only fill if current value is empty
                                if (
                                    field in compound_dict
                                    and not compound_dict[field].strip()
                                ):
                                    compound_dict[field] = value

                enriched_data.append(compound_dict)

            # Write enriched TSV
            with open(output_path, "w") as f:
                # Write header
                f.write("\t".join(header) + "\n")

                # Write enriched data
                for row_dict in enriched_data:
                    row = [row_dict.get(col, "") for col in header]
                    f.write("\t".join(row) + "\n")

        except Exception as e:
            logger.error(f"Error filling empty fields: {str(e)}")
