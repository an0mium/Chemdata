"""Patent API client.

This module provides a client for accessing patent databases with:
- Patent search by compound name, CAS, SMILES
- Patent data extraction
- Chemical structure search
- Caching and rate limiting
- LLM-based text analysis
"""

import logging
from typing import Dict, Any, List, Optional
from datetime import datetime

from ...models.compound import Compound
from ..base_client import BaseWebClient
from ..llm_utils import extract_chemical_info, analyze_patent_text
from ...pipeline.infrastructure.circuit_breaker import CircuitConfig


class PatentClient(BaseWebClient):
    """Client for patent databases."""

    def __init__(
        self,
        google_patents_key: Optional[str] = None,
        lens_api_key: Optional[str] = None,
        circuit_config: Optional[CircuitConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize patent client.

        Args:
            google_patents_key: Optional Google Patents API key
            lens_api_key: Optional Lens.org API key
            circuit_config: Optional circuit breaker config
            logger: Optional logger instance
        """
        super().__init__(circuit_config=circuit_config, logger=logger)
        self.google_patents_key = google_patents_key
        self.lens_api_key = lens_api_key

        # Base URLs
        self.google_patents_url = "https://patents.google.com/api/v1"
        self.lens_url = "https://api.lens.org/patent/search"

    def search_patents(
        self,
        query: str,
        chemical_structure: Optional[str] = None,
        from_date: Optional[str] = None,
        to_date: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Search patents by text query and/or chemical structure.

        Args:
            query: Text search query
            chemical_structure: Optional SMILES or InChI structure
            from_date: Optional start date (YYYY-MM-DD)
            to_date: Optional end date (YYYY-MM-DD)

        Returns:
            List of patent data dictionaries
        """
        results = []

        # Search Google Patents
        if self.google_patents_key:
            try:
                google_results = self._search_google_patents(
                    query=query,
                    chemical_structure=chemical_structure,
                    from_date=from_date,
                    to_date=to_date,
                )
                results.extend(google_results)
            except Exception as e:
                self.logger.error(f"Google Patents error: {str(e)}")

        # Search Lens.org
        if self.lens_api_key:
            try:
                lens_results = self._search_lens(
                    query=query,
                    chemical_structure=chemical_structure,
                    from_date=from_date,
                    to_date=to_date,
                )
                results.extend(lens_results)
            except Exception as e:
                self.logger.error(f"Lens.org error: {str(e)}")

        # Deduplicate results
        seen_ids = set()
        unique_results = []
        for result in results:
            if result["patent_number"] not in seen_ids:
                seen_ids.add(result["patent_number"])
                unique_results.append(result)

        return unique_results

    def _search_google_patents(
        self,
        query: str,
        chemical_structure: Optional[str] = None,
        from_date: Optional[str] = None,
        to_date: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Search Google Patents.

        Args:
            query: Text search query
            chemical_structure: Optional SMILES or InChI structure
            from_date: Optional start date
            to_date: Optional end date

        Returns:
            List of patent data dictionaries
        """
        # Build query
        params = {
            "q": query,
            "key": self.google_patents_key,
        }
        if from_date:
            params["from"] = from_date
        if to_date:
            params["to"] = to_date
        if chemical_structure:
            params["structure"] = chemical_structure

        # Make request
        response = self.http.get(
            f"{self.google_patents_url}/query",
            params=params,
        )
        data = response.json()

        # Process results
        results = []
        for patent in data.get("patents", []):
            # Extract chemical info using LLM
            chemical_info = extract_chemical_info(
                patent.get("description", ""),
                patent.get("claims", ""),
            )

            # Analyze patent text
            analysis = analyze_patent_text(
                title=patent.get("title", ""),
                abstract=patent.get("abstract", ""),
                description=patent.get("description", ""),
                claims=patent.get("claims", ""),
            )

            results.append(
                {
                    "patent_number": patent["publication_number"],
                    "title": patent.get("title"),
                    "abstract": patent.get("abstract"),
                    "filing_date": patent.get("filing_date"),
                    "publication_date": patent.get("publication_date"),
                    "assignee": patent.get("assignee"),
                    "inventors": patent.get("inventors", []),
                    "url": f"https://patents.google.com/patent/{patent['publication_number']}",
                    "chemical_info": chemical_info,
                    "analysis": analysis,
                }
            )

        return results

    def _search_lens(
        self,
        query: str,
        chemical_structure: Optional[str] = None,
        from_date: Optional[str] = None,
        to_date: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Search Lens.org patents.

        Args:
            query: Text search query
            chemical_structure: Optional SMILES or InChI structure
            from_date: Optional start date
            to_date: Optional end date

        Returns:
            List of patent data dictionaries
        """
        # Build query
        body = {
            "query": query,
            "token": self.lens_api_key,
        }
        if from_date:
            body["from_date"] = from_date
        if to_date:
            body["to_date"] = to_date
        if chemical_structure:
            body["structure"] = chemical_structure

        # Make request
        response = self.http.post(self.lens_url, json=body)
        data = response.json()

        # Process results
        results = []
        for patent in data.get("data", []):
            # Extract chemical info using LLM
            chemical_info = extract_chemical_info(
                patent.get("description", ""),
                patent.get("claims", ""),
            )

            # Analyze patent text
            analysis = analyze_patent_text(
                title=patent.get("title", ""),
                abstract=patent.get("abstract", ""),
                description=patent.get("description", ""),
                claims=patent.get("claims", ""),
            )

            results.append(
                {
                    "patent_number": patent["lens_id"],
                    "title": patent.get("title"),
                    "abstract": patent.get("abstract"),
                    "filing_date": patent.get("date_filed"),
                    "publication_date": patent.get("date_published"),
                    "assignee": patent.get("assignee"),
                    "inventors": patent.get("inventors", []),
                    "url": f"https://www.lens.org/{patent['lens_id']}",
                    "chemical_info": chemical_info,
                    "analysis": analysis,
                }
            )

        return results

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
            # Build search query
            query_parts = []
            if compound.name:
                query_parts.append(compound.name)
            if compound.cas_number:
                query_parts.append(compound.cas_number)
            if compound.inchikey:
                query_parts.append(compound.inchikey)

            if not query_parts:
                continue

            query = " OR ".join(query_parts)

            # Search patents
            try:
                patents = self.search_patents(
                    query=query,
                    chemical_structure=compound.smiles,
                )

                # Update compound with patent data
                compound.patents = patents
                compound.patent_count = len(patents)
                compound.patent_dates = [p["publication_date"] for p in patents if p.get("publication_date")]
                compound.patent_assignees = list(set(p["assignee"] for p in patents if p.get("assignee")))

                # Add metadata
                if not hasattr(compound, "enrichment_metadata"):
                    compound.enrichment_metadata = {}
                compound.enrichment_metadata["patent_search"] = {
                    "timestamp": datetime.now().isoformat(),
                    "query": query,
                    "total_results": len(patents),
                }

            except Exception as e:
                self.logger.error(f"Error processing patents for {compound.name}: {str(e)}")

    def get_patent_data(
        self,
        patent_number: str,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get data for a single patent.

        Args:
            patent_number: Patent number
            use_cache: Whether to use cached results

        Returns:
            Patent data dictionary or None if not found
        """
        # Try Google Patents
        if self.google_patents_key:
            try:
                response = self.http.get(
                    f"{self.google_patents_url}/patent/{patent_number}",
                    params={"key": self.google_patents_key},
                )
                if response.ok:
                    return self._process_google_patent(response.json())
            except Exception as e:
                self.logger.error(f"Google Patents error: {str(e)}")

        # Try Lens.org
        if self.lens_api_key:
            try:
                response = self.http.get(
                    f"{self.lens_url}/{patent_number}",
                    params={"token": self.lens_api_key},
                )
                if response.ok:
                    return self._process_lens_patent(response.json())
            except Exception as e:
                self.logger.error(f"Lens.org error: {str(e)}")

        return None

    def _process_google_patent(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Process Google Patents response.

        Args:
            data: Raw API response

        Returns:
            Processed patent data
        """
        # Extract chemical info using LLM
        chemical_info = extract_chemical_info(
            data.get("description", ""),
            data.get("claims", ""),
        )

        # Analyze patent text
        analysis = analyze_patent_text(
            title=data.get("title", ""),
            abstract=data.get("abstract", ""),
            description=data.get("description", ""),
            claims=data.get("claims", ""),
        )

        return {
            "patent_number": data["publication_number"],
            "title": data.get("title"),
            "abstract": data.get("abstract"),
            "filing_date": data.get("filing_date"),
            "publication_date": data.get("publication_date"),
            "assignee": data.get("assignee"),
            "inventors": data.get("inventors", []),
            "url": f"https://patents.google.com/patent/{data['publication_number']}",
            "chemical_info": chemical_info,
            "analysis": analysis,
        }

    def _process_lens_patent(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Process Lens.org response.

        Args:
            data: Raw API response

        Returns:
            Processed patent data
        """
        # Extract chemical info using LLM
        chemical_info = extract_chemical_info(
            data.get("description", ""),
            data.get("claims", ""),
        )

        # Analyze patent text
        analysis = analyze_patent_text(
            title=data.get("title", ""),
            abstract=data.get("abstract", ""),
            description=data.get("description", ""),
            claims=data.get("claims", ""),
        )

        return {
            "patent_number": data["lens_id"],
            "title": data.get("title"),
            "abstract": data.get("abstract"),
            "filing_date": data.get("date_filed"),
            "publication_date": data.get("date_published"),
            "assignee": data.get("assignee"),
            "inventors": data.get("inventors", []),
            "url": f"https://www.lens.org/{data['lens_id']}",
            "chemical_info": chemical_info,
            "analysis": analysis,
        }
