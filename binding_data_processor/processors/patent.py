"""Patent searching and compound extraction.

This module provides functionality for:
1. Searching patents using Google Patents API
2. Extracting chemical compounds from patent text
3. Analyzing patent claims and examples
4. Gathering patent metadata
5. Validating extracted structures
"""

import re
import json
import logging
from typing import Dict, List, Optional, Set, Tuple
import requests
from rdkit import Chem

from .structure import StructureProcessor


class PatentProcessor:
    """Handles patent searching and compound extraction."""

    def __init__(self, structure_processor: Optional[StructureProcessor] = None):
        """
        Initialize patent processor.

        Args:
            structure_processor: Optional StructureProcessor instance for structure validation
        """
        self.logger = logging.getLogger(__name__)
        self.structure_processor = structure_processor or StructureProcessor()

    def search_patents(
        self,
        query: str,
        api_key: str,
        max_results: int = 100,
        from_date: Optional[str] = None,
        to_date: Optional[str] = None,
    ) -> List[Dict[str, str]]:
        """
        Search patents using Google Patents API.

        Args:
            query: Search query string
            api_key: Google Patents API key
            max_results: Maximum number of results to return
            from_date: Optional start date in YYYY-MM-DD format
            to_date: Optional end date in YYYY-MM-DD format

        Returns:
            List of patent dictionaries containing:
            - number: Patent number
            - title: Patent title
            - abstract: Patent abstract
            - url: Google Patents URL
            - publication_date: Publication date
            - assignee: Patent assignee
            - inventors: List of inventors
        """
        try:
            # Build query parameters
            params = {
                "q": query,
                "key": api_key,
                "maxResults": max_results,
            }
            if from_date:
                params["publishedAfter"] = from_date
            if to_date:
                params["publishedBefore"] = to_date

            # Make API request
            response = requests.get(
                "https://patents.google.com/api/query",
                params=params,
                headers={"Accept": "application/json"},
            )
            response.raise_for_status()
            data = response.json()

            # Process results
            patents = []
            for result in data.get("results", []):
                patent = {
                    "number": result.get("patent_number", ""),
                    "title": result.get("title", ""),
                    "abstract": result.get("abstract", ""),
                    "url": f"https://patents.google.com/patent/{result.get('patent_number', '')}",
                    "publication_date": result.get("publication_date", ""),
                    "assignee": result.get("assignee", {}).get("name", ""),
                    "inventors": [
                        inv.get("name", "") for inv in result.get("inventors", [])
                    ],
                }
                patents.append(patent)

            return patents

        except requests.exceptions.RequestException as e:
            self.logger.error(f"Patent search API error: {str(e)}")
            return []
        except Exception as e:
            self.logger.error(f"Patent search error: {str(e)}")
            return []

    def extract_compounds(
        self,
        text: str,
        llm_api_key: str,
        validate: bool = True,
    ) -> List[Dict[str, str]]:
        """
        Extract chemical compounds from patent text using LLM.

        Args:
            text: Patent text to analyze
            llm_api_key: API key for LLM service
            validate: Whether to validate extracted structures

        Returns:
            List of compound dictionaries containing:
            - name: Compound name or identifier
            - smiles: SMILES string
            - inchi: InChI string (if available)
            - source: Source section in patent
            - example_id: Example number (if from example)
        """
        try:
            # Prepare extraction prompt
            prompt = self._build_extraction_prompt(text)

            # Make LLM API request
            response = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {llm_api_key}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": "gpt-4",
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.7,
                },
            )
            response.raise_for_status()

            # Parse extracted compounds
            result = response.json()
            compounds = json.loads(result["choices"][0]["message"]["content"])

            # Validate structures if requested
            if validate:
                valid_compounds = []
                for compound in compounds:
                    if "smiles" not in compound:
                        continue
                    is_valid, _ = self.structure_processor.validate_structure(
                        compound["smiles"]
                    )
                    if is_valid:
                        # Generate InChI if not provided
                        if "inchi" not in compound:
                            mol = Chem.MolFromSmiles(compound["smiles"])
                            compound["inchi"] = Chem.MolToInchi(mol)
                        valid_compounds.append(compound)
                return valid_compounds
            return compounds

        except requests.exceptions.RequestException as e:
            self.logger.error(f"LLM API error: {str(e)}")
            return []
        except Exception as e:
            self.logger.error(f"Compound extraction error: {str(e)}")
            return []

    def analyze_claims(
        self,
        text: str,
        llm_api_key: str,
    ) -> Dict[str, List[str]]:
        """
        Analyze patent claims using LLM.

        Args:
            text: Patent claims text
            llm_api_key: API key for LLM service

        Returns:
            Dictionary containing:
            - compound_claims: Claims covering compounds
            - method_claims: Claims covering methods
            - composition_claims: Claims covering compositions
            - use_claims: Claims covering uses
        """
        try:
            # Prepare claims analysis prompt
            prompt = self._build_claims_prompt(text)

            # Make LLM API request
            response = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {llm_api_key}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": "gpt-4",
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.7,
                },
            )
            response.raise_for_status()

            # Parse analyzed claims
            result = response.json()
            claims = json.loads(result["choices"][0]["message"]["content"])
            return claims

        except requests.exceptions.RequestException as e:
            self.logger.error(f"LLM API error: {str(e)}")
            return {
                "compound_claims": [],
                "method_claims": [],
                "composition_claims": [],
                "use_claims": [],
            }
        except Exception as e:
            self.logger.error(f"Claims analysis error: {str(e)}")
            return {
                "compound_claims": [],
                "method_claims": [],
                "composition_claims": [],
                "use_claims": [],
            }

    def extract_examples(self, text: str, llm_api_key: str) -> List[Dict[str, str]]:
        """
        Extract synthetic examples from patent text.

        Args:
            text: Patent text
            llm_api_key: API key for LLM service

        Returns:
            List of example dictionaries containing:
            - id: Example number/identifier
            - title: Example title/name
            - procedure: Synthetic procedure
            - compound: Product compound
            - yield: Reaction yield
            - characterization: Analytical data
        """
        try:
            # Prepare examples extraction prompt
            prompt = self._build_examples_prompt(text)

            # Make LLM API request
            response = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {llm_api_key}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": "gpt-4",
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.7,
                },
            )
            response.raise_for_status()

            # Parse extracted examples
            result = response.json()
            examples = json.loads(result["choices"][0]["message"]["content"])
            return examples

        except requests.exceptions.RequestException as e:
            self.logger.error(f"LLM API error: {str(e)}")
            return []
        except Exception as e:
            self.logger.error(f"Examples extraction error: {str(e)}")
            return []

    def _build_extraction_prompt(self, text: str) -> str:
        """Build prompt for compound extraction."""
        return f"""Extract chemical compounds and their properties from the following text. For each compound, provide:

1. Chemical names (IUPAC, common names, and any synonyms)
2. Chemical structure information:
- SMILES string
- InChI string (if available)
- Molecular formula
- Structural features (rings, functional groups)
3. Source information:
- Location in text (claims, examples, etc.)
- Example number if applicable
- Any synthetic procedures
4. Property data:
- Physical properties
- Biological activity
- Analytical data

Text: {text}

Format the response as a list of JSON objects with these fields. Include any numerical values with their units. For chemical structures, prioritize standardized identifiers (SMILES, InChI) over descriptive text."""

    def _build_claims_prompt(self, text: str) -> str:
        """Build prompt for claims analysis."""
        return f"""Analyze the following patent claims and categorize them into:

1. Compound claims (covering chemical structures)
2. Method claims (covering synthetic procedures)
3. Composition claims (covering formulations/mixtures)
4. Use claims (covering applications/uses)

For each claim, extract:
- Claim number
- Key subject matter
- Dependencies
- Scope/breadth
- Any numerical ranges/limitations

Claims text: {text}

Format the response as a JSON object with claim categories as keys and lists of analyzed claims as values."""

    def _build_examples_prompt(self, text: str) -> str:
        """Build prompt for examples extraction."""
        return f"""Extract synthetic examples from the following patent text. For each example, provide:

1. Example identification:
- Number/identifier
- Title/name of product
2. Synthetic procedure:
- Reagents and conditions
- Step-by-step process
- Workup/purification
3. Product information:
- Structure (SMILES/InChI)
- Yield
- Physical properties
4. Analytical data:
- Spectral data
- Characterization
- Purity information

Text: {text}

Format the response as a list of JSON objects with these fields. Include all numerical values and units."""
