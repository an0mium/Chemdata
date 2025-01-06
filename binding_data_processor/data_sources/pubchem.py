"""PubChem API client.

This module provides a client for accessing the PubChem API with:
- Compound lookup by name, CAS, InChIKey
- Property retrieval
- Caching
- Rate limiting
"""

from typing import Dict, Any, Optional, List
from urllib.parse import quote

from ..web_enrichment.base_client import BaseWebClient
from ..models.compound import Compound


class PubChemClient(BaseWebClient):
    """Client for PubChem API."""

    def __init__(self, **kwargs):
        """Initialize PubChem client."""
        super().__init__(**kwargs)
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

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
        response = self.http.get(f"{self.base_url}/{endpoint}")
        data = response.json()

        if not data or "IdentifierList" not in data:
            return {}

        cid = str(data["IdentifierList"].get("CID", [None])[0])
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
        response = self.http.get(f"{self.base_url}/{endpoint}")
        data = response.json()

        if not data or "IdentifierList" not in data:
            return {}

        cid = str(data["IdentifierList"].get("CID", [None])[0])
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
        response = self.http.get(f"{self.base_url}/{endpoint}")
        data = response.json()

        if not data or "IdentifierList" not in data:
            return {}

        cid = str(data["IdentifierList"].get("CID", [None])[0])
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

    def get_compound_properties(self, cid: str, properties: List[str]) -> Dict[str, Any]:
        """
        Get specific compound properties.

        Args:
            cid: PubChem CID
            properties: List of property names

        Returns:
            Property data
        """
        endpoint = f"compound/cid/{cid}/property/{','.join(properties)}/JSON"
        response = self.http.get(f"{self.base_url}/{endpoint}")
        return response.json()

    def process_compounds(
        self,
        compounds: List[Compound],
        skip_predictions: bool = False,
        use_cache: bool = True,
    ) -> None:
        """
        Process list of compounds.

        Args:
            compounds: List of compounds to process
            skip_predictions: Whether to skip ML predictions
            use_cache: Whether to use cached results
        """
        for compound in compounds:
            # Try different identifiers
            data = {}
            if compound.name:
                data = self.get_compound_by_name(compound.name)
            if not data and compound.cas_number:
                data = self.search_by_cas(compound.cas_number)
            if not data and compound.inchikey:
                data = self.search_by_inchikey(compound.inchikey)

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

    def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """
        Get data for a single compound.

        Args:
            name: Compound name
            cas_number: Optional CAS number
            use_cache: Whether to use cached results

        Returns:
            Dictionary of compound data or None if not found
        """
        # Try name first
        data = self.get_compound_by_name(name)
        if data:
            return data

        # Try CAS if provided
        if cas_number:
            data = self.search_by_cas(cas_number)
            if data:
                return data

        return None
