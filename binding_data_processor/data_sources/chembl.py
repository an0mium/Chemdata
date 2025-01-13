"""ChEMBL API client for retrieving binding and activity data.

This module provides an enhanced client for the ChEMBL API with:
- Compound searching by multiple identifiers
- Binding data retrieval
- Target data retrieval
- Psychopharmacological data enrichment
- Robust error handling and retries
- Caching support
"""

from typing import Any, Dict, List, Optional, Set
from chembl_webresource_client.new_client import new_client
import time
from functools import lru_cache

from ..web_enrichment.clients.base_client import BaseWebClient
from ..models.compound import Compound, BindingData
from ..pipeline.infrastructure.circuit_breaker import CircuitBreaker


class ChEMBLClient(BaseWebClient):
    """Enhanced client for ChEMBL API using official Python client."""

    # Activity types of interest for psychopharmacology
    ACTIVITY_TYPES = {"Ki", "IC50", "Kd", "EC50", "potency", "activity"}

    # Target types of interest
    TARGET_TYPES = {"SINGLE PROTEIN", "PROTEIN COMPLEX", "PROTEIN FAMILY", "SELECTIVITY GROUP"}

    # Receptor families of interest
    RECEPTOR_FAMILIES = {"5-HT", "NMDA", "GABA", "dopamine", "opioid", "cannabinoid", "sigma", "adrenergic"}

    def __init__(self, **kwargs):
        """Initialize ChEMBL client with caching and circuit breaker.

        Args:
            **kwargs: Additional arguments passed to BaseWebClient
        """
        super().__init__(**kwargs)
        self.circuit_breaker = CircuitBreaker(failure_threshold=5, reset_timeout=300, name="chembl")
        self._init_clients()

    def _init_clients(self):
        """Initialize ChEMBL API clients with retry and circuit breaker."""
        max_retries = 3
        for attempt in range(max_retries):
            try:
                with self.circuit_breaker:
                    self.molecule = new_client.molecule
                    self.activity = new_client.activity
                    self.target = new_client.target
                    self.mechanism = new_client.mechanism
                    return
            except Exception as e:
                self.logger.warning(f"Error initializing ChEMBL clients (attempt {attempt + 1}/{max_retries}): {str(e)}")
                if attempt < max_retries - 1:
                    time.sleep(2**attempt)  # Exponential backoff
                else:
                    raise

    @lru_cache(maxsize=1000)
    def search_compound(self, compound: Dict[str, str], use_cache: bool = True) -> Dict[str, Any]:
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
            with self.circuit_breaker:
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
                        data = self.get_compound_by_chembl_id(results[0]["molecule_chembl_id"])
                        if use_cache and self.cache_dir:
                            self.http.set_cache(cache_key, data)
                        return data

                # Try CAS number
                if compound.get("cas"):
                    results = list(self.molecule.filter(molecule_synonyms__synonym__iexact=compound["cas"]))
                    if results:
                        data = self.get_compound_by_chembl_id(results[0]["molecule_chembl_id"])
                        if use_cache and self.cache_dir:
                            self.http.set_cache(cache_key, data)
                        return data

                # Try structure search if SMILES is available
                if compound.get("smiles"):
                    results = list(self.molecule.filter(molecule_structures__canonical_smiles__flexmatch=compound["smiles"]))
                    if results:
                        data = self.get_compound_by_chembl_id(results[0]["molecule_chembl_id"])
                        if use_cache and self.cache_dir:
                            self.http.set_cache(cache_key, data)
                        return data

                # Try InChI search
                if compound.get("inchi"):
                    results = list(self.molecule.filter(molecule_structures__standard_inchi=compound["inchi"]))
                    if results:
                        data = self.get_compound_by_chembl_id(results[0]["molecule_chembl_id"])
                        if use_cache and self.cache_dir:
                            self.http.set_cache(cache_key, data)
                        return data

        except Exception as e:
            self.logger.warning(f"Error searching compound: {str(e)}")

        return {}

    def get_compound_by_chembl_id(self, chembl_id: str, use_cache: bool = True) -> Dict[str, Any]:
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
                with self.circuit_breaker:
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
                    break
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
    def search_targets(self, query: str, target_type: Optional[str] = None) -> List[Dict[str, Any]]:
        """Search for protein targets with filtering.

        Args:
            query: Search query
            target_type: Optional target type filter

        Returns:
            List of matching targets
        """
        try:
            with self.circuit_breaker:
                filters = {"pref_name__icontains": query}
                if target_type:
                    filters["target_type"] = target_type
                return list(self.target.filter(**filters))
        except Exception as e:
            self.logger.error(f"Error searching targets: {str(e)}")
            return []

    @lru_cache(maxsize=100)
    def get_target_by_chembl_id(self, target_id: str) -> Dict[str, Any]:
        """Get target data by ChEMBL target ID.

        Args:
            target_id: ChEMBL target ID

        Returns:
            Target data
        """
        try:
            with self.circuit_breaker:
                return self.target.get(target_id)
        except Exception as e:
            self.logger.error(f"Error getting target: {str(e)}")
            return {}

    def process_compounds(
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
                chembl_data = self.search_compound(
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
        # Search ChEMBL
        data = self.search_compound(
            {
                "name": name,
                "cas": cas_number,
            },
            use_cache=use_cache,
        )

        if not data:
            return None

        return data
