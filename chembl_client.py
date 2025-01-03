"""ChEMBL API client for retrieving binding data."""

from typing import Any, Dict, List
from chembl_webresource_client.new_client import new_client
from models import BindingData
from logger import LogManager


class ChEMBLClient:
    """Client for ChEMBL API using official Python client."""
    
    def __init__(self):
        """Initialize ChEMBL client."""
        self.logger = LogManager().get_logger("chembl_client")
        self._init_clients()
        
    def _init_clients(self):
        """Initialize ChEMBL API clients with retry."""
        max_retries = 3
        for attempt in range(max_retries):
            try:
                self.molecule = new_client.molecule
                self.activity = new_client.activity
                self.target = new_client.target
                return
            except Exception as e:
                self.logger.warning(f"Error initializing ChEMBL clients (attempt {attempt + 1}/{max_retries}): {str(e)}")
                if attempt < max_retries - 1:
                    import time
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    raise

    def search_compound(self, compound: Dict[str, str]) -> Dict[str, Any]:
        """
        Search for a compound using multiple identifiers.
        
        Args:
            compound: Dictionary with available identifiers (name, cas, smiles, inchi, etc.)
            
        Returns:
            Compound data
        """
        try:
            # Try exact name match first
            if compound.get('name'):
                # Try exact match
                results = list(self.molecule.filter(pref_name__iexact=compound['name']))
                if not results:
                    # Try synonym match
                    results = list(self.molecule.filter(molecule_synonyms__synonym__iexact=compound['name']))
                if not results:
                    # Try fuzzy name match
                    results = list(self.molecule.filter(pref_name__icontains=compound['name']))
                if results:
                    return self.get_compound_by_chembl_id(results[0]['molecule_chembl_id'])
                    
            # Try CAS number
            if compound.get('cas'):
                results = list(self.molecule.filter(molecule_synonyms__synonym__iexact=compound['cas']))
                if results:
                    return self.get_compound_by_chembl_id(results[0]['molecule_chembl_id'])
                    
            # Try structure search if SMILES is available
            if compound.get('smiles'):
                results = list(self.molecule.filter(molecule_structures__canonical_smiles__flexmatch=compound['smiles']))
                if results:
                    return self.get_compound_by_chembl_id(results[0]['molecule_chembl_id'])
                    
            # Try InChI search
            if compound.get('inchi'):
                results = list(self.molecule.filter(molecule_structures__standard_inchi=compound['inchi']))
                if results:
                    return self.get_compound_by_chembl_id(results[0]['molecule_chembl_id'])
        except Exception as e:
            self.logger.warning(f"Error searching compound: {str(e)}")
                
        return {}

    def get_compound_by_chembl_id(self, chembl_id: str) -> Dict[str, Any]:
        """
        Get compound data by ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            Compound data
        """
        max_retries = 3
        for attempt in range(max_retries):
            try:
                # Re-initialize clients if needed
                if attempt > 0:
                    self._init_clients()
                    
                # Get molecule details
                molecule_data = self.molecule.get(chembl_id)
                
                # Get binding data
                activities = list(self.activity.filter(
                    molecule_chembl_id=chembl_id,
                    type__in=['Ki', 'IC50', 'Kd', 'EC50'],
                    relation__in=['=', '<', '>', '<=', '>='],
                    standard_units__isnull=False
                )[:12])
                break
            except Exception as e:
                self.logger.warning(f"Error getting compound data (attempt {attempt + 1}/{max_retries}): {str(e)}")
                if attempt < max_retries - 1:
                    import time
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    return {}
        
        data = {
            'chembl_id': chembl_id,
            'binding_data': []
        }
        
        try:
            if molecule_data:
                data.update({
                    'name': molecule_data.get('pref_name', ''),
                    'synonyms': [s.get('synonym', '') for s in molecule_data.get('molecule_synonyms', []) if s.get('synonym')],
                    'chembl_url': f"https://www.ebi.ac.uk/chembl/compound/{chembl_id}"
                })
                
            for activity in activities:
                try:
                    binding = BindingData(
                        target_common_name=activity.get('target_pref_name') or 'N/A',
                        target_protein_name=(activity.get('target_components', [{}]) or [{}])[0].get('protein_name') or 'N/A',
                        target_gene_name=(activity.get('target_components', [{}]) or [{}])[0].get('gene_name') or 'N/A',
                        affinity_value=float(activity.get('value') or 0.0),
                        affinity_type=activity.get('type') or 'N/A',
                        affinity_unit=activity.get('units') or 'N/A',
                        assay_description=activity.get('assay_description') or 'N/A',
                        reference=activity.get('document_chembl_id') or 'N/A',
                        confidence_score=float(activity.get('confidence_score') or 0.0)
                    )
                    data['binding_data'].append(binding)
                except Exception as e:
                    self.logger.warning(f"Error processing activity data: {str(e)}")
                    continue
        except Exception as e:
            self.logger.warning(f"Error processing molecule data: {str(e)}")
                
        return data

    def search_targets(self, query: str) -> List[Dict[str, Any]]:
        """
        Search for protein targets.
        
        Args:
            query: Search query
            
        Returns:
            List of matching targets
        """
        return list(self.target.filter(pref_name__icontains=query))

    def get_target_by_chembl_id(self, target_id: str) -> Dict[str, Any]:
        """
        Get target data by ChEMBL target ID.
        
        Args:
            target_id: ChEMBL target ID
            
        Returns:
            Target data
        """
        return self.target.get(target_id)
