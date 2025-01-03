"""Swiss* web services integration.

This module handles:
1. SwissTargetPrediction - Target prediction
2. SwissADME - ADME properties and pharmacokinetics
3. SwissSimilarity - Chemical similarity search
"""

from typing import Dict, Any, List, Optional
import json
import time
from urllib.parse import quote

from ..http_client import HttpClient
from logger import LogManager

logger = LogManager().get_logger("web_enrichment.data_sources.swiss")


class SwissClient:
    """Client for Swiss* web services."""
    
    # Base URLs for Swiss* services
    STP_BASE_URL = "http://www.swisstargetprediction.ch/predict.php"
    ADME_BASE_URL = "http://www.swissadme.ch/predict.php"
    SIM_BASE_URL = "http://www.swisssimilarity.ch/search.php"
    
    def __init__(self, http_client: HttpClient):
        """Initialize Swiss* client."""
        self.http = http_client

    def get_target_predictions(
        self,
        smiles: str,
        organism: str = "Homo sapiens"
    ) -> Optional[Dict[str, Any]]:
        """
        Get target predictions from SwissTargetPrediction.
        
        Args:
            smiles: SMILES string of compound
            organism: Target organism (default: Homo sapiens)
            
        Returns:
            Dictionary containing target predictions or None
        """
        try:
            # Submit prediction request
            params = {
                'smiles': smiles,
                'organism': organism
            }
            
            response = self.http.make_request(
                self.STP_BASE_URL,
                params=params,
                method='POST'
            )
            
            if not response:
                return None
                
            # Parse response
            data = response.json()
            if not data.get('job_id'):
                return None
                
            job_id = data['job_id']
            
            # Poll for results
            result_url = f"{self.STP_BASE_URL}?job_id={job_id}"
            max_attempts = 10
            
            for _ in range(max_attempts):
                time.sleep(2)  # Wait between polls
                response = self.http.make_request(result_url)
                
                if response and response.json().get('status') == 'completed':
                    predictions = response.json().get('predictions', [])
                    return {
                        'predictions': [
                            {
                                'target': pred['target_name'],
                                'probability': pred['probability'],
                                'target_class': pred['target_class'],
                                'uniprot_id': pred['uniprot_id'],
                                'chembl_id': pred.get('chembl_id'),
                                'common_name': pred.get('common_name')
                            }
                            for pred in predictions
                        ],
                        'url': f"http://www.swisstargetprediction.ch/result.php?job={job_id}"
                    }
                    
        except Exception as e:
            logger.error(f"Error getting target predictions: {str(e)}")
            
        return None

    def get_adme_properties(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Get ADME properties from SwissADME.
        
        Args:
            smiles: SMILES string of compound
            
        Returns:
            Dictionary containing ADME properties or None
        """
        try:
            # Submit property calculation request
            params = {
                'smiles': smiles
            }
            
            response = self.http.make_request(
                self.ADME_BASE_URL,
                params=params,
                method='POST'
            )
            
            if not response:
                return None
                
            # Parse response
            data = response.json()
            if not data.get('job_id'):
                return None
                
            job_id = data['job_id']
            
            # Poll for results
            result_url = f"{self.ADME_BASE_URL}?job_id={job_id}"
            max_attempts = 10
            
            for _ in range(max_attempts):
                time.sleep(2)  # Wait between polls
                response = self.http.make_request(result_url)
                
                if response and response.json().get('status') == 'completed':
                    properties = response.json().get('properties', {})
                    return {
                        'physicochemical': {
                            'mw': properties.get('MW'),
                            'logp': properties.get('LogP'),
                            'hbd': properties.get('HBD'),
                            'hba': properties.get('HBA'),
                            'tpsa': properties.get('TPSA'),
                            'rotatable_bonds': properties.get('RotBonds')
                        },
                        'lipinski': {
                            'violations': properties.get('Lipinski_Violations'),
                            'is_druglike': properties.get('Druglike') == 'Yes'
                        },
                        'absorption': {
                            'gi_absorption': properties.get('GI_Absorption'),
                            'bbb_permeant': properties.get('BBB_Permeant') == 'Yes',
                            'pgp_substrate': properties.get('Pgp_Substrate') == 'Yes'
                        },
                        'metabolism': {
                            'cyp_inhibition': {
                                'cyp1a2': properties.get('CYP1A2_Inhibition') == 'Yes',
                                'cyp2c19': properties.get('CYP2C19_Inhibition') == 'Yes',
                                'cyp2c9': properties.get('CYP2C9_Inhibition') == 'Yes',
                                'cyp2d6': properties.get('CYP2D6_Inhibition') == 'Yes',
                                'cyp3a4': properties.get('CYP3A4_Inhibition') == 'Yes'
                            }
                        },
                        'url': f"http://www.swissadme.ch/result.php?job={job_id}"
                    }
                    
        except Exception as e:
            logger.error(f"Error getting ADME properties: {str(e)}")
            
        return None

    def search_similar_compounds(
        self,
        smiles: str,
        similarity_threshold: float = 0.5,
        max_results: int = 100
    ) -> Optional[Dict[str, Any]]:
        """
        Search for similar compounds using SwissSimilarity.
        
        Args:
            smiles: SMILES string of query compound
            similarity_threshold: Minimum similarity score (0-1)
            max_results: Maximum number of results to return
            
        Returns:
            Dictionary containing similar compounds or None
        """
        try:
            # Submit similarity search request
            params = {
                'smiles': smiles,
                'threshold': similarity_threshold,
                'limit': max_results
            }
            
            response = self.http.make_request(
                self.SIM_BASE_URL,
                params=params,
                method='POST'
            )
            
            if not response:
                return None
                
            # Parse response
            data = response.json()
            if not data.get('job_id'):
                return None
                
            job_id = data['job_id']
            
            # Poll for results
            result_url = f"{self.SIM_BASE_URL}?job_id={job_id}"
            max_attempts = 10
            
            for _ in range(max_attempts):
                time.sleep(2)  # Wait between polls
                response = self.http.make_request(result_url)
                
                if response and response.json().get('status') == 'completed':
                    compounds = response.json().get('compounds', [])
                    return {
                        'similar_compounds': [
                            {
                                'smiles': cmpd['smiles'],
                                'similarity': cmpd['similarity'],
                                'name': cmpd.get('name'),
                                'source': cmpd.get('source'),
                                'source_id': cmpd.get('source_id')
                            }
                            for cmpd in compounds
                        ],
                        'url': f"http://www.swisssimilarity.ch/result.php?job={job_id}"
                    }
                    
        except Exception as e:
            logger.error(f"Error searching similar compounds: {str(e)}")
            
        return None
