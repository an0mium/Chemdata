"""Swiss tools integration client.

This module provides functionality to:
1. Get target predictions from SwissTargetPrediction
2. Get ADME properties from SwissADME
3. Process and validate results
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime
import time

import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm

from .base_client import BaseWebClient
from ..models.compound import Compound


class SwissClient(BaseWebClient):
    """Client for Swiss bioinformatics tools."""

    # API endpoints
    STP_URL = "http://www.swisstargetprediction.ch/predict.php"
    ADME_URL = "http://www.swissadme.ch/predict.php"

    # Result polling intervals
    POLL_INTERVAL = 5  # seconds
    MAX_POLLS = 60  # 5 minutes total

    def __init__(
        self,
        http_client: Optional["HTTPClient"] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize Swiss tools client.
        
        Args:
            http_client: Optional HTTP client to use
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            logger: Optional logger instance
        """
        super().__init__(http_client, model_dir, cache_dir, logger)

    def process_compounds(
        self,
        compounds: List[Compound],
        skip_predictions: bool = False,
        use_cache: bool = True,
    ) -> None:
        """Process list of compounds.
        
        Args:
            compounds: List of compounds to process
            skip_predictions: Whether to skip predictions
            use_cache: Whether to use cached results
        """
        for compound in compounds:
            try:
                # Get target predictions
                if not skip_predictions:
                    targets = self._get_target_predictions(
                        compound.smiles,
                        use_cache=use_cache,
                    )
                    if targets:
                        compound.swiss_data["targets"] = targets

                # Get ADME properties
                adme = self._get_adme_properties(
                    compound.smiles,
                    use_cache=use_cache,
                )
                if adme:
                    compound.swiss_data["adme"] = adme

            except Exception as e:
                self.logger.error(f"Error processing {compound.name}: {str(e)}")

    def get_compound_data(
        self,
        name: str,
        smiles: str,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get Swiss tools data for a compound.
        
        Args:
            name: Compound name
            smiles: SMILES string
            use_cache: Whether to use cached results
            
        Returns:
            Dictionary of Swiss tools data or None if error
        """
        data = {}

        try:
            # Get target predictions
            targets = self._get_target_predictions(smiles, use_cache)
            if targets:
                data["targets"] = targets

            # Get ADME properties
            adme = self._get_adme_properties(smiles, use_cache)
            if adme:
                data["adme"] = adme

            return data if data else None

        except Exception as e:
            self.logger.error(f"Error getting Swiss data for {name}: {str(e)}")
            return None

    def _get_target_predictions(
        self,
        smiles: str,
        use_cache: bool = True,
    ) -> Optional[List[Dict[str, Any]]]:
        """Get target predictions from SwissTargetPrediction.
        
        Args:
            smiles: SMILES string
            use_cache: Whether to use cached results
            
        Returns:
            List of target predictions or None if error
        """
        try:
            # Submit prediction request
            response = self.http.post(
                self.STP_URL,
                data={"smiles": smiles},
                use_cache=use_cache,
            )
            job_id = response.json()["job_id"]

            # Poll for results
            for _ in range(self.MAX_POLLS):
                time.sleep(self.POLL_INTERVAL)
                
                response = self.http.get(
                    f"{self.STP_URL}/status/{job_id}",
                    use_cache=use_cache,
                )
                status = response.json()["status"]
                
                if status == "completed":
                    results = response.json()["results"]
                    return [
                        {
                            "target": result["target"],
                            "uniprot": result["uniprot"],
                            "gene": result["gene"],
                            "probability": float(result["probability"]),
                            "class": result["class"],
                            "known_actives": int(result["known_actives"]),
                        }
                        for result in results
                    ]
                elif status == "failed":
                    self.logger.error(
                        f"SwissTargetPrediction failed for {smiles}"
                    )
                    return None

            self.logger.error(f"SwissTargetPrediction timeout for {smiles}")
            return None

        except Exception as e:
            self.logger.error(
                f"Error getting target predictions: {str(e)}"
            )
            return None

    def _get_adme_properties(
        self,
        smiles: str,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get ADME properties from SwissADME.
        
        Args:
            smiles: SMILES string
            use_cache: Whether to use cached results
            
        Returns:
            Dictionary of ADME properties or None if error
        """
        try:
            # Submit property calculation request
            response = self.http.post(
                self.ADME_URL,
                data={"smiles": smiles},
                use_cache=use_cache,
            )
            job_id = response.json()["job_id"]

            # Poll for results
            for _ in range(self.MAX_POLLS):
                time.sleep(self.POLL_INTERVAL)
                
                response = self.http.get(
                    f"{self.ADME_URL}/status/{job_id}",
                    use_cache=use_cache,
                )
                status = response.json()["status"]
                
                if status == "completed":
                    results = response.json()["results"]
                    return {
                        # Physicochemical properties
                        "molecular_weight": float(results["MW"]),
                        "logp": float(results["LogP"]),
                        "hbd": int(results["HBD"]),
                        "hba": int(results["HBA"]),
                        "tpsa": float(results["TPSA"]),
                        "rotatable_bonds": int(results["RotBonds"]),
                        
                        # Drug-likeness
                        "lipinski": results["Lipinski"],
                        "ghose": results["Ghose"],
                        "veber": results["Veber"],
                        "egan": results["Egan"],
                        "muegge": results["Muegge"],
                        
                        # ADME predictions
                        "gi_absorption": results["GI_absorption"],
                        "bbb_permeant": results["BBB_permeant"],
                        "pgp_substrate": results["Pgp_substrate"],
                        "cyp_inhibition": {
                            "1A2": results["CYP1A2_inhibition"],
                            "2C19": results["CYP2C19_inhibition"],
                            "2C9": results["CYP2C9_inhibition"],
                            "2D6": results["CYP2D6_inhibition"],
                            "3A4": results["CYP3A4_inhibition"],
                        },
                        
                        # Medicinal chemistry
                        "pains": results["PAINS"],
                        "brenk": results["Brenk"],
                        "leadlikeness": results["Leadlikeness"],
                        "synthetic_accessibility": float(results["SA"]),
                        
                        # Metadata
                        "timestamp": datetime.now().isoformat(),
                    }
                elif status == "failed":
                    self.logger.error(
                        f"SwissADME failed for {smiles}"
                    )
                    return None

            self.logger.error(f"SwissADME timeout for {smiles}")
            return None

        except Exception as e:
            self.logger.error(
                f"Error getting ADME properties: {str(e)}"
            )
            return None
