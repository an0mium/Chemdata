"""Enhanced Swiss tools integration client.

This module provides functionality to:
1. Get target predictions from SwissTargetPrediction
2. Get ADME properties from SwissADME
3. Process and validate results

Enhanced with:
- Circuit breaker pattern for resilience
- Better error handling and recovery
- Improved metrics collection
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, TYPE_CHECKING
from datetime import datetime
import time

import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm

from .base_client import BaseWebClient
from ..models.compound import Compound
from ..pipeline.infrastructure.circuit_breaker import CircuitConfig

if TYPE_CHECKING:
    from .http_client_enhanced import HTTPClientEnhanced


class SwissClientEnhanced(BaseWebClient):
    """Enhanced client for Swiss bioinformatics tools."""

    # API endpoints
    STP_URL = "http://www.swisstargetprediction.ch/predict.php"
    ADME_URL = "http://www.swissadme.ch/predict.php"

    # Result polling intervals
    POLL_INTERVAL = 5  # seconds
    MAX_POLLS = 60  # 5 minutes total

    def __init__(
        self,
        http_client: Optional["HTTPClientEnhanced"] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
        circuit_config: Optional[CircuitConfig] = None,
    ):
        """Initialize Swiss tools client.
        
        Args:
            http_client: Optional HTTP client to use
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            logger: Optional logger instance
            circuit_config: Optional circuit breaker configuration
        """
        super().__init__(http_client, model_dir, cache_dir, logger)
        self.processed_compounds: List[str] = []
        self.failed_compounds: List[str] = []

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
                        fallback=self._get_cached_targets,
                    )
                    if targets:
                        compound.swiss_data["targets"] = targets

                # Get ADME properties
                adme = self._get_adme_properties(
                    compound.smiles,
                    use_cache=use_cache,
                    fallback=self._get_cached_adme,
                )
                if adme:
                    compound.swiss_data["adme"] = adme

                self.processed_compounds.append(compound.name)

            except Exception as e:
                self.logger.error(f"Error processing {compound.name}: {str(e)}")
                self.failed_compounds.append(compound.name)

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
            targets = self._get_target_predictions(
                smiles,
                use_cache=use_cache,
                fallback=self._get_cached_targets,
            )
            if targets:
                data["targets"] = targets

            # Get ADME properties
            adme = self._get_adme_properties(
                smiles,
                use_cache=use_cache,
                fallback=self._get_cached_adme,
            )
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
        fallback: Optional[callable] = None,
    ) -> Optional[List[Dict[str, Any]]]:
        """Get target predictions from SwissTargetPrediction.
        
        Args:
            smiles: SMILES string
            use_cache: Whether to use cached results
            fallback: Optional fallback function if service fails
            
        Returns:
            List of target predictions or None if error
        """
        try:
            # Submit prediction request
            response = self.http.post(
                self.STP_URL,
                data={"smiles": smiles},
                use_cache=use_cache,
                fallback=fallback,
            )
            job_id = response.json()["job_id"]

            # Poll for results
            for _ in range(self.MAX_POLLS):
                time.sleep(self.POLL_INTERVAL)
                
                response = self.http.get(
                    f"{self.STP_URL}/status/{job_id}",
                    use_cache=use_cache,
                    fallback=fallback,
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
        fallback: Optional[callable] = None,
    ) -> Optional[Dict[str, Any]]:
        """Get ADME properties from SwissADME.
        
        Args:
            smiles: SMILES string
            use_cache: Whether to use cached results
            fallback: Optional fallback function if service fails
            
        Returns:
            Dictionary of ADME properties or None if error
        """
        try:
            # Submit property calculation request
            response = self.http.post(
                self.ADME_URL,
                data={"smiles": smiles},
                use_cache=use_cache,
                fallback=fallback,
            )
            job_id = response.json()["job_id"]

            # Poll for results
            for _ in range(self.MAX_POLLS):
                time.sleep(self.POLL_INTERVAL)
                
                response = self.http.get(
                    f"{self.ADME_URL}/status/{job_id}",
                    use_cache=use_cache,
                    fallback=fallback,
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

    def _get_cached_targets(self, smiles: str) -> Optional[List[Dict[str, Any]]]:
        """Get cached target predictions.
        
        Args:
            smiles: SMILES string
            
        Returns:
            List of target predictions or None if not cached
        """
        if not self.cache_dir:
            return None

        try:
            cache_file = self.cache_dir / f"targets_{smiles}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                return pd.read_json(f).to_dict("records")

        except Exception as e:
            self.logger.error(f"Error reading cached targets: {str(e)}")
            return None

    def _get_cached_adme(self, smiles: str) -> Optional[Dict[str, Any]]:
        """Get cached ADME properties.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary of ADME properties or None if not cached
        """
        if not self.cache_dir:
            return None

        try:
            cache_file = self.cache_dir / f"adme_{smiles}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                return pd.read_json(f).iloc[0].to_dict()

        except Exception as e:
            self.logger.error(f"Error reading cached ADME: {str(e)}")
            return None

    def get_metrics(self) -> Dict[str, Any]:
        """Get client metrics."""
        return {
            "processed_compounds": len(self.processed_compounds),
            "failed_compounds": len(self.failed_compounds),
            "success_rate": (
                len(self.processed_compounds) /
                (len(self.processed_compounds) + len(self.failed_compounds))
                if self.processed_compounds or self.failed_compounds
                else 0
            ),
        }
