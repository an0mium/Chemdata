"""Swiss tools integration client.

This module provides functionality to:
1. Get similar compounds from SwissSimilarity
2. Get protein data from SwissProt
3. Get target predictions from SwissTargetPrediction
4. Get ADME properties from SwissADME
5. Process and validate results

Features:
- Circuit breaker pattern for resilience
- Caching with fallback
- Metrics collection
- Error tracking
- Comprehensive validation
- Batch processing
"""

import logging
import re
from pathlib import Path
from typing import Optional, Dict, Any, Union, List
from datetime import datetime

import pandas as pd

from ....pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig
from ....web_enrichment.clients.base import WebClient
from .http import HTTPClient
from ..validation.schema import ValidationError
from ...base.core import Compound


class SwissClient(WebClient):
    """Client for Swiss bioinformatics tools."""

    def __init__(
        self,
        name: str = "swiss",
        http_client: Optional[HTTPClient] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        circuit_config: Optional[CircuitBreakerConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize Swiss tools client.

        Args:
            name: Client name for circuit breaker
            http_client: Optional HTTP client to use
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            circuit_config: Optional circuit breaker configuration
            logger: Optional logger instance
        """
        super().__init__(
            name=name,
            base_url="",  # Base URLs are set per-service
            data_source="swiss",
            requests_per_second=1.0,
            max_retries=3,
            timeout=30,
            cache_ttl=3600,
            circuit_config=circuit_config,
            logger=logger,
        )

        # API endpoints
        self.swissprot_url = "https://rest.uniprot.org/uniprotkb"
        self.swissadme_url = "https://www.swissadme.ch/api"
        self.swisssimilarity_url = "https://www.swisssimilarity.ch/api"
        self.swisstargetprediction_url = "https://www.swisstargetprediction.ch/api"

        # Tracking
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
                data = self.get_all_data(
                    smiles=compound.smiles,
                    name=compound.name,
                    accession=compound.uniprot_id,
                    use_cache=use_cache,
                )

                if data:
                    compound.swiss_data.update(data)
                    self.processed_compounds.append(compound.name)
                else:
                    self.failed_compounds.append(compound.name)

            except Exception as e:
                self.logger.error(f"Error processing {compound.name}: {str(e)}")
                self.failed_compounds.append(compound.name)

    def get_swissprot_data(
        self,
        accession: str,
        name: str,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Get protein data from SwissProt.

        Args:
            accession: UniProt accession number
            name: Protein name
            use_cache: Whether to use cached results

        Returns:
            SwissProt data dictionary

        Raises:
            ValidationError: If input validation fails
        """
        self._validate_accession(accession)
        self._validate_name(name)

        response = self.request(
            method="GET",
            endpoint=f"{self.swissprot_url}/{accession}",
            params={
                "name": name,
                "format": "json",
                "include_features": True,
            },
            headers={"Accept": "application/json"},
            use_cache=use_cache,
            schema_name="swissprot",
        )

        return response

    def get_swissadme_data(
        self,
        smiles: str,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Get ADME properties from SwissADME.

        Args:
            smiles: SMILES string
            use_cache: Whether to use cached results

        Returns:
            SwissADME data dictionary

        Raises:
            ValidationError: If input validation fails
        """
        self._validate_smiles(smiles)

        response = self.request(
            method="GET",
            endpoint=f"{self.swissadme_url}/predict",
            params={
                "smiles": smiles,
                "format": "json",
            },
            headers={"Accept": "application/json"},
            use_cache=use_cache,
            schema_name="swissadme",
        )

        return response

    def get_swisssimilarity_data(
        self,
        smiles: str,
        name: str,
        threshold: float = 0.7,
        limit: int = 10,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Get similar compounds from SwissSimilarity.

        Args:
            smiles: SMILES string
            name: Compound name
            threshold: Similarity threshold (0-1)
            limit: Maximum number of results
            use_cache: Whether to use cached results

        Returns:
            SwissSimilarity data dictionary

        Raises:
            ValidationError: If input validation fails
        """
        self._validate_smiles(smiles)
        self._validate_name(name)
        self._validate_threshold(threshold)

        response = self.request(
            method="GET",
            endpoint=f"{self.swisssimilarity_url}/similar",
            params={
                "smiles": smiles,
                "name": name,
                "threshold": threshold,
                "limit": limit,
                "format": "json",
            },
            headers={"Accept": "application/json"},
            use_cache=use_cache,
            schema_name="swisssimilarity",
        )

        return response

    def get_swisstargetprediction_data(
        self,
        smiles: str,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Get target predictions from SwissTargetPrediction.

        Args:
            smiles: SMILES string
            use_cache: Whether to use cached results

        Returns:
            SwissTargetPrediction data dictionary

        Raises:
            ValidationError: If input validation fails
        """
        self._validate_smiles(smiles)

        response = self.request(
            method="GET",
            endpoint=f"{self.swisstargetprediction_url}/predict",
            params={
                "smiles": smiles,
                "format": "json",
            },
            headers={"Accept": "application/json"},
            use_cache=use_cache,
            schema_name="swisstargetprediction",
        )

        return response

    def get_all_data(
        self,
        smiles: str,
        name: str,
        accession: str,
        threshold: float = 0.7,
        limit: int = 10,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Get data from all Swiss tools.

        Args:
            smiles: SMILES string
            name: Compound name
            accession: UniProt accession
            threshold: Similarity threshold
            limit: Maximum similarity results
            use_cache: Whether to use cached results

        Returns:
            Dictionary containing all Swiss tools data
        """
        results = {}

        try:
            results["swissprot"] = self.get_swissprot_data(
                accession=accession,
                name=name,
                use_cache=use_cache,
            )
        except ValidationError as e:
            self.logger.error(f"SwissProt error: {str(e)}")
            results["swissprot"] = None

        try:
            results["swissadme"] = self.get_swissadme_data(
                smiles=smiles,
                use_cache=use_cache,
            )
        except ValidationError as e:
            self.logger.error(f"SwissADME error: {str(e)}")
            results["swissadme"] = None

        try:
            results["swisssimilarity"] = self.get_swisssimilarity_data(
                smiles=smiles,
                name=name,
                threshold=threshold,
                limit=limit,
                use_cache=use_cache,
            )
        except ValidationError as e:
            self.logger.error(f"SwissSimilarity error: {str(e)}")
            results["swisssimilarity"] = None

        try:
            results["swisstargetprediction"] = self.get_swisstargetprediction_data(
                smiles=smiles,
                use_cache=use_cache,
            )
        except ValidationError as e:
            self.logger.error(f"SwissTargetPrediction error: {str(e)}")
            results["swisstargetprediction"] = None

        return results

    def _validate_accession(self, accession: str) -> None:
        """Validate UniProt accession number.

        Args:
            accession: UniProt accession to validate

        Raises:
            TypeError: If accession is not a string
            ValueError: If accession is empty or invalid format
        """
        if not isinstance(accession, str):
            raise TypeError("Accession must be a string")
        if not accession:
            raise ValueError("Accession cannot be empty")
        # Split regex pattern to avoid line length issues
        pattern = r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|" r"[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
        if not re.match(pattern, accession):
            raise ValueError("Invalid UniProt accession format")

    def _validate_smiles(self, smiles: str) -> None:
        """Validate SMILES string.

        Args:
            smiles: SMILES string to validate

        Raises:
            TypeError: If SMILES is not a string
            ValueError: If SMILES is empty
        """
        if not isinstance(smiles, str):
            raise TypeError("SMILES must be a string")
        if not smiles:
            raise ValueError("SMILES cannot be empty")

    def _validate_name(self, name: str) -> None:
        """Validate name string.

        Args:
            name: Name to validate

        Raises:
            TypeError: If name is not a string
            ValueError: If name is empty
        """
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        if not name:
            raise ValueError("Name cannot be empty")

    def _validate_threshold(self, threshold: Union[int, float]) -> None:
        """Validate similarity threshold.

        Args:
            threshold: Threshold value to validate

        Raises:
            TypeError: If threshold is not a number
            ValueError: If threshold is not between 0 and 1
        """
        if not isinstance(threshold, (int, float)):
            raise TypeError("Threshold must be a number")
        if not 0 <= threshold <= 1:
            raise ValueError("Threshold must be between 0 and 1")

    def _get_validation_config(self) -> Dict[str, Any]:
        """Get validation configuration.

        Returns:
            Validation configuration dictionary
        """
        return {
            "swissprot": {
                "required_fields": ["accession", "name", "sequence"],
                "field_types": {
                    "accession": str,
                    "name": str,
                    "sequence": str,
                },
            },
            "swissadme": {
                "required_fields": ["smiles", "properties"],
                "field_types": {
                    "smiles": str,
                    "properties": dict,
                },
            },
            "swisssimilarity": {
                "required_fields": ["query", "results"],
                "field_types": {
                    "query": dict,
                    "results": list,
                },
            },
            "swisstargetprediction": {
                "required_fields": ["predictions"],
                "field_types": {
                    "predictions": list,
                },
            },
        }

    def _validate_swissprot_response(self, response: Dict[str, Any]) -> None:
        """Validate SwissProt API response.

        Args:
            response: Response dictionary to validate

        Raises:
            ValidationError: If response is invalid
        """
        if not isinstance(response, dict):
            raise ValidationError("Invalid response type")
        if "data" not in response or "entry" not in response["data"]:
            raise ValidationError("Missing required fields")
        entry = response["data"]["entry"]
        if not isinstance(entry.get("accession"), str):
            raise ValidationError("Invalid field types")

    def _validate_swissadme_response(self, response: Dict[str, Any]) -> None:
        """Validate SwissADME API response.

        Args:
            response: Response dictionary to validate

        Raises:
            ValidationError: If response is invalid
        """
        if not isinstance(response, dict):
            raise ValidationError("Invalid response type")
        if "data" not in response or "compound" not in response["data"]:
            raise ValidationError("Missing required fields")
        compound = response["data"]["compound"]
        if not isinstance(compound.get("smiles"), str):
            raise ValidationError("Invalid field types")

    def _validate_swisssimilarity_response(self, response: Dict[str, Any]) -> None:
        """Validate SwissSimilarity API response.

        Args:
            response: Response dictionary to validate

        Raises:
            ValidationError: If response is invalid
        """
        if not isinstance(response, dict):
            raise ValidationError("Invalid response type")
        if "data" not in response or "query" not in response["data"]:
            raise ValidationError("Missing required fields")
        query = response["data"]["query"]
        if not isinstance(query.get("smiles"), str):
            raise ValidationError("Invalid field types")

    def _validate_swisstargetprediction_response(self, response: Dict[str, Any]) -> None:
        """Validate SwissTargetPrediction API response.

        Args:
            response: Response dictionary to validate

        Raises:
            ValidationError: If response is invalid
        """
        if not isinstance(response, dict):
            raise ValidationError("Invalid response type")
        if "data" not in response or "predictions" not in response["data"]:
            raise ValidationError("Missing required fields")
        if not isinstance(response["data"]["predictions"], list):
            raise ValidationError("Invalid field types")

    def _get_cached_data(self, data_type: str, key: str) -> Optional[Dict[str, Any]]:
        """Get cached data.

        Args:
            data_type: Type of data to retrieve
            key: Cache key (SMILES or accession)

        Returns:
            Cached data or None if not found
        """
        if not self.cache_dir:
            return None

        try:
            cache_file = self.cache_dir / f"{data_type}_{key}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                data = pd.read_json(f)
                return data.to_dict("records") if len(data) > 1 else data.iloc[0].to_dict()

        except Exception as e:
            self.logger.error(f"Error reading cached {data_type}: {str(e)}")
            return None

    def get_metrics(self) -> Dict[str, Any]:
        """Get client metrics.

        Returns:
            Dictionary of metrics including:
            - Number of processed compounds
            - Number of failed compounds
            - Success rate
            - Circuit breaker metrics
        """
        metrics = super().get_metrics()
        metrics.update(
            {
                "processed_compounds": len(self.processed_compounds),
                "failed_compounds": len(self.failed_compounds),
                "success_rate": (
                    len(self.processed_compounds) / (len(self.processed_compounds) + len(self.failed_compounds)) if self.processed_compounds or self.failed_compounds else 0
                ),
            }
        )
        return metrics
