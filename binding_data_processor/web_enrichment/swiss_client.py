"""Swiss tools integration client.

This module provides functionality to:
1. Get similar compounds from SwissSimilarity
2. Get protein data from SwissProt
3. Get target predictions from SwissTargetPrediction
4. Get ADME properties from SwissADME
5. Process and validate results
"""

import logging
import re
from typing import Optional, Dict, Any, Union

from .base_client import BaseWebClient
from .base import WebClientError
from .http_client import HTTPClient
from ..pipeline.infrastructure.rate_limiter import RateLimiter, RateLimit
from ..pipeline.infrastructure.circuit_breaker import CircuitBreaker, ErrorCode, McpError

# API Rate Limits
SWISSPROT_RATE_LIMIT = RateLimit(requests=3, period=1)  # 3 requests per second
SWISSADME_RATE_LIMIT = RateLimit(requests=2, period=1)  # 2 requests per second
SWISSSIMILARITY_RATE_LIMIT = RateLimit(requests=2, period=1)  # 2 requests per second
SWISSTARGET_RATE_LIMIT = RateLimit(requests=2, period=1)  # 2 requests per second

# Circuit Breaker Config
CIRCUIT_FAILURE_THRESHOLD = 5
CIRCUIT_RESET_TIMEOUT = 300  # 5 minutes
CIRCUIT_HALF_OPEN_TIMEOUT = 60  # 1 minute


class SwissClient(BaseWebClient):
    """Client for Swiss bioinformatics tools."""

    def __init__(
        self,
        http_client: Optional[HTTPClient] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize Swiss tools client.

        Args:
            http_client: Optional HTTP client to use
            logger: Optional logger instance
        """
        super().__init__(name="swiss", http_client=http_client, logger=logger)

        # API endpoints
        self.swissprot_url = "https://rest.uniprot.org/uniprotkb"
        self.swissadme_url = "https://www.swissadme.ch/api"
        self.swisssimilarity_url = "https://www.swisssimilarity.ch/api"
        self.swisstargetprediction_url = "https://www.swisstargetprediction.ch/api"

        # Initialize rate limiters
        self.rate_limiter = RateLimiter()
        self.rate_limiter.add_limit("swissprot", SWISSPROT_RATE_LIMIT)
        self.rate_limiter.add_limit("swissadme", SWISSADME_RATE_LIMIT)
        self.rate_limiter.add_limit("swisssimilarity", SWISSSIMILARITY_RATE_LIMIT)
        self.rate_limiter.add_limit("swisstarget", SWISSTARGET_RATE_LIMIT)

        # Initialize circuit breakers
        self.swissprot_circuit = CircuitBreaker(
            "swissprot",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )
        self.swissadme_circuit = CircuitBreaker(
            "swissadme",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )
        self.swisssimilarity_circuit = CircuitBreaker(
            "swisssimilarity",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )
        self.swisstarget_circuit = CircuitBreaker(
            "swisstarget",
            failure_threshold=CIRCUIT_FAILURE_THRESHOLD,
            reset_timeout=CIRCUIT_RESET_TIMEOUT,
            half_open_timeout=CIRCUIT_HALF_OPEN_TIMEOUT,
        )

    def get_swissprot_data(
        self,
        accession: str,
        name: str,
    ) -> Dict[str, Any]:
        """Get protein data from SwissProt.

        Args:
            accession: UniProt accession number
            name: Protein name

        Returns:
            SwissProt data dictionary

        Raises:
            WebClientError: If API request fails
        """
        # Check circuit breaker
        if not self.swissprot_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "SwissProt circuit breaker open",
            )

        try:
            # Wait for rate limit
            self.rate_limiter.wait("swissprot")

            self._validate_accession(accession)
            self._validate_name(name)

            response = self.http.get(
                url=f"{self.swissprot_url}/{accession}",
                params={
                    "name": name,
                    "format": "json",
                    "include_features": True,
                },
                headers={"Accept": "application/json"},
            )

            self._validate_swissprot_response(response)
            self.swissprot_circuit.record_success()
            return response

        except Exception as e:
            self.swissprot_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"SwissProt request failed: {str(e)}",
            )

    def get_swissadme_data(
        self,
        smiles: str,
    ) -> Dict[str, Any]:
        """Get ADME properties from SwissADME.

        Args:
            smiles: SMILES string

        Returns:
            SwissADME data dictionary

        Raises:
            WebClientError: If API request fails
        """
        # Check circuit breaker
        if not self.swissadme_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "SwissADME circuit breaker open",
            )

        try:
            # Wait for rate limit
            self.rate_limiter.wait("swissadme")

            self._validate_smiles(smiles)

            response = self.http.get(
                url=f"{self.swissadme_url}/predict",
                params={
                    "smiles": smiles,
                    "format": "json",
                },
                headers={"Accept": "application/json"},
            )

            self._validate_swissadme_response(response)
            self.swissadme_circuit.record_success()
            return response

        except Exception as e:
            self.swissadme_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"SwissADME request failed: {str(e)}",
            )

    def get_swisssimilarity_data(
        self,
        smiles: str,
        name: str,
        threshold: float = 0.7,
        limit: int = 10,
    ) -> Dict[str, Any]:
        """Get similar compounds from SwissSimilarity.

        Args:
            smiles: SMILES string
            name: Compound name
            threshold: Similarity threshold (0-1)
            limit: Maximum number of results

        Returns:
            SwissSimilarity data dictionary

        Raises:
            WebClientError: If API request fails
        """
        # Check circuit breaker
        if not self.swisssimilarity_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "SwissSimilarity circuit breaker open",
            )

        try:
            # Wait for rate limit
            self.rate_limiter.wait("swisssimilarity")

            self._validate_smiles(smiles)
            self._validate_name(name)
            self._validate_threshold(threshold)

            response = self.http.get(
                url=f"{self.swisssimilarity_url}/similar",
                params={
                    "smiles": smiles,
                    "name": name,
                    "threshold": threshold,
                    "limit": limit,
                    "format": "json",
                },
                headers={"Accept": "application/json"},
            )

            self._validate_swisssimilarity_response(response)
            self.swisssimilarity_circuit.record_success()
            return response

        except Exception as e:
            self.swisssimilarity_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"SwissSimilarity request failed: {str(e)}",
            )

    def get_swisstargetprediction_data(
        self,
        smiles: str,
    ) -> Dict[str, Any]:
        """Get target predictions from SwissTargetPrediction.

        Args:
            smiles: SMILES string

        Returns:
            SwissTargetPrediction data dictionary

        Raises:
            WebClientError: If API request fails
        """
        # Check circuit breaker
        if not self.swisstarget_circuit.allow_request():
            raise McpError(
                ErrorCode.ServiceUnavailable,
                "SwissTargetPrediction circuit breaker open",
            )

        try:
            # Wait for rate limit
            self.rate_limiter.wait("swisstarget")

            self._validate_smiles(smiles)

            response = self.http.get(
                url=f"{self.swisstargetprediction_url}/predict",
                params={
                    "smiles": smiles,
                    "format": "json",
                },
                headers={"Accept": "application/json"},
            )

            self._validate_swisstargetprediction_response(response)
            self.swisstarget_circuit.record_success()
            return response

        except Exception as e:
            self.swisstarget_circuit.record_failure()
            raise McpError(
                ErrorCode.ServiceError,
                f"SwissTargetPrediction request failed: {str(e)}",
            )

    def get_all_data(
        self,
        smiles: str,
        name: str,
        accession: str,
        threshold: float = 0.7,
        limit: int = 10,
    ) -> Dict[str, Any]:
        """Get data from all Swiss tools.

        Args:
            smiles: SMILES string
            name: Compound name
            accession: UniProt accession
            threshold: Similarity threshold
            limit: Maximum similarity results

        Returns:
            Dictionary containing all Swiss tools data
        """
        results = {}

        try:
            results["swissprot"] = self.get_swissprot_data(
                accession=accession,
                name=name,
            )
        except WebClientError as e:
            self.logger.error(f"SwissProt error: {str(e)}")
            results["swissprot"] = None

        try:
            results["swissadme"] = self.get_swissadme_data(
                smiles=smiles,
            )
        except WebClientError as e:
            self.logger.error(f"SwissADME error: {str(e)}")
            results["swissadme"] = None

        try:
            results["swisssimilarity"] = self.get_swisssimilarity_data(
                smiles=smiles,
                name=name,
                threshold=threshold,
                limit=limit,
            )
        except WebClientError as e:
            self.logger.error(f"SwissSimilarity error: {str(e)}")
            results["swisssimilarity"] = None

        try:
            results["swisstargetprediction"] = self.get_swisstargetprediction_data(
                smiles=smiles,
            )
        except WebClientError as e:
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

    def _validate_swissprot_response(self, response: Dict[str, Any]) -> None:
        """Validate SwissProt API response.

        Args:
            response: Response dictionary to validate

        Raises:
            WebClientError: If response is invalid
        """
        if not isinstance(response, dict):
            raise WebClientError("Invalid response type")
        if "data" not in response or "entry" not in response["data"]:
            raise WebClientError("Missing required fields")
        entry = response["data"]["entry"]
        if not isinstance(entry.get("accession"), str):
            raise WebClientError("Invalid field types")

    def _validate_swissadme_response(self, response: Dict[str, Any]) -> None:
        """Validate SwissADME API response.

        Args:
            response: Response dictionary to validate

        Raises:
            WebClientError: If response is invalid
        """
        if not isinstance(response, dict):
            raise WebClientError("Invalid response type")
        if "data" not in response or "compound" not in response["data"]:
            raise WebClientError("Missing required fields")
        compound = response["data"]["compound"]
        if not isinstance(compound.get("smiles"), str):
            raise WebClientError("Invalid field types")

    def _validate_swisssimilarity_response(self, response: Dict[str, Any]) -> None:
        """Validate SwissSimilarity API response.

        Args:
            response: Response dictionary to validate

        Raises:
            WebClientError: If response is invalid
        """
        if not isinstance(response, dict):
            raise WebClientError("Invalid response type")
        if "data" not in response or "query" not in response["data"]:
            raise WebClientError("Missing required fields")
        query = response["data"]["query"]
        if not isinstance(query.get("smiles"), str):
            raise WebClientError("Invalid field types")

    def _validate_swisstargetprediction_response(self, response: Dict[str, Any]) -> None:
        """Validate SwissTargetPrediction API response.

        Args:
            response: Response dictionary to validate

        Raises:
            WebClientError: If response is invalid
        """
        if not isinstance(response, dict):
            raise WebClientError("Invalid response type")
        if "data" not in response or "predictions" not in response["data"]:
            raise WebClientError("Missing required fields")
        if not isinstance(response["data"]["predictions"], list):
            raise WebClientError("Invalid field types")
