"""Swiss tools client for web enrichment.

This module provides a client for:
- SwissTargetPrediction (target prediction)
- SwissADME (ADME properties)
- SwissSimilarity (structure similarity)
"""

import logging
from typing import Any, Dict, List, Optional, Set

from .base import WebClient, WebClientError, ValidationError
from ..validation.schema import DataSource, ValidationLevel
from ..validation.data import ValidationConfig
from ...pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig


class SwissValidationError(ValidationError):
    """Swiss validation error."""

    def __init__(
        self,
        message: str,
        issues: List[Any],
        tool: str,
    ):
        """Initialize error.

        Args:
            message: Error message
            issues: Validation issues
            tool: Swiss tool name
        """
        super().__init__(
            message=message,
            source=f"swiss.{tool}",
            issues=issues,
        )
        self.tool = tool


class SwissClient(WebClient):
    """Client for Swiss bioinformatics tools."""

    # Tool endpoints
    ENDPOINTS = {
        "target": "predict",
        "adme": "calculate",
        "similarity": "compare",
    }

    # Validation rules
    RULES = {
        "smiles": [
            ("required", None),
            ("pattern", r"^[A-Za-z0-9\[\]\(\)\{\}\.=#@\-\+]+$"),
        ],
        "inchi": [
            ("pattern", r"^InChI=1S?/[A-Za-z0-9/.]+$"),
        ],
        "threshold": [
            ("range", (0.0, 1.0)),
        ],
    }

    def __init__(
        self,
        base_url: str = "https://www.swisstargetprediction.ch/api",
        requests_per_second: float = 1.0,
        max_retries: int = 3,
        timeout: int = 30,
        cache_ttl: int = 3600,
        circuit_config: Optional[CircuitBreakerConfig] = None,
        validation_level: ValidationLevel = ValidationLevel.NORMAL,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize Swiss client.

        Args:
            base_url: Base URL for Swiss tools API
            requests_per_second: Maximum requests per second
            max_retries: Maximum number of retries
            timeout: Request timeout in seconds
            cache_ttl: Cache TTL in seconds
            circuit_config: Optional circuit breaker config
            validation_level: Validation strictness level
            logger: Optional logger instance
        """
        super().__init__(
            name="swiss",
            base_url=base_url,
            data_source=DataSource.SWISS,
            requests_per_second=requests_per_second,
            max_retries=max_retries,
            timeout=timeout,
            cache_ttl=cache_ttl,
            circuit_config=circuit_config,
            validation_level=validation_level,
            logger=logger,
        )

        # Track processed compounds
        self.processed_compounds: Set[str] = set()
        self.failed_compounds: Set[str] = set()

    def _get_validation_config(self) -> ValidationConfig:
        """Get validation configuration.

        Returns:
            Validation configuration
        """
        return ValidationConfig(rules=self.RULES)

    def predict_targets(
        self,
        smiles: str,
        threshold: float = 0.5,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Predict protein targets for compound.

        Args:
            smiles: Compound SMILES
            threshold: Probability threshold
            use_cache: Whether to use cache

        Returns:
            Target predictions

        Raises:
            SwissValidationError: For validation errors
            WebClientError: For other errors
        """
        try:
            # Make request
            data = self.request(
                method="POST",
                endpoint=self.ENDPOINTS["target"],
                json={
                    "smiles": smiles,
                    "threshold": threshold,
                },
                use_cache=use_cache,
                schema_name="target",
            )

            # Track success
            self.processed_compounds.add(smiles)
            return data

        except ValidationError as e:
            # Track failure
            self.failed_compounds.add(smiles)
            raise SwissValidationError(
                message=str(e),
                issues=e.issues,
                tool="target",
            )

        except WebClientError as e:
            # Track failure
            self.failed_compounds.add(smiles)
            # Re-raise with original error info
            raise WebClientError(
                message=str(e),
                source=self.name,
                status_code=e.status_code,
                details=e.details,
            )

    def calculate_adme(
        self,
        smiles: str,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Calculate ADME properties for compound.

        Args:
            smiles: Compound SMILES
            use_cache: Whether to use cache

        Returns:
            ADME properties

        Raises:
            SwissValidationError: For validation errors
            WebClientError: For other errors
        """
        try:
            # Make request
            data = self.request(
                method="POST",
                endpoint=self.ENDPOINTS["adme"],
                json={"smiles": smiles},
                use_cache=use_cache,
                schema_name="adme",
            )

            # Track success
            self.processed_compounds.add(smiles)
            return data

        except ValidationError as e:
            # Track failure
            self.failed_compounds.add(smiles)
            raise SwissValidationError(
                message=str(e),
                issues=e.issues,
                tool="adme",
            )

        except WebClientError as e:
            # Track failure
            self.failed_compounds.add(smiles)
            # Re-raise with original error info
            raise WebClientError(
                message=str(e),
                source=self.name,
                status_code=e.status_code,
                details=e.details,
            )

    def compare_structures(
        self,
        smiles: str,
        reference_smiles: str,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Compare compound structures.

        Args:
            smiles: Query compound SMILES
            reference_smiles: Reference compound SMILES
            use_cache: Whether to use cache

        Returns:
            Similarity scores

        Raises:
            SwissValidationError: For validation errors
            WebClientError: For other errors
        """
        try:
            # Make request
            data = self.request(
                method="POST",
                endpoint=self.ENDPOINTS["similarity"],
                json={
                    "smiles": smiles,
                    "reference": reference_smiles,
                },
                use_cache=use_cache,
                schema_name="similarity",
            )

            # Track success
            self.processed_compounds.add(smiles)
            return data

        except ValidationError as e:
            # Track failure
            self.failed_compounds.add(smiles)
            raise SwissValidationError(
                message=str(e),
                issues=e.issues,
                tool="similarity",
            )

        except WebClientError as e:
            # Track failure
            self.failed_compounds.add(smiles)
            # Re-raise with original error info
            raise WebClientError(
                message=str(e),
                source=self.name,
                status_code=e.status_code,
                details=e.details,
            )

    def get_metrics(self) -> Dict[str, Any]:
        """Get client metrics.

        Returns:
            Client metrics
        """
        metrics = super().get_metrics()
        metrics.update(
            {
                "processed_compounds": len(self.processed_compounds),
                "failed_compounds": len(self.failed_compounds),
            }
        )
        return metrics
