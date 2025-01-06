"""Enhanced validation for web-enriched data.

This module provides advanced validation capabilities:
1. Cross-source validation - Compare data across different sources
2. Temporal validation - Track data changes over time
3. Semantic validation - Validate meaning and relationships
4. Community consensus validation - Validate against community data
"""

from typing import Dict, List, Optional, Any, Set
from dataclasses import dataclass
from datetime import datetime, timedelta
import re

from sentence_transformers import SentenceTransformer
import torch
import numpy as np

from .schema import BaseSchema
from ..llm_utils import NootropicLLMExtractor


@dataclass
class ValidationResult:
    """Result of validation checks."""

    is_valid: bool
    confidence: float
    source_agreements: Dict[str, float]
    temporal_consistency: float
    semantic_validity: float
    community_consensus: float
    validation_date: str
    metadata: Dict[str, Any]


@dataclass
class ValidationConfig:
    """Configuration for validation checks."""

    min_source_agreement: float = 0.7
    min_temporal_consistency: float = 0.8
    min_semantic_validity: float = 0.6
    min_community_consensus: float = 0.5
    max_temporal_gap: timedelta = timedelta(days=365)
    required_sources: Set[str] = None
    semantic_threshold: float = 0.8


class EnhancedValidator:
    """Enhanced validation for web-enriched data."""

    def __init__(
        self,
        config: Optional[ValidationConfig] = None,
        llm_provider: str = "ollama/llama2",
        api_token: Optional[str] = None,
    ):
        """Initialize validator.

        Args:
            config: Optional validation configuration
            llm_provider: LLM provider for semantic validation
            api_token: Optional API token for LLM provider
        """
        self.config = config or ValidationConfig()
        self.llm_extractor = NootropicLLMExtractor(
            llm_provider=llm_provider,
            api_token=api_token,
        )
        self.semantic_model = SentenceTransformer("pritamdeka/S-PubMedBert-MS-MARCO")

    async def validate(
        self,
        data: BaseSchema,
        source_data: Dict[str, Any],
        temporal_data: List[Dict[str, Any]],
        community_data: List[Dict[str, Any]],
    ) -> ValidationResult:
        """Validate data using enhanced checks.

        Args:
            data: Data to validate
            source_data: Data from different sources
            temporal_data: Historical data points
            community_data: Community-sourced data

        Returns:
            Validation result with confidence scores
        """
        try:
            # Run all validation checks
            source_agreements = await self._validate_cross_source(data, source_data)
            temporal_consistency = self._validate_temporal(data, temporal_data)
            semantic_validity = await self._validate_semantic(data)
            community_consensus = self._validate_community(data, community_data)

            # Calculate overall validity
            scores = [
                (np.mean(list(source_agreements.values())), self.config.min_source_agreement),
                (temporal_consistency, self.config.min_temporal_consistency),
                (semantic_validity, self.config.min_semantic_validity),
                (community_consensus, self.config.min_community_consensus),
            ]

            is_valid = all(score >= threshold for score, threshold in scores)
            confidence = np.mean([score for score, _ in scores])

            return ValidationResult(
                is_valid=is_valid,
                confidence=float(confidence),
                source_agreements=source_agreements,
                temporal_consistency=float(temporal_consistency),
                semantic_validity=float(semantic_validity),
                community_consensus=float(community_consensus),
                validation_date=datetime.now().isoformat(),
                metadata={
                    "config": self.config.__dict__,
                    "validation_type": "enhanced",
                    "data_sources": list(source_data.keys()),
                    "temporal_points": len(temporal_data),
                    "community_points": len(community_data),
                },
            )

        except Exception as e:
            self.logger.error(f"Validation error: {str(e)}")
            return ValidationResult(
                is_valid=False,
                confidence=0.0,
                source_agreements={},
                temporal_consistency=0.0,
                semantic_validity=0.0,
                community_consensus=0.0,
                validation_date=datetime.now().isoformat(),
                metadata={"error": str(e)},
            )

    async def _validate_cross_source(
        self,
        data: BaseSchema,
        source_data: Dict[str, Any],
    ) -> Dict[str, float]:
        """Validate data across different sources.

        Args:
            data: Data to validate
            source_data: Data from different sources

        Returns:
            Agreement scores per source
        """
        agreements = {}

        # Check required sources
        if self.config.required_sources:
            missing = self.config.required_sources - set(source_data.keys())
            if missing:
                self.logger.warning(f"Missing required sources: {missing}")

        # Compare with each source
        for source, source_value in source_data.items():
            try:
                # Get embeddings
                data_embedding = self.semantic_model.encode(str(data))
                source_embedding = self.semantic_model.encode(str(source_value))

                # Calculate semantic similarity
                similarity = float(
                    torch.nn.functional.cosine_similarity(
                        torch.tensor(data_embedding).unsqueeze(0),
                        torch.tensor(source_embedding).unsqueeze(0),
                        dim=1,
                    )
                )

                agreements[source] = similarity

            except Exception as e:
                self.logger.error(f"Error comparing with source {source}: {str(e)}")
                agreements[source] = 0.0

        return agreements

    def _validate_temporal(
        self,
        data: BaseSchema,
        temporal_data: List[Dict[str, Any]],
    ) -> float:
        """Validate temporal consistency.

        Args:
            data: Data to validate
            temporal_data: Historical data points

        Returns:
            Temporal consistency score
        """
        if not temporal_data:
            return 1.0  # No history to compare against

        try:
            # Sort by date
            temporal_data = sorted(
                temporal_data,
                key=lambda x: datetime.fromisoformat(x["timestamp"]),
            )

            # Check temporal gaps
            current_time = datetime.now()
            last_update = datetime.fromisoformat(temporal_data[-1]["timestamp"])
            if current_time - last_update > self.config.max_temporal_gap:
                self.logger.warning(f"Data may be outdated: {current_time - last_update}")

            # Calculate consistency score
            consistencies = []
            for point in temporal_data:
                # Get embeddings
                data_embedding = self.semantic_model.encode(str(data))
                point_embedding = self.semantic_model.encode(str(point["data"]))

                # Calculate semantic similarity
                similarity = float(
                    torch.nn.functional.cosine_similarity(
                        torch.tensor(data_embedding).unsqueeze(0),
                        torch.tensor(point_embedding).unsqueeze(0),
                        dim=1,
                    )
                )
                consistencies.append(similarity)

            # Weight recent points more heavily
            weights = np.linspace(0.5, 1.0, len(consistencies))
            return float(np.average(consistencies, weights=weights))

        except Exception as e:
            self.logger.error(f"Error validating temporal consistency: {str(e)}")
            return 0.0

    async def _validate_semantic(self, data: BaseSchema) -> float:
        """Validate semantic relationships.

        Args:
            data: Data to validate

        Returns:
            Semantic validity score
        """
        try:
            # Extract semantic information
            effects = await self.llm_extractor.extract_effects(str(data))
            mechanisms = await self.llm_extractor.extract_mechanisms(str(data))
            safety = await self.llm_extractor.extract_safety(str(data))

            # Calculate validity scores
            scores = []

            # Check effects-mechanisms consistency
            if effects and mechanisms:
                effect_embedding = self.semantic_model.encode(str(effects))
                mech_embedding = self.semantic_model.encode(str(mechanisms))
                similarity = float(
                    torch.nn.functional.cosine_similarity(
                        torch.tensor(effect_embedding).unsqueeze(0),
                        torch.tensor(mech_embedding).unsqueeze(0),
                        dim=1,
                    )
                )
                scores.append(similarity)

            # Check effects-safety consistency
            if effects and safety:
                effect_embedding = self.semantic_model.encode(str(effects))
                safety_embedding = self.semantic_model.encode(str(safety))
                similarity = float(
                    torch.nn.functional.cosine_similarity(
                        torch.tensor(effect_embedding).unsqueeze(0),
                        torch.tensor(safety_embedding).unsqueeze(0),
                        dim=1,
                    )
                )
                scores.append(similarity)

            # Check mechanisms-safety consistency
            if mechanisms and safety:
                mech_embedding = self.semantic_model.encode(str(mechanisms))
                safety_embedding = self.semantic_model.encode(str(safety))
                similarity = float(
                    torch.nn.functional.cosine_similarity(
                        torch.tensor(mech_embedding).unsqueeze(0),
                        torch.tensor(safety_embedding).unsqueeze(0),
                        dim=1,
                    )
                )
                scores.append(similarity)

            return float(np.mean(scores)) if scores else 0.0

        except Exception as e:
            self.logger.error(f"Error validating semantic relationships: {str(e)}")
            return 0.0

    def _validate_community(
        self,
        data: BaseSchema,
        community_data: List[Dict[str, Any]],
    ) -> float:
        """Validate against community consensus.

        Args:
            data: Data to validate
            community_data: Community-sourced data

        Returns:
            Community consensus score
        """
        if not community_data:
            return 1.0  # No community data to validate against

        try:
            # Calculate consensus scores
            scores = []
            for point in community_data:
                # Get embeddings
                data_embedding = self.semantic_model.encode(str(data))
                point_embedding = self.semantic_model.encode(str(point["data"]))

                # Calculate semantic similarity
                similarity = float(
                    torch.nn.functional.cosine_similarity(
                        torch.tensor(data_embedding).unsqueeze(0),
                        torch.tensor(point_embedding).unsqueeze(0),
                        dim=1,
                    )
                )

                # Weight by point confidence
                confidence = point.get("confidence", 0.5)
                scores.append(similarity * confidence)

            return float(np.mean(scores)) if scores else 0.0

        except Exception as e:
            self.logger.error(f"Error validating community consensus: {str(e)}")
            return 0.0
