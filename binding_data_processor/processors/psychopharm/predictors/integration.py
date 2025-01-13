"""Integration utilities for predictors.

This module provides shared integration functionality used by multiple predictors,
helping to avoid circular dependencies between predictor modules.
"""

from typing import Optional
import logging
from pathlib import Path

from ....models.core import CompoundData
from ....models.compound.ml.predictors import PredictionResult


class BBBIntegrationMixin:
    """Mixin class providing BBB integration functionality."""

    def __init__(
        self,
        bbb_model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize BBB integration.

        Args:
            bbb_model_dir: Optional directory containing BBB models
            cache_dir: Optional directory for caching
            log_level: Logging level
        """
        self.bbb_model_dir = Path(bbb_model_dir) if bbb_model_dir else None
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.log_level = log_level

    def predict_bbb(self, compound: CompoundData) -> PredictionResult:
        """Predict BBB permeability for a compound.

        This method should be implemented by classes using this mixin.

        Args:
            compound: Compound to predict BBB permeability for

        Returns:
            PredictionResult containing BBB prediction
        """
        raise NotImplementedError

    def adjust_cns_prediction(
        self,
        base_prediction: float,
        base_confidence: float,
        bbb_result: PredictionResult,
        is_cns_related: bool = True,
    ) -> tuple[float, float]:
        """Adjust prediction based on BBB permeability.

        Args:
            base_prediction: Original prediction value
            base_confidence: Original confidence value
            bbb_result: BBB permeability prediction result
            is_cns_related: Whether the prediction is CNS-related

        Returns:
            Tuple of (adjusted_prediction, adjusted_confidence)
        """
        if base_confidence > 0.5:  # Base confidence threshold
            if is_cns_related:
                # Scale prediction by BBB permeability
                adjusted_prediction = base_prediction * bbb_result.confidence
                # Combine confidences
                adjusted_confidence = base_confidence * bbb_result.confidence
            else:
                # Non-CNS predictions less dependent on BBB
                adjusted_prediction = base_prediction
                adjusted_confidence = base_confidence

            if adjusted_confidence > 0.3:  # Adjusted confidence threshold
                return adjusted_prediction, adjusted_confidence

        return base_prediction, base_confidence
