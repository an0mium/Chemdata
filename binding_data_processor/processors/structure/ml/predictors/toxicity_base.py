"""Base classes and interfaces for toxicity prediction.

This module provides:
1. Base toxicity predictor interface
2. Common toxicity type definitions
3. Risk level calculations
4. Confidence scoring
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
from rdkit import Chem
from torch_geometric.data import Batch, Data

from ..base import PredictorBase


class ToxicityPredictorBase(PredictorBase, ABC):
    """Base class for toxicity prediction."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        uncertainty: bool = True,
        device: Optional[str] = None,
    ):
        """Initialize base predictor.

        Args:
            model_dir: Directory containing pre-trained models
            uncertainty: Whether to estimate uncertainty
            device: Device to run models on
        """
        super().__init__()
        self.model_dir = model_dir
        self.uncertainty = uncertainty
        self.device = device or "cuda" if torch.cuda.is_available() else "cpu"

    @abstractmethod
    def predict(
        self,
        compound: Union[str, Chem.Mol, Data, Batch],
        confidence_threshold: float = 0.5,
    ) -> Dict:
        """Predict compound toxicity.

        Args:
            compound: Input compound
            confidence_threshold: Minimum confidence threshold

        Returns:
            Dictionary of predictions
        """
        pass

    @abstractmethod
    def predict_batch(
        self,
        compounds: List[Union[str, Chem.Mol, Data]],
        batch_size: int = 32,
        confidence_threshold: float = 0.5,
    ) -> List[Dict]:
        """Predict toxicity for multiple compounds.

        Args:
            compounds: List of compounds
            batch_size: Batch size for predictions
            confidence_threshold: Minimum confidence threshold

        Returns:
            List of prediction dictionaries
        """
        pass

    def _calculate_confidence(
        self,
        predictions: torch.Tensor,
        uncertainties: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Calculate prediction confidence scores.

        Args:
            predictions: Model predictions
            uncertainties: Optional uncertainty estimates

        Returns:
            Confidence scores
        """
        if uncertainties is not None:
            # Use uncertainty-based confidence
            confidence = 1.0 / (1.0 + uncertainties)
        else:
            # Use prediction probability-based confidence
            confidence = torch.sigmoid(predictions)

        return confidence

    def _calculate_risk_level(self, probability: float) -> str:
        """Calculate risk level based on probability.

        Args:
            probability: Prediction probability

        Returns:
            Risk level string
        """
        if probability >= 0.8:
            return "very_high"
        elif probability >= 0.6:
            return "high"
        elif probability >= 0.4:
            return "moderate"
        elif probability >= 0.2:
            return "low"
        else:
            return "very_low"

    def _calculate_severity(self, probability: float) -> str:
        """Calculate severity level based on probability.

        Args:
            probability: Prediction probability

        Returns:
            Severity level string
        """
        if probability >= 0.8:
            return "severe"
        elif probability >= 0.6:
            return "high"
        elif probability >= 0.4:
            return "moderate"
        elif probability >= 0.2:
            return "mild"
        else:
            return "minimal"

    def _calculate_confidence_interval(
        self,
        prediction: float,
        uncertainty: float,
        confidence_level: float = 0.95,
    ) -> Tuple[float, float]:
        """Calculate confidence interval for prediction.

        Args:
            prediction: Predicted value
            uncertainty: Uncertainty estimate
            confidence_level: Confidence level (default: 95%)

        Returns:
            Tuple of (lower_bound, upper_bound)
        """
        z_score = {
            0.90: 1.645,
            0.95: 1.96,
            0.99: 2.576,
        }.get(confidence_level, 1.96)

        margin = z_score * uncertainty
        return max(0.0, prediction - margin), min(1.0, prediction + margin)
