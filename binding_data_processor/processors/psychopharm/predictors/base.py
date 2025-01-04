"""Base predictor for psychopharmacological properties.

This module provides the base PredictorBase class that defines:
1. Common ML model loading and management
2. Feature extraction and preprocessing
3. Prediction confidence scoring
4. Result validation and formatting
5. Model retraining capabilities
"""

import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
import numpy as np

from ....models.core import CompoundData
from ...structure.ml.features import EnhancedFeatureExtractor
from ...structure.ml.fingerprints import FingerprintGenerator
from ...structure.ml.descriptors import DescriptorGenerator
from ....utils.cache import CacheManager
from ..base import PredictionResult


class PredictorBase(ABC):
    """Base class for ML-based predictors."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        feature_types: Optional[List[str]] = None,
    ):
        """Initialize predictor with model and feature configuration."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.model_dir = Path(model_dir) if model_dir else None
        self.cache = CacheManager(cache_dir) if cache_dir else None
        
        # Initialize feature generators
        self.feature_types = feature_types or ["fingerprints", "descriptors"]
        self.feature_extractors = {
            "fingerprints": FingerprintGenerator(),
            "descriptors": DescriptorGenerator(),
            "enhanced": EnhancedFeatureExtractor(),
        }

        # Load models
        self.models = self._load_models()
        self.threshold = 0.5  # Default confidence threshold

    @abstractmethod
    def _load_models(self) -> Dict:
        """Load ML models from disk.
        
        Returns:
            Dict mapping model names to loaded model objects
        """
        pass

    @abstractmethod
    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate predictions for a compound.
        
        Args:
            compound: CompoundData object to generate predictions for
            
        Returns:
            PredictionResult containing prediction value and confidence
        """
        pass

    def _extract_features(
        self, compound: CompoundData, feature_type: str
    ) -> np.ndarray:
        """Extract features of specified type from compound.
        
        Args:
            compound: CompoundData object to extract features from
            feature_type: Type of features to extract ("fingerprints", "descriptors", or "enhanced")
            
        Returns:
            numpy array of extracted features
        """
        try:
            # Try to load from cache first
            if self.cache:
                cache_key = f"{compound.inchi_key}_{feature_type}"
                cached = self.cache.get(cache_key)
                if cached is not None:
                    return cached

            # Extract features
            extractor = self.feature_extractors.get(feature_type)
            if not extractor:
                raise ValueError(f"Unknown feature type: {feature_type}")

            features = extractor.generate(compound.smiles)

            # Cache features
            if self.cache:
                self.cache.set(cache_key, features)

            return features

        except Exception as e:
            self.logger.error(
                f"Feature extraction error ({feature_type}): {str(e)}"
            )
            return np.array([])

    def _combine_predictions(
        self,
        predictions: List[Tuple[Any, float]],
        weights: Optional[List[float]] = None,
    ) -> Tuple[Any, float]:
        """Combine multiple predictions with optional weights.
        
        Args:
            predictions: List of (value, confidence) tuples
            weights: Optional list of weights for each prediction
            
        Returns:
            Tuple of (combined_value, combined_confidence)
        """
        if not predictions:
            return None, 0.0

        if len(predictions) == 1:
            return predictions[0]

        # Use equal weights if none provided
        if weights is None:
            weights = [1.0 / len(predictions)] * len(predictions)

        # Separate values and confidences
        values, confidences = zip(*predictions)

        # For categorical predictions
        if isinstance(values[0], (str, bool)):
            # Weight predictions by confidence
            weighted_counts = {}
            for val, conf, weight in zip(values, confidences, weights):
                weighted_counts[val] = (
                    weighted_counts.get(val, 0) + conf * weight
                )

            # Get value with highest weighted confidence
            best_value = max(
                weighted_counts.items(), key=lambda x: x[1]
            )[0]
            
            # Calculate overall confidence
            confidence = sum(
                conf * weight
                for val, conf, weight in zip(values, confidences, weights)
                if val == best_value
            )

            return best_value, confidence

        # For numeric predictions
        else:
            # Weighted average of values
            weighted_value = sum(
                val * conf * weight
                for val, conf, weight in zip(values, confidences, weights)
            ) / sum(conf * weight for conf, weight in zip(confidences, weights))

            # Average confidence
            confidence = sum(
                conf * weight for conf, weight in zip(confidences, weights)
            )

            return weighted_value, confidence

    def _validate_prediction(
        self,
        value: Any,
        confidence: float,
        min_confidence: float = 0.5,
    ) -> Tuple[bool, str]:
        """Validate prediction value and confidence.
        
        Args:
            value: Prediction value to validate
            confidence: Confidence score to validate
            min_confidence: Minimum required confidence threshold
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if confidence < min_confidence:
            return False, f"Low confidence: {confidence:.2f}"

        if value is None:
            return False, "No prediction value"

        if isinstance(value, (int, float)):
            if not np.isfinite(value):
                return False, "Invalid numeric value"

        return True, ""

    def _format_supporting_data(
        self,
        features: Dict[str, np.ndarray],
        predictions: List[Tuple[Any, float]],
        importances: Optional[Dict[str, float]] = None,
    ) -> Dict:
        """Format supporting data for prediction result.
        
        Args:
            features: Dict mapping feature types to feature arrays
            predictions: List of (value, confidence) prediction tuples
            importances: Optional dict mapping feature names to importance scores
            
        Returns:
            Dict containing formatted supporting data
        """
        return {
            "feature_types": list(features.keys()),
            "feature_counts": {
                k: len(v) for k, v in features.items()
            },
            "model_predictions": [
                {"value": val, "confidence": conf}
                for val, conf in predictions
            ],
            "feature_importances": importances or {},
        }

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save models to disk.
        
        Args:
            save_dir: Optional directory to save models to. If not provided,
                     uses the model_dir specified during initialization.
        """
        if not self.models:
            return

        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            raise ValueError("No save directory specified")

        save_dir.mkdir(parents=True, exist_ok=True)
        
        for name, model in self.models.items():
            model_path = save_dir / f"{name}.pkl"
            np.save(model_path, model)

    def retrain(
        self,
        compounds: List[CompoundData],
        labels: List[Any],
        **kwargs
    ) -> Dict[str, float]:
        """Retrain models with new data.
        
        Args:
            compounds: List of compounds to train on
            labels: List of corresponding labels
            **kwargs: Additional training parameters
            
        Returns:
            Dict containing training metrics
            
        Raises:
            NotImplementedError: If retraining is not supported
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} does not support retraining"
        )
