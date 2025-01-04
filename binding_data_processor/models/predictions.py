"""ML prediction functionality.

This module provides:
1. PredictionsMixin - Methods for handling various ML predictions:
   - Toxicity predictions
   - Abuse potential assessment 
   - Binding affinity predictions
   - Activity predictions
   - Mechanism predictions

2. Base Classes:
   - PredictorBase - Base class for all predictors
   - FeatureExtractorBase - Base class for feature extraction
   - EnsemblePredictor - Class for combining multiple predictors

3. Data Classes:
   - PredictionResult - Container for prediction results
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional
import numpy as np
import pandas as pd


@dataclass
class PredictionResult:
    """Result of a model prediction.
    
    Attributes:
        value: The predicted value
        confidence: Confidence score between 0 and 1
        supporting_data: Additional data supporting the prediction
    """
    value: Any
    confidence: float
    supporting_data: Dict = field(default_factory=dict)


class PredictorBase:
    """Base class for ML predictors.
    
    Attributes:
        model: The underlying ML model
        feature_extractor: Component for extracting features
        scaler: Optional scaler for feature normalization
        feature_type: Type of features used by this predictor
    """

    def __init__(self):
        """Initialize predictor."""
        self.model = None
        self.feature_extractor = None
        self.scaler = None
        self.feature_type = None

    def predict(self, compound) -> PredictionResult:
        """Make prediction for compound.
        
        Args:
            compound: CompoundData instance
            
        Returns:
            PredictionResult with prediction value and confidence
        """
        # Get cached features if available
        features = compound.get_cached_features(self.feature_type)
        
        if features is None:
            # Extract features
            features = self.extract_features(compound)
            # Cache for reuse
            compound.cache_features(self.feature_type, features)
            
        # Scale features
        if self.scaler:
            features = self.scaler.transform(features.reshape(1, -1))
            
        # Make prediction
        prediction = self.model.predict(features)[0]
        confidence = self.get_confidence(features)
        
        # Get supporting data
        supporting_data = self.get_supporting_data(compound, features, prediction)
        
        return PredictionResult(
            value=prediction,
            confidence=confidence,
            supporting_data=supporting_data
        )

    def extract_features(self, compound) -> np.ndarray:
        """Extract features from compound.
        
        Args:
            compound: CompoundData instance
            
        Returns:
            Feature array
        """
        if self.feature_extractor:
            return self.feature_extractor.extract(compound)
        raise NotImplementedError
        
    def get_confidence(self, features: np.ndarray) -> float:
        """Get confidence score for prediction.
        
        Args:
            features: Feature array
            
        Returns:
            Confidence score between 0 and 1
        """
        # Default implementation returns 0.5
        return 0.5
        
    def get_supporting_data(
        self,
        compound,
        features: np.ndarray,
        prediction: Any
    ) -> Dict:
        """Get supporting data for prediction.
        
        Args:
            compound: CompoundData instance
            features: Feature array
            prediction: Model prediction
            
        Returns:
            Dictionary of supporting data
        """
        return {}


class FeatureExtractorBase:
    """Base class for feature extractors.
    
    Provides interface for extracting ML features from compounds.
    """

    def extract(self, compound) -> np.ndarray:
        """Extract features from compound.
        
        Args:
            compound: CompoundData instance
            
        Returns:
            Feature array
        """
        raise NotImplementedError


class EnsemblePredictor(PredictorBase):
    """Ensemble of multiple predictors.
    
    Combines predictions from multiple models using configurable methods.
    
    Attributes:
        predictors: List of predictors in ensemble
        weights: Optional weights for each predictor
    """

    def __init__(
        self,
        predictors: List[PredictorBase],
        weights: Optional[List[float]] = None
    ):
        """Initialize ensemble.
        
        Args:
            predictors: List of predictors to ensemble
            weights: Optional weights for each predictor
        """
        super().__init__()
        self.predictors = predictors
        self.weights = weights
        if weights and len(weights) != len(predictors):
            raise ValueError("Number of weights must match number of predictors")

    def predict(self, compound) -> PredictionResult:
        """Make ensemble prediction.
        
        Args:
            compound: CompoundData instance
            
        Returns:
            Ensemble prediction result
        """
        # Get individual predictions
        results = [
            predictor.predict(compound)
            for predictor in self.predictors
        ]
        
        # Combine predictions
        ensemble_value = self._combine_predictions(results)
        ensemble_confidence = self._combine_confidences(results)
        
        # Combine supporting data
        supporting_data = {
            f"predictor_{i}": result.supporting_data
            for i, result in enumerate(results)
        }
        supporting_data["ensemble_method"] = self._get_ensemble_method()
        supporting_data["weights"] = self.weights
        
        return PredictionResult(
            value=ensemble_value,
            confidence=ensemble_confidence,
            supporting_data=supporting_data
        )

    def _combine_predictions(
        self,
        results: List[PredictionResult]
    ) -> Any:
        """Combine individual predictions.
        
        Args:
            results: List of prediction results
            
        Returns:
            Combined prediction value
        """
        values = [result.value for result in results]
        if self.weights:
            return np.average(values, weights=self.weights)
        return np.mean(values)

    def _combine_confidences(
        self,
        results: List[PredictionResult]
    ) -> float:
        """Combine individual confidence scores.
        
        Args:
            results: List of prediction results
            
        Returns:
            Combined confidence score
        """
        confidences = [result.confidence for result in results]
        if self.weights:
            return np.average(confidences, weights=self.weights)
        return np.mean(confidences)

    def _get_ensemble_method(self) -> str:
        """Get description of ensemble method.
        
        Returns:
            String describing ensemble method
        """
        return "weighted_average" if self.weights else "mean"


class PredictionsMixin:
    """Mixin class providing ML prediction methods.
    
    Provides functionality for:
    - Storing predictions
    - Formatting predictions for output
    - Merging predictions from multiple sources
    """

    # ML predictions
    toxicity_predictions: Dict[str, Dict] = field(default_factory=dict)
    abuse_potential: Dict[str, Dict] = field(default_factory=dict)
    binding_predictions: List[Dict] = field(default_factory=list)
    activity_predictions: List[Dict] = field(default_factory=list)
    mechanism_predictions: List[Dict] = field(default_factory=list)

    # Prediction caching
    _prediction_cache: Dict[str, PredictionResult] = field(default_factory=dict)
    _feature_cache: Dict[str, np.ndarray] = field(default_factory=dict)
    _prediction_history: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(
        columns=[
            'predictor_type',
            'prediction_value',
            'confidence',
            'timestamp',
            'supporting_data',
        ]
    ))

    def get_predictions_dict(self) -> Dict:
        """Get dictionary of all ML predictions."""
        return {
            "toxicity_predictions": self._format_toxicity_predictions(),
            "abuse_potential": self._format_abuse_potential(),
            "binding_predictions": self._format_binding_predictions(),
            "activity_predictions": self._format_activity_predictions(),
            "mechanism_predictions": self._format_mechanism_predictions(),
        }

    def _format_toxicity_predictions(self) -> Dict:
        """Format toxicity predictions for output."""
        formatted = {}
        for tox_type, prediction in self.toxicity_predictions.items():
            formatted[tox_type] = {
                "probability": prediction.get("probability", 0),
                "confidence": prediction.get("confidence", 0),
                "severity": prediction.get("severity", "unknown"),
                "mechanisms": prediction.get("mechanisms", []),
                "supporting_data": prediction.get("supporting_data", {}),
            }
        return formatted

    def _format_abuse_potential(self) -> Dict:
        """Format abuse potential data for output."""
        formatted = {}
        for abuse_type, prediction in self.abuse_potential.items():
            formatted[abuse_type] = {
                "probability": prediction.get("probability", 0),
                "confidence": prediction.get("confidence", 0),
                "risk_level": prediction.get("risk_level", "unknown"),
                "mechanisms": prediction.get("mechanisms", []),
                "receptor_involvement": prediction.get("receptor_involvement", {}),
                "supporting_data": prediction.get("supporting_data", {}),
            }
        return formatted

    def _format_binding_predictions(self) -> List[Dict]:
        """Format binding predictions for output."""
        return [
            {
                "target": pred.get("target", ""),
                "probability": pred.get("probability", 0),
                "affinity_type": pred.get("affinity_type", ""),
                "predicted_value": pred.get("predicted_value"),
                "confidence": pred.get("confidence", 0),
                "supporting_features": pred.get("supporting_features", []),
            }
            for pred in self.binding_predictions
        ]

    def _format_activity_predictions(self) -> List[Dict]:
        """Format activity predictions for output."""
        return [
            {
                "activity_type": pred.get("activity_type", ""),
                "probability": pred.get("probability", 0),
                "predicted_value": pred.get("predicted_value"),
                "confidence": pred.get("confidence", 0),
                "mechanism": pred.get("mechanism", ""),
                "supporting_data": pred.get("supporting_data", {}),
            }
            for pred in self.activity_predictions
        ]

    def _format_mechanism_predictions(self) -> List[Dict]:
        """Format mechanism predictions for output."""
        return [
            {
                "mechanism": pred.get("mechanism", ""),
                "probability": pred.get("probability", 0),
                "confidence": pred.get("confidence", 0),
                "supporting_evidence": pred.get("supporting_evidence", []),
                "receptor_systems": pred.get("receptor_systems", []),
            }
            for pred in self.mechanism_predictions
        ]

    def merge_predictions(self, other: "PredictionsMixin") -> None:
        """Merge ML predictions from another instance."""
        self._merge_toxicity_predictions(other)
        self._merge_abuse_potential(other)
        self._merge_binding_predictions(other)
        self._merge_activity_predictions(other)
        self._merge_mechanism_predictions(other)

    def _merge_toxicity_predictions(self, other: "PredictionsMixin") -> None:
        """Merge toxicity predictions."""
        for tox_type, tox_data in other.toxicity_predictions.items():
            if tox_type not in self.toxicity_predictions:
                self.toxicity_predictions[tox_type] = tox_data
            else:
                # Keep prediction with higher confidence
                old_conf = self.toxicity_predictions[tox_type].get("confidence", 0)
                if tox_data.get("confidence", 0) > old_conf:
                    self.toxicity_predictions[tox_type] = tox_data

    def _merge_abuse_potential(self, other: "PredictionsMixin") -> None:
        """Merge abuse potential predictions."""
        for abuse_type, abuse_data in other.abuse_potential.items():
            if abuse_type not in self.abuse_potential:
                self.abuse_potential[abuse_type] = abuse_data
            else:
                # Keep prediction with higher confidence
                old_conf = self.abuse_potential[abuse_type].get("confidence", 0)
                if abuse_data.get("confidence", 0) > old_conf:
                    self.abuse_potential[abuse_type] = abuse_data

    def _merge_binding_predictions(self, other: "PredictionsMixin") -> None:
        """Merge binding predictions."""
        self.binding_predictions.extend(
            pred for pred in other.binding_predictions
            if pred not in self.binding_predictions
        )

    def _merge_activity_predictions(self, other: "PredictionsMixin") -> None:
        """Merge activity predictions."""
        self.activity_predictions.extend(
            pred for pred in other.activity_predictions
            if pred not in self.activity_predictions
        )

    def _merge_mechanism_predictions(self, other: "PredictionsMixin") -> None:
        """Merge mechanism predictions."""
        self.mechanism_predictions.extend(
            pred for pred in other.mechanism_predictions
            if pred not in self.mechanism_predictions
        )

    def get_cached_prediction(
        self,
        predictor_type: str,
    ) -> Optional[PredictionResult]:
        """Get cached prediction result if available."""
        return self._prediction_cache.get(predictor_type)

    def clear_prediction_cache(self) -> None:
        """Clear cached predictions."""
        self._prediction_cache.clear()

    def get_prediction_history(
        self,
        predictor_type: Optional[str] = None,
    ) -> pd.DataFrame:
        """Get prediction history, optionally filtered by type."""
        if predictor_type:
            return self._prediction_history[
                self._prediction_history['predictor_type'] == predictor_type
            ]
        return self._prediction_history

    def get_cached_features(
        self,
        feature_type: str,
    ) -> Optional[np.ndarray]:
        """Get cached features if available."""
        return self._feature_cache.get(feature_type)

    def cache_features(
        self,
        feature_type: str,
        features: np.ndarray,
    ) -> None:
        """Cache features for reuse."""
        self._feature_cache[feature_type] = features

    def clear_feature_cache(self) -> None:
        """Clear cached features."""
        self._feature_cache.clear()


class ValidationError(Exception):
    """Raised when validation fails."""
    pass
