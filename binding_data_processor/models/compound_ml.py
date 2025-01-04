"""Machine learning integration for compound data.

This module extends CompoundData with ML-specific functionality:
- Feature management and caching
- Prediction integration
- Model tracking and validation
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional
import json
import numpy as np
import pandas as pd

from .compound_base import CompoundData
from ..processors.psychopharm.base import PredictionResult
from ..processors.psychopharm.predictors.base import PredictorBase


@dataclass
class MLCompoundData(CompoundData):
    """CompoundData with ML capabilities."""
    
    # Feature Management
    _feature_cache: Dict[str, np.ndarray] = field(default_factory=dict)
    _feature_importances: Dict[str, Dict[str, float]] = field(default_factory=dict)
    _feature_scalers: Dict[str, Any] = field(default_factory=dict)

    # Prediction Integration
    _prediction_cache: Dict[str, PredictionResult] = field(default_factory=dict)
    _prediction_history: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(
        columns=[
            'predictor_type',
            'prediction_value',
            'confidence',
            'timestamp',
            'supporting_data',
        ]
    ))

    # Model Integration
    _models: Dict[str, Any] = field(default_factory=dict)
    _model_metrics: Dict[str, Dict[str, float]] = field(default_factory=dict)
    _model_versions: Dict[str, str] = field(default_factory=dict)

    def predict(
        self,
        predictor_type: str,
        predictor: PredictorBase,
        use_cache: bool = True,
    ) -> PredictionResult:
        """Run prediction using specified predictor.
        
        Args:
            predictor_type: Type of predictor to use
            predictor: Predictor instance to use
            use_cache: Whether to use cached predictions
            
        Returns:
            PredictionResult containing prediction and confidence
        """
        # Check cache
        if use_cache and predictor_type in self._prediction_cache:
            return self._prediction_cache[predictor_type]

        # Run prediction
        result = predictor.predict(self)

        # Cache result
        self._prediction_cache[predictor_type] = result

        # Update history
        self._prediction_history = pd.concat([
            self._prediction_history,
            pd.DataFrame([{
                'predictor_type': predictor_type,
                'prediction_value': result.value,
                'confidence': result.confidence,
                'timestamp': pd.Timestamp.now(),
                'supporting_data': json.dumps(result.supporting_data),
            }])
        ], ignore_index=True)

        return result

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

    def get_feature_importances(
        self,
        predictor_type: str,
    ) -> Dict[str, float]:
        """Get feature importances for a predictor."""
        return self._feature_importances.get(predictor_type, {})

    def set_feature_importances(
        self,
        predictor_type: str,
        importances: Dict[str, float],
    ) -> None:
        """Set feature importances for a predictor."""
        self._feature_importances[predictor_type] = importances

    def get_model_metrics(
        self,
        predictor_type: str,
    ) -> Dict[str, float]:
        """Get model metrics for a predictor."""
        return self._model_metrics.get(predictor_type, {})

    def set_model_metrics(
        self,
        predictor_type: str,
        metrics: Dict[str, float],
    ) -> None:
        """Set model metrics for a predictor."""
        self._model_metrics[predictor_type] = metrics

    def get_model_version(
        self,
        predictor_type: str,
    ) -> Optional[str]:
        """Get model version for a predictor."""
        return self._model_versions.get(predictor_type)

    def set_model_version(
        self,
        predictor_type: str,
        version: str,
    ) -> None:
        """Set model version for a predictor."""
        self._model_versions[predictor_type] = version

    def to_dict(self, include_predictions: bool = True) -> Dict:
        """Convert compound data to dictionary format."""
        data = super().to_dict()

        # Add ML predictions if requested
        if include_predictions:
            for pred_type, result in self._prediction_cache.items():
                data[f"{pred_type}_prediction"] = {
                    "value": result.value,
                    "confidence": result.confidence,
                    "supporting_data": result.supporting_data,
                    "feature_importances": self.get_feature_importances(pred_type),
                    "model_metrics": self.get_model_metrics(pred_type),
                    "model_version": self.get_model_version(pred_type),
                }

        return data
