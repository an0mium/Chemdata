"""Machine learning mixin for compound data.

This module provides the MLMixin class that adds ML capabilities:
- Feature management and caching
- Prediction integration and history
- Model tracking and validation
- Uncertainty estimation
- Feature importance tracking
- Model metrics and versioning

The mixin is designed to be used with CompoundData to add ML functionality
while maintaining clean separation of concerns.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, TYPE_CHECKING
import json
import numpy as np
import pandas as pd

from .types import (
    Features,
    PredictionResult,
    ModelMetrics,
    FeatureImportances,
    ModelVersion,
    PredictionDict,
)

if TYPE_CHECKING:
    from ..processors.psychopharm.predictors.base import PredictorBase


@dataclass
class MLMixin:
    """Mixin class adding ML capabilities to CompoundData."""
    
    # Feature Management
    _feature_cache: Dict[str, Features] = field(default_factory=dict)
    _feature_importances: Dict[str, FeatureImportances] = field(default_factory=dict)
    _feature_scalers: Dict[str, object] = field(default_factory=dict)

    # Prediction Integration
    _prediction_cache: PredictionDict = field(default_factory=dict)
    _prediction_history: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(
        columns=[
            'predictor_type',
            'prediction_value',
            'confidence',
            'uncertainty',  # Added uncertainty estimation
            'timestamp',
            'supporting_data',
        ]
    ))

    # Model Integration
    _model_metrics: Dict[str, ModelMetrics] = field(default_factory=dict)
    _model_versions: Dict[str, ModelVersion] = field(default_factory=dict)
    _model_params: Dict[str, Dict] = field(default_factory=dict)  # Store hyperparameters

    def predict(
        self,
        predictor_type: str,
        predictor: 'PredictorBase',
        use_cache: bool = True,
        include_uncertainty: bool = True,
    ) -> PredictionResult:
        """Run prediction using specified predictor.
        
        Args:
            predictor_type: Type of predictor to use
            predictor: Predictor instance to use
            use_cache: Whether to use cached predictions
            include_uncertainty: Whether to estimate prediction uncertainty
            
        Returns:
            PredictionResult containing prediction, confidence, and uncertainty
        """
        # Check cache
        if use_cache and predictor_type in self._prediction_cache:
            return self._prediction_cache[predictor_type]

        # Run prediction
        result = predictor.predict(self, include_uncertainty=include_uncertainty)

        # Cache result
        self._prediction_cache[predictor_type] = result

        # Update history
        self._prediction_history = pd.concat([
            self._prediction_history,
            pd.DataFrame([{
                'predictor_type': predictor_type,
                'prediction_value': result.value,
                'confidence': result.confidence,
                'uncertainty': result.uncertainty if include_uncertainty else None,
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
        include_uncertainty: bool = True,
    ) -> pd.DataFrame:
        """Get prediction history, optionally filtered by type.
        
        Args:
            predictor_type: Optional type to filter by
            include_uncertainty: Whether to include uncertainty column
            
        Returns:
            DataFrame of prediction history
        """
        history = self._prediction_history
        if predictor_type:
            history = history[history['predictor_type'] == predictor_type]
        if not include_uncertainty:
            history = history.drop('uncertainty', axis=1)
        return history

    def get_cached_features(
        self,
        feature_type: str,
    ) -> Optional[Features]:
        """Get cached features if available."""
        return self._feature_cache.get(feature_type)

    def cache_features(
        self,
        feature_type: str,
        features: Features,
        scaler: Optional[object] = None,
    ) -> None:
        """Cache features and optional scaler for reuse."""
        self._feature_cache[feature_type] = features
        if scaler is not None:
            self._feature_scalers[feature_type] = scaler

    def clear_feature_cache(self) -> None:
        """Clear cached features and scalers."""
        self._feature_cache.clear()
        self._feature_scalers.clear()

    def get_feature_importances(
        self,
        predictor_type: str,
    ) -> FeatureImportances:
        """Get feature importances for a predictor."""
        return self._feature_importances.get(predictor_type, {})

    def set_feature_importances(
        self,
        predictor_type: str,
        importances: FeatureImportances,
    ) -> None:
        """Set feature importances for a predictor."""
        self._feature_importances[predictor_type] = importances

    def get_model_metrics(
        self,
        predictor_type: str,
    ) -> ModelMetrics:
        """Get model metrics for a predictor."""
        return self._model_metrics.get(predictor_type, {})

    def set_model_metrics(
        self,
        predictor_type: str,
        metrics: ModelMetrics,
    ) -> None:
        """Set model metrics for a predictor."""
        self._model_metrics[predictor_type] = metrics

    def get_model_version(
        self,
        predictor_type: str,
    ) -> Optional[ModelVersion]:
        """Get model version for a predictor."""
        return self._model_versions.get(predictor_type)

    def set_model_version(
        self,
        predictor_type: str,
        version: ModelVersion,
    ) -> None:
        """Set model version for a predictor."""
        self._model_versions[predictor_type] = version

    def get_model_params(
        self,
        predictor_type: str,
    ) -> Dict:
        """Get model hyperparameters for a predictor."""
        return self._model_params.get(predictor_type, {})

    def set_model_params(
        self,
        predictor_type: str,
        params: Dict,
    ) -> None:
        """Set model hyperparameters for a predictor."""
        self._model_params[predictor_type] = params

    def to_dict(self, include_predictions: bool = True) -> Dict:
        """Convert ML data to dictionary format.
        
        Args:
            include_predictions: Whether to include prediction data
            
        Returns:
            Dictionary containing ML-related data
        """
        data = {}

        # Add ML predictions if requested
        if include_predictions:
            for pred_type, result in self._prediction_cache.items():
                data[f"{pred_type}_prediction"] = {
                    "value": result.value,
                    "confidence": result.confidence,
                    "uncertainty": result.uncertainty if hasattr(result, 'uncertainty') else None,
                    "supporting_data": result.supporting_data,
                    "feature_importances": self.get_feature_importances(pred_type),
                    "model_metrics": self.get_model_metrics(pred_type),
                    "model_version": self.get_model_version(pred_type),
                    "model_params": self.get_model_params(pred_type),
                }

        return data
