"""Base predictor for psychopharmacological properties.

This module provides the base classes and interfaces for implementing predictors
that analyze compounds for psychopharmacological properties including:
1. Common ML model loading and management
2. Feature extraction and preprocessing 
3. Prediction confidence scoring
4. Result validation and formatting
5. Model retraining capabilities
6. Web data enrichment
7. Ensemble methods
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import numpy as np

from ....models.compound import PsychoactiveCompound
from ....models.predictions import PredictionResult
from ...structure.ml.features import EnhancedFeatureExtractor
from ...structure.ml.fingerprints import FingerprintGenerator
from ...structure.ml.descriptors import DescriptorGenerator
from ....utils.cache import CacheManager


@dataclass
class PredictorConfig:
    """Configuration for predictors."""
    # Model parameters
    model_path: Optional[str] = None
    model_version: str = "1.0.0"
    batch_size: int = 32
    
    # Feature parameters
    feature_type: str = "fingerprint"
    feature_params: Dict[str, Any] = None
    
    # Prediction parameters
    confidence_threshold: float = 0.7
    ensemble_size: int = 5
    
    # Web enrichment
    use_web_data: bool = True
    web_sources: List[str] = None
    
    def __post_init__(self):
        """Initialize default values."""
        if self.feature_params is None:
            self.feature_params = {}
        if self.web_sources is None:
            self.web_sources = ["chembl", "pubchem", "patents"]


class PredictorBase(ABC):
    """Base class for all psychopharmacological predictors."""

    def __init__(
        self,
        config: Optional[PredictorConfig] = None,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        feature_types: Optional[List[str]] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize predictor.
        
        Args:
            config: Optional predictor configuration
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            feature_types: Optional list of feature types to use
            log_level: Logging level
        """
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        
        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / f"{self.__class__.__name__}.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

        # Initialize configuration
        self.config = config or PredictorConfig()
        self.model_dir = Path(model_dir) if model_dir else None
        self.cache = CacheManager(cache_dir) if cache_dir else None
        
        # Initialize feature generators
        self.feature_types = feature_types or ["fingerprints", "descriptors"]
        self.feature_extractors = {
            "fingerprints": FingerprintGenerator(),
            "descriptors": DescriptorGenerator(),
            "enhanced": EnhancedFeatureExtractor(),
        }

        # Load models and initialize components
        self.models = self._load_models()
        self._initialize()

        self.logger.info(f"{self.__class__.__name__} initialized successfully")

    @abstractmethod
    def _load_models(self) -> Dict:
        """Load ML models from disk.
        
        Returns:
            Dict mapping model names to loaded model objects
        """
        pass

    @abstractmethod
    def _initialize(self) -> None:
        """Initialize predictor components."""
        pass

    def _extract_features(
        self,
        compound: PsychoactiveCompound,
        feature_type: Optional[str] = None,
    ) -> np.ndarray:
        """Extract features from compound.
        
        Args:
            compound: Compound to extract features from
            feature_type: Optional specific feature type to extract
            
        Returns:
            Feature vector or dict of feature vectors if no type specified
        """
        try:
            if feature_type:
                # Extract specific feature type
                return self._extract_feature_type(compound, feature_type)
            else:
                # Extract all feature types
                features = {}
                for ft in self.feature_types:
                    features[ft] = self._extract_feature_type(compound, ft)
                return features
                
        except Exception as e:
            self.logger.error(
                f"Feature extraction error for {compound.name}: {str(e)}"
            )
            return np.array([]) if feature_type else {}

    def _extract_feature_type(
        self,
        compound: PsychoactiveCompound,
        feature_type: str,
    ) -> np.ndarray:
        """Extract specific feature type from compound.
        
        Args:
            compound: Compound to extract features from
            feature_type: Type of features to extract
            
        Returns:
            Feature vector
        """
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

    @abstractmethod
    def _predict_raw(self, features: np.ndarray) -> Tuple[Any, float]:
        """Make raw prediction from features.
        
        Args:
            features: Feature vector
            
        Returns:
            Tuple of (prediction, confidence)
        """
        pass

    @abstractmethod
    def _process_prediction(
        self,
        prediction: Any,
        confidence: float,
        compound: PsychoactiveCompound
    ) -> PredictionResult:
        """Process raw prediction into result.
        
        Args:
            prediction: Raw prediction
            confidence: Prediction confidence
            compound: Input compound
            
        Returns:
            Processed prediction result
        """
        pass

    def predict(
        self,
        compound: PsychoactiveCompound,
        use_cache: bool = True
    ) -> PredictionResult:
        """Make prediction for compound.
        
        Args:
            compound: Compound to make prediction for
            use_cache: Whether to use cached predictions
            
        Returns:
            Prediction result
        """
        try:
            # Check cache
            if use_cache:
                cached = self._check_cache(compound)
                if cached is not None:
                    return cached

            # Extract features
            features = self._extract_features(compound)

            # Make prediction
            prediction, confidence = self._predict_raw(features)

            # Validate prediction
            is_valid, error = self._validate_prediction(
                prediction,
                confidence,
                self.config.confidence_threshold
            )
            if not is_valid:
                self.logger.warning(
                    f"Invalid prediction for {compound.name}: {error}"
                )
                return PredictionResult(
                    value=None,
                    confidence=0.0,
                    metadata={"error": error}
                )

            # Process result
            result = self._process_prediction(prediction, confidence, compound)

            # Cache result
            if use_cache:
                self._cache_prediction(compound, result)

            return result
            
        except Exception as e:
            self.logger.error(
                f"Prediction error for {compound.name}: {str(e)}"
            )
            return PredictionResult(
                value=None,
                confidence=0.0,
                metadata={"error": str(e)}
            )

    def predict_batch(
        self,
        compounds: List[PsychoactiveCompound],
        use_cache: bool = True
    ) -> List[PredictionResult]:
        """Make predictions for multiple compounds.
        
        Args:
            compounds: List of compounds
            use_cache: Whether to use cached predictions
            
        Returns:
            List of prediction results
        """
        results = []
        batch = []
        batch_compounds = []

        for compound in compounds:
            # Check cache
            if use_cache:
                cached = self._check_cache(compound)
                if cached is not None:
                    results.append(cached)
                    continue

            # Add to batch
            features = self._extract_features(compound)
            batch.append(features)
            batch_compounds.append(compound)

            # Process batch
            if len(batch) >= self.config.batch_size:
                batch_results = self._process_batch(batch, batch_compounds)
                results.extend(batch_results)
                batch = []
                batch_compounds = []

        # Process remaining compounds
        if batch:
            batch_results = self._process_batch(batch, batch_compounds)
            results.extend(batch_results)

        return results

    def _process_batch(
        self,
        batch: List[np.ndarray],
        compounds: List[PsychoactiveCompound]
    ) -> List[PredictionResult]:
        """Process a batch of compounds.
        
        Args:
            batch: List of feature vectors
            compounds: Corresponding compounds
            
        Returns:
            List of prediction results
        """
        # Stack features
        features = np.stack(batch)

        # Make predictions
        predictions = []
        confidences = []
        for i in range(len(features)):
            pred, conf = self._predict_raw(features[i])
            predictions.append(pred)
            confidences.append(conf)

        # Process results
        results = []
        for pred, conf, compound in zip(predictions, confidences, compounds):
            result = self._process_prediction(pred, conf, compound)
            self._cache_prediction(compound, result)
            results.append(result)

        return results

    def _check_cache(
        self,
        compound: PsychoactiveCompound
    ) -> Optional[PredictionResult]:
        """Check for cached prediction.
        
        Args:
            compound: Compound to check cache for
            
        Returns:
            Cached prediction if available, else None
        """
        return compound.get_cached_prediction(self.__class__.__name__)

    def _cache_prediction(
        self,
        compound: PsychoactiveCompound,
        result: PredictionResult
    ) -> None:
        """Cache prediction result.
        
        Args:
            compound: Compound to cache prediction for
            result: Prediction result to cache
        """
        compound.cache_prediction(self.__class__.__name__, result)

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
        """Save models to disk with versioning.
        
        Args:
            save_dir: Optional directory to save models to
        """
        if not self.models:
            return

        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            raise ValueError("No save directory specified")

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving models to {save_dir}")

        # Save with versioning
        timestamp = np.datetime64('now').astype(str)
        version_dir = save_dir / f"version_{timestamp}"
        version_dir.mkdir(exist_ok=True)

        try:
            for name, model in self.models.items():
                # Save current version
                model_path = save_dir / f"{name}.pkl"
                np.save(model_path, model)
                
                # Save versioned copy
                version_path = version_dir / f"{name}.pkl"
                np.save(version_path, model)
                
                self.logger.debug(f"Saved model {name}")

        except Exception as e:
            self.logger.error(f"Error saving models: {str(e)}")

    def retrain(
        self,
        compounds: List[PsychoactiveCompound],
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


class EnsemblePredictor(PredictorBase):
    """Base class for ensemble predictors."""

    def __init__(self, config: Optional[PredictorConfig] = None, **kwargs):
        """Initialize ensemble predictor.
        
        Args:
            config: Optional predictor configuration
            **kwargs: Additional arguments passed to PredictorBase
        """
        super().__init__(config, **kwargs)
        self._models = []
        self._initialize_ensemble()

    @abstractmethod
    def _initialize_ensemble(self) -> None:
        """Initialize ensemble models."""
        pass

    def _predict_raw(self, features: np.ndarray) -> Tuple[Any, float]:
        """Make ensemble prediction.
        
        Args:
            features: Feature vector
            
        Returns:
            Tuple of (prediction, confidence)
        """
        # Get predictions from each model
        predictions = []
        confidences = []
        for model in self._models:
            pred, conf = self._predict_with_model(model, features)
            predictions.append(pred)
            confidences.append(conf)

        # Aggregate predictions
        final_pred = self._aggregate_predictions(predictions)
        final_conf = self._aggregate_confidences(confidences)

        return final_pred, final_conf

    @abstractmethod
    def _predict_with_model(
        self,
        model: Any,
        features: np.ndarray
    ) -> Tuple[Any, float]:
        """Make prediction with single model.
        
        Args:
            model: Model to use
            features: Feature vector
            
        Returns:
            Tuple of (prediction, confidence)
        """
        pass

    @abstractmethod
    def _aggregate_predictions(self, predictions: List[Any]) -> Any:
        """Aggregate predictions from ensemble.
        
        Args:
            predictions: List of predictions
            
        Returns:
            Aggregated prediction
        """
        pass

    @abstractmethod
    def _aggregate_confidences(self, confidences: List[float]) -> float:
        """Aggregate confidences from ensemble.
        
        Args:
            confidences: List of confidence values
            
        Returns:
            Aggregated confidence
        """
        pass


class WebEnrichedPredictor(PredictorBase):
    """Base class for predictors that use web data."""

    def predict(
        self,
        compound: PsychoactiveCompound,
        use_cache: bool = True
    ) -> PredictionResult:
        """Make prediction using both ML and web data.
        
        Args:
            compound: Compound to make prediction for
            use_cache: Whether to use cached predictions
            
        Returns:
            Prediction result
        """
        # Get ML prediction
        ml_result = super().predict(compound, use_cache)

        # Skip web enrichment if disabled
        if not self.config.use_web_data:
            return ml_result

        # Get web data
        web_data = self._get_web_data(compound)
        if not web_data:
            return ml_result

        # Combine predictions
        return self._combine_predictions(ml_result, web_data, compound)

    @abstractmethod
    def _get_web_data(self, compound: PsychoactiveCompound) -> Optional[Dict]:
        """Get relevant web data for compound.
        
        Args:
            compound: Compound to get data for
            
        Returns:
            Dictionary of web data if available
        """
        pass

    @abstractmethod
    def _combine_predictions(
        self,
        ml_result: PredictionResult,
        web_data: Dict,
        compound: PsychoactiveCompound
    ) -> PredictionResult:
        """Combine ML and web-based predictions.
        
        Args:
            ml_result: ML prediction result
            web_data: Web data
            compound: Input compound
            
        Returns:
            Combined prediction result
        """
        pass
