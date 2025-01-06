"""Base predictor for ML models.

This module provides the base predictor class that all ML models inherit from:
- Feature extraction
- Model management
- Uncertainty estimation
- Validation
- Metadata tracking
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from rdkit import Chem

from ....models.core import CompoundData
from ....models.predictions import PredictionResult
from ..features.fingerprints import FingerprintGenerator
from .uncertainty import UncertaintyEstimator
from ..training.validation import ModelValidator


class BasePredictor:
    """Base class for all predictors."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        feature_types: Optional[List[str]] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize predictor.

        Args:
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            feature_types: Optional list of feature types to use
            log_level: Logging level
        """
        self.model_dir = Path(model_dir) if model_dir else None
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.feature_types = feature_types or ["morgan", "maccs", "rdkit"]

        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)

        # Initialize components
        self.fingerprints = FingerprintGenerator()
        self.uncertainty = UncertaintyEstimator()
        self.validator = ModelValidator()

        # Initialize model
        self.model = None
        self._load_models()

    def _load_models(self) -> Dict:
        """Load trained models.

        Returns:
            Dictionary of loaded models
        """
        if not self.model_dir or not self.model_dir.exists():
            self.logger.info("No models found, initializing new models")
            return self._initialize_models()

        try:
            self.logger.info(f"Loading models from {self.model_dir}")
            models = {}
            for model_path in self.model_dir.glob("*.pkl"):
                models[model_path.stem] = self._load_model(model_path)
            return models

        except Exception as e:
            self.logger.error(f"Error loading models: {str(e)}")
            self.logger.info("Initializing new models")
            return self._initialize_models()

    def _initialize_models(self) -> Dict:
        """Initialize new models.

        Returns:
            Dictionary of initialized models
        """
        raise NotImplementedError("Subclasses must implement _initialize_models")

    def _load_model(self, path: Path) -> Any:
        """Load model from file.

        Args:
            path: Path to model file

        Returns:
            Loaded model
        """
        raise NotImplementedError("Subclasses must implement _load_model")

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save models to directory.

        Args:
            save_dir: Optional directory to save models to
        """
        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            raise ValueError("No save directory specified")

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving models to {save_dir}")

        try:
            for name, model in self.models.items():
                path = save_dir / f"{name}.pkl"
                self._save_model(model, path)

        except Exception as e:
            self.logger.error(f"Error saving models: {str(e)}")

    def _save_model(self, model: Any, path: Path) -> None:
        """Save model to file.

        Args:
            model: Model to save
            path: Path to save model to
        """
        raise NotImplementedError("Subclasses must implement _save_model")

    def _extract_features(
        self,
        compound: CompoundData,
        feature_types: Optional[List[str]] = None,
    ) -> np.ndarray:
        """Extract features from compound.

        Args:
            compound: Compound to extract features from
            feature_types: Optional list of feature types to use

        Returns:
            Feature array
        """
        feature_types = feature_types or self.feature_types
        features = []

        # Convert SMILES to RDKit mol
        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {compound.smiles}")

        # Generate fingerprints
        for fp_type in feature_types:
            fp = self.fingerprints.generate(mol, fp_type=fp_type)
            features.append(fp)

        return np.concatenate(features)

    def predict(
        self,
        compound: CompoundData,
        **kwargs,
    ) -> PredictionResult:
        """Make prediction with uncertainty.

        Args:
            compound: Compound to make prediction for
            **kwargs: Additional keyword arguments

        Returns:
            Prediction result with uncertainty
        """
        try:
            # Extract features
            features = self._extract_features(compound)

            # Make prediction
            prediction = self._predict(features, **kwargs)

            # Estimate uncertainty
            uncertainty = self.uncertainty.estimate(
                model=self.model,
                features=features,
                prediction=prediction,
            )

            # Get metadata
            metadata = self._get_metadata(compound, features, prediction)

            return PredictionResult(
                value=prediction,
                confidence=1.0 - uncertainty,
                supporting_data=metadata,
            )

        except Exception as e:
            self.logger.error(f"Error making prediction: {str(e)}")
            return PredictionResult(
                value=None,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _predict(
        self,
        features: np.ndarray,
        **kwargs,
    ) -> Any:
        """Make raw prediction.

        Args:
            features: Feature array
            **kwargs: Additional keyword arguments

        Returns:
            Raw prediction
        """
        raise NotImplementedError("Subclasses must implement _predict")

    def _get_metadata(
        self,
        compound: CompoundData,
        features: np.ndarray,
        prediction: Any,
    ) -> Dict[str, Any]:
        """Get prediction metadata.

        Args:
            compound: Input compound
            features: Extracted features
            prediction: Raw prediction

        Returns:
            Prediction metadata
        """
        return {
            "compound": {
                "name": compound.name,
                "smiles": compound.smiles,
                "cas_number": compound.cas_number,
            },
            "features": {
                "types": self.feature_types,
                "shape": features.shape,
            },
            "model": {
                "type": self.__class__.__name__,
                "version": self._get_version(),
            },
        }

    def _get_version(self) -> str:
        """Get model version.

        Returns:
            Model version string
        """
        return "0.1.0"  # Default version

    def validate(
        self,
        compounds: List[CompoundData],
        labels: List[Any],
        metrics: Optional[List[str]] = None,
    ) -> Dict[str, float]:
        """Validate model on test data.

        Args:
            compounds: Test compounds
            labels: True labels
            metrics: Optional list of metrics to calculate

        Returns:
            Dictionary of validation metrics
        """
        # Extract features
        features = []
        for compound in compounds:
            try:
                feature = self._extract_features(compound)
                features.append(feature)
            except Exception as e:
                self.logger.warning(f"Error extracting features for {compound.name}: {str(e)}")

        if not features:
            raise ValueError("No valid features extracted")

        # Convert to arrays
        X = np.array(features)
        y = np.array(labels)

        # Validate
        return self.validator.validate(
            model=self.model,
            X=X,
            y=y,
            metrics=metrics,
        )

    def train(
        self,
        compounds: List[CompoundData],
        labels: List[Any],
        **kwargs,
    ) -> Dict[str, float]:
        """Train model on data.

        Args:
            compounds: Training compounds
            labels: Training labels
            **kwargs: Additional training parameters

        Returns:
            Dictionary of training metrics
        """
        raise NotImplementedError("Subclasses must implement train")
