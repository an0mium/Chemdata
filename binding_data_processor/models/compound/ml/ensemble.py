"""Ensemble model management for compound predictions.

This module provides:
1. Ensemble model base classes
2. Model combination strategies
3. Uncertainty estimation
4. Cross-validation support
5. Model selection
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin
from sklearn.ensemble import (
    RandomForestClassifier,
    RandomForestRegressor,
    GradientBoostingClassifier,
    GradientBoostingRegressor,
    VotingClassifier,
    VotingRegressor,
    StackingClassifier,
    StackingRegressor,
)
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler

from ....models.core import CompoundData
from ....processors.structure.ml.features import EnhancedFeatureExtractor


class EnsembleBase:
    """Base class for ensemble models."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        feature_types: Optional[List[str]] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize ensemble base.

        Args:
            model_dir: Optional directory for model storage
            cache_dir: Optional directory for feature cache
            feature_types: Optional list of feature types to use
            log_level: Logging level
        """
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)

        # Initialize paths
        self.model_dir = Path(model_dir) if model_dir else None
        self.cache_dir = Path(cache_dir) if cache_dir else None

        # Initialize feature extraction
        self.feature_types = feature_types or ["fingerprints", "descriptors"]
        self.feature_extractor = EnhancedFeatureExtractor(
            feature_types=self.feature_types,
        )

        # Initialize scalers
        self.scalers: Dict[str, StandardScaler] = {}

        # Initialize models
        self.models: Dict[str, Union[BaseEstimator, List[BaseEstimator]]] = {}

        # Initialize metrics
        self.metrics: Dict[str, Dict[str, float]] = {}

    def _extract_features(
        self,
        compound: CompoundData,
        feature_type: str,
    ) -> np.ndarray:
        """Extract features for a compound.

        Args:
            compound: Compound to extract features for
            feature_type: Type of features to extract

        Returns:
            Feature array
        """
        return self.feature_extractor.extract_features([compound.mol], [feature_type])[feature_type]

    def _scale_features(
        self,
        features: Dict[str, np.ndarray],
        model_name: str,
    ) -> np.ndarray:
        """Scale features for a model.

        Args:
            features: Dictionary of feature arrays
            model_name: Name of model to scale for

        Returns:
            Scaled feature array
        """
        # Combine features
        X = np.hstack([features[ft] for ft in self.feature_types])

        # Scale if scaler exists
        if model_name in self.scalers:
            X = self.scalers[model_name].transform(X)

        return X

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save models to disk."""
        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            return

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving models to {save_dir}")

        # Save with versioning
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
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

                self.logger.debug(f"Saved model: {name}")

            # Save scalers
            for name, scaler in self.scalers.items():
                # Save current version
                scaler_path = save_dir / f"{name}_scaler.pkl"
                np.save(scaler_path, scaler)

                # Save versioned copy
                version_path = version_dir / f"{name}_scaler.pkl"
                np.save(version_path, scaler)

                self.logger.debug(f"Saved scaler: {name}")

            # Save metrics
            metrics_path = save_dir / "metrics.json"
            pd.DataFrame(self.metrics).to_json(metrics_path)

        except Exception as e:
            self.logger.error(f"Error saving models: {str(e)}")

    def load_models(self, load_dir: Optional[str] = None) -> None:
        """Load models from disk."""
        load_dir = Path(load_dir) if load_dir else self.model_dir
        if not load_dir or not load_dir.exists():
            return

        self.logger.info(f"Loading models from {load_dir}")

        try:
            # Load models
            for model_path in load_dir.glob("*.pkl"):
                if model_path.stem.endswith("_scaler"):
                    continue
                name = model_path.stem
                self.models[name] = np.load(model_path, allow_pickle=True)
                self.logger.debug(f"Loaded model: {name}")

            # Load scalers
            for scaler_path in load_dir.glob("*_scaler.pkl"):
                name = scaler_path.stem.replace("_scaler", "")
                self.scalers[name] = np.load(scaler_path, allow_pickle=True)
                self.logger.debug(f"Loaded scaler: {name}")

            # Load metrics
            metrics_path = load_dir / "metrics.json"
            if metrics_path.exists():
                self.metrics = pd.read_json(metrics_path).to_dict()

        except Exception as e:
            self.logger.error(f"Error loading models: {str(e)}")


class ClassifierEnsemble(EnsembleBase):
    """Ensemble of classification models."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        feature_types: Optional[List[str]] = None,
        log_level: int = logging.INFO,
        voting: str = "soft",
        n_jobs: int = -1,
    ):
        """Initialize classifier ensemble.

        Args:
            model_dir: Optional directory for model storage
            cache_dir: Optional directory for feature cache
            feature_types: Optional list of feature types to use
            log_level: Logging level
            voting: Voting strategy ('hard' or 'soft')
            n_jobs: Number of jobs for parallel processing
        """
        super().__init__(model_dir, cache_dir, feature_types, log_level)
        self.voting = voting
        self.n_jobs = n_jobs

    def create_ensemble(
        self,
        name: str,
        estimators: Optional[List[Tuple[str, BaseEstimator]]] = None,
        ensemble_type: str = "voting",
        **kwargs: Any,
    ) -> None:
        """Create a new ensemble model.

        Args:
            name: Name for the ensemble
            estimators: Optional list of (name, estimator) tuples
            ensemble_type: Type of ensemble ('voting' or 'stacking')
            **kwargs: Additional arguments for ensemble
        """
        if estimators is None:
            # Default estimators
            estimators = [
                (
                    "rf",
                    RandomForestClassifier(
                        n_estimators=100,
                        max_depth=10,
                        n_jobs=self.n_jobs,
                    ),
                ),
                (
                    "gb",
                    GradientBoostingClassifier(
                        n_estimators=100,
                        max_depth=5,
                    ),
                ),
            ]

        if ensemble_type == "voting":
            model = VotingClassifier(
                estimators=estimators,
                voting=self.voting,
                n_jobs=self.n_jobs,
                **kwargs,
            )
        elif ensemble_type == "stacking":
            model = StackingClassifier(
                estimators=estimators,
                n_jobs=self.n_jobs,
                **kwargs,
            )
        else:
            raise ValueError(f"Unknown ensemble type: {ensemble_type}")

        self.models[name] = model
        self.logger.info(f"Created {ensemble_type} ensemble: {name}")

    def fit(
        self,
        name: str,
        X: np.ndarray,
        y: np.ndarray,
        **kwargs: Any,
    ) -> Dict[str, float]:
        """Fit an ensemble model.

        Args:
            name: Name of ensemble to fit
            X: Feature matrix
            y: Target values
            **kwargs: Additional fit arguments

        Returns:
            Dictionary of metrics
        """
        if name not in self.models:
            raise ValueError(f"Unknown model: {name}")

        # Create and fit scaler
        self.scalers[name] = StandardScaler()
        X_scaled = self.scalers[name].fit_transform(X)

        # Fit model
        self.models[name].fit(X_scaled, y, **kwargs)

        # Calculate metrics
        cv = KFold(n_splits=5, shuffle=True, random_state=42)
        scores = cross_val_score(
            self.models[name],
            X_scaled,
            y,
            cv=cv,
            n_jobs=self.n_jobs,
        )

        metrics = {
            "accuracy_mean": float(scores.mean()),
            "accuracy_std": float(scores.std()),
        }
        self.metrics[name] = metrics

        self.logger.info(f"Fitted {name}: accuracy={metrics['accuracy_mean']:.3f} " f"(±{metrics['accuracy_std']:.3f})")

        return metrics

    def predict(
        self,
        name: str,
        compound: CompoundData,
    ) -> Tuple[Any, float]:
        """Generate predictions for a compound.

        Args:
            name: Name of ensemble to use
            compound: Compound to predict for

        Returns:
            Tuple of (prediction, confidence)
        """
        if name not in self.models:
            raise ValueError(f"Unknown model: {name}")

        # Extract features
        features = {}
        for feature_type in self.feature_types:
            features[feature_type] = self._extract_features(
                compound,
                feature_type,
            )

        # Scale features
        X = self._scale_features(features, name)

        # Get prediction and confidence
        model = self.models[name]
        if self.voting == "soft":
            probs = model.predict_proba(X)[0]
            pred_idx = np.argmax(probs)
            prediction = model.classes_[pred_idx]
            confidence = float(probs[pred_idx])
        else:
            prediction = model.predict(X)[0]
            # Use proportion of agreeing estimators as confidence
            predictions = np.array([est.predict(X)[0] for name, est in model.estimators_])
            confidence = float(np.mean(predictions == prediction))

        return prediction, confidence


class RegressorEnsemble(EnsembleBase):
    """Ensemble of regression models."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        feature_types: Optional[List[str]] = None,
        log_level: int = logging.INFO,
        n_jobs: int = -1,
    ):
        """Initialize regressor ensemble.

        Args:
            model_dir: Optional directory for model storage
            cache_dir: Optional directory for feature cache
            feature_types: Optional list of feature types to use
            log_level: Logging level
            n_jobs: Number of jobs for parallel processing
        """
        super().__init__(model_dir, cache_dir, feature_types, log_level)
        self.n_jobs = n_jobs

    def create_ensemble(
        self,
        name: str,
        estimators: Optional[List[Tuple[str, BaseEstimator]]] = None,
        ensemble_type: str = "voting",
        **kwargs: Any,
    ) -> None:
        """Create a new ensemble model.

        Args:
            name: Name for the ensemble
            estimators: Optional list of (name, estimator) tuples
            ensemble_type: Type of ensemble ('voting' or 'stacking')
            **kwargs: Additional arguments for ensemble
        """
        if estimators is None:
            # Default estimators
            estimators = [
                (
                    "rf",
                    RandomForestRegressor(
                        n_estimators=100,
                        max_depth=10,
                        n_jobs=self.n_jobs,
                    ),
                ),
                (
                    "gb",
                    GradientBoostingRegressor(
                        n_estimators=100,
                        max_depth=5,
                    ),
                ),
            ]

        if ensemble_type == "voting":
            model = VotingRegressor(
                estimators=estimators,
                n_jobs=self.n_jobs,
                **kwargs,
            )
        elif ensemble_type == "stacking":
            model = StackingRegressor(
                estimators=estimators,
                n_jobs=self.n_jobs,
                **kwargs,
            )
        else:
            raise ValueError(f"Unknown ensemble type: {ensemble_type}")

        self.models[name] = model
        self.logger.info(f"Created {ensemble_type} ensemble: {name}")

    def fit(
        self,
        name: str,
        X: np.ndarray,
        y: np.ndarray,
        **kwargs: Any,
    ) -> Dict[str, float]:
        """Fit an ensemble model.

        Args:
            name: Name of ensemble to fit
            X: Feature matrix
            y: Target values
            **kwargs: Additional fit arguments

        Returns:
            Dictionary of metrics
        """
        if name not in self.models:
            raise ValueError(f"Unknown model: {name}")

        # Create and fit scaler
        self.scalers[name] = StandardScaler()
        X_scaled = self.scalers[name].fit_transform(X)

        # Fit model
        self.models[name].fit(X_scaled, y, **kwargs)

        # Calculate metrics
        cv = KFold(n_splits=5, shuffle=True, random_state=42)
        scores = cross_val_score(
            self.models[name],
            X_scaled,
            y,
            cv=cv,
            scoring="r2",
            n_jobs=self.n_jobs,
        )

        metrics = {
            "r2_mean": float(scores.mean()),
            "r2_std": float(scores.std()),
        }
        self.metrics[name] = metrics

        self.logger.info(f"Fitted {name}: R²={metrics['r2_mean']:.3f} " f"(±{metrics['r2_std']:.3f})")

        return metrics

    def predict(
        self,
        name: str,
        compound: CompoundData,
    ) -> Tuple[float, float]:
        """Generate predictions for a compound.

        Args:
            name: Name of ensemble to use
            compound: Compound to predict for

        Returns:
            Tuple of (prediction, confidence)
        """
        if name not in self.models:
            raise ValueError(f"Unknown model: {name}")

        # Extract features
        features = {}
        for feature_type in self.feature_types:
            features[feature_type] = self._extract_features(
                compound,
                feature_type,
            )

        # Scale features
        X = self._scale_features(features, name)

        # Get predictions from all estimators
        predictions = np.array([est.predict(X)[0] for name, est in self.models[name].estimators_])

        # Calculate prediction and confidence
        prediction = float(np.mean(predictions))
        confidence = 1.0 - float(np.std(predictions))

        return prediction, confidence
