"""Base functionality for BBB permeability prediction.

This module provides the BBBPredictorBase class that implements core functionality for:
1. Feature extraction and scaling
2. Model management and versioning
3. Prediction history tracking
4. Basic statistics
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler

from ....models.core import CompoundData
from ....models.psychopharm import BBBPermeability
from .base import PredictorBase


class BBBPredictorBase(PredictorBase):
    """Base class for BBB permeability prediction."""

    # Default model paths
    DEFAULT_MODEL_DIR = Path("models/bbb")
    MODEL_FILENAMES = {
        "rf_classifier": "bbb_rf_classifier.pkl",
        "gb_classifier": "bbb_gb_classifier.pkl",
        "rf_regressor": "bbb_rf_regressor.pkl",
        "gb_regressor": "bbb_gb_regressor.pkl",
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize BBB predictor base."""
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        
        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / "bbb_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            feature_types=["fingerprints", "descriptors", "enhanced"],
        )

        # Initialize feature scalers
        self.scalers = self._initialize_scalers()

        # Configure model weights
        self.model_weights = {
            "rf_classifier": 0.4,
            "gb_classifier": 0.3,
            "rf_regressor": 0.3,
        }

        # Initialize prediction history
        self.prediction_history = pd.DataFrame(
            columns=[
                'compound_name',
                'permeability_class',
                'confidence',
                'permeability_score',
                'timestamp',
            ]
        )

        self.logger.info("BBBPredictorBase initialized successfully")

    def _initialize_scalers(self) -> Dict[str, StandardScaler]:
        """Initialize feature scalers for each feature type."""
        self.logger.debug("Initializing feature scalers")
        return {
            feature_type: StandardScaler()
            for feature_type in self.feature_types
        }

    def _load_models(self) -> Dict:
        """Load BBB prediction models from disk."""
        models = {}
        model_dir = Path(self.model_dir) if self.model_dir else self.DEFAULT_MODEL_DIR

        if model_dir.exists():
            self.logger.info(f"Loading models from {model_dir}")
            try:
                for model_name, filename in self.MODEL_FILENAMES.items():
                    model_path = model_dir / filename
                    if model_path.exists():
                        self.logger.debug(f"Loading {model_name} from {model_path}")
                        models[model_name] = np.load(model_path, allow_pickle=True)
                    else:
                        self.logger.warning(
                            f"Model file not found: {model_path}, initializing new model"
                        )
                        models[model_name] = self._initialize_model(model_name)

            except Exception as e:
                self.logger.error(f"Error loading models: {str(e)}")
                self.logger.info("Initializing new models")
                models = self._initialize_models()
        else:
            self.logger.info(f"Model directory not found: {model_dir}")
            self.logger.info("Initializing new models")
            models = self._initialize_models()

        return models

    def _initialize_model(self, model_name: str):
        """Initialize a single model with appropriate configuration."""
        self.logger.debug(f"Initializing new {model_name} model")
        
        if model_name in ["rf_classifier"]:
            return RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                random_state=42,
                n_jobs=-1,
                verbose=1,
            )
        elif model_name in ["gb_classifier", "gb_regressor"]:
            return GradientBoostingRegressor(
                n_estimators=100,
                max_depth=5,
                random_state=42,
                verbose=1,
            )
        elif model_name == "rf_regressor":
            return RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                random_state=42,
                n_jobs=-1,
                verbose=1,
            )
        else:
            raise ValueError(f"Unknown model type: {model_name}")

    def _initialize_models(self) -> Dict:
        """Initialize new BBB prediction models."""
        self.logger.info("Initializing new models")
        return {
            name: self._initialize_model(name)
            for name in self.MODEL_FILENAMES
        }

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save models to disk with versioning."""
        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            save_dir = self.DEFAULT_MODEL_DIR

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving models to {save_dir}")

        # Save models with versioning
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        version_dir = save_dir / f"version_{timestamp}"
        version_dir.mkdir(exist_ok=True)

        for model_name, model in self.models.items():
            try:
                # Save current version
                model_path = save_dir / self.MODEL_FILENAMES[model_name]
                np.save(model_path, model)
                
                # Save versioned copy
                version_path = version_dir / self.MODEL_FILENAMES[model_name]
                np.save(version_path, model)
                
                self.logger.debug(f"Saved {model_name} to {model_path}")
            except Exception as e:
                self.logger.error(f"Error saving {model_name}: {str(e)}")

        # Save scalers
        try:
            scaler_path = save_dir / "scalers.pkl"
            np.save(scaler_path, self.scalers)
            self.logger.debug(f"Saved scalers to {scaler_path}")
        except Exception as e:
            self.logger.error(f"Error saving scalers: {str(e)}")

        # Save prediction history
        try:
            history_path = save_dir / "prediction_history.csv"
            self.prediction_history.to_csv(history_path, index=False)
            self.logger.debug(f"Saved prediction history to {history_path}")
        except Exception as e:
            self.logger.error(f"Error saving prediction history: {str(e)}")

    def get_prediction_statistics(self) -> pd.DataFrame:
        """Get statistics about predictions made so far."""
        stats = pd.DataFrame()
        
        # Class distribution
        stats['class_dist'] = (
            self.prediction_history['permeability_class'].value_counts()
        )
        
        # Average confidence by class
        stats['avg_confidence'] = (
            self.prediction_history.groupby('permeability_class')['confidence'].mean()
        )
        
        # Score statistics
        stats['score_mean'] = self.prediction_history['permeability_score'].mean()
        stats['score_std'] = self.prediction_history['permeability_score'].std()
        
        return stats

    def _predict_permeability_class(
        self,
        model: RandomForestClassifier,
        features: Dict[str, np.ndarray],
        model_name: str,
    ) -> Tuple[str, float]:
        """Predict BBB permeability class."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[model_name].transform(X)

        # Get class probabilities
        probs = model.predict_proba(X_scaled)[0]
        pred_idx = np.argmax(probs)
        
        # Map to BBBPermeability class
        permeability_class = BBBPermeability(
            model.classes_[pred_idx]
        ).value
        
        return permeability_class, float(probs[pred_idx])

    def _predict_permeability_score(
        self,
        model: GradientBoostingRegressor,
        features: Dict[str, np.ndarray],
        model_name: str,
    ) -> Tuple[float, float]:
        """Predict BBB permeability score."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[model_name].transform(X)

        # Get prediction and confidence
        score = model.predict(X_scaled)[0]
        confidence = 1.0 - np.std(
            [est.predict(X_scaled)[0] for est in model.estimators_]
        )

        return float(score), float(confidence)

    def _extract_compound_features(
        self, compound: CompoundData
    ) -> Dict[str, np.ndarray]:
        """Extract features from compound."""
        features = {}
        for feature_type in self.feature_types:
            features[feature_type] = self._extract_features(compound, feature_type)
            self.logger.debug(
                f"Extracted {feature_type} features: "
                f"shape={features[feature_type].shape}"
            )
        return features

    def _extract_training_features(
        self, compounds: List[CompoundData]
    ) -> np.ndarray:
        """Extract features for training."""
        X = []
        for compound in compounds:
            features = []
            for feature_type in self.feature_types:
                feat = self._extract_features(compound, feature_type)
                features.append(feat)
            X.append(np.hstack(features))
        return np.vstack(X)

    def _train_classification_models(
        self, X: np.ndarray, labels: List[BBBPermeability]
    ) -> Dict[str, float]:
        """Train classification models."""
        metrics = {}
        self.logger.debug("Training classification models")
        
        for name, model in self.models.items():
            if name in ["rf_classifier", "gb_classifier"]:
                model.fit(X, labels)
                score = model.score(X, labels)
                metrics[f"{name}_score"] = score
                self.logger.debug(f"{name} training score: {score:.3f}")
        
        return metrics

    def _train_regression_models(
        self, X: np.ndarray, scores: List[float]
    ) -> Dict[str, float]:
        """Train regression models."""
        metrics = {}
        self.logger.debug("Training regression models")
        
        for name, model in self.models.items():
            if name in ["rf_regressor", "gb_regressor"]:
                model.fit(X, scores)
                score = model.score(X, scores)
                metrics[f"{name}_score"] = score
                self.logger.debug(f"{name} training score: {score:.3f}")
        
        return metrics

    def export_predictions(
        self,
        output_path: str,
        columns: Optional[List[str]] = None,
    ) -> None:
        """Export prediction history to TSV file."""
        if columns is None:
            columns = self.prediction_history.columns
        
        self.prediction_history[columns].to_csv(
            output_path,
            sep='\t',
            index=False,
        )
        self.logger.info(f"Exported predictions to {output_path}")
