"""Receptor binding profile predictor.

This module provides the ReceptorProfilePredictor class that:
1. Predicts binding affinities for multiple receptor types
2. Predicts activity types (agonist/antagonist/etc.)
3. Uses ensemble of ML models for robust predictions
4. Provides confidence scores and supporting data
5. Handles model loading and feature extraction
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler

from ....models.core import CompoundData, TargetData
from ..base import PredictionResult
from .base import PredictorBase


class ReceptorProfilePredictor(PredictorBase):
    """Predict receptor binding profiles."""

    # Key receptor families to predict
    RECEPTOR_FAMILIES = {
        "serotonin": {
            "5-HT1A", "5-HT1B", "5-HT1D", "5-HT2A", "5-HT2B", "5-HT2C",
            "5-HT3", "5-HT4", "5-HT5A", "5-HT6", "5-HT7"
        },
        "dopamine": {
            "D1", "D2", "D3", "D4", "D5"
        },
        "norepinephrine": {
            "α1A", "α1B", "α1D", "α2A", "α2B", "α2C",
            "β1", "β2", "β3"
        },
        "glutamate": {
            "NMDA", "AMPA", "Kainate", "mGluR1", "mGluR2", "mGluR3",
            "mGluR4", "mGluR5", "mGluR6", "mGluR7", "mGluR8"
        },
        "gaba": {
            "GABA-A", "GABA-B", "GABA-C"
        },
        "opioid": {
            "μ", "κ", "δ", "NOP"
        },
        "cannabinoid": {
            "CB1", "CB2"
        },
        "histamine": {
            "H1", "H2", "H3", "H4"
        },
        "sigma": {
            "σ1", "σ2"
        },
    }

    # Activity types to predict
    ACTIVITY_TYPES = [
        "full_agonist",
        "partial_agonist",
        "antagonist",
        "inverse_agonist",
        "allosteric_modulator",
    ]

    # Default model paths
    DEFAULT_MODEL_DIR = Path("models/receptors")
    MODEL_FILENAMES = {
        "affinity": "affinity_predictor.pkl",
        "activity": "activity_predictor.pkl",
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        receptor_families: Optional[Dict[str, Set[str]]] = None,
    ):
        """Initialize receptor profile predictor."""
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        
        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / "receptor_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

        # Use custom receptor families if provided
        self.receptor_families = receptor_families or self.RECEPTOR_FAMILIES
        
        # Initialize base class
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            feature_types=["fingerprints", "descriptors", "enhanced"],
        )

        # Initialize feature scalers
        self.scalers = self._initialize_scalers()

        # Initialize prediction history
        self.prediction_history = pd.DataFrame(
            columns=[
                'compound_name',
                'receptor',
                'affinity_value',
                'affinity_confidence',
                'activity_type',
                'activity_confidence',
                'timestamp',
            ]
        )

        self.logger.info("ReceptorProfilePredictor initialized successfully")

    def _initialize_scalers(self) -> Dict[str, StandardScaler]:
        """Initialize feature scalers for each feature type."""
        self.logger.debug("Initializing feature scalers")
        return {
            feature_type: StandardScaler()
            for feature_type in self.feature_types
        }

    def _load_models(self) -> Dict:
        """Load receptor prediction models from disk."""
        models = {}
        model_dir = Path(self.model_dir) if self.model_dir else self.DEFAULT_MODEL_DIR

        if model_dir.exists():
            self.logger.info(f"Loading models from {model_dir}")
            try:
                # Load affinity prediction models
                affinity_models = {}
                for family, receptors in self.receptor_families.items():
                    family_models = {}
                    for receptor in receptors:
                        model_path = model_dir / f"affinity_{receptor}.pkl"
                        if model_path.exists():
                            self.logger.debug(f"Loading affinity model for {receptor}")
                            family_models[receptor] = np.load(
                                model_path, allow_pickle=True
                            )
                        else:
                            self.logger.warning(
                                f"Model not found for {receptor}, initializing new model"
                            )
                            family_models[receptor] = self._initialize_affinity_model()
                    affinity_models[family] = family_models
                models["affinity"] = affinity_models

                # Load activity type prediction models
                activity_models = {}
                for family, receptors in self.receptor_families.items():
                    family_models = {}
                    for receptor in receptors:
                        model_path = model_dir / f"activity_{receptor}.pkl"
                        if model_path.exists():
                            self.logger.debug(f"Loading activity model for {receptor}")
                            family_models[receptor] = np.load(
                                model_path, allow_pickle=True
                            )
                        else:
                            self.logger.warning(
                                f"Model not found for {receptor}, initializing new model"
                            )
                            family_models[receptor] = self._initialize_activity_model()
                    activity_models[family] = family_models
                models["activity"] = activity_models

            except Exception as e:
                self.logger.error(f"Error loading models: {str(e)}")
                self.logger.info("Initializing new models")
                models = self._initialize_models()
        else:
            self.logger.info(f"Model directory not found: {model_dir}")
            self.logger.info("Initializing new models")
            models = self._initialize_models()

        return models

    def _initialize_affinity_model(self) -> RandomForestRegressor:
        """Initialize new affinity prediction model."""
        return RandomForestRegressor(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1,
            verbose=1,
        )

    def _initialize_activity_model(self) -> GradientBoostingRegressor:
        """Initialize new activity type prediction model."""
        return GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=1,
        )

    def _initialize_models(self) -> Dict:
        """Initialize all receptor prediction models."""
        self.logger.info("Initializing new models")
        
        models = {
            "affinity": {},
            "activity": {},
        }

        for family, receptors in self.receptor_families.items():
            # Initialize affinity models
            affinity_models = {
                receptor: self._initialize_affinity_model()
                for receptor in receptors
            }
            models["affinity"][family] = affinity_models

            # Initialize activity models
            activity_models = {
                receptor: self._initialize_activity_model()
                for receptor in receptors
            }
            models["activity"][family] = activity_models

        return models

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

        try:
            # Save affinity models
            for family, receptors in self.models["affinity"].items():
                for receptor, model in receptors.items():
                    # Save current version
                    model_path = save_dir / f"affinity_{receptor}.pkl"
                    np.save(model_path, model)
                    
                    # Save versioned copy
                    version_path = version_dir / f"affinity_{receptor}.pkl"
                    np.save(version_path, model)
                    
                    self.logger.debug(f"Saved affinity model for {receptor}")

            # Save activity models
            for family, receptors in self.models["activity"].items():
                for receptor, model in receptors.items():
                    # Save current version
                    model_path = save_dir / f"activity_{receptor}.pkl"
                    np.save(model_path, model)
                    
                    # Save versioned copy
                    version_path = version_dir / f"activity_{receptor}.pkl"
                    np.save(version_path, model)
                    
                    self.logger.debug(f"Saved activity model for {receptor}")

            # Save scalers
            scaler_path = save_dir / "scalers.pkl"
            np.save(scaler_path, self.scalers)
            self.logger.debug("Saved feature scalers")

            # Save prediction history
            history_path = save_dir / "prediction_history.csv"
            self.prediction_history.to_csv(history_path, index=False)
            self.logger.debug("Saved prediction history")

        except Exception as e:
            self.logger.error(f"Error saving models: {str(e)}")

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate receptor binding predictions for a compound."""
        self.logger.debug(f"Generating predictions for {compound.name}")
        try:
            # Extract features
            features = {}
            for feature_type in self.feature_types:
                features[feature_type] = self._extract_features(
                    compound, feature_type
                )
                self.logger.debug(
                    f"Extracted {feature_type} features: shape={features[feature_type].shape}"
                )

            # Generate predictions for each receptor
            predictions = {}
            for family, receptors in self.receptor_families.items():
                family_predictions = {}
                for receptor in receptors:
                    # Predict binding affinity
                    affinity_value, affinity_conf = self._predict_affinity(
                        self.models["affinity"][family][receptor],
                        features,
                        receptor,
                    )
                    
                    # Predict activity type
                    activity_type, activity_conf = self._predict_activity_type(
                        self.models["activity"][family][receptor],
                        features,
                        receptor,
                    )

                    # Store predictions
                    family_predictions[receptor] = {
                        "affinity_value": affinity_value,
                        "affinity_confidence": affinity_conf,
                        "activity_type": activity_type,
                        "activity_confidence": activity_conf,
                    }

                    # Update prediction history
                    self.prediction_history = pd.concat([
                        self.prediction_history,
                        pd.DataFrame([{
                            'compound_name': compound.name,
                            'receptor': receptor,
                            'affinity_value': affinity_value,
                            'affinity_confidence': affinity_conf,
                            'activity_type': activity_type,
                            'activity_confidence': activity_conf,
                            'timestamp': pd.Timestamp.now(),
                        }])
                    ], ignore_index=True)

                predictions[family] = family_predictions

            # Create target data objects
            targets = []
            for family, receptors in predictions.items():
                for receptor, pred in receptors.items():
                    if pred["affinity_confidence"] > 0.5:  # Confidence threshold
                        target = TargetData(
                            common_name=receptor,
                            affinity_value=pred["affinity_value"],
                            affinity_type="Ki",  # Assuming Ki predictions
                            affinity_unit="nM",
                            activity_type=pred["activity_type"],
                            confidence=min(
                                pred["affinity_confidence"],
                                pred["activity_confidence"]
                            ),
                        )
                        targets.append(target)

            # Format result
            result = PredictionResult(
                value=targets,
                confidence=np.mean([
                    t.confidence for t in targets
                ]) if targets else 0.0,
                supporting_data={
                    "predictions": predictions,
                    **self._format_supporting_data(
                        features,
                        [(t.common_name, t.confidence) for t in targets],
                    ),
                },
            )
            
            self.logger.info(
                f"Generated predictions for {compound.name}: "
                f"{len(targets)} targets above confidence threshold"
            )
            
            return result

        except Exception as e:
            self.logger.error(
                f"Error predicting receptor profiles for {compound.name}: {str(e)}",
                exc_info=True
            )
            return PredictionResult(
                value=[],
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _predict_affinity(
        self,
        model: RandomForestRegressor,
        features: Dict[str, np.ndarray],
        receptor: str,
    ) -> Tuple[float, float]:
        """Predict binding affinity for a receptor."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"affinity_{receptor}"].transform(X)

        # Get prediction and confidence
        affinity = model.predict(X_scaled)[0]
        confidence = 1.0 - np.std(
            [tree.predict(X_scaled)[0] for tree in model.estimators_]
        )

        return float(affinity), float(confidence)

    def _predict_activity_type(
        self,
        model: GradientBoostingRegressor,
        features: Dict[str, np.ndarray],
        receptor: str,
    ) -> Tuple[str, float]:
        """Predict activity type for a receptor."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"activity_{receptor}"].transform(X)

        # Get prediction probabilities for each activity type
        probs = model.predict_proba(X_scaled)[0]
        pred_idx = np.argmax(probs)
        
        return self.ACTIVITY_TYPES[pred_idx], float(probs[pred_idx])

    def get_prediction_statistics(self) -> pd.DataFrame:
        """Get statistics about predictions made so far."""
        stats = pd.DataFrame()
        
        # Affinity value distribution
        stats['affinity_mean'] = (
            self.prediction_history.groupby('receptor')['affinity_value'].mean()
        )
        stats['affinity_std'] = (
            self.prediction_history.groupby('receptor')['affinity_value'].std()
        )
        
        # Activity type distribution
        stats['activity_dist'] = (
            self.prediction_history.groupby('receptor')['activity_type']
            .value_counts(normalize=True)
        )
        
        # Average confidence scores
        stats['affinity_confidence'] = (
            self.prediction_history.groupby('receptor')['affinity_confidence'].mean()
        )
        stats['activity_confidence'] = (
            self.prediction_history.groupby('receptor')['activity_confidence'].mean()
        )
        
        return stats

    def retrain(
        self,
        compounds: List[CompoundData],
        affinity_data: Dict[str, List[float]],
        activity_data: Dict[str, List[str]],
        **kwargs,
    ) -> Dict[str, float]:
        """Retrain models with new data."""
        self.logger.info(
            f"Retraining models with {len(compounds)} compounds"
        )
        
        # Extract features
        X = []
        for compound in compounds:
            features = []
            for feature_type in self.feature_types:
                feat = self._extract_features(compound, feature_type)
                features.append(feat)
            X.append(np.hstack(features))
        X = np.vstack(X)

        metrics = {}

        # Train models for each receptor
        for family, receptors in self.receptor_families.items():
            for receptor in receptors:
                # Train affinity model if data available
                if receptor in affinity_data:
                    self.logger.debug(f"Training affinity model for {receptor}")
                    model = self.models["affinity"][family][receptor]
                    model.fit(X, affinity_data[receptor])
                    score = model.score(X, affinity_data[receptor])
                    metrics[f"affinity_{receptor}_score"] = score
                    self.logger.debug(
                        f"Affinity model for {receptor} training score: {score:.3f}"
                    )

                # Train activity model if data available
                if receptor in activity_data:
                    self.logger.debug(f"Training activity model for {receptor}")
                    model = self.models["activity"][family][receptor]
                    model.fit(X, activity_data[receptor])
                    score = model.score(X, activity_data[receptor])
                    metrics[f"activity_{receptor}_score"] = score
                    self.logger.debug(
                        f"Activity model for {receptor} training score: {score:.3f}"
                    )

        # Save updated models
        self.save_models()
        
        self.logger.info("Model retraining completed successfully")
        return metrics
