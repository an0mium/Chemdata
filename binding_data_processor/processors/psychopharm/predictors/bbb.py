"""Blood-brain barrier permeability predictor.

This module provides the BBBPredictor class that:
1. Predicts BBB permeability class and score
2. Predicts transporter interactions (P-gp, BCRP, etc.)
3. Uses ensemble of ML models for robust predictions
4. Provides confidence scores and supporting data
5. Handles model loading and feature extraction
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler

from ....models.core import CompoundData
from ....models.psychopharm import (
    BBBPermeability, PsychopharmMixin, NootropicMechanism, PsychoactiveClass
)
from ..base import PredictionResult
from .base import PredictorBase


class BBBPredictor(PredictorBase):
    """Predict blood-brain barrier permeability properties."""

    # BBB-related transporters
    TRANSPORTERS = {
        "efflux": {
            "p_glycoprotein", "bcrp", "mrp1", "mrp2", "mrp4",
        },
        "uptake": {
            "lat1", "mct1", "glut1", "oatp1a2", "oat3",
        },
        "specialized": {
            "organic_cation_transporter", "amino_acid_transporter",
            "peptide_transporter", "fatty_acid_transporter",
        },
    }

    # Receptor-mediated transport mechanisms
    RECEPTOR_TRANSPORTERS = {
        "transferrin_receptor": {"iron_transport", "antibody_transport"},
        "insulin_receptor": {"insulin_transport", "protein_transport"},
        "leptin_receptor": {"peptide_transport"},
        "ldl_receptor": {"lipid_transport"},
    }

    # Duration-based features
    DURATION_FEATURES = {
        "rapid": {"onset": "<15min", "half_life": "<2h"},
        "intermediate": {"onset": "15-60min", "half_life": "2-12h"},
        "slow": {"onset": ">60min", "half_life": ">12h"},
    }

    # Default model paths
    DEFAULT_MODEL_DIR = Path("models/bbb")
    MODEL_FILENAMES = {
        "rf_classifier": "bbb_rf_classifier.pkl",
        "gb_classifier": "bbb_gb_classifier.pkl",
        "rf_regressor": "bbb_rf_regressor.pkl",
        "gb_regressor": "bbb_gb_regressor.pkl",
        "pgp_classifier": "pgp_classifier.pkl",
        "transporter_classifier": "transporter_classifier.pkl",
        "ensemble_classifier": "ensemble_classifier.pkl",
        "receptor_classifier": "receptor_classifier.pkl",
        "duration_classifier": "duration_classifier.pkl",
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        transporters: Optional[Dict[str, Set[str]]] = None,
        receptor_transporters: Optional[Dict[str, Set[str]]] = None,
    ):
        """Initialize BBB predictor with models and configuration."""
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

        # Use custom transporters if provided
        self.transporters = transporters or self.TRANSPORTERS
        self.receptor_transporters = receptor_transporters or self.RECEPTOR_TRANSPORTERS

        # Configure model weights for ensemble
        self.model_weights = {
            "rf_classifier": 0.25,
            "gb_classifier": 0.15,
            "rf_regressor": 0.15,
            "ensemble_classifier": 0.25,
            "receptor_classifier": 0.1,
            "duration_classifier": 0.1,
        }

        # Initialize prediction history with enhanced data
        self.prediction_history = pd.DataFrame(
            columns=[
                'compound_name',
                'permeability_class',
                'confidence',
                'permeability_score',
                'p_gp_substrate',
                'transporter_category',
                'transporter',
                'is_substrate',
                'transporter_confidence',
                'receptor_mediated',
                'receptor_type',
                'duration_class',
                'half_life',
                'cross_tolerance',
                'timestamp',
            ]
        )

        self.logger.info("BBBPredictor initialized successfully")

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
        
        if model_name in [
            "rf_classifier", "pgp_classifier", "transporter_classifier",
            "receptor_classifier", "duration_classifier"
        ]:
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
        elif model_name == "ensemble_classifier":
            return RandomForestClassifier(
                n_estimators=200,
                max_depth=15,
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

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate BBB permeability predictions for a compound."""
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

            # Get permeability predictions
            permeability_preds = []
            for model_name, model in self.models.items():
                if model_name in ["rf_classifier", "gb_classifier", "ensemble_classifier"]:
                    pred, conf = self._predict_permeability_class(
                        model, features, model_name
                    )
                    permeability_preds.append((pred, conf))
                    self.logger.debug(
                        f"{model_name} prediction: {pred} (confidence: {conf:.3f})"
                    )

            # Get score predictions
            score_preds = []
            for model_name, model in self.models.items():
                if model_name in ["rf_regressor", "gb_regressor"]:
                    pred, conf = self._predict_permeability_score(
                        model, features, model_name
                    )
                    score_preds.append((pred, conf))
                    self.logger.debug(
                        f"{model_name} prediction: {pred:.3f} (confidence: {conf:.3f})"
                    )

            # Predict transporter interactions
            transporter_predictions = {}
            for category, transporters in self.transporters.items():
                category_predictions = {}
                for transporter in transporters:
                    # Predict substrate status
                    is_substrate, conf = self._predict_transporter_substrate(
                        self.models["transporter_classifier"],
                        features,
                        transporter,
                    )
                    
                    if conf > 0.5:  # Confidence threshold
                        category_predictions[transporter] = {
                            "is_substrate": is_substrate,
                            "confidence": conf,
                        }

                        # Update prediction history
                        self.prediction_history = pd.concat([
                            self.prediction_history,
                            pd.DataFrame([{
                                'compound_name': compound.name,
                                'permeability_class': None,
                                'confidence': None,
                                'permeability_score': None,
                                'p_gp_substrate': None,
                                'transporter_category': category,
                                'transporter': transporter,
                                'is_substrate': is_substrate,
                                'transporter_confidence': conf,
                                'receptor_mediated': False,
                                'receptor_type': None,
                                'duration_class': None,
                                'half_life': None,
                                'cross_tolerance': None,
                                'timestamp': pd.Timestamp.now(),
                            }])
                        ], ignore_index=True)

                if category_predictions:
                    transporter_predictions[category] = category_predictions

            # Predict receptor-mediated transport
            receptor_predictions = {}
            for receptor, mechanisms in self.receptor_transporters.items():
                is_substrate, conf = self._predict_receptor_transport(
                    self.models["receptor_classifier"],
                    features,
                    receptor,
                )
                
                if conf > 0.5:  # Confidence threshold
                    receptor_predictions[receptor] = {
                        "is_substrate": is_substrate,
                        "confidence": conf,
                        "mechanisms": mechanisms,
                    }

            # Predict duration class and half-life
            duration_class, half_life = self._predict_duration_metrics(
                self.models["duration_classifier"],
                features,
            )

            # Special handling for P-gp as it's particularly important
            is_pgp, pgp_confidence = self._predict_pgp_substrate(
                self.models["pgp_classifier"], features
            )
            self.logger.debug(
                f"P-gp substrate prediction: {is_pgp} "
                f"(confidence: {pgp_confidence:.3f})"
            )

            # Combine predictions
            permeability_class, class_confidence = self._combine_predictions(
                permeability_preds,
                weights=[
                    self.model_weights["rf_classifier"],
                    self.model_weights["gb_classifier"],
                    self.model_weights["ensemble_classifier"],
                ],
            )
            self.logger.debug(
                f"Combined permeability class: {permeability_class} "
                f"(confidence: {class_confidence:.3f})"
            )

            score, score_confidence = self._combine_predictions(
                score_preds,
                weights=[
                    self.model_weights["rf_regressor"],
                    self.model_weights["gb_regressor"],
                ],
            )
            self.logger.debug(
                f"Combined permeability score: {score:.3f} "
                f"(confidence: {score_confidence:.3f})"
            )

            # Get cross-tolerance data if compound has PsychopharmMixin
            cross_tolerance = []
            if isinstance(compound, PsychopharmMixin):
                cross_tolerance = list(compound.cross_tolerance)

            # Format result
            result = PredictionResult(
                value=BBBPermeability(permeability_class),
                confidence=class_confidence,
                supporting_data={
                    "permeability_score": float(score),
                    "score_confidence": float(score_confidence),
                    "p_gp_substrate": bool(is_pgp),
                    "p_gp_confidence": float(pgp_confidence),
                    "transporters": transporter_predictions,
                    "receptor_transport": receptor_predictions,
                    "duration_metrics": {
                        "class": duration_class,
                        "half_life": half_life,
                    },
                    "cross_tolerance": cross_tolerance,
                    **self._format_supporting_data(
                        features,
                        permeability_preds + score_preds,
                    ),
                },
            )
            
            # Update prediction history
            self.prediction_history = pd.concat([
                self.prediction_history,
                pd.DataFrame([{
                    'compound_name': compound.name,
                    'permeability_class': permeability_class,
                    'confidence': class_confidence,
                    'permeability_score': score,
                    'p_gp_substrate': is_pgp,
                    'transporter_category': None,
                    'transporter': None,
                    'is_substrate': None,
                    'transporter_confidence': None,
                    'receptor_mediated': bool(receptor_predictions),
                    'receptor_type': next(iter(receptor_predictions), None),
                    'duration_class': duration_class,
                    'half_life': half_life,
                    'cross_tolerance': ','.join(cross_tolerance),
                    'timestamp': pd.Timestamp.now(),
                }])
            ], ignore_index=True)
            
            self.logger.info(
                f"Generated predictions for {compound.name}: "
                f"{result.value.value} (confidence: {result.confidence:.3f})"
            )
            
            return result

        except Exception as e:
            self.logger.error(
                "Error predicting BBB properties",
                f"Compound {compound.name}: {str(e)}",
                exc_info=True
            )
            return PredictionResult(
                value=BBBPermeability.UNKNOWN,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

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
        
        # P-gp substrate frequency
        stats['pgp_substrate_freq'] = (
            self.prediction_history['p_gp_substrate'].value_counts(normalize=True)
        )
        
        # Score statistics
        stats['score_mean'] = self.prediction_history['permeability_score'].mean()
        stats['score_std'] = self.prediction_history['permeability_score'].std()
        
        # Duration class distribution
        stats['duration_dist'] = (
            self.prediction_history['duration_class'].value_counts(normalize=True)
        )
        
        # Receptor-mediated transport frequency
        stats['receptor_mediated_freq'] = (
            self.prediction_history['receptor_mediated'].value_counts(normalize=True)
        )
        
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

    def _predict_transporter_substrate(
        self,
        model: RandomForestClassifier,
        features: Dict[str, np.ndarray],
        transporter: str,
    ) -> Tuple[bool, float]:
        """Predict transporter substrate status."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"transporter_{transporter}"].transform(X)

        # Get prediction probabilities
        probs = model.predict_proba(X_scaled)[0]
        is_substrate = bool(np.argmax(probs))
        confidence = float(probs[int(is_substrate)])

        return is_substrate, confidence

    def _predict_receptor_transport(
        self,
        model: RandomForestClassifier,
        features: Dict[str, np.ndarray],
        receptor: str,
    ) -> Tuple[bool, float]:
        """Predict receptor-mediated transport."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"receptor_{receptor}"].transform(X)

        # Get prediction probabilities
        probs = model.predict_proba(X_scaled)[0]
        is_substrate = bool(np.argmax(probs))
        confidence = float(probs[int(is_substrate)])

        return is_substrate, confidence

    def _predict_duration_metrics(
        self,
        model: RandomForestClassifier,
        features: Dict[str, np.ndarray],
    ) -> Tuple[str, Optional[str]]:
        """Predict duration class and half-life."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers["duration"].transform(X)

        # Get prediction probabilities
        probs = model.predict_proba(X_scaled)[0]
        pred_idx = np.argmax(probs)
        
        # Map to duration class
        duration_class = list(self.DURATION_FEATURES.keys())[pred_idx]
        half_life = self.DURATION_FEATURES[duration_class]["half_life"]

        return duration_class, half_life

    def _predict_pgp_substrate(
        self,
        model: RandomForestClassifier,
        features: Dict[str, np.ndarray],
    ) -> Tuple[bool, float]:
        """Predict P-glycoprotein substrate status."""
        # Combine features
        X = np.hstack([features[ft] for ft in self.feature_types])
        
        # Get prediction probabilities
        probs = model.predict_proba(X)[0]
        is_substrate = bool(np.argmax(probs))
        confidence = float(probs[int(is_substrate)])

        return is_substrate, confidence

    def retrain(
        self,
        compounds: List[CompoundData],
        labels: List[BBBPermeability],
        scores: Optional[List[float]] = None,
        pgp_labels: Optional[List[bool]] = None,
        transporter_data: Optional[Dict[str, List[bool]]] = None,
        receptor_data: Optional[Dict[str, List[bool]]] = None,
        duration_labels: Optional[List[str]] = None,
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

        # Train classification models
        if labels:
            self.logger.debug("Training classification models")
            for name, model in self.models.items():
                if name in ["rf_classifier", "gb_classifier", "ensemble_classifier"]:
                    model.fit(X, labels)
                    score = model.score(X, labels)
                    metrics[f"{name}_score"] = score
                    self.logger.debug(f"{name} training score: {score:.3f}")

        # Train regression models
        if scores:
            self.logger.debug("Training regression models")
            for name, model in self.models.items():
                if name in ["rf_regressor", "gb_regressor"]:
                    model.fit(X, scores)
                    score = model.score(X, scores)
                    metrics[f"{name}_score"] = score
                    self.logger.debug(f"{name} training score: {score:.3f}")

        # Train P-gp model
        if pgp_labels:
            self.logger.debug("Training P-gp substrate model")
            self.models["pgp_classifier"].fit(X, pgp_labels)
            score = self.models["pgp_classifier"].score(X, pgp_labels)
            metrics["pgp_classifier_score"] = score
            self.logger.debug(f"P-gp classifier training score: {score:.3f}")

        # Train transporter models
        if transporter_data:
            self.logger.debug("Training transporter models")
            for transporter, labels in transporter_data.items():
                self.models["transporter_classifier"].fit(X, labels)
                score = self.models["transporter_classifier"].score(X, labels)
                metrics[f"transporter_{transporter}_score"] = score
                self.logger.debug(
                    f"Transporter classifier for {transporter} "
                    f"training score: {score:.3f}"
                )

        # Train receptor models
        if receptor_data:
            self.logger.debug("Training receptor models")
            for receptor, labels in receptor_data.items():
                self.models["receptor_classifier"].fit(X, labels)
                score = self.models["receptor_classifier"].score(X, labels)
                metrics[f"receptor_{receptor}_score"] = score
                self.logger.debug(
                    f"Receptor classifier for {receptor} "
                    f"training score: {score:.3f}"
                )

        # Train duration model
        if duration_labels:
            self.logger.debug("Training duration model")
            self.models["duration_classifier"].fit(X, duration_labels)
            score = self.models["duration_classifier"].score(X, duration_labels)
            metrics["duration_classifier_score"] = score
            self.logger.debug(f"Duration classifier training score: {score:.3f}")

        # Save updated models
        self.save_models()
        
        self.logger.info("Model retraining completed successfully")
        return metrics
