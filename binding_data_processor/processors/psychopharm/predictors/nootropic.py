"""Nootropic effects predictor.

This module provides the NootropicPredictor class that:
1. Predicts nootropic mechanisms of action
2. Predicts cognitive enhancement effects
3. Predicts potential side effects
4. Uses ensemble of ML models for robust predictions
5. Provides confidence scores and supporting data
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler

from ....models.core import CompoundData
from ....models.psychopharm import NootropicMechanism
from ..base import PredictionResult
from .base import PredictorBase


class NootropicPredictor(PredictorBase):
    """Predict nootropic effects and mechanisms."""

    # Cognitive domains and specific effects
    COGNITIVE_DOMAINS = {
        "memory": {
            "working_memory", "long_term_memory", "memory_formation",
            "memory_recall", "memory_consolidation", "spatial_memory",
            "verbal_memory", "episodic_memory", "procedural_memory",
        },
        "attention": {
            "sustained_attention", "divided_attention", "selective_attention",
            "attention_switching", "focus", "concentration", "alertness",
            "vigilance", "mental_clarity",
        },
        "learning": {
            "skill_acquisition", "pattern_recognition", "associative_learning",
            "cognitive_flexibility", "neuroplasticity", "learning_speed",
            "error_correction", "habit_formation",
        },
        "executive_function": {
            "planning", "decision_making", "problem_solving",
            "abstract_thinking", "cognitive_control", "inhibition",
            "task_switching", "working_memory_updating",
        },
        "processing": {
            "processing_speed", "mental_processing", "cognitive_throughput",
            "information_processing", "reaction_time", "mental_speed",
            "cognitive_efficiency",
        },
    }

    # Side effect categories
    SIDE_EFFECTS = {
        "physical": {
            "headache", "insomnia", "anxiety", "jitters", "nausea",
            "appetite_changes", "blood_pressure_changes", "heart_rate_changes",
        },
        "cognitive": {
            "brain_fog", "confusion", "memory_issues", "attention_problems",
            "mood_changes", "irritability", "mental_fatigue",
        },
        "tolerance": {
            "acute_tolerance", "chronic_tolerance", "withdrawal_effects",
            "dependence_risk", "rebound_effects",
        },
    }

    # Default model paths
    DEFAULT_MODEL_DIR = Path("models/nootropic")
    MODEL_FILENAMES = {
        "mechanism": "mechanism_predictor.pkl",
        "effects": "effects_predictor.pkl",
        "side_effects": "side_effects_predictor.pkl",
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        cognitive_domains: Optional[Dict[str, Set[str]]] = None,
        side_effects: Optional[Dict[str, Set[str]]] = None,
    ):
        """Initialize nootropic predictor."""
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        
        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / "nootropic_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

        # Use custom categories if provided
        self.cognitive_domains = cognitive_domains or self.COGNITIVE_DOMAINS
        self.side_effects = side_effects or self.SIDE_EFFECTS
        
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
                'mechanism',
                'mechanism_confidence',
                'domain',
                'effect',
                'effect_score',
                'effect_confidence',
                'side_effect_category',
                'side_effect',
                'risk_score',
                'risk_confidence',
                'timestamp',
            ]
        )

        self.logger.info("NootropicPredictor initialized successfully")

    def _initialize_scalers(self) -> Dict[str, StandardScaler]:
        """Initialize feature scalers for each feature type."""
        self.logger.debug("Initializing feature scalers")
        return {
            feature_type: StandardScaler()
            for feature_type in self.feature_types
        }

    def _load_models(self) -> Dict:
        """Load nootropic prediction models from disk."""
        models = {}
        model_dir = Path(self.model_dir) if self.model_dir else self.DEFAULT_MODEL_DIR

        if model_dir.exists():
            self.logger.info(f"Loading models from {model_dir}")
            try:
                # Load mechanism prediction model
                mechanism_path = model_dir / "mechanism_predictor.pkl"
                if mechanism_path.exists():
                    self.logger.debug("Loading mechanism prediction model")
                    models["mechanism"] = np.load(mechanism_path, allow_pickle=True)
                else:
                    self.logger.warning(
                        "Mechanism predictor not found, initializing new model"
                    )
                    models["mechanism"] = self._initialize_mechanism_model()

                # Load effect prediction models
                effect_models = {}
                for domain, effects in self.cognitive_domains.items():
                    domain_models = {}
                    for effect in effects:
                        model_path = model_dir / f"effect_{effect}.pkl"
                        if model_path.exists():
                            self.logger.debug(f"Loading effect model for {effect}")
                            domain_models[effect] = np.load(
                                model_path, allow_pickle=True
                            )
                        else:
                            self.logger.warning(
                                f"Model not found for {effect}, initializing new model"
                            )
                            domain_models[effect] = self._initialize_effect_model()
                    effect_models[domain] = domain_models
                models["effects"] = effect_models

                # Load side effect prediction models
                side_effect_models = {}
                for category, effects in self.side_effects.items():
                    category_models = {}
                    for effect in effects:
                        model_path = model_dir / f"side_effect_{effect}.pkl"
                        if model_path.exists():
                            self.logger.debug(f"Loading side effect model for {effect}")
                            category_models[effect] = np.load(
                                model_path, allow_pickle=True
                            )
                        else:
                            self.logger.warning(
                                f"Model not found for {effect}, initializing new model"
                            )
                            category_models[effect] = self._initialize_side_effect_model()
                    side_effect_models[category] = category_models
                models["side_effects"] = side_effect_models

            except Exception as e:
                self.logger.error(f"Error loading models: {str(e)}")
                self.logger.info("Initializing new models")
                models = self._initialize_models()
        else:
            self.logger.info(f"Model directory not found: {model_dir}")
            self.logger.info("Initializing new models")
            models = self._initialize_models()

        return models

    def _initialize_mechanism_model(self) -> RandomForestClassifier:
        """Initialize new mechanism prediction model."""
        return RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1,
            verbose=1,
        )

    def _initialize_effect_model(self) -> GradientBoostingRegressor:
        """Initialize new effect prediction model."""
        return GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=1,
        )

    def _initialize_side_effect_model(self) -> GradientBoostingRegressor:
        """Initialize new side effect prediction model."""
        return GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=1,
        )

    def _initialize_models(self) -> Dict:
        """Initialize all nootropic prediction models."""
        self.logger.info("Initializing new models")
        
        models = {
            "mechanism": self._initialize_mechanism_model(),
            "effects": {},
            "side_effects": {},
        }

        # Initialize effect models
        for domain, effects in self.cognitive_domains.items():
            domain_models = {
                effect: self._initialize_effect_model()
                for effect in effects
            }
            models["effects"][domain] = domain_models

        # Initialize side effect models
        for category, effects in self.side_effects.items():
            category_models = {
                effect: self._initialize_side_effect_model()
                for effect in effects
            }
            models["side_effects"][category] = category_models

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
            # Save mechanism prediction model
            mechanism_path = save_dir / "mechanism_predictor.pkl"
            np.save(mechanism_path, self.models["mechanism"])
            version_path = version_dir / "mechanism_predictor.pkl"
            np.save(version_path, self.models["mechanism"])
            self.logger.debug("Saved mechanism prediction model")

            # Save effect prediction models
            for domain, effects in self.models["effects"].items():
                for effect, model in effects.items():
                    # Save current version
                    model_path = save_dir / f"effect_{effect}.pkl"
                    np.save(model_path, model)
                    
                    # Save versioned copy
                    version_path = version_dir / f"effect_{effect}.pkl"
                    np.save(version_path, model)
                    
                    self.logger.debug(f"Saved effect model for {effect}")

            # Save side effect prediction models
            for category, effects in self.models["side_effects"].items():
                for effect, model in effects.items():
                    # Save current version
                    model_path = save_dir / f"side_effect_{effect}.pkl"
                    np.save(model_path, model)
                    
                    # Save versioned copy
                    version_path = version_dir / f"side_effect_{effect}.pkl"
                    np.save(version_path, model)
                    
                    self.logger.debug(f"Saved side effect model for {effect}")

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
        """Generate nootropic predictions for a compound."""
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

            # Predict mechanism of action
            mechanism, mechanism_confidence = self._predict_mechanism(
                self.models["mechanism"],
                features,
            )
            self.logger.debug(
                f"Mechanism prediction: {mechanism} (confidence: {mechanism_confidence:.3f})"
            )

            # Predict cognitive effects
            effect_predictions = {}
            for domain, effects in self.cognitive_domains.items():
                domain_predictions = {}
                for effect in effects:
                    # Predict effect score
                    score, conf = self._predict_effect(
                        self.models["effects"][domain][effect],
                        features,
                        effect,
                    )
                    
                    if conf > 0.5:  # Confidence threshold
                        domain_predictions[effect] = {
                            "score": score,
                            "confidence": conf,
                        }

                        # Update prediction history
                        self.prediction_history = pd.concat([
                            self.prediction_history,
                            pd.DataFrame([{
                                'compound_name': compound.name,
                                'mechanism': mechanism,
                                'mechanism_confidence': mechanism_confidence,
                                'domain': domain,
                                'effect': effect,
                                'effect_score': score,
                                'effect_confidence': conf,
                                'side_effect_category': None,
                                'side_effect': None,
                                'risk_score': None,
                                'risk_confidence': None,
                                'timestamp': pd.Timestamp.now(),
                            }])
                        ], ignore_index=True)

                if domain_predictions:
                    effect_predictions[domain] = domain_predictions

            # Predict side effects
            side_effect_predictions = {}
            for category, effects in self.side_effects.items():
                category_predictions = {}
                for effect in effects:
                    # Predict risk score
                    risk, conf = self._predict_side_effect(
                        self.models["side_effects"][category][effect],
                        features,
                        effect,
                    )
                    
                    if conf > 0.5:  # Confidence threshold
                        category_predictions[effect] = {
                            "risk": risk,
                            "confidence": conf,
                        }

                        # Update prediction history
                        self.prediction_history = pd.concat([
                            self.prediction_history,
                            pd.DataFrame([{
                                'compound_name': compound.name,
                                'mechanism': mechanism,
                                'mechanism_confidence': mechanism_confidence,
                                'domain': None,
                                'effect': None,
                                'effect_score': None,
                                'effect_confidence': None,
                                'side_effect_category': category,
                                'side_effect': effect,
                                'risk_score': risk,
                                'risk_confidence': conf,
                                'timestamp': pd.Timestamp.now(),
                            }])
                        ], ignore_index=True)

                if category_predictions:
                    side_effect_predictions[category] = category_predictions

            # Format result
            result = PredictionResult(
                value=NootropicMechanism(mechanism),
                confidence=mechanism_confidence,
                supporting_data={
                    "effects": effect_predictions,
                    "side_effects": side_effect_predictions,
                    **self._format_supporting_data(
                        features,
                        [(effect, pred["confidence"])
                         for preds in effect_predictions.values()
                         for effect, pred in preds.items()],
                    ),
                },
            )
            
            self.logger.info(
                f"Generated predictions for {compound.name}: "
                f"{result.value.value} (confidence: {result.confidence:.3f})"
            )
            
            return result

        except Exception as e:
            self.logger.error(
                "Error predicting nootropic properties",
                f"Compound {compound.name}: {str(e)}",
                exc_info=True
            )
            return PredictionResult(
                value=NootropicMechanism.UNKNOWN,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _predict_mechanism(
        self,
        model: RandomForestClassifier,
        features: Dict[str, np.ndarray],
    ) -> Tuple[str, float]:
        """Predict mechanism of action."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers["mechanism"].transform(X)

        # Get class probabilities
        probs = model.predict_proba(X_scaled)[0]
        pred_idx = np.argmax(probs)
        
        # Map to NootropicMechanism
        mechanism = NootropicMechanism(
            model.classes_[pred_idx]
        ).value
        
        return mechanism, float(probs[pred_idx])

    def _predict_effect(
        self,
        model: GradientBoostingRegressor,
        features: Dict[str, np.ndarray],
        effect: str,
    ) -> Tuple[float, float]:
        """Predict cognitive effect score."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"effect_{effect}"].transform(X)

        # Get prediction and confidence
        score = model.predict(X_scaled)[0]
        confidence = 1.0 - np.std(
            [est.predict(X_scaled)[0] for est in model.estimators_]
        )

        return float(score), float(confidence)

    def _predict_side_effect(
        self,
        model: GradientBoostingRegressor,
        features: Dict[str, np.ndarray],
        effect: str,
    ) -> Tuple[float, float]:
        """Predict side effect risk score."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"side_effect_{effect}"].transform(X)

        # Get prediction and confidence
        risk = model.predict(X_scaled)[0]
        confidence = 1.0 - np.std(
            [est.predict(X_scaled)[0] for est in model.estimators_]
        )

        return float(risk), float(confidence)

    def get_prediction_statistics(self) -> pd.DataFrame:
        """Get statistics about predictions made so far."""
        stats = pd.DataFrame()
        
        # Mechanism distribution
        stats['mechanism_dist'] = (
            self.prediction_history['mechanism'].value_counts(normalize=True)
        )
        
        # Average mechanism confidence
        stats['mechanism_confidence'] = (
            self.prediction_history.groupby('mechanism')['mechanism_confidence'].mean()
        )
        
        # Effect score distribution
        stats['effect_score'] = (
            self.prediction_history.groupby(['domain', 'effect'])['effect_score'].mean()
        )
        
        # Effect confidence distribution
        stats['effect_confidence'] = (
            self.prediction_history.groupby(['domain', 'effect'])['effect_confidence'].mean()
        )
        
        # Side effect risk distribution
        stats['risk_score'] = (
            self.prediction_history.groupby(['side_effect_category', 'side_effect'])['risk_score'].mean()
        )
        
        # Side effect confidence distribution
        stats['risk_confidence'] = (
            self.prediction_history.groupby(['side_effect_category', 'side_effect'])['risk_confidence'].mean()
        )
        
        return stats

    def retrain(
        self,
        compounds: List[CompoundData],
        mechanism_labels: List[NootropicMechanism],
        effect_data: Dict[str, Dict[str, List[float]]],
        side_effect_data: Dict[str, Dict[str, List[float]]],
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

        # Train mechanism prediction model
        self.logger.debug("Training mechanism prediction model")
        self.models["mechanism"].fit(X, mechanism_labels)
        score = self.models["mechanism"].score(X, mechanism_labels)
        metrics["mechanism_score"] = score
        self.logger.debug(f"Mechanism prediction score: {score:.3f}")

        # Train effect prediction models
        for domain, effects in effect_data.items():
            for effect, scores in effects.items():
                if effect in self.cognitive_domains.get(domain, []):
                    self.logger.debug(f"Training effect model for {effect}")
                    model = self.models["effects"][domain][effect]
                    model.fit(X, scores)
                    score = model.score(X, scores)
                    metrics[f"effect_{effect}_score"] = score
                    self.logger.debug(
                        f"Effect model for {effect} training score: {score:.3f}"
                    )

        # Train side effect prediction models
        for category, effects in side_effect_data.items():
            for effect, risks in effects.items():
                if effect in self.side_effects.get(category, []):
                    self.logger.debug(f"Training side effect model for {effect}")
                    model = self.models["side_effects"][category][effect]
                    model.fit(X, risks)
                    score = model.score(X, risks)
                    metrics[f"side_effect_{effect}_score"] = score
                    self.logger.debug(
                        f"Side effect model for {effect} training score: {score:.3f}"
                    )

        # Save updated models
        self.save_models()
        
        self.logger.info("Model retraining completed successfully")
        return metrics
