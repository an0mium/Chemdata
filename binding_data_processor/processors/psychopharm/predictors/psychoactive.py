"""Psychoactive class predictor.

This module provides the PsychoactiveClassPredictor class that:
1. Predicts primary and secondary psychoactive classes
2. Predicts specific effects and their intensities
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
from ....models.psychopharm import PsychoactiveClass
from ..base import PredictionResult
from .base import PredictorBase


class PsychoactiveClassPredictor(PredictorBase):
    """Predict psychoactive classifications and effects."""

    # Effect categories and specific effects
    EFFECT_CATEGORIES = {
        "cognitive": {
            "attention", "memory", "learning", "decision_making",
            "cognitive_enhancement", "cognitive_suppression",
            "thought_acceleration", "thought_deceleration",
            "conceptual_thinking", "abstract_thinking",
            "creativity", "analysis", "mental_clarity",
        },
        "emotional": {
            "euphoria", "dysphoria", "anxiety_relief", "anxiety_increase",
            "mood_lift", "mood_suppression", "emotional_enhancement",
            "emotional_suppression", "empathy", "sociability",
            "motivation", "apathy",
        },
        "perceptual": {
            "visual_enhancement", "visual_suppression", "visual_distortion",
            "visual_hallucination", "auditory_enhancement", "auditory_suppression",
            "auditory_distortion", "auditory_hallucination", "tactile_enhancement",
            "tactile_suppression", "time_distortion", "dissociation",
        },
        "physical": {
            "stimulation", "sedation", "analgesia", "numbness",
            "muscle_relaxation", "muscle_tension", "motor_enhancement",
            "motor_suppression", "appetite_enhancement", "appetite_suppression",
            "nausea", "dizziness",
        },
        "consciousness": {
            "consciousness_expansion", "consciousness_suppression",
            "spiritual_enhancement", "ego_dissolution", "dream_enhancement",
            "delirium", "confusion", "clarity",
        },
    }

    # Default model paths
    DEFAULT_MODEL_DIR = Path("models/psychoactive")
    MODEL_FILENAMES = {
        "class": "class_predictor.pkl",
        "effects": "effects_predictor.pkl",
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        effect_categories: Optional[Dict[str, Set[str]]] = None,
    ):
        """Initialize psychoactive class predictor."""
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        
        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / "psychoactive_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

        # Use custom effect categories if provided
        self.effect_categories = effect_categories or self.EFFECT_CATEGORIES
        
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
                'primary_class',
                'class_confidence',
                'effect_category',
                'effect',
                'intensity',
                'effect_confidence',
                'timestamp',
            ]
        )

        self.logger.info("PsychoactiveClassPredictor initialized successfully")

    def _initialize_scalers(self) -> Dict[str, StandardScaler]:
        """Initialize feature scalers for each feature type."""
        self.logger.debug("Initializing feature scalers")
        return {
            feature_type: StandardScaler()
            for feature_type in self.feature_types
        }

    def _load_models(self) -> Dict:
        """Load psychoactive prediction models from disk."""
        models = {}
        model_dir = Path(self.model_dir) if self.model_dir else self.DEFAULT_MODEL_DIR

        if model_dir.exists():
            self.logger.info(f"Loading models from {model_dir}")
            try:
                # Load class prediction model
                class_path = model_dir / "class_predictor.pkl"
                if class_path.exists():
                    self.logger.debug("Loading class prediction model")
                    models["class"] = np.load(class_path, allow_pickle=True)
                else:
                    self.logger.warning(
                        "Class predictor not found, initializing new model"
                    )
                    models["class"] = self._initialize_class_model()

                # Load effect prediction models
                effect_models = {}
                for category, effects in self.effect_categories.items():
                    category_models = {}
                    for effect in effects:
                        model_path = model_dir / f"effect_{effect}.pkl"
                        if model_path.exists():
                            self.logger.debug(f"Loading effect model for {effect}")
                            category_models[effect] = np.load(
                                model_path, allow_pickle=True
                            )
                        else:
                            self.logger.warning(
                                f"Model not found for {effect}, initializing new model"
                            )
                            category_models[effect] = self._initialize_effect_model()
                    effect_models[category] = category_models
                models["effects"] = effect_models

            except Exception as e:
                self.logger.error(f"Error loading models: {str(e)}")
                self.logger.info("Initializing new models")
                models = self._initialize_models()
        else:
            self.logger.info(f"Model directory not found: {model_dir}")
            self.logger.info("Initializing new models")
            models = self._initialize_models()

        return models

    def _initialize_class_model(self) -> RandomForestClassifier:
        """Initialize new class prediction model."""
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

    def _initialize_models(self) -> Dict:
        """Initialize all psychoactive prediction models."""
        self.logger.info("Initializing new models")
        
        models = {
            "class": self._initialize_class_model(),
            "effects": {},
        }

        # Initialize effect models
        for category, effects in self.effect_categories.items():
            category_models = {
                effect: self._initialize_effect_model()
                for effect in effects
            }
            models["effects"][category] = category_models

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
            # Save class prediction model
            class_path = save_dir / "class_predictor.pkl"
            np.save(class_path, self.models["class"])
            version_path = version_dir / "class_predictor.pkl"
            np.save(version_path, self.models["class"])
            self.logger.debug("Saved class prediction model")

            # Save effect prediction models
            for category, effects in self.models["effects"].items():
                for effect, model in effects.items():
                    # Save current version
                    model_path = save_dir / f"effect_{effect}.pkl"
                    np.save(model_path, model)
                    
                    # Save versioned copy
                    version_path = version_dir / f"effect_{effect}.pkl"
                    np.save(version_path, model)
                    
                    self.logger.debug(f"Saved effect model for {effect}")

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
        """Generate psychoactive predictions for a compound."""
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

            # Predict psychoactive class
            primary_class, class_confidence = self._predict_class(
                self.models["class"],
                features,
            )
            self.logger.debug(
                f"Class prediction: {primary_class} (confidence: {class_confidence:.3f})"
            )

            # Predict effects for each category
            effect_predictions = {}
            for category, effects in self.effect_categories.items():
                category_predictions = {}
                for effect in effects:
                    # Predict effect intensity
                    intensity, conf = self._predict_effect(
                        self.models["effects"][category][effect],
                        features,
                        effect,
                    )
                    
                    if conf > 0.5:  # Confidence threshold
                        category_predictions[effect] = {
                            "intensity": intensity,
                            "confidence": conf,
                        }

                        # Update prediction history
                        self.prediction_history = pd.concat([
                            self.prediction_history,
                            pd.DataFrame([{
                                'compound_name': compound.name,
                                'primary_class': primary_class,
                                'class_confidence': class_confidence,
                                'effect_category': category,
                                'effect': effect,
                                'intensity': intensity,
                                'effect_confidence': conf,
                                'timestamp': pd.Timestamp.now(),
                            }])
                        ], ignore_index=True)

                if category_predictions:
                    effect_predictions[category] = category_predictions

            # Format result
            result = PredictionResult(
                value=PsychoactiveClass(primary_class),
                confidence=class_confidence,
                supporting_data={
                    "effects": effect_predictions,
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
                "Error predicting psychoactive properties",
                f"Compound {compound.name}: {str(e)}",
                exc_info=True
            )
            return PredictionResult(
                value=PsychoactiveClass.UNKNOWN,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _predict_class(
        self,
        model: RandomForestClassifier,
        features: Dict[str, np.ndarray],
    ) -> Tuple[str, float]:
        """Predict psychoactive class."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers["class"].transform(X)

        # Get class probabilities
        probs = model.predict_proba(X_scaled)[0]
        pred_idx = np.argmax(probs)
        
        # Map to PsychoactiveClass
        psychoactive_class = PsychoactiveClass(
            model.classes_[pred_idx]
        ).value
        
        return psychoactive_class, float(probs[pred_idx])

    def _predict_effect(
        self,
        model: GradientBoostingRegressor,
        features: Dict[str, np.ndarray],
        effect: str,
    ) -> Tuple[float, float]:
        """Predict effect intensity."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"effect_{effect}"].transform(X)

        # Get prediction and confidence
        intensity = model.predict(X_scaled)[0]
        confidence = 1.0 - np.std(
            [est.predict(X_scaled)[0] for est in model.estimators_]
        )

        return float(intensity), float(confidence)

    def get_prediction_statistics(self) -> pd.DataFrame:
        """Get statistics about predictions made so far."""
        stats = pd.DataFrame()
        
        # Class distribution
        stats['class_dist'] = (
            self.prediction_history['primary_class'].value_counts(normalize=True)
        )
        
        # Average class confidence
        stats['class_confidence'] = (
            self.prediction_history.groupby('primary_class')['class_confidence'].mean()
        )
        
        # Effect intensity distribution
        stats['effect_intensity'] = (
            self.prediction_history.groupby(['effect_category', 'effect'])['intensity'].mean()
        )
        
        # Effect confidence distribution
        stats['effect_confidence'] = (
            self.prediction_history.groupby(['effect_category', 'effect'])['effect_confidence'].mean()
        )
        
        return stats

    def retrain(
        self,
        compounds: List[CompoundData],
        class_labels: List[PsychoactiveClass],
        effect_data: Dict[str, Dict[str, List[float]]],
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

        # Train class prediction model
        self.logger.debug("Training class prediction model")
        self.models["class"].fit(X, class_labels)
        score = self.models["class"].score(X, class_labels)
        metrics["class_score"] = score
        self.logger.debug(f"Class prediction score: {score:.3f}")

        # Train effect prediction models
        for category, effects in effect_data.items():
            for effect, intensities in effects.items():
                if effect in self.effect_categories.get(category, []):
                    self.logger.debug(f"Training effect model for {effect}")
                    model = self.models["effects"][category][effect]
                    model.fit(X, intensities)
                    score = model.score(X, intensities)
                    metrics[f"effect_{effect}_score"] = score
                    self.logger.debug(
                        f"Effect model for {effect} training score: {score:.3f}"
                    )

        # Save updated models
        self.save_models()
        
        self.logger.info("Model retraining completed successfully")
        return metrics
