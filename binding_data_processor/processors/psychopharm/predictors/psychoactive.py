"""Psychoactive effects and classification predictor.

This module provides predictors for analyzing psychoactive properties:
1. Primary and secondary psychoactive class prediction
2. Effect prediction and intensity estimation
3. Activity type classification
4. Web data enrichment from community sources
5. Uncertainty estimation and confidence scoring
6. Model versioning and prediction history
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler

from .base import WebEnrichedPredictor, PredictorConfig
from ....models.compound import PsychoactiveCompound
from ....models.psychopharm import PsychoactiveClass
from ....models.compound.ml.predictors import PredictionResult
from ....web_enrichment.manager import WebEnrichmentManager


class PsychoactiveConfig(PredictorConfig):
    """Configuration for psychoactive prediction."""
    
    # Effect categories to predict
    EFFECT_CATEGORIES = {
        "cognitive": {
            "attention", "memory", "learning", "decision_making",
            "cognitive_enhancement", "cognitive_suppression",
            "thought_acceleration", "thought_deceleration",
            "conceptual_thinking", "abstract_thinking",
            "creativity", "analysis", "mental_clarity",
            "focus_enhancement", "confusion",
        },
        "emotional": {
            "euphoria", "dysphoria", "anxiety_relief", "anxiety_increase",
            "mood_lift", "mood_suppression", "emotional_enhancement",
            "emotional_suppression", "empathy", "sociability",
            "motivation", "apathy", "irritability",
        },
        "perceptual": {
            "visual_enhancement", "visual_suppression", "visual_distortion",
            "visual_hallucination", "auditory_enhancement", "auditory_suppression",
            "auditory_distortion", "auditory_hallucination", "tactile_enhancement",
            "tactile_suppression", "time_distortion", "dissociation",
            "taste_enhancement", "smell_enhancement",
        },
        "physical": {
            "stimulation", "sedation", "analgesia", "numbness",
            "muscle_relaxation", "muscle_tension", "motor_enhancement",
            "motor_suppression", "appetite_enhancement", "appetite_suppression",
            "nausea", "dizziness", "respiratory_depression",
        },
        "consciousness": {
            "consciousness_expansion", "consciousness_suppression",
            "spiritual_enhancement", "ego_dissolution", "dream_enhancement",
            "delirium", "confusion", "clarity", "unity_experience",
        },
    }

    # Activity types to classify
    ACTIVITY_TYPES = [
        "hallucinogen",
        "stimulant", 
        "depressant",
        "dissociative",
        "deliriant",
        "entactogen",
        "nootropic",
        "oneirogen",
        "psychedelic",
        "antipsychotic",
    ]

    # Default model paths
    DEFAULT_MODEL_DIR = Path("models/psychoactive")

    def __init__(
        self,
        effect_categories: Optional[Dict[str, Set[str]]] = None,
        activity_types: Optional[List[str]] = None,
        **kwargs
    ):
        """Initialize config.
        
        Args:
            effect_categories: Optional custom effect categories
            activity_types: Optional custom activity types
            **kwargs: Additional config parameters
        """
        super().__init__(**kwargs)
        self.effect_categories = effect_categories or self.EFFECT_CATEGORIES
        self.activity_types = activity_types or self.ACTIVITY_TYPES
        
        # Flatten effect list
        self.effects = []
        for effects in self.effect_categories.values():
            self.effects.extend(effects)


class PsychoactivePredictor(WebEnrichedPredictor):
    """Predicts psychoactive properties of compounds."""

    def __init__(
        self,
        config: Optional[PsychoactiveConfig] = None,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize predictor.
        
        Args:
            config: Optional predictor configuration
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            log_level: Logging level
        """
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

        super().__init__(
            config or PsychoactiveConfig(),
            model_dir=model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
        )
        
        # Initialize feature scalers
        self.scalers = self._initialize_scalers()
        
        # Web enrichment
        self.web_manager = WebEnrichmentManager()

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
                'activity_type',
                'activity_confidence',
                'timestamp',
            ]
        )

        self.logger.info("PsychoactivePredictor initialized successfully")

    def _initialize_scalers(self) -> Dict[str, StandardScaler]:
        """Initialize feature scalers for each feature type."""
        self.logger.debug("Initializing feature scalers")
        return {
            feature_type: StandardScaler()
            for feature_type in self.feature_types
        }

    def _load_models(self) -> Dict:
        """Load psychoactive prediction models."""
        models = {
            "class": None,
            "effects": {},
            "activity": None,
        }

        if self.model_dir:
            try:
                # Load class prediction model
                class_path = self.model_dir / "class_predictor.pkl"
                if class_path.exists():
                    self.logger.debug("Loading class prediction model")
                    models["class"] = np.load(class_path, allow_pickle=True)
                else:
                    self.logger.warning(
                        "Class predictor not found, initializing new model"
                    )
                    models["class"] = self._initialize_class_model()

                # Load effect prediction models
                for effect in self.config.effects:
                    model_path = self.model_dir / f"effect_{effect}.pkl"
                    if model_path.exists():
                        self.logger.debug(f"Loading effect model for {effect}")
                        models["effects"][effect] = np.load(model_path, allow_pickle=True)
                    else:
                        self.logger.warning(
                            f"Model not found for {effect}, initializing new model"
                        )
                        models["effects"][effect] = self._initialize_effect_model()

                # Load activity type model
                model_path = self.model_dir / "activity_type.pkl"
                if model_path.exists():
                    self.logger.debug("Loading activity type model")
                    models["activity"] = np.load(model_path, allow_pickle=True)
                else:
                    self.logger.warning(
                        "Activity type model not found, initializing new model"
                    )
                    models["activity"] = self._initialize_activity_model()

            except Exception as e:
                self.logger.error(f"Error loading models: {str(e)}")
                self.logger.info("Initializing new models")
                models = self._initialize_models()
        else:
            self.logger.info("No model directory provided, initializing new models")
            models = self._initialize_models()

        return models

    def _initialize_models(self) -> Dict:
        """Initialize new prediction models."""
        models = {
            "class": self._initialize_class_model(),
            "effects": {},
            "activity": self._initialize_activity_model(),
        }
        
        # Initialize effect models
        for effect in self.config.effects:
            models["effects"][effect] = self._initialize_effect_model()
            
        return models

    def _initialize_class_model(self) -> RandomForestClassifier:
        """Initialize new class prediction model."""
        return RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1,
            verbose=0
        )

    def _initialize_effect_model(self) -> GradientBoostingRegressor:
        """Initialize new effect prediction model."""
        return GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=0
        )

    def _initialize_activity_model(self) -> RandomForestClassifier:
        """Initialize new activity type model."""
        return RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            verbose=0
        )

    def _predict_raw(self, features: np.ndarray) -> Tuple[Dict, float]:
        """Generate psychoactive predictions."""
        predictions = {
            "class": None,
            "effects": {},
            "activity_type": None,
        }
        confidences = []
        
        try:
            # Predict psychoactive class
            class_model = self.models["class"]
            class_probs = class_model.predict_proba([features])[0]
            class_idx = np.argmax(class_probs)
            primary_class = PsychoactiveClass(
                class_model.classes_[class_idx]
            ).value
            class_conf = float(class_probs[class_idx])
            
            predictions["class"] = {
                "value": primary_class,
                "confidence": class_conf,
            }
            confidences.append(class_conf)
            
            # Predict effects
            for effect, model in self.models["effects"].items():
                # Predict effect intensity
                intensity = float(model.predict([features])[0])
                
                # Get prediction confidence
                intensity_conf = 1.0 - np.std([
                    tree.predict([features])[0]
                    for tree in model.estimators_
                ])
                
                predictions["effects"][effect] = {
                    "intensity": intensity,
                    "confidence": intensity_conf,
                }
                confidences.append(intensity_conf)

            # Predict activity type
            activity_model = self.models["activity"]
            activity_probs = activity_model.predict_proba([features])[0]
            activity_type = self.config.activity_types[np.argmax(activity_probs)]
            activity_conf = float(np.max(activity_probs))
            
            predictions["activity_type"] = {
                "type": activity_type,
                "confidence": activity_conf,
            }
            confidences.append(activity_conf)
            
            # Calculate overall confidence
            confidence = float(np.mean(confidences))
            
            return predictions, confidence
            
        except Exception as e:
            self.logger.error(f"Prediction error: {str(e)}")
            return {}, 0.0

    def _process_prediction(
        self,
        prediction: Dict,
        confidence: float,
        compound: PsychoactiveCompound
    ) -> PredictionResult:
        """Process psychoactive predictions."""
        try:
            # Format predictions by category
            results = {
                "class": prediction["class"],
                "effects": {},
                "activity_type": prediction["activity_type"],
            }
            
            for category, effects in self.config.effect_categories.items():
                category_results = {}
                for effect in effects:
                    if effect in prediction["effects"]:
                        category_results[effect] = prediction["effects"][effect]
                        
                        # Update prediction history
                        self.prediction_history = pd.concat([
                            self.prediction_history,
                            pd.DataFrame([{
                                'compound_name': compound.name,
                                'primary_class': prediction["class"]["value"],
                                'class_confidence': prediction["class"]["confidence"],
                                'effect_category': category,
                                'effect': effect,
                                'intensity': prediction["effects"][effect]["intensity"],
                                'effect_confidence': prediction["effects"][effect]["confidence"],
                                'activity_type': prediction["activity_type"]["type"],
                                'activity_confidence': prediction["activity_type"]["confidence"],
                                'timestamp': pd.Timestamp.now(),
                            }])
                        ], ignore_index=True)
                        
                results["effects"][category] = category_results
                
            # Add metadata
            metadata = {
                "effect_categories": self.config.effect_categories,
                "activity_types": self.config.activity_types,
                "model_version": self.config.model_version,
            }
                
            return PredictionResult(
                value=results,
                confidence=confidence,
                metadata=metadata
            )
            
        except Exception as e:
            self.logger.error(f"Error processing predictions: {str(e)}")
            return PredictionResult(
                value={},
                confidence=0.0,
                metadata={"error": str(e)}
            )

    def _get_web_data(self, compound: PsychoactiveCompound) -> Optional[Dict]:
        """Get psychoactive data from web sources."""
        if not self.config.use_web_data:
            return None
            
        try:
            web_data = {
                "class": None,
                "effects": {},
                "activity_type": None,
            }
            
            # Get community reports
            community_data = self.web_manager.get_community_data(
                compound.name,
                sources=[
                    "psychonautwiki",
                    "erowid",
                    "tripsit",
                ]
            )
            if community_data:
                web_data["community"] = community_data
                
            # Get social media mentions
            social_data = self.web_manager.get_social_data(
                compound.name,
                sources=[
                    "reddit",
                    "twitter",
                    "bluesky",
                ]
            )
            if social_data:
                web_data["social"] = social_data
                
            return web_data if web_data["community"] or web_data["social"] else None
            
        except Exception as e:
            self.logger.error(f"Error getting web data: {str(e)}")
            return None

    def _get_web_effect_data(
        self,
        web_data: Dict,
        effect: str
    ) -> Tuple[List[float], List[float]]:
        """Extract effect data from web sources.
        
        Args:
            web_data: Web data dictionary
            effect: Effect name to extract data for
            
        Returns:
            Tuple of (intensities, confidences)
        """
        intensities = []
        confidences = []
        
        # Check community data
        if "community" in web_data:
            for source, data in web_data["community"].items():
                if effect in data:
                    intensities.append(data[effect].get("intensity"))
                    confidences.append(data[effect].get("confidence", 0.5))
                    
        # Check social data
        if "social" in web_data:
            for source, data in web_data["social"].items():
                if effect in data:
                    intensities.append(data[effect].get("intensity"))
                    confidences.append(data[effect].get("confidence", 0.5))
                    
        return intensities, confidences

    def _combine_effect_predictions(
        self,
        ml_data: Dict,
        web_intensities: List[float],
        web_confidences: List[float]
    ) -> Dict:
        """Combine ML and web predictions for an effect.
        
        Args:
            ml_data: ML prediction data
            web_intensities: Web-based intensity values
            web_confidences: Web-based confidence values
            
        Returns:
            Combined prediction dictionary
        """
        if not web_intensities:
            return ml_data
            
        # Weight ML and web predictions
        ml_weight = ml_data.get("confidence", 0.0)
        web_weight = np.mean(web_confidences)
        
        # Calculate weighted average
        combined_intensity = (
            ml_data["intensity"] * ml_weight +
            np.mean(web_intensities) * web_weight
        ) / (ml_weight + web_weight)
        
        combined_conf = (ml_weight + web_weight) / 2
        
        return {
            "intensity": float(combined_intensity),
            "confidence": float(combined_conf),
        }

    def _combine_predictions(
        self,
        ml_result: PredictionResult,
        web_data: Dict,
        compound: PsychoactiveCompound
    ) -> PredictionResult:
        """Combine ML and web-based predictions."""
        try:
            combined_results = {
                "class": ml_result.value["class"],
                "effects": {},
                "activity_type": ml_result.value["activity_type"],
            }
            
            # Process each effect category
            for category, effects in self.config.effect_categories.items():
                category_results = {}
                
                for effect in effects:
                    ml_data = ml_result.value["effects"].get(category, {}).get(effect, {})
                    if not ml_data:
                        continue
                        
                    # Get web data for effect
                    web_intensities, web_confs = self._get_web_effect_data(web_data, effect)
                    
                    # Combine predictions
                    category_results[effect] = self._combine_effect_predictions(
                        ml_data,
                        web_intensities,
                        web_confs
                    )
                        
                if category_results:
                    combined_results["effects"][category] = category_results
                    
            # Calculate overall confidence
            effect_confidences = [
                effect["confidence"]
                for category in combined_results["effects"].values()
                for effect in category.values()
            ]
            confidence = float(np.mean(effect_confidences)) if effect_confidences else 0.0
                    
            return PredictionResult(
                value=combined_results,
                confidence=confidence,
                metadata=ml_result.metadata
            )
            
        except Exception as e:
            self.logger.error(f"Error combining predictions: {str(e)}")
            return ml_result

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save models to disk with versioning."""
        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            save_dir = self.config.DEFAULT_MODEL_DIR

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
            for effect, model in self.models["effects"].items():
                # Save current version
                model_path = save_dir / f"effect_{effect}.pkl"
                np.save(model_path, model)
                
                # Save versioned copy
                version_path = version_dir / f"effect_{effect}.pkl"
                np.save(version_path, model)
                
                self.logger.debug(f"Saved effect model for {effect}")

            # Save activity type model
            activity_path = save_dir / "activity_type.pkl"
            np.save(activity_path, self.models["activity"])
            version_path = version_dir / "activity_type.pkl"
            np.save(version_path, self.models["activity"])
            self.logger.debug("Saved activity type model")

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
        
        # Activity type distribution
        stats['activity_dist'] = (
            self.prediction_history['activity_type'].value_counts(normalize=True)
        )
        
        # Average activity confidence
        stats['activity_confidence'] = (
            self.prediction_history.groupby('activity_type')['activity_confidence'].mean()
        )
        
        # Effect intensity distribution
        stats['effect_intensity'] = (
            self.prediction_history.groupby(['effect_category', 'effect'])['intensity'].mean()
        )
        
        # Effect confidence distribution
        stats['effect_confidence'] = (
            self.prediction_history.groupby(['effect_category', 'effect'])
            ['effect_confidence'].mean()
        )
        
        return stats

    def retrain(
        self,
        compounds: List[PsychoactiveCompound],
        class_labels: List[PsychoactiveClass],
        activity_labels: List[str],
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

        # Train activity type model
        self.logger.debug("Training activity type model")
        self.models["activity"].fit(X, activity_labels)
        score = self.models["activity"].score(X, activity_labels)
        metrics["activity_score"] = score
        self.logger.debug(f"Activity type score: {score:.3f}")

        # Train effect prediction models
        for category, effects in effect_data.items():
            for effect, intensities in effects.items():
                if effect in self.config.effect_categories.get(category, []):
                    self.logger.debug(f"Training effect model for {effect}")
                    model = self.models["effects"][effect]
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
