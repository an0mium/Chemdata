"""Base nootropic predictor.

This module provides the core functionality for predicting nootropic properties
of chemical compounds. It handles:
1. Model loading and management
2. Feature extraction and scaling
3. Basic nootropic predictions
4. Prediction history tracking
5. BBB permeability integration
6. Ensemble model support
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Set

import pandas as pd

from .....models.core import CompoundData
from .....models.psychopharm import NootropicMechanism
from ...base import PredictorBase
from . import model_loading, prediction, features


class NootropicPredictorBase(PredictorBase):
    """Base class for nootropic prediction."""

    # Cognitive domains and specific effects
    COGNITIVE_DOMAINS = {
        "memory": {
            "working_memory",
            "long_term_memory",
            "memory_formation",
            "memory_recall",
            "memory_consolidation",
            "spatial_memory",
            "verbal_memory",
            "episodic_memory",
            "procedural_memory",
        },
        "attention": {
            "sustained_attention",
            "divided_attention",
            "selective_attention",
            "attention_switching",
            "focus",
            "concentration",
            "alertness",
            "vigilance",
            "mental_clarity",
        },
        "learning": {
            "skill_acquisition",
            "pattern_recognition",
            "associative_learning",
            "cognitive_flexibility",
            "neuroplasticity",
            "learning_speed",
            "error_correction",
            "habit_formation",
        },
        "executive_function": {
            "planning",
            "decision_making",
            "problem_solving",
            "abstract_thinking",
            "cognitive_control",
            "inhibition",
            "task_switching",
            "working_memory_updating",
        },
        "processing": {
            "processing_speed",
            "mental_processing",
            "cognitive_throughput",
            "information_processing",
            "reaction_time",
            "mental_speed",
            "cognitive_efficiency",
        },
    }

    # Side effect categories
    SIDE_EFFECTS = {
        "physical": {
            "headache",
            "insomnia",
            "anxiety",
            "jitters",
            "nausea",
            "appetite_changes",
            "blood_pressure_changes",
            "heart_rate_changes",
        },
        "cognitive": {
            "brain_fog",
            "confusion",
            "memory_issues",
            "attention_problems",
            "mood_changes",
            "irritability",
            "mental_fatigue",
        },
        "tolerance": {
            "acute_tolerance",
            "chronic_tolerance",
            "withdrawal_effects",
            "dependence_risk",
            "rebound_effects",
        },
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        cognitive_domains: Optional[Dict[str, Set[str]]] = None,
        side_effects: Optional[Dict[str, Set[str]]] = None,
    ):
        """Initialize nootropic predictor.

        Args:
            model_dir: Directory containing model files
            cache_dir: Directory for caching predictions
            log_level: Logging level
            cognitive_domains: Custom cognitive domain categories
            side_effects: Custom side effect categories
        """
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)

        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / "nootropic_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
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

        # Initialize models and scalers
        self.models = model_loading.load_models(
            model_dir=self.model_dir,
            cognitive_domains=self.cognitive_domains,
            side_effects=self.side_effects,
            logger=self.logger,
        )
        self.scalers = model_loading.initialize_scalers(
            feature_types=self.feature_types,
            cognitive_domains=self.cognitive_domains,
            side_effects=self.side_effects,
            logger=self.logger,
        )

        # Initialize prediction history
        self.prediction_history = pd.DataFrame(
            columns=[
                "compound_name",
                "nootropic_class",
                "class_confidence",
                "domain",
                "effect",
                "effect_score",
                "effect_confidence",
                "side_effect_category",
                "side_effect",
                "risk_score",
                "risk_confidence",
                "bbb_permeability",
                "bbb_confidence",
                "timestamp",
            ]
        )

        self.logger.info("NootropicPredictorBase initialized successfully")

    def predict(self, compound: CompoundData) -> Dict:
        """Generate nootropic predictions for a compound.

        Args:
            compound: Compound to predict nootropic effects for

        Returns:
            Dictionary containing:
            - nootropic_class: Overall nootropic classification
            - class_confidence: Confidence in classification
            - cognitive: Predicted cognitive effects
            - side_effects: Predicted side effects
            - bbb_prediction: BBB permeability prediction
        """
        self.logger.debug(f"Generating predictions for {compound.name}")

        try:
            # Extract features
            compound_features = {}
            for feature_type in self.feature_types:
                compound_features[feature_type] = features.extract_features(compound, feature_type)
                shape = compound_features[feature_type].shape
                feat_msg = f"Extracted {feature_type} features: shape={shape}"
                self.logger.debug(feat_msg)

            # Predict nootropic class
            noot_class, class_confidence = prediction.predict_class(
                self.models["class"],
                compound_features,
                self.feature_types,
                self.scalers,
                self.logger,
            )
            class_msg = f"Nootropic class: {noot_class} (confidence: {class_confidence:.3f})"
            self.logger.debug(class_msg)

            # Get BBB prediction
            bbb_result = prediction.predict_bbb(compound, self.logger)
            pred_part = f"BBB prediction: {bbb_result.value.value}"
            conf_part = f"(confidence: {bbb_result.confidence:.3f})"
            self.logger.debug(f"{pred_part} {conf_part}")

            # Predict cognitive effects with BBB integration
            effect_predictions = prediction.predict_all_effects(
                compound,
                self.models,
                self.cognitive_domains,
                self.feature_types,
                self.scalers,
                noot_class,
                class_confidence,
                bbb_result,
                self.logger,
            )

            # Predict side effects with BBB integration
            side_effect_predictions = prediction.predict_all_side_effects(
                compound,
                self.models,
                self.side_effects,
                self.feature_types,
                self.scalers,
                noot_class,
                class_confidence,
                bbb_result,
                self.logger,
            )

            # Update prediction history
            self._update_history(
                compound,
                noot_class,
                class_confidence,
                effect_predictions,
                side_effect_predictions,
                bbb_result,
            )

            return {
                "nootropic_class": NootropicMechanism(noot_class),
                "class_confidence": class_confidence,
                "cognitive": effect_predictions,
                "side_effects": side_effect_predictions,
                "bbb_prediction": {
                    "value": bbb_result.value.value,
                    "confidence": bbb_result.confidence,
                    **bbb_result.supporting_data,
                },
            }

        except Exception as e:
            error_msg = f"Error predicting nootropic effects for {compound.name}: {str(e)}"
            self.logger.error(error_msg, exc_info=True)
            return {
                "nootropic_class": NootropicMechanism.UNKNOWN,
                "class_confidence": 0.0,
                "error": str(e),
            }

    def _update_history(
        self,
        compound: CompoundData,
        noot_class: str,
        class_confidence: float,
        effect_predictions: Dict,
        side_effect_predictions: Dict,
        bbb_result: Dict,
    ) -> None:
        """Update prediction history with new predictions."""
        # Add cognitive effect predictions
        for domain, effects in effect_predictions.items():
            for effect, data in effects.items():
                self.prediction_history = pd.concat(
                    [
                        self.prediction_history,
                        pd.DataFrame(
                            [
                                {
                                    "compound_name": compound.name,
                                    "nootropic_class": noot_class,
                                    "class_confidence": class_confidence,
                                    "domain": domain,
                                    "effect": effect,
                                    "effect_score": data["score"],
                                    "effect_confidence": data["confidence"],
                                    "side_effect_category": None,
                                    "side_effect": None,
                                    "risk_score": None,
                                    "risk_confidence": None,
                                    "bbb_permeability": bbb_result.value.value,
                                    "bbb_confidence": bbb_result.confidence,
                                    "timestamp": pd.Timestamp.now(),
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )

        # Add side effect predictions
        for category, effects in side_effect_predictions.items():
            for effect, data in effects.items():
                self.prediction_history = pd.concat(
                    [
                        self.prediction_history,
                        pd.DataFrame(
                            [
                                {
                                    "compound_name": compound.name,
                                    "nootropic_class": noot_class,
                                    "class_confidence": class_confidence,
                                    "domain": None,
                                    "effect": None,
                                    "effect_score": None,
                                    "effect_confidence": None,
                                    "side_effect_category": category,
                                    "side_effect": effect,
                                    "risk_score": data["risk"],
                                    "risk_confidence": data["confidence"],
                                    "bbb_permeability": bbb_result.value.value,
                                    "bbb_confidence": bbb_result.confidence,
                                    "timestamp": pd.Timestamp.now(),
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )

    def get_prediction_statistics(self) -> pd.DataFrame:
        """Get statistics about predictions made so far."""
        stats = pd.DataFrame()

        # Nootropic class distribution
        nootropic_counts = self.prediction_history["nootropic_class"].value_counts(normalize=True)
        stats["class_dist"] = nootropic_counts

        # Average class confidence
        class_conf = self.prediction_history.groupby("nootropic_class")["class_confidence"].mean()
        stats["class_confidence"] = class_conf

        # Effect score distribution
        effect_groups = ["domain", "effect"]
        effect_scores = self.prediction_history.groupby(effect_groups)["effect_score"].mean()
        stats["effect_score"] = effect_scores

        # Effect confidence distribution
        effect_conf = self.prediction_history.groupby(effect_groups)["effect_confidence"].mean()
        stats["effect_confidence"] = effect_conf

        # Side effect risk distribution
        side_effect_groups = ["side_effect_category", "side_effect"]
        risk_scores = self.prediction_history.groupby(side_effect_groups)["risk_score"].mean()
        stats["risk_score"] = risk_scores

        # Side effect confidence distribution
        risk_conf = self.prediction_history.groupby(side_effect_groups)["risk_confidence"].mean()
        stats["risk_confidence"] = risk_conf

        # BBB permeability distribution
        bbb_counts = self.prediction_history["bbb_permeability"].value_counts(normalize=True)
        stats["bbb_permeability"] = bbb_counts

        # Average BBB confidence
        bbb_conf = self.prediction_history.groupby("bbb_permeability")["bbb_confidence"].mean()
        stats["bbb_confidence"] = bbb_conf

        return stats
