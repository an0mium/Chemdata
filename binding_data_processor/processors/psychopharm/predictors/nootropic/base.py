"""Nootropic effects prediction.

This module provides prediction of nootropic effects for chemical compounds,
including:
1. Cognitive enhancement prediction
2. Memory effects prediction
3. Focus/attention effects
4. Neuroprotective effects
5. Side effect profiling
6. Mechanism of action prediction
7. BBB permeability integration
8. Literature evidence integration
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from .....models.core import CompoundData, PredictionResult
from ..base import BasePredictor
from .features import extract_nootropic_features
from .model_loading import load_nootropic_models


class NootropicPredictor(BasePredictor):
    """Predictor for nootropic effects of compounds."""

    def __init__(
        self,
        model_dir: str = "../models/nootropic",
        cache_dir: Optional[str] = None,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
        log_level: int = logging.INFO,
    ):
        """Initialize nootropic predictor.

        Args:
            model_dir: Directory containing trained models
            cache_dir: Optional directory for caching results
            device: Device to run models on
            log_level: Logging level
        """
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            device=device,
            log_level=log_level,
        )

        # Load models
        self.models = load_nootropic_models(
            model_dir=model_dir,
            device=device,
        )

        # Initialize feature extraction
        self.feature_types = ["molecular", "pharmacophore", "binding", "literature", "community"]

        # Effect categories
        self.effect_categories = {
            "cognitive": ["memory", "focus", "learning", "reasoning"],
            "neuroprotective": ["antioxidant", "anti_inflammatory", "neuroplasticity", "neurogenesis"],
            "side_effects": ["headache", "insomnia", "anxiety", "tolerance"],
        }

        # Mechanism categories
        self.mechanism_types = ["cholinergic", "glutamatergic", "dopaminergic", "serotonergic", "nootropic_unknown"]

    def predict(
        self,
        compound: CompoundData,
        include_features: bool = False,
    ) -> PredictionResult:
        """Predict nootropic effects for a compound.

        Args:
            compound: Compound to make predictions for
            include_features: Whether to include extracted features in result

        Returns:
            Prediction results including effects and confidence scores
        """
        try:
            # Extract features
            features = extract_nootropic_features(compound, feature_types=self.feature_types, models=self.models, device=self.device)

            # Make predictions
            predictions = {}
            confidences = {}

            # Predict effects
            for category, effects in self.effect_categories.items():
                predictions[category] = {}
                confidences[category] = {}

                for effect in effects:
                    if category in self.models and effect in self.models[category]:
                        model = self.models[category][effect]
                        pred, conf = model.predict(features)
                        predictions[category][effect] = float(pred)
                        confidences[category][effect] = float(conf)

            # Predict mechanism
            mechanism_preds = []
            mechanism_confs = []
            for mtype in self.mechanism_types:
                if "mechanism" in self.models and mtype in self.models["mechanism"]:
                    model = self.models["mechanism"][mtype]
                    pred, conf = model.predict(features)
                    mechanism_preds.append((mtype, float(pred)))
                    mechanism_confs.append(float(conf))

            # Get top mechanism prediction
            if mechanism_preds:
                top_mechanism = max(mechanism_preds, key=lambda x: x[1])
                mechanism_confidence = np.mean(mechanism_confs)
            else:
                top_mechanism = ("nootropic_unknown", 0.0)
                mechanism_confidence = 0.0

            # Construct result
            result = PredictionResult(
                value=top_mechanism[0],
                confidence=float(mechanism_confidence),
                supporting_data={
                    "effects": predictions,
                    "effect_confidences": confidences,
                    "mechanism_probabilities": dict(mechanism_preds),
                },
            )

            if include_features:
                result.supporting_data["features"] = features

            return result

        except Exception as e:
            self.logger.error(f"Error predicting nootropic effects: {str(e)}")
            return PredictionResult(value="ERROR", confidence=0.0, supporting_data={"error": str(e)})

    def get_feature_importance(
        self,
        effect: Optional[str] = None,
        mechanism: Optional[str] = None,
    ) -> Dict[str, float]:
        """Get feature importance scores.

        Args:
            effect: Optional specific effect to get importances for
            mechanism: Optional specific mechanism to get importances for

        Returns:
            Dictionary of feature names to importance scores
        """
        try:
            importances = {}

            if effect:
                # Get importance for specific effect
                category = None
                for cat, effects in self.effect_categories.items():
                    if effect in effects:
                        category = cat
                        break

                if category and category in self.models and effect in self.models[category]:
                    model = self.models[category][effect]
                    importances = model.feature_importance()

            elif mechanism:
                # Get importance for specific mechanism
                if mechanism in self.mechanism_types:
                    if "mechanism" in self.models and mechanism in self.models["mechanism"]:
                        model = self.models["mechanism"][mechanism]
                        importances = model.feature_importance()

            else:
                # Get average importance across all models
                all_importances = []

                # Effect models
                for category, effects in self.effect_categories.items():
                    if category in self.models:
                        for effect in effects:
                            if effect in self.models[category]:
                                model = self.models[category][effect]
                                all_importances.append(model.feature_importance())

                # Mechanism models
                if "mechanism" in self.models:
                    for mtype in self.mechanism_types:
                        if mtype in self.models["mechanism"]:
                            model = self.models["mechanism"][mtype]
                            all_importances.append(model.feature_importance())

                # Average importances
                if all_importances:
                    features = set().union(*[set(imp.keys()) for imp in all_importances])
                    importances = {}
                    for feature in features:
                        scores = [imp.get(feature, 0.0) for imp in all_importances]
                        importances[feature] = float(np.mean(scores))

            return importances

        except Exception as e:
            self.logger.error(f"Error getting feature importance: {str(e)}")
            return {}

    def retrain(
        self,
        compounds: List[CompoundData],
        labels: List[str],
        effects: Dict[str, Dict[str, List[float]]],
        mechanisms: Optional[Dict[str, List[float]]] = None,
    ) -> Dict[str, float]:
        """Retrain models with new data.

        Args:
            compounds: List of compounds to train on
            labels: List of mechanism labels
            effects: Dictionary of effect categories to effects to scores
            mechanisms: Optional dictionary of mechanism types to scores

        Returns:
            Dictionary of metric names to scores
        """
        try:
            metrics = {}

            # Extract features
            features = []
            for compound in compounds:
                feat = extract_nootropic_features(compound, feature_types=self.feature_types, models=self.models, device=self.device)
                features.append(feat)

            # Train effect models
            for category, effect_dict in effects.items():
                if category not in self.models:
                    self.models[category] = {}

                for effect, scores in effect_dict.items():
                    if effect not in self.models[category]:
                        self.models[category][effect] = self._create_model()

                    model = self.models[category][effect]
                    score = model.train(features, scores)
                    metrics[f"{category}_{effect}_score"] = float(score)

            # Train mechanism models
            if mechanisms:
                if "mechanism" not in self.models:
                    self.models["mechanism"] = {}

                for mtype, scores in mechanisms.items():
                    if mtype not in self.models["mechanism"]:
                        self.models["mechanism"][mtype] = self._create_model()

                    model = self.models["mechanism"][mtype]
                    score = model.train(features, scores)
                    metrics[f"mechanism_{mtype}_score"] = float(score)

            # Save updated models
            self._save_models()

            return metrics

        except Exception as e:
            self.logger.error(f"Error retraining models: {str(e)}")
            return {"error": str(e)}

    def _create_model(self):
        """Create a new model instance."""
        raise NotImplementedError

    def _save_models(self):
        """Save all models to disk."""
        try:
            os.makedirs(self.model_dir, exist_ok=True)

            for category, models in self.models.items():
                category_dir = os.path.join(self.model_dir, category)
                os.makedirs(category_dir, exist_ok=True)

                for name, model in models.items():
                    model_path = os.path.join(category_dir, f"{name}.pt")
                    torch.save(model.state_dict(), model_path)

        except Exception as e:
            self.logger.error(f"Error saving models: {str(e)}")
