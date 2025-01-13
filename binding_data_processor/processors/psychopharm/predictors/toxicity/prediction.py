"""Prediction module for toxicity prediction.

This module provides functionality for making toxicity predictions using
trained models, including:
- Mechanism prediction with uncertainty estimation
- Organ toxicity prediction with BBB integration
- Safety concern prediction with confidence scoring
"""

import logging
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor

from .....models.core import CompoundData
from .....models.psychopharm import ToxicityClass
from .....models.compound.ml.ensemble import ClassifierEnsemble, RegressorEnsemble
from ..base import EnsemblePredictor
from .....models.compound.ml.predictors import PredictionResult
from .features import extract_features, extract_all_features, extract_enhanced_features
from . import model_loading
from ..integration import BBBIntegrationMixin

logger = logging.getLogger(__name__)


def predict_mechanism(
    compound: CompoundData,
    models: Dict,
    feature_types: List[str],
    scalers: Dict[str, object],
) -> Tuple[str, float]:
    """Predict toxicity mechanism with uncertainty estimation.

    Args:
        compound: Compound to predict mechanism for
        models: Dictionary of trained models
        feature_types: Types of features to use
        scalers: Feature scalers for each prediction type

    Returns:
        Tuple of (predicted mechanism, confidence score)
    """
    # Extract features with enhanced descriptors
    features = extract_all_features(compound, enhanced=True)

    # Get ensemble predictions
    mechanism_ensemble = models.get("mechanism")
    if not mechanism_ensemble:
        logger.error("Mechanism ensemble not found")
        return ToxicityClass.UNKNOWN.value, 0.0

    # Get predictions from individual models
    predictions = []
    probabilities = []
    for model in mechanism_ensemble.estimators_:
        pred = model.predict([features])[0]
        prob = model.predict_proba([features])[0]
        predictions.append(pred)
        probabilities.append(prob)

    predictions = np.array(predictions)
    probabilities = np.array(probabilities)

    # Get ensemble prediction
    mechanism = mechanism_ensemble.predict([features])[0]
    raw_confidence = np.max(mechanism_ensemble.predict_proba([features])[0])

    # Calculate uncertainty metrics
    disagreement = 1 - (np.sum(predictions == mechanism) / len(predictions))
    prob_variance = np.var(probabilities, axis=0)

    # Combine metrics for final confidence
    confidence = raw_confidence * (1 - disagreement) * (1 - np.mean(prob_variance))

    return mechanism, confidence


def predict_organ_effects(
    compound: CompoundData,
    models: Dict,
    organ_toxicity: Dict[str, List[str]],
    feature_types: List[str],
    scalers: Dict[str, object],
    bbb_result: PredictionResult,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Predict organ-specific toxicity effects with BBB integration.

    Args:
        compound: Compound to predict effects for
        models: Dictionary of trained models
        organ_toxicity: Dictionary mapping organs to toxicity types
        feature_types: Types of features to use
        scalers: Feature scalers for each prediction type
        bbb_result: BBB prediction result for CNS effects

    Returns:
        Nested dictionary of organ toxicity predictions
    """
    # Extract features with enhanced descriptors
    features = extract_all_features(compound, enhanced=True)

    organ_predictions = {}
    for organ, toxicities in organ_toxicity.items():
        organ_predictions_inner = {}
        for toxicity in toxicities:
            ensemble = models.get(f"organ_{organ}_{toxicity}")
            if not ensemble:
                continue

            # Get base prediction
            severity = ensemble.predict([features])[0]
            conf = 1.0 - np.std([est.predict([features])[0] for est in ensemble.estimators_])

            # Adjust based on BBB permeability for CNS effects
            if conf > 0.5:  # Base confidence threshold
                if organ == "brain":
                    # Scale severity by BBB permeability
                    adjusted_severity = severity * bbb_result.confidence
                    # Combine confidences
                    adjusted_conf = conf * bbb_result.confidence
                else:
                    # Non-CNS effects less dependent on BBB
                    adjusted_severity = severity
                    adjusted_conf = conf

                if adjusted_conf > 0.3:  # Adjusted confidence threshold
                    organ_predictions_inner[toxicity] = {
                        "severity": adjusted_severity,
                        "confidence": adjusted_conf,
                        "bbb_factor": bbb_result.confidence,
                    }

        if organ_predictions_inner:
            organ_predictions[organ] = organ_predictions_inner

    return organ_predictions


def predict_safety_concerns(
    compound: CompoundData,
    models: Dict,
    safety_concerns: Dict[str, List[str]],
    feature_types: List[str],
    scalers: Dict[str, object],
    bbb_result: PredictionResult,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Predict safety concerns with BBB integration.

    Args:
        compound: Compound to predict concerns for
        models: Dictionary of trained models
        safety_concerns: Dictionary mapping categories to concerns
        feature_types: Types of features to use
        scalers: Feature scalers for each prediction type
        bbb_result: BBB prediction result for CNS-related concerns

    Returns:
        Nested dictionary of safety concern predictions
    """
    # Extract features with enhanced descriptors
    features = extract_all_features(compound, enhanced=True)

    concern_predictions = {}
    for category, concerns in safety_concerns.items():
        category_predictions = {}
        for concern in concerns:
            ensemble = models.get(f"safety_{category}_{concern}")
            if not ensemble:
                continue

            # Get base prediction
            risk = ensemble.predict([features])[0]
            conf = 1.0 - np.std([est.predict([features])[0] for est in ensemble.estimators_])

            # Adjust based on BBB permeability for CNS-related concerns
            if conf > 0.5:  # Base confidence threshold
                if category in ["neurological", "cognitive", "behavioral"]:
                    # Scale risk by BBB permeability
                    adjusted_risk = risk * bbb_result.confidence
                    # Combine confidences
                    adjusted_conf = conf * bbb_result.confidence
                else:
                    # Non-CNS concerns less dependent on BBB
                    adjusted_risk = risk
                    adjusted_conf = conf

                if adjusted_conf > 0.3:  # Adjusted confidence threshold
                    category_predictions[concern] = {
                        "risk": adjusted_risk,
                        "confidence": adjusted_conf,
                        "bbb_factor": bbb_result.confidence,
                    }

        if category_predictions:
            concern_predictions[category] = category_predictions

    return concern_predictions


def predict_all_effects(
    compound: CompoundData,
    models: Dict,
    toxicity_mechanisms: Dict[str, List[str]],
    organ_toxicity: Dict[str, List[str]],
    safety_concerns: Dict[str, List[str]],
    feature_types: List[str],
    scalers: Dict[str, object],
    bbb_result: Optional[PredictionResult] = None,
) -> PredictionResult:
    """Generate comprehensive toxicity predictions.

    Args:
        compound: Compound to predict effects for
        models: Dictionary of trained models
        toxicity_mechanisms: Dictionary mapping categories to mechanisms
        organ_toxicity: Dictionary mapping organs to toxicity types
        safety_concerns: Dictionary mapping categories to concerns
        feature_types: Types of features to use
        scalers: Feature scalers for each prediction type
        bbb_result: Optional BBB prediction result

    Returns:
        PredictionResult containing all toxicity predictions
    """
    try:
        # Get BBB prediction if not provided
        if not bbb_result:
            from ..bbb import BBBPredictorEnhanced

            bbb_predictor = BBBPredictorEnhanced()
            bbb_result = bbb_predictor.predict(compound)

        # Predict mechanism
        mechanism, mechanism_confidence = predict_mechanism(
            compound,
            models,
            feature_types,
            scalers,
        )

        # Predict organ toxicity
        organ_predictions = predict_organ_effects(
            compound,
            models,
            organ_toxicity,
            feature_types,
            scalers,
            bbb_result,
        )

        # Predict safety concerns
        concern_predictions = predict_safety_concerns(
            compound,
            models,
            safety_concerns,
            feature_types,
            scalers,
            bbb_result,
        )

        # Format result
        result = PredictionResult(
            value=ToxicityClass(mechanism),
            confidence=mechanism_confidence,
            supporting_data={
                "bbb_prediction": {
                    "value": bbb_result.value,
                    "confidence": bbb_result.confidence,
                    **bbb_result.supporting_data,
                },
                "organs": organ_predictions,
                "concerns": concern_predictions,
            },
        )

        return result

    except Exception as e:
        logger.error(f"Error generating predictions: {str(e)}")
        return PredictionResult(
            value=ToxicityClass.UNKNOWN,
            confidence=0.0,
            supporting_data={"error": str(e)},
        )
