"""Prediction utilities for nootropic effects."""

import logging
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor

from .....models.core import CompoundData
from .....models.psychopharm import NootropicMechanism
from ...base import PredictionResult


def predict_mechanism(
    model: RandomForestClassifier,
    features: Dict[str, np.ndarray],
    feature_types: List[str],
    scalers: Dict[str, object],
    logger: logging.Logger,
) -> Tuple[str, float]:
    """Predict mechanism of action."""
    # Combine and scale features
    X = np.hstack([features[ft] for ft in feature_types])
    X_scaled = scalers["mechanism"].transform(X)

    # Get class probabilities
    probs = model.predict_proba(X_scaled)[0]
    pred_idx = np.argmax(probs)

    # Map to NootropicMechanism
    mechanism = NootropicMechanism(model.classes_[pred_idx]).value

    return mechanism, float(probs[pred_idx])


def predict_effect(
    model: GradientBoostingRegressor,
    features: Dict[str, np.ndarray],
    feature_types: List[str],
    scalers: Dict[str, object],
    effect: str,
) -> Tuple[float, float]:
    """Predict cognitive effect score."""
    # Combine and scale features
    X = np.hstack([features[ft] for ft in feature_types])
    X_scaled = scalers[f"effect_{effect}"].transform(X)

    # Get prediction and confidence
    score = model.predict(X_scaled)[0]
    confidence = 1.0 - np.std([est.predict(X_scaled)[0] for est in model.estimators_])

    return float(score), float(confidence)


def predict_side_effect(
    model: GradientBoostingRegressor,
    features: Dict[str, np.ndarray],
    feature_types: List[str],
    scalers: Dict[str, object],
    effect: str,
) -> Tuple[float, float]:
    """Predict side effect risk score."""
    # Combine and scale features
    X = np.hstack([features[ft] for ft in feature_types])
    X_scaled = scalers[f"side_effect_{effect}"].transform(X)

    # Get prediction and confidence
    risk = model.predict(X_scaled)[0]
    confidence = 1.0 - np.std([est.predict(X_scaled)[0] for est in model.estimators_])

    return float(risk), float(confidence)


def predict_all_effects(
    compound: CompoundData,
    models: Dict,
    cognitive_domains: Dict,
    feature_types: List[str],
    scalers: Dict[str, object],
    logger: logging.Logger,
) -> Tuple[Dict, pd.DataFrame]:
    """Predict all cognitive effects for a compound."""
    effect_predictions = {}
    history_records = []

    # Extract features
    features = {}
    for feature_type in feature_types:
        features[feature_type] = extract_features(compound, feature_type)
        logger.debug(f"Extracted {feature_type} features: shape={features[feature_type].shape}")

    # Predict mechanism
    mechanism, mechanism_confidence = predict_mechanism(
        models["mechanism"],
        features,
        feature_types,
        scalers,
        logger,
    )
    logger.debug(f"Mechanism prediction: {mechanism} (confidence: {mechanism_confidence:.3f})")

    # Predict cognitive effects
    for domain, effects in cognitive_domains.items():
        domain_predictions = {}
        for effect in effects:
            # Predict effect score
            score, conf = predict_effect(
                models["effects"][domain][effect],
                features,
                feature_types,
                scalers,
                effect,
            )

            if conf > 0.5:  # Confidence threshold
                domain_predictions[effect] = {
                    "score": score,
                    "confidence": conf,
                }

                # Add to history records
                history_records.append(
                    {
                        "compound_name": compound.name,
                        "mechanism": mechanism,
                        "mechanism_confidence": mechanism_confidence,
                        "domain": domain,
                        "effect": effect,
                        "effect_score": score,
                        "effect_confidence": conf,
                        "side_effect_category": None,
                        "side_effect": None,
                        "risk_score": None,
                        "risk_confidence": None,
                        "timestamp": pd.Timestamp.now(),
                    }
                )

        if domain_predictions:
            effect_predictions[domain] = domain_predictions

    return effect_predictions, pd.DataFrame(history_records)


def predict_all_side_effects(
    compound: CompoundData,
    models: Dict,
    side_effects: Dict,
    feature_types: List[str],
    scalers: Dict[str, object],
    mechanism: str,
    mechanism_confidence: float,
    logger: logging.Logger,
) -> Tuple[Dict, pd.DataFrame]:
    """Predict all side effects for a compound."""
    side_effect_predictions = {}
    history_records = []

    # Extract features
    features = {}
    for feature_type in feature_types:
        features[feature_type] = extract_features(compound, feature_type)

    # Predict side effects
    for category, effects in side_effects.items():
        category_predictions = {}
        for effect in effects:
            # Predict risk score
            risk, conf = predict_side_effect(
                models["side_effects"][category][effect],
                features,
                feature_types,
                scalers,
                effect,
            )

            if conf > 0.5:  # Confidence threshold
                category_predictions[effect] = {
                    "risk": risk,
                    "confidence": conf,
                }

                # Add to history records
                history_records.append(
                    {
                        "compound_name": compound.name,
                        "mechanism": mechanism,
                        "mechanism_confidence": mechanism_confidence,
                        "domain": None,
                        "effect": None,
                        "effect_score": None,
                        "effect_confidence": None,
                        "side_effect_category": category,
                        "side_effect": effect,
                        "risk_score": risk,
                        "risk_confidence": conf,
                        "timestamp": pd.Timestamp.now(),
                    }
                )

        if category_predictions:
            side_effect_predictions[category] = category_predictions

    return side_effect_predictions, pd.DataFrame(history_records)


def predict_bbb(compound: CompoundData, logger: logging.Logger) -> PredictionResult:
    """Predict BBB permeability using consolidated BBB predictor.

    Args:
        compound: Compound to predict BBB permeability for
        logger: Logger instance for tracking prediction

    Returns:
        PredictionResult containing:
        - value: BBB permeability prediction
        - confidence: Prediction confidence
        - supporting_data: Additional prediction data
    """
    from ...bbb import BBBPredictorEnhanced

    # Initialize BBB predictor with web enrichment
    bbb_predictor = BBBPredictorEnhanced(model_dir="models/bbb", cache_dir="cache/bbb", logger=logger)

    # Get BBB prediction with supporting data
    result = bbb_predictor.predict(compound)

    return result


def extract_features(compound: CompoundData, feature_type: str) -> np.ndarray:
    """Extract features from compound."""
    # This is a placeholder - actual implementation would be in the base class
    pass
