"""Model loading utilities for toxicity prediction.

This module handles:
1. Loading pre-trained models from disk
2. Initializing new models when needed
3. Managing model versioning
4. Initializing feature scalers
"""

import logging
from pathlib import Path
from typing import Dict, Set

import numpy as np
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler


def _load_class_model(model_dir: Path, logger: logging.Logger) -> RandomForestClassifier:
    """Load toxicity class prediction model."""
    class_path = model_dir / "class_predictor.pkl"
    if class_path.exists():
        logger.debug("Loading toxicity class model")
        return np.load(class_path, allow_pickle=True)
    else:
        logger.warning("Class predictor not found, initializing new model")
        return initialize_class_model()


def _load_mechanism_models(
    model_dir: Path,
    toxicity_mechanisms: Dict[str, Set[str]],
    logger: logging.Logger,
) -> Dict:
    """Load mechanism prediction models."""
    mechanism_models = {}
    for category, mechanisms in toxicity_mechanisms.items():
        category_models = {}
        for mechanism in mechanisms:
            model_path = model_dir / f"mechanism_{mechanism}.pkl"
            if model_path.exists():
                logger.debug(f"Loading mechanism model for {mechanism}")
                category_models[mechanism] = np.load(model_path, allow_pickle=True)
            else:
                msg = f"Model not found for {mechanism}, initializing new model"
                logger.warning(msg)
                category_models[mechanism] = initialize_mechanism_model()
        mechanism_models[category] = category_models
    return mechanism_models


def _load_organ_models(
    model_dir: Path,
    organ_toxicity: Dict[str, Set[str]],
    logger: logging.Logger,
) -> Dict:
    """Load organ toxicity prediction models."""
    organ_models = {}
    for organ, toxicities in organ_toxicity.items():
        organ_type_models = {}
        for toxicity in toxicities:
            model_path = model_dir / f"organ_{organ}_{toxicity}.pkl"
            if model_path.exists():
                logger.debug(f"Loading organ model for {organ}_{toxicity}")
                organ_type_models[toxicity] = np.load(model_path, allow_pickle=True)
            else:
                msg = f"Model not found for {organ}_{toxicity}, initializing new model"
                logger.warning(msg)
                organ_type_models[toxicity] = initialize_organ_model()
        organ_models[organ] = organ_type_models
    return organ_models


def _load_concern_models(
    model_dir: Path,
    safety_concerns: Dict[str, Set[str]],
    logger: logging.Logger,
) -> Dict:
    """Load safety concern prediction models."""
    concern_models = {}
    for category, concerns in safety_concerns.items():
        category_models = {}
        for concern in concerns:
            model_path = model_dir / f"concern_{concern}.pkl"
            if model_path.exists():
                logger.debug(f"Loading concern model for {concern}")
                category_models[concern] = np.load(model_path, allow_pickle=True)
            else:
                msg = f"Model not found for {concern}, initializing new model"
                logger.warning(msg)
                category_models[concern] = initialize_concern_model()
        concern_models[category] = category_models
    return concern_models


def load_models(
    model_dir: str,
    toxicity_mechanisms: Dict[str, Set[str]],
    organ_toxicity: Dict[str, Set[str]],
    safety_concerns: Dict[str, Set[str]],
    logger: logging.Logger,
) -> Dict:
    """Load toxicity prediction models.

    Args:
        model_dir: Directory containing model files
        toxicity_mechanisms: Toxicity mechanism categories
        organ_toxicity: Organ toxicity categories
        safety_concerns: Safety concern categories
        logger: Logger instance

    Returns:
        Dictionary containing loaded models
    """
    models = {}
    model_dir = Path(model_dir) if model_dir else Path("models/toxicity")

    if model_dir.exists():
        logger.info(f"Loading models from {model_dir}")
        try:
            # Load all model types
            models["class"] = _load_class_model(model_dir, logger)
            models["mechanisms"] = _load_mechanism_models(model_dir, toxicity_mechanisms, logger)
            models["organs"] = _load_organ_models(model_dir, organ_toxicity, logger)
            models["concerns"] = _load_concern_models(model_dir, safety_concerns, logger)

        except Exception as e:
            logger.error(f"Error loading models: {str(e)}")
            logger.info("Initializing new models")
            models = initialize_models(
                toxicity_mechanisms,
                organ_toxicity,
                safety_concerns,
            )
    else:
        logger.info(f"Model directory not found: {model_dir}")
        logger.info("Initializing new models")
        models = initialize_models(
            toxicity_mechanisms,
            organ_toxicity,
            safety_concerns,
        )

    return models


def initialize_models(
    toxicity_mechanisms: Dict[str, Set[str]],
    organ_toxicity: Dict[str, Set[str]],
    safety_concerns: Dict[str, Set[str]],
) -> Dict:
    """Initialize new toxicity prediction models.

    Args:
        toxicity_mechanisms: Toxicity mechanism categories
        organ_toxicity: Organ toxicity categories
        safety_concerns: Safety concern categories

    Returns:
        Dictionary containing initialized models
    """
    models = {
        "class": initialize_class_model(),
        "mechanisms": {},
        "organs": {},
        "concerns": {},
    }

    # Initialize mechanism models
    for category, mechanisms in toxicity_mechanisms.items():
        category_models = {mechanism: initialize_mechanism_model() for mechanism in mechanisms}
        models["mechanisms"][category] = category_models

    # Initialize organ toxicity models
    for organ, toxicities in organ_toxicity.items():
        organ_type_models = {toxicity: initialize_organ_model() for toxicity in toxicities}
        models["organs"][organ] = organ_type_models

    # Initialize safety concern models
    for category, concerns in safety_concerns.items():
        category_models = {concern: initialize_concern_model() for concern in concerns}
        models["concerns"][category] = category_models

    return models


def initialize_class_model() -> RandomForestClassifier:
    """Initialize new toxicity class prediction model."""
    return RandomForestClassifier(
        n_estimators=100,
        max_depth=10,
        random_state=42,
        n_jobs=-1,
        verbose=1,
    )


def initialize_mechanism_model() -> GradientBoostingRegressor:
    """Initialize new mechanism prediction model."""
    return GradientBoostingRegressor(
        n_estimators=100,
        max_depth=5,
        random_state=42,
        verbose=1,
    )


def initialize_organ_model() -> GradientBoostingRegressor:
    """Initialize new organ toxicity prediction model."""
    return GradientBoostingRegressor(
        n_estimators=100,
        max_depth=5,
        random_state=42,
        verbose=1,
    )


def initialize_concern_model() -> GradientBoostingRegressor:
    """Initialize new safety concern prediction model."""
    return GradientBoostingRegressor(
        n_estimators=100,
        max_depth=5,
        random_state=42,
        verbose=1,
    )


def initialize_scalers(
    feature_types: list,
    toxicity_mechanisms: Dict[str, Set[str]],
    organ_toxicity: Dict[str, Set[str]],
    safety_concerns: Dict[str, Set[str]],
    logger: logging.Logger,
) -> Dict[str, StandardScaler]:
    """Initialize feature scalers.

    Args:
        feature_types: List of feature types to scale
        toxicity_mechanisms: Toxicity mechanism categories
        organ_toxicity: Organ toxicity categories
        safety_concerns: Safety concern categories
        logger: Logger instance

    Returns:
        Dictionary mapping feature types to scalers
    """
    logger.debug("Initializing feature scalers")
    scalers = {}

    # Base feature scalers
    for feature_type in feature_types:
        scalers[feature_type] = StandardScaler()

    # Mechanism scalers
    for category, mechanisms in toxicity_mechanisms.items():
        for mechanism in mechanisms:
            scalers[f"mechanism_{mechanism}"] = StandardScaler()

    # Organ toxicity scalers
    for organ, toxicities in organ_toxicity.items():
        for toxicity in toxicities:
            scalers[f"organ_{organ}_{toxicity}"] = StandardScaler()

    # Safety concern scalers
    for category, concerns in safety_concerns.items():
        for concern in concerns:
            scalers[f"concern_{concern}"] = StandardScaler()

    return scalers
