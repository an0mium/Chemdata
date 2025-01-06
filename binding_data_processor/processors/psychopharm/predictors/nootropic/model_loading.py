"""Model loading utilities for nootropic prediction."""

import logging
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor


def load_mechanism_model(
    model_dir: Path,
    logger: logging.Logger,
    initialize_fn: callable,
) -> RandomForestClassifier:
    """Load mechanism prediction model."""
    mechanism_path = model_dir / "mechanism_predictor.pkl"
    if mechanism_path.exists():
        logger.debug("Loading mechanism prediction model")
        return np.load(mechanism_path, allow_pickle=True)

    logger.warning("Mechanism predictor not found, initializing new model")
    return initialize_fn()


def load_effect_models(
    model_dir: Path,
    logger: logging.Logger,
    cognitive_domains: Dict,
    initialize_fn: callable,
) -> Dict:
    """Load effect prediction models."""
    effect_models = {}
    for domain, effects in cognitive_domains.items():
        domain_models = {}
        for effect in effects:
            model_path = model_dir / f"effect_{effect}.pkl"
            if model_path.exists():
                logger.debug(f"Loading effect model for {effect}")
                domain_models[effect] = np.load(model_path, allow_pickle=True)
            else:
                logger.warning(f"Model not found for {effect}, initializing new model")
                domain_models[effect] = initialize_fn()
        effect_models[domain] = domain_models
    return effect_models


def load_side_effect_models(
    model_dir: Path,
    logger: logging.Logger,
    side_effects: Dict,
    initialize_fn: callable,
) -> Dict:
    """Load side effect prediction models."""
    side_effect_models = {}
    for category, effects in side_effects.items():
        category_models = {}
        for effect in effects:
            model_path = model_dir / f"side_effect_{effect}.pkl"
            if model_path.exists():
                logger.debug(f"Loading side effect model for {effect}")
                category_models[effect] = np.load(model_path, allow_pickle=True)
            else:
                logger.warning(f"Model not found for {effect}, initializing new model")
                category_models[effect] = initialize_fn()
        side_effect_models[category] = category_models
    return side_effect_models


def load_all_models(
    model_dir: Optional[Path],
    default_dir: Path,
    logger: logging.Logger,
    cognitive_domains: Dict,
    side_effects: Dict,
    initialize_mechanism_fn: callable,
    initialize_effect_fn: callable,
    initialize_side_effect_fn: callable,
) -> Dict:
    """Load all nootropic prediction models."""
    model_dir = model_dir if model_dir else default_dir

    if not model_dir.exists():
        logger.info(f"Model directory not found: {model_dir}")
        logger.info("Initializing new models")
        return {
            "mechanism": initialize_mechanism_fn(),
            "effects": {
                domain: {effect: initialize_effect_fn() for effect in effects}
                for domain, effects in cognitive_domains.items()
            },
            "side_effects": {
                category: {effect: initialize_side_effect_fn() for effect in effects}
                for category, effects in side_effects.items()
            },
        }

    logger.info(f"Loading models from {model_dir}")
    try:
        return {
            "mechanism": load_mechanism_model(model_dir, logger, initialize_mechanism_fn),
            "effects": load_effect_models(model_dir, logger, cognitive_domains, initialize_effect_fn),
            "side_effects": load_side_effect_models(model_dir, logger, side_effects, initialize_side_effect_fn),
        }
    except Exception as e:
        logger.error(f"Error loading models: {str(e)}")
        logger.info("Initializing new models")
        return {
            "mechanism": initialize_mechanism_fn(),
            "effects": {
                domain: {effect: initialize_effect_fn() for effect in effects}
                for domain, effects in cognitive_domains.items()
            },
            "side_effects": {
                category: {effect: initialize_side_effect_fn() for effect in effects}
                for category, effects in side_effects.items()
            },
        }
