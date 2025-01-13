"""Model loading utilities for nootropic prediction.

This module handles loading and initialization of ML models for predicting
nootropic effects, including:
1. Loading pretrained models from disk
2. Model configuration and setup
3. Device placement (CPU/GPU)
4. Model validation
5. Support for both PyTorch and scikit-learn models
6. Fallback model initialization
7. Comprehensive error handling
"""

import logging
import os
from pathlib import Path
from typing import Dict, Any, Optional, Union, Callable

import numpy as np
import torch
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor


def load_nootropic_models(
    model_dir: str,
    device: str = "cuda" if torch.cuda.is_available() else "cpu",
    cognitive_domains: Optional[Dict[str, list]] = None,
    side_effects: Optional[Dict[str, list]] = None,
    initialize_mechanism_fn: Optional[Callable] = None,
    initialize_effect_fn: Optional[Callable] = None,
    initialize_side_effect_fn: Optional[Callable] = None,
) -> Dict[str, Dict[str, Any]]:
    """Load nootropic prediction models.

    This function supports loading both PyTorch and scikit-learn models.
    For PyTorch models, it handles device placement.
    For scikit-learn models, it uses numpy load/save.

    Args:
        model_dir: Directory containing model files
        device: Device to load PyTorch models to
        cognitive_domains: Optional mapping of cognitive domains to effects
        side_effects: Optional mapping of side effect categories to effects
        initialize_mechanism_fn: Optional function to initialize new mechanism models
        initialize_effect_fn: Optional function to initialize new effect models
        initialize_side_effect_fn: Optional function to initialize new side effect models

    Returns:
        Dictionary of model categories to models
    """
    logger = logging.getLogger(__name__)
    models = {}

    try:
        model_dir = Path(model_dir)
        if not model_dir.exists():
            logger.warning(f"Model directory {model_dir} does not exist")
            return _initialize_new_models(
                cognitive_domains,
                side_effects,
                initialize_mechanism_fn,
                initialize_effect_fn,
                initialize_side_effect_fn,
                logger,
            )

        # Load mechanism model
        mechanism_path = model_dir / "mechanism_predictor"
        if (mechanism_path.with_suffix(".pt")).exists():
            # PyTorch model
            models["mechanism"] = torch.load(
                mechanism_path.with_suffix(".pt"),
                map_location=device,
            )
            logger.debug("Loaded PyTorch mechanism model")
        elif (mechanism_path.with_suffix(".pkl")).exists():
            # scikit-learn model
            models["mechanism"] = np.load(
                mechanism_path.with_suffix(".pkl"),
                allow_pickle=True,
            )
            logger.debug("Loaded scikit-learn mechanism model")
        elif initialize_mechanism_fn:
            logger.warning("Mechanism model not found, initializing new model")
            models["mechanism"] = initialize_mechanism_fn()

        # Load effect models
        models["effects"] = {}
        if cognitive_domains:
            for domain, effects in cognitive_domains.items():
                domain_dir = model_dir / domain
                if not domain_dir.exists():
                    if initialize_effect_fn:
                        logger.warning(f"Domain directory {domain} not found, initializing new models")
                        models["effects"][domain] = {effect: initialize_effect_fn() for effect in effects}
                    continue

                models["effects"][domain] = {}
                for effect in effects:
                    effect_path = domain_dir / f"effect_{effect}"
                    if (effect_path.with_suffix(".pt")).exists():
                        # PyTorch model
                        models["effects"][domain][effect] = torch.load(
                            effect_path.with_suffix(".pt"),
                            map_location=device,
                        )
                        logger.debug(f"Loaded PyTorch effect model for {effect}")
                    elif (effect_path.with_suffix(".pkl")).exists():
                        # scikit-learn model
                        models["effects"][domain][effect] = np.load(
                            effect_path.with_suffix(".pkl"),
                            allow_pickle=True,
                        )
                        logger.debug(f"Loaded scikit-learn effect model for {effect}")
                    elif initialize_effect_fn:
                        logger.warning(f"Effect model for {effect} not found, initializing new model")
                        models["effects"][domain][effect] = initialize_effect_fn()

        # Load side effect models
        models["side_effects"] = {}
        if side_effects:
            for category, effects in side_effects.items():
                category_dir = model_dir / category
                if not category_dir.exists():
                    if initialize_side_effect_fn:
                        logger.warning(f"Category directory {category} not found, initializing new models")
                        models["side_effects"][category] = {effect: initialize_side_effect_fn() for effect in effects}
                    continue

                models["side_effects"][category] = {}
                for effect in effects:
                    effect_path = category_dir / f"side_effect_{effect}"
                    if (effect_path.with_suffix(".pt")).exists():
                        # PyTorch model
                        models["side_effects"][category][effect] = torch.load(
                            effect_path.with_suffix(".pt"),
                            map_location=device,
                        )
                        logger.debug(f"Loaded PyTorch side effect model for {effect}")
                    elif (effect_path.with_suffix(".pkl")).exists():
                        # scikit-learn model
                        models["side_effects"][category][effect] = np.load(
                            effect_path.with_suffix(".pkl"),
                            allow_pickle=True,
                        )
                        logger.debug(f"Loaded scikit-learn side effect model for {effect}")
                    elif initialize_side_effect_fn:
                        logger.warning(f"Side effect model for {effect} not found, initializing new model")
                        models["side_effects"][category][effect] = initialize_side_effect_fn()

    except Exception as e:
        logger.error(f"Error loading nootropic models: {str(e)}")
        if any([cognitive_domains, side_effects, initialize_mechanism_fn, initialize_effect_fn, initialize_side_effect_fn]):
            logger.info("Initializing new models due to loading error")
            return _initialize_new_models(
                cognitive_domains,
                side_effects,
                initialize_mechanism_fn,
                initialize_effect_fn,
                initialize_side_effect_fn,
                logger,
            )

    return models


def save_nootropic_models(
    models: Dict[str, Dict[str, Any]],
    model_dir: Union[str, Path],
    use_pytorch: bool = True,
) -> None:
    """Save nootropic prediction models.

    Args:
        models: Dictionary of models to save
        model_dir: Directory to save models to
        use_pytorch: Whether to save as PyTorch models (True) or scikit-learn models (False)
    """
    logger = logging.getLogger(__name__)
    model_dir = Path(model_dir)

    try:
        model_dir.mkdir(parents=True, exist_ok=True)

        # Save mechanism model
        if "mechanism" in models:
            mechanism_path = model_dir / "mechanism_predictor"
            if use_pytorch:
                torch.save(models["mechanism"], mechanism_path.with_suffix(".pt"))
            else:
                np.save(mechanism_path.with_suffix(".pkl"), models["mechanism"])

        # Save effect models
        if "effects" in models:
            for domain, domain_models in models["effects"].items():
                domain_dir = model_dir / domain
                domain_dir.mkdir(exist_ok=True)
                for effect, model in domain_models.items():
                    effect_path = domain_dir / f"effect_{effect}"
                    if use_pytorch:
                        torch.save(model, effect_path.with_suffix(".pt"))
                    else:
                        np.save(effect_path.with_suffix(".pkl"), model)

        # Save side effect models
        if "side_effects" in models:
            for category, category_models in models["side_effects"].items():
                category_dir = model_dir / category
                category_dir.mkdir(exist_ok=True)
                for effect, model in category_models.items():
                    effect_path = category_dir / f"side_effect_{effect}"
                    if use_pytorch:
                        torch.save(model, effect_path.with_suffix(".pt"))
                    else:
                        np.save(effect_path.with_suffix(".pkl"), model)

    except Exception as e:
        logger.error(f"Error saving nootropic models: {str(e)}")
        raise


def _initialize_new_models(
    cognitive_domains: Optional[Dict[str, list]],
    side_effects: Optional[Dict[str, list]],
    initialize_mechanism_fn: Optional[Callable],
    initialize_effect_fn: Optional[Callable],
    initialize_side_effect_fn: Optional[Callable],
    logger: logging.Logger,
) -> Dict[str, Dict[str, Any]]:
    """Initialize new models when loading fails or directory doesn't exist."""
    models = {}

    if initialize_mechanism_fn:
        models["mechanism"] = initialize_mechanism_fn()

    if cognitive_domains and initialize_effect_fn:
        models["effects"] = {domain: {effect: initialize_effect_fn() for effect in effects} for domain, effects in cognitive_domains.items()}

    if side_effects and initialize_side_effect_fn:
        models["side_effects"] = {category: {effect: initialize_side_effect_fn() for effect in effects} for category, effects in side_effects.items()}

    return models
