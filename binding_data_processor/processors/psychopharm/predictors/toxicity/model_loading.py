"""Model loading utilities for toxicity prediction.

This module handles:
1. Loading pre-trained models from disk with version control
2. Model validation and verification
3. Ensemble model management
4. Feature scaler initialization and management
5. Model initialization with optimized defaults
"""

import logging
import os
from pathlib import Path
from typing import Dict, Optional, List, Set

import numpy as np
import joblib
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler


logger = logging.getLogger(__name__)


def load_models(
    model_dir: Optional[str] = None,
    toxicity_mechanisms: Optional[Dict[str, Set[str]]] = None,
    organ_toxicity: Optional[Dict[str, Set[str]]] = None,
    safety_concerns: Optional[Dict[str, Set[str]]] = None,
    feature_types: Optional[List[str]] = None,
) -> Dict:
    """Load trained toxicity prediction models.

    Args:
        model_dir: Optional directory containing trained models
        toxicity_mechanisms: Optional toxicity mechanism categories
        organ_toxicity: Optional organ toxicity categories
        safety_concerns: Optional safety concern categories
        feature_types: Optional list of feature types to load models for

    Returns:
        Dictionary of loaded models
    """
    if not model_dir:
        model_dir = os.getenv("TOXICITY_MODEL_DIR", "models/toxicity")
    model_dir = Path(model_dir)

    if not feature_types:
        feature_types = ["fingerprints", "descriptors", "enhanced"]

    models = {}
    try:
        # Load mechanism classifier
        mechanism_path = model_dir / "mechanism" / "ensemble.joblib"
        if mechanism_path.exists():
            models["mechanism"] = joblib.load(mechanism_path)
            logger.debug("Loaded mechanism ensemble model")
        else:
            logger.warning("Mechanism ensemble not found, initializing new model")
            models["mechanism"] = initialize_class_model()

        # Load mechanism-specific models
        mechanism_dir = model_dir / "mechanisms"
        if mechanism_dir.exists():
            for category in mechanism_dir.iterdir():
                if category.is_dir():
                    category_models = {}
                    for mechanism in category.iterdir():
                        if mechanism.suffix == ".joblib":
                            model_name = mechanism.stem
                            category_models[model_name] = joblib.load(mechanism)
                            logger.debug(f"Loaded mechanism model: {category.name}/{model_name}")
                    if category_models:
                        models[f"mechanism_{category.name}"] = category_models

        # Initialize missing mechanism models if categories provided
        if toxicity_mechanisms:
            for category, mechanisms in toxicity_mechanisms.items():
                if f"mechanism_{category}" not in models:
                    models[f"mechanism_{category}"] = {m: initialize_mechanism_model() for m in mechanisms}
                    logger.debug(f"Initialized new models for mechanism category: {category}")

        # Load organ toxicity models
        organ_dir = model_dir / "organs"
        if organ_dir.exists():
            for organ in organ_dir.iterdir():
                if organ.is_dir():
                    organ_models = {}
                    for toxicity in organ.iterdir():
                        if toxicity.suffix == ".joblib":
                            model_name = toxicity.stem
                            organ_models[model_name] = joblib.load(toxicity)
                            logger.debug(f"Loaded organ model: {organ.name}/{model_name}")
                    if organ_models:
                        models[f"organ_{organ.name}"] = organ_models

        # Initialize missing organ models if categories provided
        if organ_toxicity:
            for organ, toxicities in organ_toxicity.items():
                if f"organ_{organ}" not in models:
                    models[f"organ_{organ}"] = {t: initialize_organ_model() for t in toxicities}
                    logger.debug(f"Initialized new models for organ: {organ}")

        # Load safety concern models
        safety_dir = model_dir / "safety"
        if safety_dir.exists():
            for category in safety_dir.iterdir():
                if category.is_dir():
                    category_models = {}
                    for concern in category.iterdir():
                        if concern.suffix == ".joblib":
                            model_name = concern.stem
                            category_models[model_name] = joblib.load(concern)
                            logger.debug(f"Loaded safety model: {category.name}/{model_name}")
                    if category_models:
                        models[f"safety_{category.name}"] = category_models

        # Initialize missing safety models if categories provided
        if safety_concerns:
            for category, concerns in safety_concerns.items():
                if f"safety_{category}" not in models:
                    models[f"safety_{category}"] = {c: initialize_concern_model() for c in concerns}
                    logger.debug(f"Initialized new models for safety category: {category}")

        # Load feature scalers
        scaler_dir = model_dir / "scalers"
        if scaler_dir.exists():
            scalers = {}
            for scaler_file in scaler_dir.iterdir():
                if scaler_file.suffix == ".joblib":
                    scaler_name = scaler_file.stem
                    scalers[scaler_name] = joblib.load(scaler_file)
                    logger.debug(f"Loaded scaler: {scaler_name}")
            if scalers:
                models["scalers"] = scalers
        else:
            # Initialize new scalers if not found
            models["scalers"] = initialize_scalers(
                feature_types=feature_types,
                toxicity_mechanisms=toxicity_mechanisms,
                organ_toxicity=organ_toxicity,
                safety_concerns=safety_concerns,
            )
            logger.debug("Initialized new feature scalers")

        if not models:
            logger.warning(f"No models found in {model_dir}")
            return {}

        logger.info(f"Successfully loaded {len(models)} model groups")
        return models

    except Exception as e:
        logger.error(f"Error loading models: {str(e)}")
        return {}


def validate_models(models: Dict) -> bool:
    """Validate loaded models.

    Args:
        models: Dictionary of loaded models

    Returns:
        True if all models are valid, False otherwise
    """
    try:
        # Check mechanism ensemble
        if "mechanism" not in models:
            logger.error("Missing mechanism ensemble model")
            return False
        if not isinstance(models["mechanism"], RandomForestClassifier):
            logger.error("Invalid mechanism ensemble model type")
            return False

        # Check mechanism models
        for key, value in models.items():
            if key.startswith("mechanism_"):
                for model in value.values():
                    if not isinstance(model, GradientBoostingRegressor):
                        logger.error(f"Invalid mechanism model type in {key}")
                        return False

        # Check organ models
        for key, value in models.items():
            if key.startswith("organ_"):
                for model in value.values():
                    if not isinstance(model, GradientBoostingRegressor):
                        logger.error(f"Invalid organ model type in {key}")
                        return False

        # Check safety models
        for key, value in models.items():
            if key.startswith("safety_"):
                for model in value.values():
                    if not isinstance(model, GradientBoostingRegressor):
                        logger.error(f"Invalid safety model type in {key}")
                        return False

        # Check scalers
        if "scalers" in models:
            for scaler in models["scalers"].values():
                if not hasattr(scaler, "transform"):
                    logger.error("Invalid scaler type")
                    return False

        return True

    except Exception as e:
        logger.error(f"Error validating models: {str(e)}")
        return False


def save_models(
    models: Dict,
    save_dir: str,
    feature_types: Optional[List[str]] = None,
) -> bool:
    """Save trained toxicity prediction models.

    Args:
        models: Dictionary of models to save
        save_dir: Directory to save models in
        feature_types: Optional list of feature types models were trained on

    Returns:
        True if models were saved successfully, False otherwise
    """
    try:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)

        # Save mechanism ensemble
        if "mechanism" in models:
            mechanism_dir = save_dir / "mechanism"
            mechanism_dir.mkdir(exist_ok=True)
            joblib.dump(models["mechanism"], mechanism_dir / "ensemble.joblib")
            logger.debug("Saved mechanism ensemble model")

        # Save mechanism models
        for key, value in models.items():
            if key.startswith("mechanism_"):
                category = key.split("_")[1]
                category_dir = save_dir / "mechanisms" / category
                category_dir.mkdir(parents=True, exist_ok=True)
                for model_name, model in value.items():
                    joblib.dump(model, category_dir / f"{model_name}.joblib")
                    logger.debug(f"Saved mechanism model: {category}/{model_name}")

        # Save organ models
        for key, value in models.items():
            if key.startswith("organ_"):
                organ = key.split("_")[1]
                organ_dir = save_dir / "organs" / organ
                organ_dir.mkdir(parents=True, exist_ok=True)
                for model_name, model in value.items():
                    joblib.dump(model, organ_dir / f"{model_name}.joblib")
                    logger.debug(f"Saved organ model: {organ}/{model_name}")

        # Save safety models
        for key, value in models.items():
            if key.startswith("safety_"):
                category = key.split("_")[1]
                category_dir = save_dir / "safety" / category
                category_dir.mkdir(parents=True, exist_ok=True)
                for model_name, model in value.items():
                    joblib.dump(model, category_dir / f"{model_name}.joblib")
                    logger.debug(f"Saved safety model: {category}/{model_name}")

        # Save scalers
        if "scalers" in models:
            scaler_dir = save_dir / "scalers"
            scaler_dir.mkdir(exist_ok=True)
            for scaler_name, scaler in models["scalers"].items():
                joblib.dump(scaler, scaler_dir / f"{scaler_name}.joblib")
                logger.debug(f"Saved scaler: {scaler_name}")

        # Save feature types
        if feature_types:
            with open(save_dir / "feature_types.txt", "w") as f:
                f.write("\n".join(feature_types))
            logger.debug("Saved feature types")

        logger.info(f"Successfully saved models to {save_dir}")
        return True

    except Exception as e:
        logger.error(f"Error saving models: {str(e)}")
        return False


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


def initialize_models(
    toxicity_mechanisms: Optional[Dict[str, Set[str]]] = None,
    organ_toxicity: Optional[Dict[str, Set[str]]] = None,
    safety_concerns: Optional[Dict[str, Set[str]]] = None,
) -> Dict:
    """Initialize all toxicity prediction models.

    Args:
        toxicity_mechanisms: Optional toxicity mechanism categories
        organ_toxicity: Optional organ toxicity categories
        safety_concerns: Optional safety concern categories

    Returns:
        Dictionary of initialized models
    """
    models = {}

    # Initialize mechanism classifier
    models["mechanism"] = initialize_class_model()
    logger.debug("Initialized mechanism classifier")

    # Initialize mechanism-specific models
    if toxicity_mechanisms:
        for category, mechanisms in toxicity_mechanisms.items():
            category_models = {m: initialize_mechanism_model() for m in mechanisms}
            models[f"mechanism_{category}"] = category_models
            logger.debug(f"Initialized models for mechanism category: {category}")

    # Initialize organ toxicity models
    if organ_toxicity:
        for organ, toxicities in organ_toxicity.items():
            organ_models = {t: initialize_organ_model() for t in toxicities}
            models[f"organ_{organ}"] = organ_models
            logger.debug(f"Initialized models for organ: {organ}")

    # Initialize safety concern models
    if safety_concerns:
        for category, concerns in safety_concerns.items():
            concern_models = {c: initialize_concern_model() for c in concerns}
            models[f"safety_{category}"] = concern_models
            logger.debug(f"Initialized models for safety category: {category}")

    logger.info(f"Successfully initialized {len(models)} model groups")
    return models


def initialize_scalers(
    feature_types: List[str],
    toxicity_mechanisms: Optional[Dict[str, Set[str]]] = None,
    organ_toxicity: Optional[Dict[str, Set[str]]] = None,
    safety_concerns: Optional[Dict[str, Set[str]]] = None,
) -> Dict[str, StandardScaler]:
    """Initialize feature scalers.

    Args:
        feature_types: List of feature types to scale
        toxicity_mechanisms: Optional toxicity mechanism categories
        organ_toxicity: Optional organ toxicity categories
        safety_concerns: Optional safety concern categories

    Returns:
        Dictionary mapping feature types to scalers
    """
    scalers = {}

    # Base feature scalers
    for feature_type in feature_types:
        scalers[feature_type] = StandardScaler()

    # Mechanism scalers
    if toxicity_mechanisms:
        for category, mechanisms in toxicity_mechanisms.items():
            for mechanism in mechanisms:
                scalers[f"mechanism_{mechanism}"] = StandardScaler()

    # Organ toxicity scalers
    if organ_toxicity:
        for organ, toxicities in organ_toxicity.items():
            for toxicity in toxicities:
                scalers[f"organ_{organ}_{toxicity}"] = StandardScaler()

    # Safety concern scalers
    if safety_concerns:
        for category, concerns in safety_concerns.items():
            for concern in concerns:
                scalers[f"concern_{concern}"] = StandardScaler()

    return scalers
