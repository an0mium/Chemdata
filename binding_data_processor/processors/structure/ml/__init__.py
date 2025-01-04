"""Machine learning models for chemical structure analysis.

This module provides ML models for:
1. Activity prediction
- Binding affinity prediction
- Target prediction
- Activity classification
- Structure-activity relationships

2. Toxicity prediction
- Toxicity classification
- Side effect prediction
- Drug-drug interactions
- Metabolism prediction

3. Property prediction
- Physicochemical properties
- ADME properties
- Drug-likeness scores
- Synthetic accessibility

4. Deep learning
- Graph neural networks
- Molecular fingerprints
- Attention mechanisms
- Transfer learning

5. Model interpretation
- Feature importance
- Attention visualization
- Uncertainty estimation
- Model explanation
"""

import logging
from typing import Dict, List, Optional, Union

from .activity import ActivityPredictor
from .toxicity import ToxicityPredictor
from .property import PropertyPredictor
from .graph import GraphNeuralNetwork
from .deepchem import DeepChemModel
from .base import BaseMLModel
from .utils import (
    ModelRegistry,
    DataPreprocessor,
    FeatureExtractor,
    ModelEvaluator,
    UncertaintyEstimator,
)

# Configure logging
logger = logging.getLogger(__name__)

# Export public interface
__all__ = [
    # Models
    "ActivityPredictor",
    "ToxicityPredictor",
    "PropertyPredictor",
    "GraphNeuralNetwork",
    "DeepChemModel",
    "BaseMLModel",
    # Utilities
    "ModelRegistry",
    "DataPreprocessor",
    "FeatureExtractor",
    "ModelEvaluator",
    "UncertaintyEstimator",
]

# Default model configurations
DEFAULT_MODEL_CONFIGS = {
    "activity": {
        "model_type": "graph",
        "hidden_size": 128,
        "num_layers": 3,
        "dropout": 0.1,
        "learning_rate": 0.001,
        "batch_size": 32,
        "num_epochs": 100,
        "early_stopping": True,
        "patience": 10,
    },
    "toxicity": {
        "model_type": "ensemble",
        "base_models": ["rf", "xgb", "lgb"],
        "meta_model": "lr",
        "cv_folds": 5,
        "use_probabilities": True,
    },
    "property": {
        "model_type": "deepchem",
        "featurizer": "graph_conv",
        "splitter": "random",
        "transformers": ["normalization"],
        "model_dir": "models/property",
    },
}


def get_model(
    model_type: str, config: Optional[Dict] = None, pretrained: bool = True, **kwargs
) -> Union[
    ActivityPredictor,
    ToxicityPredictor,
    PropertyPredictor,
    GraphNeuralNetwork,
    DeepChemModel,
]:
    """
    Get ML model instance.

    Args:
        model_type: Type of model to create
            - activity: Activity prediction
            - toxicity: Toxicity prediction
            - property: Property prediction
            - graph: Graph neural networks
            - deepchem: DeepChem models
        config: Model configuration
        pretrained: Whether to load pretrained weights
        **kwargs: Additional model arguments

    Returns:
        ML model instance

    Raises:
        ValueError: If model_type is unknown
    """
    # Get default config
    model_config = DEFAULT_MODEL_CONFIGS.get(model_type, {}).copy()
    if config:
        model_config.update(config)

    # Add config to kwargs
    kwargs["config"] = model_config

    # Map model types to classes
    models = {
        "activity": ActivityPredictor,
        "toxicity": ToxicityPredictor,
        "property": PropertyPredictor,
        "graph": GraphNeuralNetwork,
        "deepchem": DeepChemModel,
    }

    if model_type not in models:
        raise ValueError(
            f"Unknown model type: {model_type}. "
            f"Available types: {list(models.keys())}"
        )

    try:
        model = models[model_type](**kwargs)
        if pretrained:
            model.load_pretrained()
        return model
    except Exception as e:
        logger.error(f"Error creating {model_type} model: {str(e)}")
        raise


def get_all_models(
    config: Optional[Dict] = None, pretrained: bool = True, **kwargs
) -> Dict[
    str,
    Union[
        ActivityPredictor,
        ToxicityPredictor,
        PropertyPredictor,
        GraphNeuralNetwork,
        DeepChemModel,
    ],
]:
    """
    Get instances of all available models.

    Args:
        config: Model configuration
        pretrained: Whether to load pretrained weights
        **kwargs: Additional model arguments

    Returns:
        Dictionary mapping model names to instances
    """
    return {
        name: get_model(name, config, pretrained, **kwargs)
        for name in [
            "activity",
            "toxicity",
            "property",
            "graph",
            "deepchem",
        ]
    }


# Model registry instance
registry = ModelRegistry()


def register_model(
    name: str, model_class: type, config: Optional[Dict] = None, **kwargs
) -> None:
    """
    Register a new model type.

    Args:
        name: Model name
        model_class: Model class
        config: Default configuration
        **kwargs: Additional registration arguments
    """
    registry.register(name, model_class, config, **kwargs)


def get_registered_model(
    name: str, config: Optional[Dict] = None, **kwargs
) -> BaseMLModel:
    """
    Get registered model instance.

    Args:
        name: Model name
        config: Model configuration
        **kwargs: Additional model arguments

    Returns:
        Model instance
    """
    return registry.get(name, config, **kwargs)
