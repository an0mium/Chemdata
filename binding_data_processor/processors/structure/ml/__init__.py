"""Machine learning module for structure processing.

This module provides comprehensive ML capabilities for:

1. Activity prediction
- Binding affinity prediction 
- Target prediction
- Activity classification
- Structure-activity relationships

2. Structure analysis
- Molecular fingerprints
- Graph representations
- Pharmacophore detection
- Structural similarity

3. Property prediction
- Physicochemical properties
- ADME properties
- Drug-likeness scores
- Synthetic accessibility

4. Deep learning
- Graph neural networks
- Attention mechanisms
- Transfer learning
- Uncertainty estimation

5. Model interpretation
- Feature importance
- Attention visualization 
- Uncertainty quantification
- Model explanation
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union, Any
import numpy as np

# Configure logging
logger = logging.getLogger(__name__)

# Selectively import DeepChem components to avoid DGL dependencies
from deepchem.models.sklearn_models import SklearnModel
from deepchem.models.torch_models.layers import get_activation
from deepchem.data import Dataset
from deepchem.feat import Featurizer
from deepchem.models import Model as DCModel
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from xgboost import XGBRegressor


# Try importing DGL
try:
    import dgl
    import dgl.graphbolt

    HAS_DGL = True
except ImportError:
    logger.warning("DGL not available. Graph-based models will not be available.")
    HAS_DGL = False


class DeepChemModel:
    """Wrapper around DeepChem model functionality.

    This class provides a standardized interface for using DeepChem models
    within the binding_data_processor framework.
    """

    def __init__(
        self,
        model_type: str,
        model_params: Optional[Dict] = None,
        featurizer: Optional[dc.feat.Featurizer] = None,
    ):
        """Initialize DeepChem model wrapper.

        Args:
            model_type: Type of DeepChem model to use
            model_params: Optional model parameters
            featurizer: Optional DeepChem featurizer
        """
        self.model_type = model_type
        self.model_params = model_params or {}
        self.featurizer = featurizer
        self.model = self._init_model()

    def _init_model(self) -> DCModel:
        """Initialize the underlying DeepChem model."""
        if self.model_type in ["graph_conv", "mpnn", "attentive_fp"]:
            logger.warning(f"Graph-based models not available. Falling back to RandomForest for {self.model_type}")
            return SklearnModel(model=RandomForestRegressor(**self.model_params))
        elif self.model_type == "rf":
            return SklearnModel(model=RandomForestRegressor(**self.model_params))
        elif self.model_type == "xgb":
            return SklearnModel(model=XGBRegressor(**self.model_params))
        elif self.model_type == "svm":
            return SklearnModel(model=SVR(**self.model_params))
        else:
            raise ValueError(f"Unsupported model type: {self.model_type}. Available types: rf, xgb, svm")

    def fit(self, dataset: dc.data.Dataset, **kwargs) -> None:
        """Train the model on a dataset."""
        self.model.fit(dataset, **kwargs)

    def predict(self, dataset: dc.data.Dataset, **kwargs) -> np.ndarray:
        """Generate predictions for a dataset."""
        return self.model.predict(dataset, **kwargs)

    def save(self, model_dir: str) -> None:
        """Save model to directory."""
        self.model.save_checkpoint(model_dir=model_dir)

    def load(self, model_dir: str) -> None:
        """Load model from directory."""
        self.model.restore(model_dir=model_dir)


# Core functionality
from .base import MLProcessor, MLPredictorConfig

# Activity prediction
from .activity.binding import BindingPredictor, BindingSitePredictor, BindingModePredictor
from .activity.protein import ProteinPredictor, ProteinStructureAnalyzer
from .activity.interaction import InteractionPredictor
from .activity.base import ActivityPredictor, ActivityPredictorConfig

# Property prediction
from .property import PropertyPredictor

# Toxicity prediction
from .toxicity import ToxicityPredictor

# Graph Neural Networks
from .models.gnn import EnhancedGNN as GraphNeuralNetwork
from .models.gnn import EnhancedGNN  # Also export under original name

# Default configurations for various model types
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
    "property": {
        "model_type": "ensemble",
        "base_models": ["rf", "xgb", "lgb"],
        "meta_model": "lr",
        "cv_folds": 5,
        "use_probabilities": True,
    },
    "toxicity": {
        "model_type": "ensemble",
        "base_models": ["rf", "xgb", "svm"],
        "meta_model": "lr",
        "cv_folds": 5,
        "use_probabilities": True,
        "endpoints": [
            "acute_toxicity",
            "carcinogenicity",
            "mutagenicity",
            "reproductive_toxicity",
            "hepatotoxicity",
            "cardiotoxicity",
            "neurotoxicity",
        ],
    },
}

__all__ = [
    # Base classes
    "MLPredictor",
    "MLPredictorConfig",
    "MLProcessor",
    # Activity prediction
    "ActivityPredictor",
    "ActivityPredictorConfig",
    "BindingPredictor",
    "BindingSitePredictor",
    "BindingModePredictor",
    "InteractionPredictor",
    # Protein prediction and analysis
    "ProteinPredictor",
    "ProteinStructureAnalyzer",
    # Property prediction
    "PropertyPredictor",
    # Toxicity prediction
    "ToxicityPredictor",
    # Graph Neural Networks
    "GraphNeuralNetwork",
    "EnhancedGNN",
    # DeepChem integration
    "DeepChemModel",
]
