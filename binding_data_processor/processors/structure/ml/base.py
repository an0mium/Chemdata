"""Base classes for machine learning processors with hybrid PyTorch/sklearn support.

This module provides base classes for ML-enabled molecular processing with:
1. Neural network model support via PyTorch
2. Traditional ML model support via scikit-learn 
3. Feature preprocessing and validation
4. Model persistence and evaluation
5. Comprehensive error handling and logging
6. Cross-validation and metrics
7. Model registry and management
8. Device handling and optimization
"""

import abc
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union
from dataclasses import dataclass

import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem
from sklearn.base import BaseEstimator
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, Dataset

from ..base import BaseStructureProcessor
from .core import ModelBase as BaseMLModel
from .neural_net import BaseNeuralNet
from .utils import MolecularFeaturizer

logger = logging.getLogger(__name__)


@dataclass
class MLPredictorConfig:
    """Configuration for ML predictors."""

    model_type: str
    hidden_size: int = 128
    num_layers: int = 3
    dropout: float = 0.1
    learning_rate: float = 0.001
    batch_size: int = 32
    num_epochs: int = 100
    early_stopping: bool = True
    patience: int = 10
    model_dir: str = "models"
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    seed: int = 42


class MLProcessor(BaseStructureProcessor):
    """Base class for ML-enabled molecular processors."""

    def __init__(
        self,
        model_dir: Optional[Union[str, Path]] = None,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
        config: Optional[Dict] = None,
    ):
        """Initialize processor.

        Args:
            model_dir: Directory containing model files
            device: Device to run models on
            config: Optional configuration dictionary
        """
        super().__init__(config or {})
        self.model_dir = Path(model_dir) if model_dir else None
        self.device = torch.device(device)
        self.logger = logging.getLogger(self.__class__.__name__)

        # Feature processing
        self.featurizer = MolecularFeaturizer()
        self.scaler = StandardScaler()

        # Model management
        self.nn_models: Dict[str, BaseNeuralNet] = {}
        self.ml_models: Dict[str, BaseEstimator] = {}
        self.feature_importance: Dict[str, np.ndarray] = {}
        self.cross_val_scores: Dict[str, Dict[str, List[float]]] = {}

    def _validate_input(self, compound: Union[str, Chem.Mol]) -> Optional[Chem.Mol]:
        """Validate and convert input to RDKit molecule."""
        try:
            if isinstance(compound, str):
                mol = Chem.MolFromSmiles(compound)
            else:
                mol = compound

            if mol is None:
                raise ValueError("Invalid compound input")

            return mol
        except Exception as e:
            self.logger.error(f"Error validating input: {str(e)}")
            return None

    def load_model(self, name: str, model: Union[BaseNeuralNet, BaseEstimator], path: Optional[str] = None) -> bool:
        """Load ML model.

        Args:
            name: Model name
            model: Model instance
            path: Optional path to model file

        Returns:
            Success status
        """
        try:
            if isinstance(model, nn.Module):
                self.nn_models[name] = model
                if path:
                    return model.load(path)
                elif self.model_dir:
                    model_path = self.model_dir / f"{name}.pt"
                    if model_path.exists():
                        return model.load(model_path)
            else:
                self.ml_models[name] = model
                if path:
                    import joblib

                    model_data = joblib.load(path)
                    self.ml_models[name] = model_data["model"]
                    self.scaler = model_data.get("scaler", self.scaler)
                    if "feature_importance" in model_data:
                        self.feature_importance[name] = model_data["feature_importance"]
                    if "cross_val_scores" in model_data:
                        self.cross_val_scores[name] = model_data["cross_val_scores"]
                    return True
            return False
        except Exception as e:
            self.logger.error(f"Error loading model {name}: {str(e)}")
            return False

    def save_model(self, name: str, path: str) -> bool:
        """Save trained model.

        Args:
            name: Model name
            path: Save path

        Returns:
            Success status
        """
        try:
            if name in self.nn_models:
                return self.nn_models[name].save(path)
            elif name in self.ml_models:
                import joblib

                model_data = {
                    "model": self.ml_models[name],
                    "scaler": self.scaler,
                    "feature_importance": self.feature_importance.get(name),
                    "cross_val_scores": self.cross_val_scores.get(name),
                }
                joblib.dump(model_data, path)
                return True
            return False
        except Exception as e:
            self.logger.error(f"Error saving model {name}: {str(e)}")
            return False

    def validate_model(
        self,
        model_name: str,
        X: np.ndarray,
        y: np.ndarray,
        cv: int = 5,
        metrics: Optional[List[str]] = None,
    ) -> Dict[str, Dict[str, float]]:
        """Validate model using cross-validation.

        Args:
            model_name: Model name
            X: Feature matrix
            y: Target values
            cv: Number of CV folds
            metrics: Metrics to compute

        Returns:
            Dictionary of validation metrics
        """
        try:
            if metrics is None:
                metrics = ["accuracy", "roc_auc", "precision", "recall"]

            if model_name in self.ml_models:
                model = self.ml_models[model_name]
                scores = cross_val_score(model, X, y, cv=cv, scoring=metrics)

                results = {}
                for metric in metrics:
                    metric_scores = scores[f"test_{metric}"]
                    results[metric] = {
                        "mean": float(np.mean(metric_scores)),
                        "std": float(np.std(metric_scores)),
                    }

                self.cross_val_scores[model_name] = results
                return results

            elif model_name in self.nn_models:
                # TODO: Implement neural net validation
                pass

            raise ValueError(f"Model {model_name} not found")

        except Exception as e:
            self.logger.error(f"Error validating model: {str(e)}")
            return {}

    def get_feature_importance(
        self,
        model_name: str,
        feature_names: Optional[List[str]] = None,
    ) -> Dict[str, float]:
        """Get feature importance scores."""
        try:
            if model_name not in self.feature_importance:
                raise ValueError(f"No feature importance for model {model_name}")

            importance = self.feature_importance[model_name]
            if feature_names is None:
                feature_names = [f"feature_{i}" for i in range(len(importance))]

            return dict(zip(feature_names, importance.tolist()))

        except Exception as e:
            self.logger.error(f"Error getting feature importance: {str(e)}")
            return {}

    def get_model_info(self, model_name: str) -> Dict:
        """Get model information."""
        try:
            info = {
                "type": None,
                "parameters": None,
                "feature_importance": None,
                "cross_val_scores": None,
            }

            if model_name in self.ml_models:
                model = self.ml_models[model_name]
                info["type"] = type(model).__name__
                info["parameters"] = model.get_params()

            elif model_name in self.nn_models:
                model = self.nn_models[model_name]
                info["type"] = type(model).__name__
                info["parameters"] = {name: param.shape for name, param in model.named_parameters()}

            else:
                raise ValueError(f"Model {model_name} not found")

            if model_name in self.feature_importance:
                info["feature_importance"] = self.feature_importance[model_name].tolist()

            if model_name in self.cross_val_scores:
                info["cross_val_scores"] = self.cross_val_scores[model_name]

            return info

        except Exception as e:
            self.logger.error(f"Error getting model info: {str(e)}")
            return {}

    @abc.abstractmethod
    def train(self, *args, **kwargs):
        """Train model(s)."""
        pass

    @abc.abstractmethod
    def predict(self, *args, **kwargs):
        """Make predictions."""
        pass


class MLPredictor(MLProcessor):
    """Base class for machine learning predictors."""

    def __init__(
        self,
        config: Optional[Union[Dict, MLPredictorConfig]] = None,
        **kwargs,
    ):
        """Initialize ML predictor.

        Args:
            config: Predictor configuration
            **kwargs: Additional arguments passed to parent classes
        """
        # Convert dict config to dataclass
        if isinstance(config, dict):
            self.config = MLPredictorConfig(**config)
        elif isinstance(config, MLPredictorConfig):
            self.config = config
        else:
            self.config = MLPredictorConfig()

        # Initialize parent class with config
        super().__init__(
            model_dir=self.config.model_dir,
            device=self.config.device,
            config=kwargs.get("config", {}),
        )

        self.model = None
        self.is_trained = False

    def preprocess(self, data: Any) -> Any:
        """Preprocess input data.

        Args:
            data: Input data

        Returns:
            Preprocessed data
        """
        raise NotImplementedError

    def train(
        self,
        train_data: Any,
        val_data: Optional[Any] = None,
        **kwargs,
    ) -> Dict[str, float]:
        """Train the model.

        Args:
            train_data: Training data
            val_data: Validation data
            **kwargs: Additional training arguments

        Returns:
            Dictionary of training metrics
        """
        raise NotImplementedError

    def predict(self, data: Any, **kwargs) -> Any:
        """Make predictions.

        Args:
            data: Input data
            **kwargs: Additional prediction arguments

        Returns:
            Model predictions
        """
        raise NotImplementedError

    def evaluate(
        self,
        test_data: Any,
        metrics: Optional[List[str]] = None,
        **kwargs,
    ) -> Dict[str, float]:
        """Evaluate model performance.

        Args:
            test_data: Test data
            metrics: List of metrics to compute
            **kwargs: Additional evaluation arguments

        Returns:
            Dictionary of evaluation metrics
        """
        raise NotImplementedError


class MoleculeDataset(Dataset):
    """Dataset class for molecular data."""

    def __init__(
        self,
        features: np.ndarray,
        labels: Optional[np.ndarray] = None,
        transform: Optional[callable] = None,
    ):
        """Initialize dataset.

        Args:
            features: Feature matrix
            labels: Optional label array
            transform: Optional transform function
        """
        self.features = torch.FloatTensor(features)
        self.labels = torch.FloatTensor(labels) if labels is not None else None
        self.transform = transform

    def __len__(self) -> int:
        """Get dataset length."""
        return len(self.features)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        """Get dataset item."""
        x = self.features[idx]
        if self.transform:
            x = self.transform(x)
        if self.labels is not None:
            return x, self.labels[idx]
        return x, None


class ModelRegistry:
    """Registry for ML models."""

    def __init__(self):
        """Initialize registry."""
        self._models: Dict[str, Tuple[type, Optional[Dict]]] = {}
        self.logger = logging.getLogger(self.__class__.__name__)

    def register(self, name: str, model_class: type, config: Optional[Dict] = None, **kwargs) -> None:
        """Register model."""
        if name in self._models:
            self.logger.warning(f"Overwriting existing model registration: {name}")
        self._models[name] = (model_class, config)

    def get_model(self, name: str, **kwargs) -> Union[BaseNeuralNet, BaseEstimator]:
        """Get model instance."""
        if name not in self._models:
            raise ValueError(f"Model not found: {name}")
        model_class, config = self._models[name]
        if config:
            kwargs = {**config, **kwargs}
        return model_class(**kwargs)

    def list_models(self) -> List[str]:
        """List registered models."""
        return list(self._models.keys())
