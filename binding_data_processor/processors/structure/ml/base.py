"""Base classes for machine learning processors with hybrid PyTorch/sklearn support.

This module provides base classes for ML-enabled molecular processing with:
1. Neural network model support via PyTorch
2. Traditional ML model support via scikit-learn
3. Feature preprocessing and validation
4. Model persistence and evaluation
5. Comprehensive error handling and logging
"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem
from sklearn.base import BaseEstimator
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, Dataset

from ..base import StructureProcessor


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
        """Get dataset item.

        Args:
            idx: Item index

        Returns:
            Tuple of (features, labels)
        """
        x = self.features[idx]
        if self.transform:
            x = self.transform(x)
        if self.labels is not None:
            return x, self.labels[idx]
        return x, None


class BaseNeuralNet(nn.Module):
    """Base neural network model."""

    def __init__(self):
        """Initialize model."""
        super().__init__()
        self.logger = logging.getLogger(self.__class__.__name__)

    def save(self, path: Union[str, Path]) -> bool:
        """Save model weights.

        Args:
            path: Path to save model

        Returns:
            True if successful
        """
        try:
            torch.save(self.state_dict(), path)
            self.logger.info(f"Saved model to {path}")
            return True
        except Exception as e:
            self.logger.error(f"Error saving model: {str(e)}")
            return False

    def load(self, path: Union[str, Path]) -> bool:
        """Load model weights.

        Args:
            path: Path to load model from

        Returns:
            True if successful
        """
        try:
            self.load_state_dict(torch.load(path))
            self.logger.info(f"Loaded model from {path}")
            return True
        except Exception as e:
            self.logger.error(f"Error loading model: {str(e)}")
            return False


class MLProcessor(StructureProcessor):
    """Base class for ML-enabled molecular processors."""

    def __init__(self):
        """Initialize processor."""
        super().__init__()
        self.logger = logging.getLogger(self.__class__.__name__)

        # PyTorch setup
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.nn_models: Dict[str, BaseNeuralNet] = {}

        # Sklearn setup
        self.scaler = StandardScaler()
        self.ml_models: Dict[str, BaseEstimator] = {}
        self.feature_importance: Dict[str, np.ndarray] = {}
        self.cross_val_scores: Dict[str, Dict[str, List[float]]] = {}

    def _validate_input(self, compound: Union[str, Chem.Mol]) -> Optional[Chem.Mol]:
        """Validate and convert input to RDKit molecule.

        Args:
            compound: Input compound

        Returns:
            RDKit molecule or None
        """
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

    def _to_device(self, tensor: torch.Tensor) -> torch.Tensor:
        """Move tensor to correct device."""
        return tensor.to(self.device)

    def _numpy_to_tensor(self, array: np.ndarray, dtype=torch.float32) -> torch.Tensor:
        """Convert numpy array to tensor."""
        return torch.tensor(array, dtype=dtype, device=self.device)

    def _tensor_to_numpy(self, tensor: torch.Tensor) -> np.ndarray:
        """Convert tensor to numpy array."""
        return tensor.detach().cpu().numpy()

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
                # Implement neural net validation
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
        """Get feature importance scores.

        Args:
            model_name: Model name
            feature_names: Feature names

        Returns:
            Feature importance scores
        """
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

    def save_model(self, model_name: str, path: str) -> bool:
        """Save trained model.

        Args:
            model_name: Model name
            path: Save path

        Returns:
            True if successful
        """
        try:
            if model_name in self.ml_models:
                import joblib

                model_data = {
                    "model": self.ml_models[model_name],
                    "scaler": self.scaler,
                    "feature_importance": self.feature_importance.get(model_name),
                    "cross_val_scores": self.cross_val_scores.get(model_name),
                }
                joblib.dump(model_data, path)
                return True

            elif model_name in self.nn_models:
                return self.nn_models[model_name].save(path)

            raise ValueError(f"Model {model_name} not found")

        except Exception as e:
            self.logger.error(f"Error saving model: {str(e)}")
            return False

    def load_model(self, model_name: str, path: str) -> bool:
        """Load trained model.

        Args:
            model_name: Model name
            path: Load path

        Returns:
            True if successful
        """
        try:
            if path.endswith(".pt"):
                # Load PyTorch model
                model = BaseNeuralNet()
                if model.load(path):
                    self.nn_models[model_name] = model
                    return True
                return False

            else:
                # Load sklearn model
                import joblib

                model_data = joblib.load(path)
                self.ml_models[model_name] = model_data["model"]
                self.scaler = model_data["scaler"]
                if "feature_importance" in model_data:
                    self.feature_importance[model_name] = model_data[
                        "feature_importance"
                    ]
                if "cross_val_scores" in model_data:
                    self.cross_val_scores[model_name] = model_data["cross_val_scores"]
                return True

        except Exception as e:
            self.logger.error(f"Error loading model: {str(e)}")
            return False

    def get_model_info(self, model_name: str) -> Dict:
        """Get model information.

        Args:
            model_name: Model name

        Returns:
            Model information
        """
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
                info["parameters"] = {
                    name: param.shape for name, param in model.named_parameters()
                }

            else:
                raise ValueError(f"Model {model_name} not found")

            if model_name in self.feature_importance:
                info["feature_importance"] = self.feature_importance[
                    model_name
                ].tolist()

            if model_name in self.cross_val_scores:
                info["cross_val_scores"] = self.cross_val_scores[model_name]

            return info

        except Exception as e:
            self.logger.error(f"Error getting model info: {str(e)}")
            return {}

    def preprocess_features(
        self,
        features: np.ndarray,
        fit: bool = False,
    ) -> np.ndarray:
        """Preprocess feature matrix.

        Args:
            features: Feature matrix
            fit: Whether to fit scaler

        Returns:
            Preprocessed features
        """
        try:
            if fit:
                return self.scaler.fit_transform(features)
            return self.scaler.transform(features)

        except Exception as e:
            self.logger.error(f"Error preprocessing features: {str(e)}")
            return features

    @abstractmethod
    def train(self, *args, **kwargs):
        """Train model(s)."""
        pass

    @abstractmethod
    def predict(self, *args, **kwargs):
        """Make predictions."""
        pass
