"""Core ML model functionality.

This module provides:
1. Base model interfaces and abstract classes
2. Common model functionality and utilities
3. Model configuration and validation
4. Model persistence and serialization
5. Device and resource management
6. Training and evaluation utilities
7. Error handling and logging
"""

import abc
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import torch
import torch.nn as nn
from sklearn.base import BaseEstimator


class ModelBase(abc.ABC):
    """Base class for all ML models."""

    def __init__(
        self,
        model_dir: Optional[Union[str, Path]] = None,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
        **kwargs,
    ):
        """Initialize model.

        Args:
            model_dir: Directory for model files
            device: Device to run model on
            **kwargs: Additional model parameters
        """
        self.model_dir = Path(model_dir) if model_dir else None
        self.device = torch.device(device)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.config = kwargs
        self._model: Optional[Union[nn.Module, BaseEstimator]] = None
        self._is_trained = False

    @property
    def is_trained(self) -> bool:
        """Whether model is trained."""
        return self._is_trained

    @property
    def model(self) -> Optional[Union[nn.Module, BaseEstimator]]:
        """Get underlying model."""
        return self._model

    @abc.abstractmethod
    def forward(self, *args, **kwargs):
        """Forward pass."""
        pass

    @abc.abstractmethod
    def train_step(self, *args, **kwargs):
        """Single training step."""
        pass

    @abc.abstractmethod
    def validate(self, *args, **kwargs):
        """Validate model."""
        pass

    def save(self, path: Union[str, Path]) -> bool:
        """Save model.

        Args:
            path: Save path

        Returns:
            Success status
        """
        try:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)

            if isinstance(self._model, nn.Module):
                state = {
                    "model_state": self._model.state_dict(),
                    "config": self.config,
                    "is_trained": self._is_trained,
                }
                torch.save(state, path)
            else:
                import joblib

                joblib.dump(
                    {
                        "model": self._model,
                        "config": self.config,
                        "is_trained": self._is_trained,
                    },
                    path,
                )

            return True

        except Exception as e:
            self.logger.error(f"Error saving model: {str(e)}")
            return False

    def load(self, path: Union[str, Path]) -> bool:
        """Load model.

        Args:
            path: Load path

        Returns:
            Success status
        """
        try:
            path = Path(path)
            if not path.exists():
                raise FileNotFoundError(f"Model file not found: {path}")

            if isinstance(self._model, nn.Module):
                state = torch.load(path, map_location=self.device)
                self._model.load_state_dict(state["model_state"])
                self.config.update(state.get("config", {}))
                self._is_trained = state.get("is_trained", False)
            else:
                import joblib

                data = joblib.load(path)
                self._model = data["model"]
                self.config.update(data.get("config", {}))
                self._is_trained = data.get("is_trained", False)

            return True

        except Exception as e:
            self.logger.error(f"Error loading model: {str(e)}")
            return False

    def to_device(self, device: Optional[Union[str, torch.device]] = None) -> "ModelBase":
        """Move model to device.

        Args:
            device: Target device

        Returns:
            Self for chaining
        """
        if device is not None:
            self.device = torch.device(device)
            if isinstance(self._model, nn.Module):
                self._model.to(self.device)
        return self

    def get_config(self) -> Dict[str, Any]:
        """Get model configuration.

        Returns:
            Configuration dictionary
        """
        return self.config.copy()

    def update_config(self, config: Dict[str, Any]) -> None:
        """Update model configuration.

        Args:
            config: New configuration values
        """
        self.config.update(config)

    def get_parameter_count(self) -> int:
        """Get number of model parameters.

        Returns:
            Parameter count
        """
        if isinstance(self._model, nn.Module):
            return sum(p.numel() for p in self._model.parameters())
        return 0

    def get_device(self) -> torch.device:
        """Get current device.

        Returns:
            Current device
        """
        return self.device

    def get_state(self) -> Dict[str, Any]:
        """Get model state.

        Returns:
            State dictionary
        """
        state = {
            "config": self.config,
            "device": str(self.device),
            "is_trained": self._is_trained,
        }
        if isinstance(self._model, nn.Module):
            state["model_state"] = self._model.state_dict()
        return state

    def set_state(self, state: Dict[str, Any]) -> None:
        """Set model state.

        Args:
            state: State dictionary
        """
        self.config = state.get("config", {})
        if "device" in state:
            self.device = torch.device(state["device"])
        self._is_trained = state.get("is_trained", False)
        if isinstance(self._model, nn.Module) and "model_state" in state:
            self._model.load_state_dict(state["model_state"])


class ModelConfig:
    """Model configuration."""

    def __init__(self, **kwargs):
        """Initialize config.

        Args:
            **kwargs: Configuration parameters
        """
        self.params = kwargs

    def __getattr__(self, name: str) -> Any:
        """Get config parameter.

        Args:
            name: Parameter name

        Returns:
            Parameter value
        """
        if name not in self.params:
            raise AttributeError(f"Config parameter not found: {name}")
        return self.params[name]

    def __getitem__(self, key: str) -> Any:
        """Get config parameter.

        Args:
            key: Parameter name

        Returns:
            Parameter value
        """
        return self.params[key]

    def get(self, key: str, default: Any = None) -> Any:
        """Get config parameter with default.

        Args:
            key: Parameter name
            default: Default value

        Returns:
            Parameter value
        """
        return self.params.get(key, default)

    def update(self, **kwargs) -> None:
        """Update config parameters.

        Args:
            **kwargs: New parameters
        """
        self.params.update(kwargs)

    def to_dict(self) -> Dict:
        """Convert config to dictionary.

        Returns:
            Config dictionary
        """
        return self.params.copy()


class PredictionStats:
    """Model prediction statistics."""

    def __init__(self):
        """Initialize stats."""
        self.metrics: Dict[str, float] = {}
        self.predictions: List[Any] = []
        self.targets: List[Any] = []
        self.errors: List[str] = []

    def add_metric(self, name: str, value: float) -> None:
        """Add metric value.

        Args:
            name: Metric name
            value: Metric value
        """
        self.metrics[name] = value

    def add_prediction(self, pred: Any, target: Any) -> None:
        """Add prediction and target.

        Args:
            pred: Prediction
            target: Target value
        """
        self.predictions.append(pred)
        self.targets.append(target)

    def add_error(self, error: str) -> None:
        """Add error message.

        Args:
            error: Error message
        """
        self.errors.append(error)

    def get_metrics(self) -> Dict[str, float]:
        """Get all metrics.

        Returns:
            Dictionary of metrics
        """
        return self.metrics.copy()

    def get_predictions(self) -> List[Any]:
        """Get all predictions.

        Returns:
            List of predictions
        """
        return self.predictions.copy()

    def get_targets(self) -> List[Any]:
        """Get all targets.

        Returns:
            List of targets
        """
        return self.targets.copy()

    def get_errors(self) -> List[str]:
        """Get all errors.

        Returns:
            List of errors
        """
        return self.errors.copy()

    def clear(self) -> None:
        """Clear all stats."""
        self.metrics.clear()
        self.predictions.clear()
        self.targets.clear()
        self.errors.clear()
