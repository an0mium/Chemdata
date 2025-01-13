"""Model registry for managing ML models.

This module provides:
1. Model registration and retrieval
2. Model configuration management
3. Model versioning
4. Model metadata tracking
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Type, Any

from .neural_net import BaseNeuralNet


class ModelRegistry:
    """Registry for managing ML models."""

    def __init__(self):
        """Initialize registry."""
        self.logger = logging.getLogger(__name__)
        self._models: Dict[str, Type[BaseNeuralNet]] = {}
        self._configs: Dict[str, Dict] = {}
        self._metadata: Dict[str, Dict] = {}

    def register(self, name: str, model_class: Type[BaseNeuralNet], config: Optional[Dict] = None, **metadata) -> None:
        """Register a model.

        Args:
            name: Model name
            model_class: Model class
            config: Default configuration
            **metadata: Additional metadata
        """
        if name in self._models:
            self.logger.warning(f"Model {name} already registered, overwriting")

        self._models[name] = model_class
        self._configs[name] = config or {}
        self._metadata[name] = metadata

        self.logger.info(f"Registered model {name}")

    def get_model(self, name: str, config: Optional[Dict] = None, **kwargs) -> BaseNeuralNet:
        """Get model instance.

        Args:
            name: Model name
            config: Model configuration
            **kwargs: Additional model arguments

        Returns:
            Model instance

        Raises:
            KeyError: If model not found
        """
        if name not in self._models:
            raise KeyError(f"Model {name} not found")

        # Get model class and merge configs
        model_class = self._models[name]
        model_config = self._configs[name].copy()
        if config:
            model_config.update(config)

        # Create and return model instance
        try:
            return model_class(config=model_config, **kwargs)
        except Exception as e:
            self.logger.error(f"Error creating model {name}: {str(e)}")
            raise

    def list_models(self) -> Dict[str, Dict[str, Any]]:
        """List registered models.

        Returns:
            Dictionary of model info
        """
        return {
            name: {
                "class": model_class.__name__,
                "config": self._configs[name],
                **self._metadata[name],
            }
            for name, model_class in self._models.items()
        }

    def get_config(self, name: str) -> Dict:
        """Get model configuration.

        Args:
            name: Model name

        Returns:
            Model configuration

        Raises:
            KeyError: If model not found
        """
        if name not in self._configs:
            raise KeyError(f"Model {name} not found")
        return self._configs[name].copy()

    def get_metadata(self, name: str) -> Dict:
        """Get model metadata.

        Args:
            name: Model name

        Returns:
            Model metadata

        Raises:
            KeyError: If model not found
        """
        if name not in self._metadata:
            raise KeyError(f"Model {name} not found")
        return self._metadata[name].copy()

    def unregister(self, name: str) -> None:
        """Unregister a model.

        Args:
            name: Model name

        Raises:
            KeyError: If model not found
        """
        if name not in self._models:
            raise KeyError(f"Model {name} not found")

        del self._models[name]
        del self._configs[name]
        del self._metadata[name]

        self.logger.info(f"Unregistered model {name}")

    def clear(self) -> None:
        """Clear all registered models."""
        self._models.clear()
        self._configs.clear()
        self._metadata.clear()
        self.logger.info("Cleared model registry")
