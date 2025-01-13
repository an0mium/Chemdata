"""Base classes for neural network models.

This module provides base classes for neural network models with:
1. Common model architecture components
2. Training and evaluation utilities
3. Model persistence
4. Device management
5. Optimization utilities
"""

import abc
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

logger = logging.getLogger(__name__)


class BaseNeuralNet(nn.Module, abc.ABC):
    """Base class for neural network models."""

    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_dims: Optional[List[int]] = None,
        dropout: float = 0.1,
        device: Optional[str] = None,
    ):
        """Initialize neural network.

        Args:
            input_dim: Input dimension
            output_dim: Output dimension
            hidden_dims: Hidden layer dimensions
            dropout: Dropout probability
            device: Device to run model on
        """
        super().__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.hidden_dims = hidden_dims or [64, 32]
        self.dropout = dropout
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.logger = logging.getLogger(self.__class__.__name__)

        # Build model architecture
        self.layers = self._build_layers()
        self.to(self.device)

    @abc.abstractmethod
    def _build_layers(self) -> nn.ModuleList:
        """Build model layers.

        Returns:
            List of model layers
        """
        pass

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass.

        Args:
            x: Input tensor

        Returns:
            Output tensor
        """
        for layer in self.layers:
            x = layer(x)
        return x

    def save(self, path: Union[str, Path]) -> bool:
        """Save model.

        Args:
            path: Save path

        Returns:
            Success status
        """
        try:
            torch.save(
                {
                    "state_dict": self.state_dict(),
                    "input_dim": self.input_dim,
                    "output_dim": self.output_dim,
                    "hidden_dims": self.hidden_dims,
                    "dropout": self.dropout,
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
            checkpoint = torch.load(path, map_location=self.device)
            self.load_state_dict(checkpoint["state_dict"])
            self.input_dim = checkpoint["input_dim"]
            self.output_dim = checkpoint["output_dim"]
            self.hidden_dims = checkpoint["hidden_dims"]
            self.dropout = checkpoint["dropout"]
            return True
        except Exception as e:
            self.logger.error(f"Error loading model: {str(e)}")
            return False

    def train_step(
        self,
        batch: Dict[str, torch.Tensor],
        optimizer: torch.optim.Optimizer,
        criterion: nn.Module,
    ) -> Dict[str, float]:
        """Training step.

        Args:
            batch: Batch of data
            optimizer: Optimizer
            criterion: Loss criterion

        Returns:
            Dictionary of metrics
        """
        self.train()
        optimizer.zero_grad()

        x = batch["features"].to(self.device)
        y = batch["labels"].to(self.device)

        outputs = self(x)
        loss = criterion(outputs, y)
        loss.backward()
        optimizer.step()

        return {"loss": loss.item()}

    def validate_step(
        self,
        batch: Dict[str, torch.Tensor],
        criterion: nn.Module,
    ) -> Dict[str, float]:
        """Validation step.

        Args:
            batch: Batch of data
            criterion: Loss criterion

        Returns:
            Dictionary of metrics
        """
        self.eval()
        with torch.no_grad():
            x = batch["features"].to(self.device)
            y = batch["labels"].to(self.device)

            outputs = self(x)
            loss = criterion(outputs, y)

            return {"loss": loss.item()}
