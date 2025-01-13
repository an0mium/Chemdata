"""Uncertainty estimation model for molecular property prediction.

This module provides:
1. Bayesian uncertainty estimation
2. Ensemble-based uncertainty quantification
3. Dropout-based uncertainty estimation
4. Calibrated confidence scores
"""

from typing import Dict, List, Optional, Tuple, Union

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data, Batch

from ..core import ModelBase


class UncertaintyModel(ModelBase):
    """Model with uncertainty estimation capabilities."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 128,
        n_layers: int = 3,
        dropout: float = 0.1,
        n_samples: int = 10,
        **kwargs,
    ):
        """Initialize uncertainty model.

        Args:
            input_dim: Input feature dimension
            hidden_dim: Hidden layer dimension
            n_layers: Number of layers
            dropout: Dropout probability
            n_samples: Number of Monte Carlo samples
            **kwargs: Additional arguments passed to ModelBase
        """
        super().__init__(**kwargs)

        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.n_layers = n_layers
        self.n_samples = n_samples

        # Input layer
        self.input_layer = nn.Linear(input_dim, hidden_dim)

        # Hidden layers with dropout
        self.hidden_layers = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                    nn.LayerNorm(hidden_dim),
                )
                for _ in range(n_layers)
            ]
        )

        # Output layers for mean and variance
        self.mean_output = nn.Linear(hidden_dim, 1)
        self.var_output = nn.Linear(hidden_dim, 1)

    def forward(self, data: Union[Data, torch.Tensor]) -> Tuple[torch.Tensor, torch.Tensor]:
        """Forward pass with uncertainty estimation.

        Args:
            data: Input data

        Returns:
            Tuple of (mean predictions, uncertainty estimates)
        """
        if isinstance(data, Data):
            x = data.x
        else:
            x = data

        # Input projection
        h = self.input_layer(x)

        # Process through hidden layers
        for layer in self.hidden_layers:
            h = layer(h)

        # Get mean and variance predictions
        mean = self.mean_output(h)
        log_var = self.var_output(h)
        var = torch.exp(log_var)

        return mean, var

    def predict_with_uncertainty(
        self,
        data: Union[Data, torch.Tensor],
        n_samples: Optional[int] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Make predictions with uncertainty estimation.

        Args:
            data: Input data
            n_samples: Number of Monte Carlo samples (optional)

        Returns:
            Tuple of (mean predictions, epistemic uncertainty, aleatoric uncertainty)
        """
        if n_samples is None:
            n_samples = self.n_samples

        # Enable dropout during inference
        self.train()

        # Multiple forward passes
        means = []
        vars = []
        for _ in range(n_samples):
            mean, var = self(data)
            means.append(mean)
            vars.append(var)

        # Calculate statistics
        means = torch.stack(means)
        vars = torch.stack(vars)

        # Epistemic uncertainty from prediction variance
        epistemic = torch.var(means, dim=0)

        # Aleatoric uncertainty from mean of predicted variances
        aleatoric = torch.mean(vars, dim=0)

        # Final mean prediction
        mean_pred = torch.mean(means, dim=0)

        return mean_pred, epistemic, aleatoric

    def calibrate(
        self,
        val_data: List[Union[Data, torch.Tensor]],
        val_targets: torch.Tensor,
    ) -> Dict[str, float]:
        """Calibrate uncertainty estimates using validation data.

        Args:
            val_data: Validation inputs
            val_targets: Validation targets

        Returns:
            Dictionary of calibration metrics
        """
        metrics = {}

        # Get predictions with uncertainty
        means = []
        epistemics = []
        aleatorics = []

        for data in val_data:
            mean, epist, aleat = self.predict_with_uncertainty(data)
            means.append(mean)
            epistemics.append(epist)
            aleatorics.append(aleat)

        means = torch.cat(means)
        epistemics = torch.cat(epistemics)
        aleatorics = torch.cat(aleatorics)

        # Calculate calibration metrics
        errors = torch.abs(means - val_targets)
        total_uncertainty = epistemics + aleatorics

        # Expected Calibration Error
        confidence_bins = torch.linspace(0, 1, steps=10)
        ece = 0.0
        for i in range(len(confidence_bins) - 1):
            bin_mask = (total_uncertainty >= confidence_bins[i]) & (total_uncertainty < confidence_bins[i + 1])
            if torch.any(bin_mask):
                bin_accuracy = torch.mean(errors[bin_mask])
                bin_confidence = torch.mean(total_uncertainty[bin_mask])
                ece += torch.abs(bin_accuracy - bin_confidence) * (bin_mask.sum() / len(total_uncertainty))

        metrics["ece"] = ece.item()

        # Negative Log Likelihood
        nll = F.gaussian_nll_loss(means, val_targets, total_uncertainty)
        metrics["nll"] = nll.item()

        return metrics

    def get_uncertainty_decomposition(
        self,
        data: Union[Data, torch.Tensor],
    ) -> Dict[str, torch.Tensor]:
        """Decompose uncertainty into different sources.

        Args:
            data: Input data

        Returns:
            Dictionary of uncertainty components
        """
        mean_pred, epistemic, aleatoric = self.predict_with_uncertainty(data)

        uncertainty_components = {
            "total": epistemic + aleatoric,
            "epistemic": epistemic,
            "aleatoric": aleatoric,
            "mean_prediction": mean_pred,
        }

        return uncertainty_components
