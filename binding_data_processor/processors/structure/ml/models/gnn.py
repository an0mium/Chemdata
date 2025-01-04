"""Enhanced Graph Neural Network models for molecular property prediction.

This module provides:
1. Graph attention network with residual connections
2. Multi-head attention for better feature learning
3. Batch normalization and dropout for regularization
4. Uncertainty estimation capabilities
5. Feature extraction and interpretation
"""

from typing import Dict, List, Optional, Tuple, Union

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Batch, Data
from torch_geometric.nn import GATConv, global_mean_pool, global_add_pool

from ..base import ModelBase


class EnhancedGNN(ModelBase):
    """Enhanced Graph Neural Network with attention and uncertainty."""

    def __init__(
        self,
        input_dim: int = 74,  # RDKit atom features
        hidden_dim: int = 128,
        output_dim: int = 1,
        num_layers: int = 3,
        heads: int = 4,
        dropout: float = 0.1,
        residual: bool = True,
        uncertainty: bool = True,
        readout: str = "mean",
    ):
        """Initialize model.

        Args:
            input_dim: Input feature dimension
            hidden_dim: Hidden layer dimension
            output_dim: Output dimension
            num_layers: Number of GNN layers
            heads: Number of attention heads
            dropout: Dropout probability
            residual: Whether to use residual connections
            uncertainty: Whether to estimate uncertainty
            readout: Readout function ("mean" or "add")
        """
        super().__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim
        self.uncertainty = uncertainty
        self.readout = readout

        # Input layer
        self.conv1 = GATConv(input_dim, hidden_dim // heads, heads=heads)

        # Hidden layers
        self.convs = nn.ModuleList(
            [
                GATConv(hidden_dim, hidden_dim // heads, heads=heads)
                for _ in range(num_layers - 2)
            ]
        )

        # Output layers
        self.conv_out = GATConv(hidden_dim, hidden_dim // heads, heads=heads)
        self.fc1 = nn.Linear(hidden_dim, hidden_dim // 2)
        self.fc2 = nn.Linear(hidden_dim // 2, output_dim)

        # Uncertainty estimation
        if uncertainty:
            self.log_var = nn.Linear(hidden_dim // 2, output_dim)

        # Additional components
        self.batch_norms = nn.ModuleList(
            [nn.BatchNorm1d(hidden_dim) for _ in range(num_layers)]
        )
        self.dropout = nn.Dropout(dropout)
        self.residual = residual

    def forward(
        self, data: Union[Data, Batch]
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass.

        Args:
            data: Input graph data

        Returns:
            Predictions and optionally uncertainties
        """
        x, edge_index, batch = data.x, data.edge_index, data.batch

        # Initial convolution
        h = self.conv1(x, edge_index)
        h = self.batch_norms[0](h)
        h = F.elu(h)
        h = self.dropout(h)
        h_prev = h if self.residual else None

        # Hidden convolutions
        for i, conv in enumerate(self.convs):
            h_new = conv(h, edge_index)
            h_new = self.batch_norms[i + 1](h_new)
            h_new = F.elu(h_new)
            h_new = self.dropout(h_new)
            h = h_new + h_prev if self.residual else h_new
            h_prev = h

        # Final convolution
        h = self.conv_out(h, edge_index)

        # Readout
        if self.readout == "mean":
            h = global_mean_pool(h, batch)
        else:
            h = global_add_pool(h, batch)

        # MLP head
        h = F.elu(self.fc1(h))
        h = self.dropout(h)

        # Predictions
        pred = self.fc2(h)

        # Return with uncertainty if enabled
        if self.uncertainty:
            log_var = self.log_var(h)
            return pred, torch.exp(log_var)
        return pred

    def get_attention_weights(self, data: Union[Data, Batch]) -> List[torch.Tensor]:
        """Get attention weights for interpretation.

        Args:
            data: Input graph data

        Returns:
            List of attention weight tensors
        """
        weights = []
        x = data.x
        edge_index = data.edge_index

        # Get weights from each layer
        _, attn = self.conv1(x, edge_index, return_attention_weights=True)
        weights.append(attn)

        for conv in self.convs:
            _, attn = conv(x, edge_index, return_attention_weights=True)
            weights.append(attn)

        _, attn = self.conv_out(x, edge_index, return_attention_weights=True)
        weights.append(attn)

        return weights

    def get_node_embeddings(self, data: Union[Data, Batch]) -> torch.Tensor:
        """Get node embeddings for visualization.

        Args:
            data: Input graph data

        Returns:
            Node embedding tensor
        """
        x, edge_index = data.x, data.edge_index

        # Get embeddings before readout
        h = self.conv1(x, edge_index)
        for conv, bn in zip(self.convs, self.batch_norms):
            h = conv(h, edge_index)
            h = bn(h)
            h = F.elu(h)
        h = self.conv_out(h, edge_index)

        return h

    def predict_with_uncertainty(
        self, data: Union[Data, Batch], n_samples: int = 10
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Make predictions with Monte Carlo dropout uncertainty.

        Args:
            data: Input graph data
            n_samples: Number of MC samples

        Returns:
            Mean predictions and uncertainties
        """
        self.train()  # Enable dropout
        preds = []

        # Multiple forward passes
        for _ in range(n_samples):
            if self.uncertainty:
                pred, _ = self(data)
            else:
                pred = self(data)
            preds.append(pred)

        # Calculate statistics
        preds = torch.stack(preds)
        mean = preds.mean(dim=0)
        std = preds.std(dim=0)

        self.eval()  # Restore eval mode
        return mean, std
