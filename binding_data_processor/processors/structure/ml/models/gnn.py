"""Enhanced Graph Neural Network models for molecular property prediction.

This module provides:
1. Graph attention network with residual connections
2. Multi-head attention for better feature learning
3. Batch normalization and dropout for regularization
4. Uncertainty estimation capabilities
5. Feature extraction and interpretation
6. Multiple pooling strategies
7. Monte Carlo dropout for uncertainty
"""

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Batch, Data
from torch_geometric.nn import GATConv, global_mean_pool, global_add_pool, global_max_pool

from ..neural_net import BaseNeuralNet


class EnhancedGNN(BaseNeuralNet):
    """Enhanced Graph Neural Network with attention and uncertainty."""

    def __init__(
        self,
        input_dim: int = 74,  # RDKit atom features
        hidden_dim: int = 256,
        output_dim: int = 1,
        num_layers: int = 4,
        heads: int = 8,
        dropout: float = 0.2,
        residual: bool = True,
        uncertainty: bool = True,
        readout: Union[str, List[str]] = ["mean", "max", "add"],
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
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
            readout: Readout functions to use
            device: Device to run model on
        """
        super().__init__(input_dim, output_dim, hidden_dims=[hidden_dim] * num_layers, dropout=dropout, device=device)
        self.hidden_dim = hidden_dim
        self.uncertainty = uncertainty
        self.readout = readout if isinstance(readout, list) else [readout]
        self.residual = residual
        self.num_layers = num_layers
        self.heads = heads

        # Input layer
        self.input = GATConv(
            input_dim,
            hidden_dim // heads,
            heads=heads,
            dropout=dropout,
        )

        # Hidden layers
        self.convs = nn.ModuleList()
        self.batch_norms = nn.ModuleList()

        for _ in range(num_layers - 1):
            self.convs.append(
                GATConv(
                    hidden_dim,
                    hidden_dim // heads,
                    heads=heads,
                    dropout=dropout,
                )
            )
            self.batch_norms.append(nn.BatchNorm1d(hidden_dim))

        # Calculate output dimension based on readout functions
        readout_dim = hidden_dim * len(self.readout)

        # Output layers
        self.output = nn.Sequential(
            nn.Linear(readout_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
        )

        # Prediction head
        self.pred_head = nn.Linear(hidden_dim // 2, output_dim)

        # Uncertainty head
        if uncertainty:
            self.uncert_head = nn.Linear(hidden_dim // 2, output_dim)

        self.to(device)

    def _build_layers(self) -> nn.ModuleList:
        """Build model layers.

        Returns:
            List of model layers
        """
        layers = nn.ModuleList([self.input, *self.convs, *self.batch_norms, self.output, self.pred_head])
        if self.uncertainty:
            layers.append(self.uncert_head)
        return layers

    def forward(self, data: Union[Data, Batch]) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass.

        Args:
            data: Input graph data

        Returns:
            Predictions and optionally uncertainties
        """
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = x.to(self.device)
        edge_index = edge_index.to(self.device)
        batch = batch.to(self.device)

        # Initial convolution
        h = self.input(x, edge_index)
        h_prev = h if self.residual else None

        # Hidden convolutions with batch norm and residual connections
        for i, (conv, bn) in enumerate(zip(self.convs, self.batch_norms)):
            h_new = conv(h, edge_index)
            h_new = bn(h_new)
            h_new = F.elu(h_new)
            h = h_new + h_prev if self.residual else h_new
            h_prev = h

        # Multiple readout functions
        pooled = []
        if "mean" in self.readout:
            pooled.append(global_mean_pool(h, batch))
        if "max" in self.readout:
            pooled.append(global_max_pool(h, batch))
        if "add" in self.readout:
            pooled.append(global_add_pool(h, batch))

        # Combine pooled features
        h = torch.cat(pooled, dim=1)

        # Output layers
        h = self.output(h)

        # Predictions
        pred = self.pred_head(h)

        # Return with uncertainty if enabled
        if self.uncertainty:
            log_var = self.uncert_head(h)
            return pred, torch.exp(log_var)
        return pred

    def predict_with_uncertainty(self, data: Union[Data, Batch], n_samples: int = 30) -> Tuple[torch.Tensor, torch.Tensor]:
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
        uncertainty = preds.std(dim=0)

        self.eval()  # Restore eval mode
        return mean, uncertainty

    def get_attention_weights(self, data: Union[Data, Batch]) -> List[torch.Tensor]:
        """Get attention weights for interpretation.

        Args:
            data: Input graph data

        Returns:
            List of attention weight tensors
        """
        weights = []
        x = data.x.to(self.device)
        edge_index = data.edge_index.to(self.device)

        # Get weights from each layer
        _, attn = self.input(x, edge_index, return_attention_weights=True)
        weights.append(attn)

        for conv in self.convs:
            _, attn = conv(x, edge_index, return_attention_weights=True)
            weights.append(attn)

        return weights

    def get_node_embeddings(self, data: Union[Data, Batch]) -> torch.Tensor:
        """Get node embeddings before readout.

        Args:
            data: Input graph data

        Returns:
            Node embedding tensor
        """
        x = data.x.to(self.device)
        edge_index = data.edge_index.to(self.device)

        # Get embeddings before readout
        h = self.input(x, edge_index)

        for conv, bn in zip(self.convs, self.batch_norms):
            h_new = conv(h, edge_index)
            h = bn(h_new)
            h = F.elu(h)

        return h

    def get_feature_importance(self, data: Union[Data, Batch]) -> Dict[str, np.ndarray]:
        """Calculate feature importance scores.

        Args:
            data: PyG Data object

        Returns:
            Dictionary of feature importance scores
        """
        x = data.x.to(self.device)
        edge_index = data.edge_index.to(self.device)

        importance_scores = {}

        # Get attention weights from all layers
        weights = self.get_attention_weights(data)

        # Node importance from attention
        node_importance = torch.zeros(x.size(0)).to(self.device)
        for layer_weights in weights:
            node_importance += layer_weights.mean(dim=1)
        importance_scores["node_importance"] = node_importance.cpu().detach().numpy() / len(weights)

        # Edge importance from attention
        edge_importance = torch.zeros(edge_index.size(1)).to(self.device)
        for layer_weights in weights:
            edge_importance += layer_weights
        importance_scores["edge_importance"] = edge_importance.cpu().detach().numpy() / len(weights)

        # Feature importance from input layer
        feature_importance = self.input.weight.mean(dim=0)
        importance_scores["feature_importance"] = feature_importance.cpu().detach().numpy()

        return importance_scores
