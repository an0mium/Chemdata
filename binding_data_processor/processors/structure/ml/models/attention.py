"""Attention model implementation for molecular property prediction.

This module provides:
1. Self-attention mechanism for molecular graphs
2. Multi-head attention layers
3. Position-wise feed-forward networks
4. Attention-based feature extraction
"""

from typing import Dict, List, Optional, Tuple, Union

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data, Batch

from ..core import ModelBase


class AttentionModel(ModelBase):
    """Attention-based model for molecular property prediction."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 128,
        n_heads: int = 4,
        n_layers: int = 3,
        dropout: float = 0.1,
        **kwargs,
    ):
        """Initialize attention model.

        Args:
            input_dim: Input feature dimension
            hidden_dim: Hidden layer dimension
            n_heads: Number of attention heads
            n_layers: Number of attention layers
            dropout: Dropout probability
            **kwargs: Additional arguments passed to ModelBase
        """
        super().__init__(**kwargs)

        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.n_heads = n_heads
        self.n_layers = n_layers

        # Input projection
        self.input_proj = nn.Linear(input_dim, hidden_dim)

        # Multi-head attention layers
        self.attention_layers = nn.ModuleList([MultiHeadAttention(hidden_dim, n_heads, dropout) for _ in range(n_layers)])

        # Feed-forward layers
        self.ff_layers = nn.ModuleList([PositionwiseFeedForward(hidden_dim, dropout) for _ in range(n_layers)])

        # Layer normalization
        self.layer_norms1 = nn.ModuleList([nn.LayerNorm(hidden_dim) for _ in range(n_layers)])
        self.layer_norms2 = nn.ModuleList([nn.LayerNorm(hidden_dim) for _ in range(n_layers)])

        # Output layers
        self.output_norm = nn.LayerNorm(hidden_dim)
        self.output_proj = nn.Linear(hidden_dim, 1)

        self.dropout = nn.Dropout(dropout)

    def forward(self, data: Data) -> torch.Tensor:
        """Forward pass through attention model.

        Args:
            data: Input graph data

        Returns:
            Predictions tensor
        """
        # Get node features and edge index
        x = data.x
        edge_index = data.edge_index

        # Initial projection
        h = self.input_proj(x)

        # Process through attention layers
        for i in range(self.n_layers):
            # Multi-head attention
            attended = self.attention_layers[i](h, edge_index)
            h = self.layer_norms1[i](h + self.dropout(attended))

            # Feed-forward
            ff_out = self.ff_layers[i](h)
            h = self.layer_norms2[i](h + self.dropout(ff_out))

        # Global pooling
        h = torch.mean(h, dim=0)

        # Output projection
        h = self.output_norm(h)
        out = self.output_proj(h)

        return out

    def get_attention_weights(self, data: Data) -> List[torch.Tensor]:
        """Get attention weights for interpretability.

        Args:
            data: Input graph data

        Returns:
            List of attention weight tensors
        """
        weights = []
        h = self.input_proj(data.x)

        for layer in self.attention_layers:
            _, attn_weights = layer(h, data.edge_index, return_weights=True)
            weights.append(attn_weights)

        return weights


class MultiHeadAttention(nn.Module):
    """Multi-head attention layer."""

    def __init__(self, hidden_dim: int, n_heads: int, dropout: float = 0.1):
        """Initialize multi-head attention.

        Args:
            hidden_dim: Hidden dimension
            n_heads: Number of attention heads
            dropout: Dropout probability
        """
        super().__init__()

        assert hidden_dim % n_heads == 0

        self.hidden_dim = hidden_dim
        self.n_heads = n_heads
        self.head_dim = hidden_dim // n_heads

        self.q_linear = nn.Linear(hidden_dim, hidden_dim)
        self.k_linear = nn.Linear(hidden_dim, hidden_dim)
        self.v_linear = nn.Linear(hidden_dim, hidden_dim)

        self.output_linear = nn.Linear(hidden_dim, hidden_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        return_weights: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass.

        Args:
            x: Input tensor
            edge_index: Edge indices
            return_weights: Whether to return attention weights

        Returns:
            Output tensor and optionally attention weights
        """
        batch_size, n_nodes = x.size(0), x.size(1)

        # Linear transformations
        q = self.q_linear(x).view(batch_size, n_nodes, self.n_heads, self.head_dim)
        k = self.k_linear(x).view(batch_size, n_nodes, self.n_heads, self.head_dim)
        v = self.v_linear(x).view(batch_size, n_nodes, self.n_heads, self.head_dim)

        # Compute attention scores
        scores = torch.matmul(q, k.transpose(-2, -1)) / torch.sqrt(torch.tensor(self.head_dim, dtype=torch.float))

        # Mask attention for graph structure
        mask = torch.zeros_like(scores, dtype=torch.bool)
        mask[edge_index[0], edge_index[1]] = True
        scores = scores.masked_fill(~mask, float("-inf"))

        # Attention weights
        attn_weights = F.softmax(scores, dim=-1)
        attn_weights = self.dropout(attn_weights)

        # Apply attention
        out = torch.matmul(attn_weights, v)
        out = out.view(batch_size, n_nodes, self.hidden_dim)
        out = self.output_linear(out)

        if return_weights:
            return out, attn_weights
        return out


class PositionwiseFeedForward(nn.Module):
    """Position-wise feed-forward network."""

    def __init__(self, hidden_dim: int, dropout: float = 0.1):
        """Initialize feed-forward network.

        Args:
            hidden_dim: Hidden dimension
            dropout: Dropout probability
        """
        super().__init__()

        self.ff = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass.

        Args:
            x: Input tensor

        Returns:
            Output tensor
        """
        return self.ff(x)
