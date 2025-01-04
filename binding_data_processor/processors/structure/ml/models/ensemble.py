"""Ensemble models for improved prediction accuracy and uncertainty estimation.

This module provides:
1. Model ensembling with different architectures
2. Uncertainty estimation through ensemble variance
3. Feature importance analysis
4. Model selection and weighting
5. Ensemble diversity metrics
"""

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from torch_geometric.data import Batch, Data

from .gnn import EnhancedGNN
from ..base import ModelBase
from ..utils import mol_to_graph, compute_fingerprints


class EnsembleModel(ModelBase):
    """Ensemble of multiple models for robust prediction."""

    def __init__(
        self,
        model_configs: List[Dict],
        device: Optional[str] = None,
        weights: Optional[List[float]] = None,
        uncertainty: bool = True,
    ):
        """Initialize ensemble.

        Args:
            model_configs: List of model configurations
            device: Device to run models on
            weights: Optional model weights
            uncertainty: Whether to estimate uncertainty
        """
        super().__init__()
        self.device = device or "cuda" if torch.cuda.is_available() else "cpu"
        self.uncertainty = uncertainty

        # Initialize models
        self.models = []
        self.model_types = []
        for config in model_configs:
            model_type = config.pop("type", "gnn")
            self.model_types.append(model_type)

            if model_type == "gnn":
                model = EnhancedGNN(
                    uncertainty=uncertainty, device=self.device, **config
                ).to(self.device)
                self.models.append(model)
            elif model_type == "rf_regressor":
                model = RandomForestRegressor(**config)
                self.models.append(model)
            elif model_type == "rf_classifier":
                model = RandomForestClassifier(**config)
                self.models.append(model)
            else:
                raise ValueError(f"Unknown model type: {model_type}")

        # Set model weights
        if weights is None:
            self.weights = [1.0 / len(self.models)] * len(self.models)
        else:
            if len(weights) != len(self.models):
                raise ValueError("Number of weights must match number of models")
            weight_sum = sum(weights)
            self.weights = [w / weight_sum for w in weights]

    def forward(
        self, data: Union[Data, Batch, np.ndarray]
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass through ensemble.

        Args:
            data: Input data (graph or features)

        Returns:
            Predictions and optionally uncertainties
        """
        predictions = []
        uncertainties = []

        # Get predictions from each model
        for model, model_type, weight in zip(
            self.models, self.model_types, self.weights
        ):
            if model_type == "gnn":
                if self.uncertainty:
                    pred, uncert = model(data)
                    predictions.append(weight * pred)
                    uncertainties.append(weight * uncert)
                else:
                    pred = model(data)
                    predictions.append(weight * pred)
            else:
                # Convert graph to features for classical models
                if isinstance(data, (Data, Batch)):
                    features = compute_fingerprints(data)
                else:
                    features = data

                if hasattr(model, "predict_proba"):
                    pred = torch.tensor(
                        model.predict_proba([features])[0],
                        device=self.device,
                    )
                else:
                    pred = torch.tensor(
                        model.predict([features]),
                        device=self.device,
                    )
                predictions.append(weight * pred)

        # Combine predictions
        ensemble_pred = torch.stack(predictions).sum(dim=0)

        if self.uncertainty:
            # Combine uncertainties
            ensemble_uncert = torch.stack(uncertainties).sum(dim=0)
            # Add model variance uncertainty
            pred_variance = torch.var(torch.stack(predictions), dim=0)
            ensemble_uncert = ensemble_uncert + pred_variance
            return ensemble_pred, ensemble_uncert

        return ensemble_pred

    def predict_with_uncertainty(
        self, data: Union[Data, Batch, np.ndarray], n_samples: int = 10
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Make predictions with uncertainty estimation.

        Args:
            data: Input data
            n_samples: Number of MC samples

        Returns:
            Mean predictions and uncertainties
        """
        predictions = []

        # Multiple forward passes
        for _ in range(n_samples):
            if self.uncertainty:
                pred, _ = self(data)
            else:
                pred = self(data)
            predictions.append(pred)

        # Calculate statistics
        predictions = torch.stack(predictions)
        mean = predictions.mean(dim=0)
        std = predictions.std(dim=0)

        return mean, std

    def get_feature_importance(
        self, data: Union[Data, Batch, np.ndarray]
    ) -> Dict[str, np.ndarray]:
        """Get feature importance scores.

        Args:
            data: Input data

        Returns:
            Dictionary of importance scores per model
        """
        importance_scores = {}

        for model, model_type in zip(self.models, self.model_types):
            if model_type == "gnn":
                # Use attention weights for GNN
                if isinstance(data, (Data, Batch)):
                    attn_weights = model.get_attention_weights(data)
                    importance = (
                        torch.cat([w.mean(dim=0) for w in attn_weights]).cpu().numpy()
                    )
                else:
                    continue
            else:
                # Use feature importance for classical models
                if hasattr(model, "feature_importances_"):
                    importance = model.feature_importances_
                else:
                    continue

            importance_scores[f"{model_type}_importance"] = importance

        return importance_scores

    def get_ensemble_diversity(self) -> Dict[str, float]:
        """Calculate diversity metrics for ensemble.

        Returns:
            Dictionary of diversity metrics
        """
        metrics = {}

        # Pairwise disagreement
        n_models = len(self.models)
        if n_models > 1:
            disagreement = 0
            count = 0
            for i in range(n_models):
                for j in range(i + 1, n_models):
                    if self.model_types[i] == self.model_types[j]:
                        # Compare model parameters/weights
                        if isinstance(self.models[i], nn.Module):
                            params_i = dict(self.models[i].named_parameters())
                            params_j = dict(self.models[j].named_parameters())
                            diff = 0
                            for name in params_i:
                                if name in params_j:
                                    diff += torch.mean(
                                        torch.abs(params_i[name] - params_j[name])
                                    ).item()
                            disagreement += diff
                            count += 1

            if count > 0:
                metrics["param_diversity"] = disagreement / count

        # Architecture diversity
        type_counts = {}
        for model_type in self.model_types:
            type_counts[model_type] = type_counts.get(model_type, 0) + 1
        metrics["architecture_diversity"] = len(type_counts) / len(self.models)

        return metrics
