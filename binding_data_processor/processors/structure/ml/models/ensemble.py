"""Ensemble models for improved prediction accuracy and uncertainty estimation.

This module provides:
1. Model ensembling with different architectures
2. Uncertainty estimation through ensemble variance
3. Feature importance analysis
4. Model selection and weighting
5. Ensemble diversity metrics
6. Model persistence and serialization
7. Device management and optimization
8. Support for both neural networks and classical ML models
"""

import logging
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from torch_geometric.data import Batch, Data

from .gnn import EnhancedGNN
from ..core import ModelBase
from ..utils import mol_to_graph, compute_fingerprints


class EnsembleModel(ModelBase):
    """Ensemble of multiple models for robust prediction."""

    def __init__(
        self,
        model_configs: Optional[List[Dict]] = None,
        models: Optional[List[Union[nn.Module, BaseEstimator]]] = None,
        device: Optional[str] = None,
        weights: Optional[List[float]] = None,
        uncertainty: bool = True,
        **kwargs,
    ):
        """Initialize ensemble.

        Args:
            model_configs: List of model configurations
            models: Pre-initialized models to use
            device: Device to run models on
            weights: Optional model weights
            uncertainty: Whether to estimate uncertainty
            **kwargs: Additional arguments passed to ModelBase
        """
        super().__init__(device=device or "cuda" if torch.cuda.is_available() else "cpu", **kwargs)
        self.uncertainty = uncertainty
        self.logger = logging.getLogger(self.__class__.__name__)

        # Initialize models
        self.models = []
        self.model_types = []

        if models is not None:
            # Use pre-initialized models
            self.models = models
            self.model_types = ["custom"] * len(models)
        elif model_configs is not None:
            # Initialize models from configs
            for config in model_configs:
                model_type = config.pop("type", "gnn")
                self.model_types.append(model_type)

                if model_type == "gnn":
                    model = EnhancedGNN(uncertainty=uncertainty, device=self.device, **config).to(self.device)
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
            self.weights = [1.0 / len(self.models)] * len(self.models) if self.models else []
        else:
            if len(weights) != len(self.models):
                raise ValueError("Number of weights must match number of models")
            weight_sum = sum(weights)
            self.weights = [w / weight_sum for w in weights]

    def add_model(
        self,
        model: Union[nn.Module, BaseEstimator],
        weight: float = 1.0,
        model_type: str = "custom",
    ) -> None:
        """Add model to ensemble.

        Args:
            model: Model to add
            weight: Weight for model predictions
            model_type: Type identifier for the model
        """
        self.models.append(model)
        self.model_types.append(model_type)

        # Renormalize weights
        total = sum(self.weights) + weight
        self.weights = [w / total for w in self.weights]
        self.weights.append(weight / total)

    def remove_model(self, index: int) -> None:
        """Remove model from ensemble.

        Args:
            index: Index of model to remove
        """
        if 0 <= index < len(self.models):
            self.models.pop(index)
            self.model_types.pop(index)
            self.weights.pop(index)

            # Renormalize weights
            if self.weights:
                total = sum(self.weights)
                self.weights = [w / total for w in self.weights]
        else:
            raise IndexError(f"Invalid model index: {index}")

    def forward(self, data: Union[Data, Batch, torch.Tensor, np.ndarray]) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass through ensemble.

        Args:
            data: Input data (graph, tensor, or features)

        Returns:
            Predictions and optionally uncertainties
        """
        if not self.models:
            raise RuntimeError("No models in ensemble")

        predictions = []
        uncertainties = []

        # Get predictions from each model
        for model, model_type, weight in zip(self.models, self.model_types, self.weights):
            if isinstance(model, nn.Module):
                if self.uncertainty and hasattr(model, "predict_with_uncertainty"):
                    pred, uncert = model(data)
                    predictions.append(weight * pred)
                    uncertainties.append(weight * uncert)
                else:
                    pred = model(data)
                    predictions.append(weight * pred)
            else:
                # Convert input for classical models
                if isinstance(data, (Data, Batch)):
                    features = compute_fingerprints(data)
                elif isinstance(data, torch.Tensor):
                    features = data.cpu().numpy()
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

        if self.uncertainty and uncertainties:
            # Combine uncertainties
            ensemble_uncert = torch.stack(uncertainties).sum(dim=0)
            # Add model variance uncertainty
            pred_variance = torch.var(torch.stack(predictions), dim=0)
            ensemble_uncert = ensemble_uncert + pred_variance
            return ensemble_pred, ensemble_uncert

        return ensemble_pred

    def predict_with_uncertainty(self, data: Union[Data, Batch, torch.Tensor, np.ndarray], n_samples: int = 10) -> Tuple[torch.Tensor, torch.Tensor]:
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

    def get_feature_importance(self, data: Union[Data, Batch, torch.Tensor, np.ndarray]) -> Dict[str, np.ndarray]:
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
                    importance = torch.cat([w.mean(dim=0) for w in attn_weights]).cpu().numpy()
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
                                    diff += torch.mean(torch.abs(params_i[name] - params_j[name])).item()
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

    def train_step(
        self,
        data: Union[Data, Batch, torch.Tensor],
        target: torch.Tensor,
        optimizer: torch.optim.Optimizer,
    ) -> float:
        """Training step for ensemble.

        Args:
            data: Input batch
            target: Target batch
            optimizer: Optimizer instance

        Returns:
            Loss value
        """
        optimizer.zero_grad()
        output = self(data)
        if isinstance(output, tuple):
            output = output[0]  # Use predictions only, not uncertainty
        loss = nn.functional.mse_loss(output, target)
        loss.backward()
        optimizer.step()
        return loss.item()

    def validate(
        self,
        data: Union[Data, Batch, torch.Tensor],
        target: torch.Tensor,
    ) -> Dict[str, float]:
        """Validate ensemble.

        Args:
            data: Validation inputs
            target: Validation targets

        Returns:
            Dictionary of validation metrics
        """
        with torch.no_grad():
            output = self(data)
            if isinstance(output, tuple):
                output, uncertainty = output
                metrics = {
                    "val_uncertainty": uncertainty.mean().item(),
                }
            else:
                metrics = {}

            metrics.update(
                {
                    "val_loss": nn.functional.mse_loss(output, target).item(),
                    "val_mae": nn.functional.l1_loss(output, target).item(),
                }
            )

        return metrics

    def to_device(self, device: Optional[Union[str, torch.device]] = None) -> "EnsembleModel":
        """Move ensemble to device.

        Args:
            device: Target device

        Returns:
            Self for chaining
        """
        if device is not None:
            self.device = torch.device(device)
            for i, model in enumerate(self.models):
                if isinstance(model, nn.Module):
                    self.models[i] = model.to(self.device)
        return self

    def save(self, path: str) -> bool:
        """Save ensemble model.

        Args:
            path: Save path

        Returns:
            Success status
        """
        try:
            state = {
                "weights": self.weights,
                "model_types": self.model_types,
                "uncertainty": self.uncertainty,
                "models": [],
            }
            for model in self.models:
                if isinstance(model, nn.Module):
                    state["models"].append(
                        {
                            "type": "torch",
                            "state": model.state_dict(),
                        }
                    )
                else:
                    import joblib

                    state["models"].append(
                        {
                            "type": "sklearn",
                            "state": joblib.dumps(model),
                        }
                    )
            torch.save(state, path)
            return True
        except Exception as e:
            self.logger.error(f"Error saving ensemble: {str(e)}")
            return False

    def load(self, path: str) -> bool:
        """Load ensemble model.

        Args:
            path: Load path

        Returns:
            Success status
        """
        try:
            state = torch.load(path, map_location=self.device)
            self.weights = state["weights"]
            self.model_types = state["model_types"]
            self.uncertainty = state.get("uncertainty", False)
            self.models = []

            for model_state in state["models"]:
                if model_state["type"] == "torch":
                    model = nn.Module()  # Create appropriate model instance
                    model.load_state_dict(model_state["state"])
                    model.to(self.device)
                    self.models.append(model)
                else:
                    import joblib

                    model = joblib.loads(model_state["state"])
                    self.models.append(model)
            return True
        except Exception as e:
            self.logger.error(f"Error loading ensemble: {str(e)}")
            return False
