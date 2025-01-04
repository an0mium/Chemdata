"""Ensemble models for activity prediction.

This module provides ensemble learning capabilities:
1. Model ensembling strategies
2. Uncertainty estimation
3. Consensus predictions
4. Cross-validation
5. Stacking and blending
"""

import logging
from typing import Dict, List, Optional, Union, Tuple
import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem

from .activity import ActivityPredictor
from .base import MLProcessor, ModelBase
from .features import MolecularFeaturizer
from .utils import mol_to_graph


class EnsemblePredictor(MLProcessor):
    """Ensemble of multiple activity prediction models."""

    def __init__(
        self,
        model_types: Optional[List[str]] = None,
        device: Optional[str] = None,
        n_models: int = 5,
        dropout_rate: float = 0.1,
    ):
        """Initialize ensemble predictor.

        Args:
            model_types: List of model types to include
            device: Device to run models on
            n_models: Number of models in ensemble
            dropout_rate: Dropout rate for uncertainty estimation
        """
        super().__init__()
        self.device = device or "cuda" if torch.cuda.is_available() else "cpu"
        self.n_models = n_models
        self.dropout_rate = dropout_rate

        # Default model types if none provided
        if model_types is None:
            model_types = ["gcn", "gat", "transformer"]
        self.model_types = model_types

        # Initialize base predictors
        self.predictors = []
        for _ in range(n_models):
            for model_type in model_types:
                predictor = ActivityPredictor(
                    device=self.device, model_type=model_type, dropout=dropout_rate
                )
                self.predictors.append(predictor)

    def predict_with_uncertainty(
        self, compound: Union[str, Chem.Mol], n_samples: int = 30
    ) -> Dict:
        """Get predictions with uncertainty estimates.

        Args:
            compound: Input compound
            n_samples: Number of MC dropout samples

        Returns:
            Dictionary containing:
            - mean predictions
            - standard deviations
            - confidence intervals
        """
        predictions = {"binding": [], "activity": [], "toxicity": []}

        # Get multiple predictions with dropout
        for _ in range(n_samples):
            for predictor in self.predictors:
                result = predictor.predict_compound(compound)
                predictions["binding"].append(result["binding_affinity"])
                predictions["activity"].append(result["activity_probs"])
                predictions["toxicity"].append(result["toxicity_probs"])

        # Calculate statistics
        results = {}
        for pred_type in predictions:
            values = np.array(predictions[pred_type])
            results[pred_type] = {
                "mean": np.mean(values, axis=0),
                "std": np.std(values, axis=0),
                "ci_lower": np.percentile(values, 2.5, axis=0),
                "ci_upper": np.percentile(values, 97.5, axis=0),
            }

        return results

    def consensus_predict(
        self, compound: Union[str, Chem.Mol], threshold: float = 0.5
    ) -> Dict:
        """Get consensus predictions across models.

        Args:
            compound: Input compound
            threshold: Consensus threshold

        Returns:
            Consensus predictions
        """
        predictions = []
        for predictor in self.predictors:
            pred = predictor.predict_compound(compound)
            predictions.append(pred)

        # Get consensus predictions
        consensus = {}

        # Binding affinity consensus
        binding_values = [p["binding_affinity"] for p in predictions]
        consensus["binding_affinity"] = float(np.median(binding_values))

        # Activity consensus
        activity_probs = np.array([p["activity_probs"] for p in predictions])
        mean_probs = np.mean(activity_probs, axis=0)
        consensus["activities"] = [
            (activity, prob)
            for activity, prob in zip(ActivityPredictor.ACTIVITY_TYPES, mean_probs)
            if prob > threshold
        ]

        # Toxicity consensus
        toxicity_probs = np.array([p["toxicity_probs"] for p in predictions])
        mean_probs = np.mean(toxicity_probs, axis=0)
        consensus["toxicities"] = [
            (tox, prob)
            for tox, prob in zip(ActivityPredictor.TOXICITY_TYPES, mean_probs)
            if prob > threshold
        ]

        return consensus

    def cross_validate(
        self,
        compounds: List[Union[str, Chem.Mol]],
        labels: np.ndarray,
        n_folds: int = 5,
    ) -> Dict:
        """Perform cross-validation.

        Args:
            compounds: List of compounds
            labels: Ground truth labels
            n_folds: Number of CV folds

        Returns:
            Cross-validation metrics
        """
        from sklearn.model_selection import KFold

        metrics = {"binding_rmse": [], "activity_auc": [], "toxicity_auc": []}

        kf = KFold(n_splits=n_folds, shuffle=True)
        for train_idx, val_idx in kf.split(compounds):
            # Train models
            train_compounds = [compounds[i] for i in train_idx]
            train_labels = labels[train_idx]

            for predictor in self.predictors:
                predictor.train(train_compounds, train_labels)

            # Evaluate
            val_compounds = [compounds[i] for i in val_idx]
            val_labels = labels[val_idx]

            val_metrics = self.evaluate(val_compounds, val_labels)
            for metric in metrics:
                metrics[metric].append(val_metrics[metric])

        # Calculate mean and std of metrics
        results = {}
        for metric in metrics:
            results[metric] = {
                "mean": float(np.mean(metrics[metric])),
                "std": float(np.std(metrics[metric])),
            }

        return results

    def evaluate(
        self, compounds: List[Union[str, Chem.Mol]], labels: np.ndarray
    ) -> Dict:
        """Evaluate models on test data.

        Args:
            compounds: Test compounds
            labels: Ground truth labels

        Returns:
            Evaluation metrics
        """
        from sklearn.metrics import roc_auc_score, mean_squared_error

        predictions = []
        for compound in compounds:
            pred = self.consensus_predict(compound)
            predictions.append(pred)

        # Calculate metrics
        binding_pred = [p["binding_affinity"] for p in predictions]
        binding_rmse = np.sqrt(mean_squared_error(labels[:, 0], binding_pred))

        activity_pred = np.array([p["activity_probs"] for p in predictions])
        activity_auc = roc_auc_score(labels[:, 1:7], activity_pred)

        toxicity_pred = np.array([p["toxicity_probs"] for p in predictions])
        toxicity_auc = roc_auc_score(labels[:, 7:], toxicity_pred)

        return {
            "binding_rmse": float(binding_rmse),
            "activity_auc": float(activity_auc),
            "toxicity_auc": float(toxicity_auc),
        }

    def calibrate(
        self, compounds: List[Union[str, Chem.Mol]], labels: np.ndarray
    ) -> None:
        """Calibrate model probabilities.

        Args:
            compounds: Calibration compounds
            labels: Ground truth labels
        """
        from sklearn.calibration import CalibratedClassifierCV

        # Get predictions
        predictions = []
        for compound in compounds:
            pred = self.predict_with_uncertainty(compound)
            predictions.append(pred)

        # Calibrate each predictor
        for predictor in self.predictors:
            predictor.calibrate(compounds, labels)

        self.logger.info("Models calibrated successfully")
