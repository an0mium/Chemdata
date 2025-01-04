"""Binding affinity prediction using ensemble models.

This module provides:
1. Binding affinity prediction for multiple targets
2. Uncertainty estimation
3. Feature importance analysis
4. Confidence scoring
5. Integration with web data sources
"""

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
from rdkit import Chem
from torch_geometric.data import Batch, Data

from ..models.ensemble import EnsembleModel
from ..models.gnn import EnhancedGNN
from ..utils import mol_to_graph, compute_fingerprints
from ...base import PredictorBase


class AffinityPredictor(PredictorBase):
    """Binding affinity prediction with uncertainty estimation."""

    # Default model configurations
    DEFAULT_GNN_CONFIG = {
        "type": "gnn",
        "input_dim": 74,  # RDKit atom features
        "hidden_dim": 128,
        "output_dim": 1,
        "num_layers": 3,
        "heads": 4,
        "dropout": 0.1,
        "residual": True,
        "uncertainty": True,
    }

    DEFAULT_RF_CONFIG = {
        "type": "rf_regressor",
        "n_estimators": 100,
        "max_depth": 10,
        "n_jobs": -1,
        "random_state": 42,
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        use_ensemble: bool = True,
        uncertainty: bool = True,
        device: Optional[str] = None,
        model_configs: Optional[List[Dict]] = None,
    ):
        """Initialize predictor.

        Args:
            model_dir: Directory containing pre-trained models
            use_ensemble: Whether to use ensemble models
            uncertainty: Whether to estimate uncertainty
            device: Device to run models on
            model_configs: Optional model configurations
        """
        super().__init__()
        self.model_dir = model_dir
        self.use_ensemble = use_ensemble
        self.uncertainty = uncertainty
        self.device = device or "cuda" if torch.cuda.is_available() else "cpu"

        # Set up model configurations
        if model_configs is None:
            if use_ensemble:
                model_configs = [
                    self.DEFAULT_GNN_CONFIG,
                    self.DEFAULT_GNN_CONFIG.copy(),  # Different random init
                    self.DEFAULT_RF_CONFIG,
                ]
            else:
                model_configs = [self.DEFAULT_GNN_CONFIG]

        # Initialize models
        if use_ensemble:
            self.model = EnsembleModel(
                model_configs=model_configs,
                device=self.device,
                uncertainty=uncertainty,
            )
        else:
            self.model = EnhancedGNN(
                **model_configs[0],
                device=self.device,
                uncertainty=uncertainty,
            ).to(self.device)

        # Load pre-trained models if available
        if model_dir:
            self._load_models()

    def predict(
        self,
        compound: Union[str, Chem.Mol, Data, Batch],
        confidence_threshold: float = 0.5,
    ) -> Dict:
        """Predict binding affinity.

        Args:
            compound: Input compound
            confidence_threshold: Minimum confidence threshold

        Returns:
            Dictionary containing:
            - binding_affinity: Predicted binding affinity (Ki in nM)
            - uncertainty: Prediction uncertainty (if enabled)
            - confidence: Prediction confidence score
            - feature_importance: Feature importance scores
        """
        try:
            # Convert input to appropriate format
            if isinstance(compound, str):
                mol = Chem.MolFromSmiles(compound)
                if mol is None:
                    raise ValueError("Invalid SMILES string")
                data = mol_to_graph(mol)
            elif isinstance(compound, Chem.Mol):
                data = mol_to_graph(compound)
            else:
                data = compound

            # Make prediction
            if self.uncertainty:
                pred, uncert = self.model(data)
                confidence = self._calculate_confidence(pred, uncert)
            else:
                pred = self.model(data)
                uncert = None
                confidence = self._calculate_confidence(pred)

            # Get feature importance
            importance = self.model.get_feature_importance(data)

            # Prepare results
            results = {
                "binding_affinity": float(pred.item()),
                "confidence": float(confidence),
                "feature_importance": importance,
            }

            if uncert is not None:
                results["uncertainty"] = float(uncert.item())

            # Filter by confidence
            if confidence < confidence_threshold:
                results["warning"] = "Low confidence prediction"

            return results

        except Exception as e:
            self.logger.error(f"Error in affinity prediction: {str(e)}")
            return {}

    def predict_batch(
        self,
        compounds: List[Union[str, Chem.Mol, Data]],
        batch_size: int = 32,
        confidence_threshold: float = 0.5,
    ) -> List[Dict]:
        """Predict binding affinities for multiple compounds.

        Args:
            compounds: List of compounds
            batch_size: Batch size for predictions
            confidence_threshold: Minimum confidence threshold

        Returns:
            List of prediction dictionaries
        """
        try:
            results = []

            # Process in batches
            for i in range(0, len(compounds), batch_size):
                batch = compounds[i : i + batch_size]

                # Convert batch to graphs
                graphs = []
                for comp in batch:
                    if isinstance(comp, str):
                        mol = Chem.MolFromSmiles(comp)
                        if mol is not None:
                            graphs.append(mol_to_graph(mol))
                    elif isinstance(comp, Chem.Mol):
                        graphs.append(mol_to_graph(comp))
                    else:
                        graphs.append(comp)

                if not graphs:
                    continue

                # Create batch
                batch_data = Batch.from_data_list(graphs)

                # Get predictions
                if self.uncertainty:
                    preds, uncerts = self.model(batch_data)
                    confidences = self._calculate_confidence(preds, uncerts)
                else:
                    preds = self.model(batch_data)
                    uncerts = None
                    confidences = self._calculate_confidence(preds)

                # Get feature importance
                importance = self.model.get_feature_importance(batch_data)

                # Process results
                for j in range(len(graphs)):
                    result = {
                        "binding_affinity": float(preds[j].item()),
                        "confidence": float(confidences[j]),
                        "feature_importance": {k: v[j] for k, v in importance.items()},
                    }

                    if uncerts is not None:
                        result["uncertainty"] = float(uncerts[j].item())

                    if result["confidence"] < confidence_threshold:
                        result["warning"] = "Low confidence prediction"

                    results.append(result)

            return results

        except Exception as e:
            self.logger.error(f"Error in batch prediction: {str(e)}")
            return []

    def _calculate_confidence(
        self,
        predictions: torch.Tensor,
        uncertainties: Optional[torch.Tensor] = None,
    ) -> Union[float, torch.Tensor]:
        """Calculate prediction confidence scores.

        Args:
            predictions: Model predictions
            uncertainties: Optional uncertainty estimates

        Returns:
            Confidence scores
        """
        if uncertainties is not None:
            # Use uncertainty-based confidence
            confidence = 1.0 / (1.0 + uncertainties)
        else:
            # Use prediction magnitude-based confidence
            confidence = torch.sigmoid(predictions)

        return confidence

    def _load_models(self) -> None:
        """Load pre-trained models."""
        try:
            if self.use_ensemble:
                model_path = f"{self.model_dir}/affinity_ensemble.pt"
            else:
                model_path = f"{self.model_dir}/affinity_gnn.pt"

            self.model.load_state_dict(torch.load(model_path, map_location=self.device))
            self.logger.info(f"Loaded model from {model_path}")

        except Exception as e:
            self.logger.error(f"Error loading models: {str(e)}")
