"""Activity prediction for psychopharmacological compounds.

This module provides:
1. Activity type prediction (agonist, antagonist, etc.)
2. Psychopharmacological effect prediction
3. Receptor interaction prediction
4. Multi-label classification
5. Confidence scoring and uncertainty estimation
6. Web data integration
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import torch
from rdkit import Chem
from torch_geometric.data import Batch, Data

from ....models.compound.ml.ensemble import EnsembleModel
from .models import EnhancedGNN
from .models.graph_utils import mol_to_graph
from ...base import PredictorBase


logger = logging.getLogger(__name__)


class ActivityPredictor(PredictorBase):
    """Activity prediction with uncertainty estimation."""

    # Activity type categories
    ACTIVITY_TYPES = [
        "agonist",
        "antagonist",
        "partial_agonist",
        "inverse_agonist",
        "positive_modulator",
        "negative_modulator",
    ]

    # Psychopharmacological effects
    PSYCHO_EFFECTS = [
        "anxiolytic",
        "antidepressant",
        "antipsychotic",
        "sedative",
        "stimulant",
        "psychedelic",
        "dissociative",
        "entactogenic",
        "nootropic",
        "euphoriant",
    ]

    # Receptor interactions
    RECEPTOR_TYPES = [
        "5-HT2A",
        "5-HT2B",
        "5-HT2C",
        "5-HT1A",
        "D2",
        "NMDA",
        "GABA-A",
        "mu-opioid",
        "kappa-opioid",
        "sigma",
        "CB1",
        "CB2",
    ]

    # Default model configurations
    DEFAULT_GNN_CONFIG = {
        "type": "gnn",
        "input_dim": 74,  # RDKit atom features
        "hidden_dim": 256,
        "output_dim": len(ACTIVITY_TYPES) + len(PSYCHO_EFFECTS) + len(RECEPTOR_TYPES),
        "num_layers": 4,
        "heads": 8,
        "dropout": 0.2,
        "residual": True,
        "uncertainty": True,
    }

    DEFAULT_RF_CONFIG = {
        "type": "rf_classifier",
        "n_estimators": 200,
        "max_depth": 15,
        "class_weight": "balanced",
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
        """Predict compound activities.

        Args:
            compound: Input compound
            confidence_threshold: Minimum confidence threshold

        Returns:
            Dictionary containing:
            - activity_types: Predicted activity types with probabilities
            - psycho_effects: Predicted psychopharmacological effects
            - receptor_interactions: Predicted receptor interactions
            - uncertainties: Prediction uncertainties (if enabled)
            - confidence_scores: Prediction confidence scores
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
                preds, uncerts = self.model(data)
                confidences = self._calculate_confidence(preds, uncerts)
            else:
                preds = self.model(data)
                uncerts = None
                confidences = self._calculate_confidence(preds)

            # Split predictions by category
            results = self._process_predictions(preds, confidences, uncerts, confidence_threshold)

            # Get feature importance
            results["feature_importance"] = self.model.get_feature_importance(data)

            return results

        except Exception as e:
            logger.error(f"Error in activity prediction: {str(e)}")
            return {}

    def predict_batch(
        self,
        compounds: List[Union[str, Chem.Mol, Data]],
        batch_size: int = 32,
        confidence_threshold: float = 0.5,
    ) -> List[Dict]:
        """Predict activities for multiple compounds.

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

                # Process results for each compound
                for j in range(len(graphs)):
                    pred = preds[j]
                    conf = confidences[j]
                    unc = uncerts[j] if uncerts is not None else None

                    # Process predictions
                    result = self._process_predictions(pred, conf, unc, confidence_threshold)

                    # Add feature importance
                    result["feature_importance"] = {k: v[j] for k, v in importance.items()}

                    results.append(result)

            return results

        except Exception as e:
            logger.error(f"Error in batch prediction: {str(e)}")
            return []

    def _process_predictions(
        self,
        predictions: torch.Tensor,
        confidences: torch.Tensor,
        uncertainties: Optional[torch.Tensor] = None,
        confidence_threshold: float = 0.5,
    ) -> Dict:
        """Process raw predictions into structured results.

        Args:
            predictions: Raw model predictions
            confidences: Confidence scores
            uncertainties: Optional uncertainty estimates
            confidence_threshold: Minimum confidence threshold

        Returns:
            Processed prediction dictionary
        """
        # Convert tensors to numpy
        preds = predictions.cpu().numpy()
        confs = confidences.cpu().numpy()
        uncs = uncertainties.cpu().numpy() if uncertainties is not None else None

        # Get indices for each category
        n_activities = len(self.ACTIVITY_TYPES)
        n_effects = len(self.PSYCHO_EFFECTS)

        activity_slice = slice(0, n_activities)
        effects_slice = slice(n_activities, n_activities + n_effects)
        receptor_slice = slice(n_activities + n_effects, None)

        # Process activity types
        activities = []
        for i, (act_type, prob, conf) in enumerate(zip(self.ACTIVITY_TYPES, preds[activity_slice], confs[activity_slice])):
            if conf >= confidence_threshold:
                activities.append(
                    {
                        "type": act_type,
                        "probability": float(prob),
                        "confidence": float(conf),
                    }
                )
                if uncs is not None:
                    activities[-1]["uncertainty"] = float(uncs[i])

        # Process psychopharmacological effects
        effects = []
        for i, (effect, prob, conf) in enumerate(zip(self.PSYCHO_EFFECTS, preds[effects_slice], confs[effects_slice])):
            if conf >= confidence_threshold:
                effects.append(
                    {
                        "effect": effect,
                        "probability": float(prob),
                        "confidence": float(conf),
                    }
                )
                if uncs is not None:
                    effects[-1]["uncertainty"] = float(uncs[i + n_activities])

        # Process receptor interactions
        receptors = []
        for i, (receptor, prob, conf) in enumerate(zip(self.RECEPTOR_TYPES, preds[receptor_slice], confs[receptor_slice])):
            if conf >= confidence_threshold:
                receptors.append(
                    {
                        "receptor": receptor,
                        "probability": float(prob),
                        "confidence": float(conf),
                    }
                )
                if uncs is not None:
                    receptors[-1]["uncertainty"] = float(uncs[i + n_activities + n_effects])

        return {
            "activity_types": activities,
            "psycho_effects": effects,
            "receptor_interactions": receptors,
        }

    def _calculate_confidence(
        self,
        predictions: torch.Tensor,
        uncertainties: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
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
            # Use prediction probability-based confidence
            confidence = torch.sigmoid(predictions)

        return confidence

    def _load_models(self) -> None:
        """Load pre-trained models."""
        try:
            if self.use_ensemble:
                model_path = f"{self.model_dir}/activity_ensemble.pt"
            else:
                model_path = f"{self.model_dir}/activity_gnn.pt"

            self.model.load_state_dict(torch.load(model_path, map_location=self.device))
            logger.info(f"Loaded model from {model_path}")

        except Exception as e:
            logger.error(f"Error loading models: {str(e)}")
