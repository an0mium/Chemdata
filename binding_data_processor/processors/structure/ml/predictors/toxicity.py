"""Enhanced toxicity prediction for psychopharmacological compounds.

This module provides comprehensive functionality for:
1. Toxicity type prediction and risk assessment
2. Severity level estimation and mechanism prediction
3. Adverse effect prediction and drug interaction analysis
4. Safety profiling and risk assessment
5. Uncertainty estimation and confidence scoring
6. Web data integration and enrichment
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import torch
import torch.nn.functional as F
from rdkit import Chem
from torch_geometric.data import Batch, Data
from sklearn.ensemble import RandomForestClassifier

from ..models.ensemble import EnsembleModel
from ..models.gnn import EnhancedGNN
from ..utils import mol_to_graph, compute_fingerprints
from ...base import PredictorBase
from .toxicity_types import (
    ALL_TOXICITY_TYPES,
    SEVERITY_LEVELS,
    RISK_LEVELS,
    ToxicityType,
)


class ToxicityPredictor(PredictorBase):
    """Enhanced toxicity prediction with uncertainty estimation."""

    # Additional toxicity categories specific to psychopharmacology
    PSYCHO_TOXICITY = {
        "serotonin_syndrome": {
            "severity_levels": ["mild", "moderate", "severe", "life_threatening"],
            "mechanisms": [
                "serotonin_excess",
                "receptor_overstimulation",
                "metabolic_inhibition",
            ],
            "related_effects": [
                "hyperthermia",
                "neuromuscular_effects",
                "autonomic_instability",
            ],
        },
        "anticholinergic_toxicity": {
            "severity_levels": ["mild", "moderate", "severe"],
            "mechanisms": [
                "receptor_antagonism",
                "cognitive_impairment",
                "peripheral_effects",
            ],
            "related_effects": [
                "confusion",
                "tachycardia",
                "urinary_retention",
            ],
        },
        "nmda_toxicity": {
            "severity_levels": ["mild", "moderate", "severe"],
            "mechanisms": [
                "receptor_antagonism",
                "glutamate_disruption",
                "neurotoxicity",
            ],
            "related_effects": [
                "dissociation",
                "cognitive_impairment",
                "neurotoxicity",
            ],
        },
    }

    # Drug interaction categories
    INTERACTION_TYPES = {
        "pharmacodynamic": [
            "serotonergic",
            "dopaminergic",
            "gabaergic",
            "glutamatergic",
            "cholinergic",
        ],
        "pharmacokinetic": [
            "cyp_inhibition",
            "cyp_induction",
            "p_glycoprotein",
            "plasma_binding",
        ],
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        use_ensemble: bool = True,
        uncertainty: bool = True,
        device: Optional[str] = None,
        model_configs: Optional[List[Dict]] = None,
        web_enrichment: bool = True,
    ):
        """Initialize predictor.

        Args:
            model_dir: Directory containing pre-trained models
            use_ensemble: Whether to use ensemble models
            uncertainty: Whether to estimate uncertainty
            device: Device to run models on
            model_configs: Optional model configurations
            web_enrichment: Whether to include web data enrichment
        """
        super().__init__()
        self.model_dir = Path(model_dir) if model_dir else None
        self.use_ensemble = use_ensemble
        self.uncertainty = uncertainty
        self.device = device or "cuda" if torch.cuda.is_available() else "cpu"
        self.web_enrichment = web_enrichment

        # Combine toxicity types
        self.toxicity_types = {
            **ALL_TOXICITY_TYPES,
            **self.PSYCHO_TOXICITY,
        }

        # Set up model configurations
        if model_configs is None:
            output_dim = (
                len(self.toxicity_types)
                + len(self.INTERACTION_TYPES["pharmacodynamic"])
                + len(self.INTERACTION_TYPES["pharmacokinetic"])
            )

            if use_ensemble:
                model_configs = [
                    {
                        "type": "gnn",
                        "input_dim": 74,  # RDKit atom features
                        "hidden_dim": 256,
                        "output_dim": output_dim,
                        "num_layers": 4,
                        "heads": 8,
                        "dropout": 0.2,
                        "residual": True,
                    },
                    {
                        "type": "gnn",
                        "input_dim": 74,
                        "hidden_dim": 256,
                        "output_dim": output_dim,
                        "num_layers": 4,
                        "heads": 8,
                        "dropout": 0.2,
                        "residual": True,
                    },
                    {
                        "type": "rf_classifier",
                        "n_estimators": 200,
                        "max_depth": 15,
                        "class_weight": "balanced",
                        "n_jobs": -1,
                        "random_state": 42,
                    },
                ]
            else:
                model_configs = [
                    {
                        "type": "gnn",
                        "input_dim": 74,
                        "hidden_dim": 256,
                        "output_dim": output_dim,
                        "num_layers": 4,
                        "heads": 8,
                        "dropout": 0.2,
                        "residual": True,
                    }
                ]

        # Initialize models
        if use_ensemble:
            self.model = EnsembleModel(
                model_configs=model_configs,
                device=self.device,
                uncertainty=uncertainty,
            )
        else:
            self.model = EnhancedGNN(**model_configs[0]).to(self.device)

        # Initialize web enrichment if enabled
        if web_enrichment:
            from web_enrichment import WebDataEnricher

            self.web_enricher = WebDataEnricher()

        # Load pre-trained models if available
        if model_dir:
            self._load_models()

    def predict(
        self,
        compound: Union[str, Chem.Mol, Data, Batch],
        confidence_threshold: float = 0.5,
        include_web_data: bool = True,
    ) -> Dict:
        """Predict compound toxicity.

        Args:
            compound: Input compound
            confidence_threshold: Minimum confidence threshold
            include_web_data: Whether to include web data

        Returns:
            Dictionary containing:
            - toxicity_types: Predicted toxicity types with probabilities
            - severity_levels: Predicted severity levels
            - mechanisms: Predicted toxicity mechanisms
            - drug_interactions: Predicted drug interactions
            - risks: Risk assessments
            - uncertainties: Prediction uncertainties (if enabled)
            - confidence_scores: Prediction confidence scores
            - web_data: Enriched web data (if enabled)
        """
        try:
            # Convert input to appropriate format
            if isinstance(compound, str):
                mol = Chem.MolFromSmiles(compound)
                if mol is None:
                    raise ValueError("Invalid SMILES string")
                data = mol_to_graph(mol)
                smiles = compound
            elif isinstance(compound, Chem.Mol):
                data = mol_to_graph(compound)
                smiles = Chem.MolToSmiles(compound)
            else:
                data = compound
                smiles = None

            # Make prediction
            if self.uncertainty:
                preds, uncerts = self.model(data)
                confidences = self._calculate_confidence(preds, uncerts)
            else:
                preds = self.model(data)
                uncerts = None
                confidences = self._calculate_confidence(preds)

            # Process predictions
            results = self._process_predictions(
                preds, confidences, uncerts, confidence_threshold
            )

            # Add web data if enabled and SMILES available
            if self.web_enrichment and include_web_data and smiles:
                web_data = self.web_enricher.enrich_compound(smiles)
                results["web_data"] = web_data

                # Update predictions with web data
                self._update_predictions_with_web_data(results, web_data)

            return results

        except Exception as e:
            self.logger.error(f"Error in toxicity prediction: {str(e)}")
            return {}

    def _process_predictions(
        self,
        predictions: torch.Tensor,
        confidences: torch.Tensor,
        uncertainties: Optional[torch.Tensor] = None,
        confidence_threshold: float = 0.5,
    ) -> Dict:
        """Process raw predictions into structured results."""
        # Convert tensors to numpy
        preds = predictions.cpu().numpy()
        confs = confidences.cpu().numpy()
        uncs = uncertainties.cpu().numpy() if uncertainties is not None else None

        results = {
            "toxicity_types": [],
            "severity_levels": {},
            "mechanisms": {},
            "risks": {},
            "drug_interactions": {
                "pharmacodynamic": [],
                "pharmacokinetic": [],
            },
        }

        # Process toxicity types
        n_tox = len(self.toxicity_types)
        n_pd = len(self.INTERACTION_TYPES["pharmacodynamic"])
        n_pk = len(self.INTERACTION_TYPES["pharmacokinetic"])

        # Toxicity predictions
        for i, (tox_id, tox_type) in enumerate(self.toxicity_types.items()):
            if confs[i] >= confidence_threshold:
                tox_result = self._create_toxicity_result(
                    tox_id,
                    tox_type,
                    preds[i],
                    confs[i],
                    uncs[i] if uncs is not None else None,
                )
                results["toxicity_types"].append(tox_result)

                # Update related data
                self._update_severity_levels(results, tox_id, tox_result)
                self._update_mechanisms(results, tox_id, tox_type)
                self._update_risks(results, tox_id, preds[i], confs[i])

        # Pharmacodynamic interactions
        pd_slice = slice(n_tox, n_tox + n_pd)
        for i, interaction in enumerate(self.INTERACTION_TYPES["pharmacodynamic"]):
            if confs[pd_slice][i] >= confidence_threshold:
                results["drug_interactions"]["pharmacodynamic"].append(
                    self._create_interaction_result(
                        interaction,
                        preds[pd_slice][i],
                        confs[pd_slice][i],
                        uncs[pd_slice][i] if uncs is not None else None,
                    )
                )

        # Pharmacokinetic interactions
        pk_slice = slice(n_tox + n_pd, n_tox + n_pd + n_pk)
        for i, interaction in enumerate(self.INTERACTION_TYPES["pharmacokinetic"]):
            if confs[pk_slice][i] >= confidence_threshold:
                results["drug_interactions"]["pharmacokinetic"].append(
                    self._create_interaction_result(
                        interaction,
                        preds[pk_slice][i],
                        confs[pk_slice][i],
                        uncs[pk_slice][i] if uncs is not None else None,
                    )
                )

        return results

    def _create_toxicity_result(
        self,
        tox_id: str,
        tox_type: Dict,
        pred: float,
        conf: float,
        unc: Optional[float] = None,
    ) -> Dict:
        """Create toxicity prediction result."""
        result = {
            "id": tox_id,
            "name": tox_type.get("name", tox_id),
            "probability": float(pred),
            "confidence": float(conf),
        }

        if unc is not None:
            result["uncertainty"] = float(unc)

        # Add severity level
        severity_levels = tox_type.get("severity_levels", [])
        if severity_levels:
            severity_idx = int(pred * (len(severity_levels) - 1))
            result["severity"] = severity_levels[severity_idx]

        # Add mechanisms
        mechanisms = tox_type.get("mechanisms", [])
        if mechanisms:
            result["mechanisms"] = mechanisms

        # Add related effects
        related_effects = tox_type.get("related_effects", [])
        if related_effects:
            result["related_effects"] = related_effects

        return result

    def _create_interaction_result(
        self,
        interaction: str,
        pred: float,
        conf: float,
        unc: Optional[float] = None,
    ) -> Dict:
        """Create interaction prediction result."""
        result = {
            "type": interaction,
            "probability": float(pred),
            "confidence": float(conf),
            "severity": self._calculate_severity(pred),
        }

        if unc is not None:
            result["uncertainty"] = float(unc)

        return result

    def _update_severity_levels(
        self, results: Dict, tox_id: str, tox_result: Dict
    ) -> None:
        """Update severity levels in results."""
        if "severity" in tox_result:
            results["severity_levels"][tox_id] = {
                "level": tox_result["severity"],
                "description": SEVERITY_LEVELS.get(tox_result["severity"], "Unknown"),
            }

    def _update_mechanisms(self, results: Dict, tox_id: str, tox_type: Dict) -> None:
        """Update mechanisms in results."""
        mechanisms = tox_type.get("mechanisms", [])
        if mechanisms:
            results["mechanisms"][tox_id] = mechanisms

    def _update_risks(
        self, results: Dict, tox_id: str, pred: float, conf: float
    ) -> None:
        """Update risk assessments in results."""
        risk_score = float(pred * conf)
        if risk_score < 0.2:
            risk_level = "very_low"
        elif risk_score < 0.4:
            risk_level = "low"
        elif risk_score < 0.6:
            risk_level = "moderate"
        elif risk_score < 0.8:
            risk_level = "high"
        else:
            risk_level = "very_high"

        results["risks"][tox_id] = {
            "level": risk_level,
            "score": risk_score,
            "description": RISK_LEVELS[risk_level],
        }

    def _update_predictions_with_web_data(self, results: Dict, web_data: Dict) -> None:
        """Update predictions with web data."""
        if not web_data:
            return

        # Update toxicity information
        if "toxicity" in web_data:
            for tox_type in results["toxicity_types"]:
                web_tox = web_data["toxicity"].get(tox_type["id"], {})
                if web_tox:
                    tox_type["web_references"] = web_tox.get("references", [])
                    tox_type["reported_effects"] = web_tox.get("effects", [])

        # Update drug interactions
        if "drug_interactions" in web_data:
            web_interactions = web_data["drug_interactions"]
            for interaction_type in ["pharmacodynamic", "pharmacokinetic"]:
                for interaction in results["drug_interactions"][interaction_type]:
                    web_int = web_interactions.get(interaction["type"], {})
                    if web_int:
                        interaction["web_references"] = web_int.get("references", [])
                        interaction["reported_cases"] = web_int.get("cases", [])

    def _calculate_confidence(
        self,
        predictions: torch.Tensor,
        uncertainties: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Calculate prediction confidence scores."""
        if uncertainties is not None:
            # Use uncertainty-based confidence
            confidence = 1.0 / (1.0 + uncertainties)
        else:
            # Use prediction probability-based confidence
            confidence = torch.sigmoid(predictions)

        return confidence

    def _calculate_severity(self, probability: float) -> str:
        """Calculate severity level based on probability."""
        if probability >= 0.8:
            return "high"
        elif probability >= 0.5:
            return "moderate"
        else:
            return "low"

    def _load_models(self) -> None:
        """Load pre-trained models."""
        try:
            if self.use_ensemble:
                model_path = f"{self.model_dir}/toxicity_ensemble.pt"
            else:
                model_path = f"{self.model_dir}/toxicity_gnn.pt"

            self.model.load_state_dict(torch.load(model_path, map_location=self.device))
            self.logger.info(f"Loaded model from {model_path}")

        except Exception as e:
            self.logger.error(f"Error loading models: {str(e)}")
