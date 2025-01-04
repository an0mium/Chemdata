"""Enhanced abuse potential prediction for psychopharmacological compounds.

This module provides comprehensive functionality for:
1. Recreational abuse potential prediction
2. Addiction risk assessment 
3. Dependence liability prediction
4. Withdrawal risk estimation
5. Tolerance development prediction
6. Reward pathway analysis
7. Uncertainty estimation and confidence scoring
8. Web data integration and enrichment
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
from .abuse_types import (
    ABUSE_CATEGORIES,
    ABUSE_EFFECTS,
    RECEPTOR_SYSTEMS,
    RISK_LEVELS,
    WITHDRAWAL_SEVERITY,
    get_abuse_category_info,
    get_receptor_info,
    get_risk_level_info,
    get_withdrawal_info,
)


class AbusePotentialPredictor(PredictorBase):
    """Enhanced abuse potential prediction with uncertainty estimation."""

    # Default model configurations
    DEFAULT_GNN_CONFIG = {
        "type": "gnn",
        "input_dim": 74,  # RDKit atom features
        "hidden_dim": 256,
        "output_dim": (
            len(ABUSE_CATEGORIES)
            + len(ABUSE_EFFECTS)
            + sum(len(receptors) for receptors in RECEPTOR_SYSTEMS.values())
        ),
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
        """Predict compound abuse potential.

        Args:
            compound: Input compound
            confidence_threshold: Minimum confidence threshold
            include_web_data: Whether to include web data

        Returns:
            Dictionary containing:
            - abuse_potential: Predicted abuse categories with probabilities
            - abuse_effects: Predicted abuse effects
            - receptor_affinities: Predicted receptor system interactions
            - mechanisms: Predicted abuse mechanisms
            - risk_levels: Risk level assessments
            - withdrawal: Withdrawal risk assessment
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
            self.logger.error(f"Error in abuse potential prediction: {str(e)}")
            return {}

    def predict_batch(
        self,
        compounds: List[Union[str, Chem.Mol, Data]],
        batch_size: int = 32,
        confidence_threshold: float = 0.5,
        include_web_data: bool = True,
    ) -> List[Dict]:
        """Predict abuse potential for multiple compounds.

        Args:
            compounds: List of compounds
            batch_size: Batch size for predictions
            confidence_threshold: Minimum confidence threshold
            include_web_data: Whether to include web data

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
                smiles_list = []
                for comp in batch:
                    if isinstance(comp, str):
                        mol = Chem.MolFromSmiles(comp)
                        if mol is not None:
                            graphs.append(mol_to_graph(mol))
                            smiles_list.append(comp)
                    elif isinstance(comp, Chem.Mol):
                        graphs.append(mol_to_graph(comp))
                        smiles_list.append(Chem.MolToSmiles(comp))
                    else:
                        graphs.append(comp)
                        smiles_list.append(None)

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

                # Process results for each compound
                for j, smiles in enumerate(smiles_list):
                    pred = preds[j]
                    conf = confidences[j]
                    unc = uncerts[j] if uncerts is not None else None

                    # Process predictions
                    result = self._process_predictions(
                        pred, conf, unc, confidence_threshold
                    )

                    # Add web data if enabled
                    if self.web_enrichment and include_web_data and smiles:
                        web_data = self.web_enricher.enrich_compound(smiles)
                        result["web_data"] = web_data
                        self._update_predictions_with_web_data(result, web_data)

                    results.append(result)

            return results

        except Exception as e:
            self.logger.error(f"Error in batch prediction: {str(e)}")
            return []

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
            "abuse_potential": [],
            "abuse_effects": [],
            "receptor_affinities": [],
            "mechanisms": {},
            "risk_levels": {},
            "withdrawal": {},
        }

        # Get indices for each category
        n_categories = len(ABUSE_CATEGORIES)
        n_effects = len(ABUSE_EFFECTS)
        n_receptors = sum(len(receptors) for receptors in RECEPTOR_SYSTEMS.values())

        # Process abuse categories
        for i, (category, info) in enumerate(ABUSE_CATEGORIES.items()):
            if confs[i] >= confidence_threshold:
                category_result = self._create_category_result(
                    category,
                    info,
                    preds[i],
                    confs[i],
                    uncs[i] if uncs is not None else None,
                )
                results["abuse_potential"].append(category_result)

                # Update mechanisms and risk levels
                self._update_mechanisms(results, category, info, preds[i])
                self._update_risk_levels(results, category, category_result)

                # Special handling for withdrawal risk
                if category == "withdrawal_risk":
                    self._update_withdrawal_info(results, category_result)

        # Process abuse effects
        effects_slice = slice(n_categories, n_categories + n_effects)
        for i, (effect, subtypes) in enumerate(ABUSE_EFFECTS.items()):
            if confs[effects_slice][i] >= confidence_threshold:
                results["abuse_effects"].append(
                    self._create_effect_result(
                        effect,
                        subtypes,
                        preds[effects_slice][i],
                        confs[effects_slice][i],
                        uncs[effects_slice][i] if uncs is not None else None,
                    )
                )

        # Process receptor affinities
        receptor_slice = slice(n_categories + n_effects, None)
        receptor_idx = 0
        for system, receptors in RECEPTOR_SYSTEMS.items():
            for receptor in receptors:
                if confs[receptor_slice][receptor_idx] >= confidence_threshold:
                    results["receptor_affinities"].append(
                        self._create_receptor_result(
                            system,
                            receptor,
                            preds[receptor_slice][receptor_idx],
                            confs[receptor_slice][receptor_idx],
                            (
                                uncs[receptor_slice][receptor_idx]
                                if uncs is not None
                                else None
                            ),
                        )
                    )
                receptor_idx += 1

        return results

    def _create_category_result(
        self,
        category: str,
        info: Dict,
        pred: float,
        conf: float,
        unc: Optional[float] = None,
    ) -> Dict:
        """Create abuse category prediction result."""
        result = {
            "category": category,
            "probability": float(pred),
            "confidence": float(conf),
            "description": info["description"],
        }

        if unc is not None:
            result["uncertainty"] = float(unc)

        # Add risk level
        risk_levels = info["risk_levels"]
        risk_idx = int(pred * (len(risk_levels) - 1))
        result["risk_level"] = risk_levels[risk_idx]
        result["risk_description"] = get_risk_level_info(risk_levels[risk_idx])

        # Add mechanisms
        result["mechanisms"] = info["mechanisms"]

        return result

    def _create_effect_result(
        self,
        effect: str,
        subtypes: List[str],
        pred: float,
        conf: float,
        unc: Optional[float] = None,
    ) -> Dict:
        """Create abuse effect prediction result."""
        result = {
            "effect": effect,
            "subtypes": subtypes,
            "probability": float(pred),
            "confidence": float(conf),
            "severity": self._calculate_severity(pred),
        }

        if unc is not None:
            result["uncertainty"] = float(unc)

        return result

    def _create_receptor_result(
        self,
        system: str,
        receptor: str,
        pred: float,
        conf: float,
        unc: Optional[float] = None,
    ) -> Dict:
        """Create receptor affinity prediction result."""
        result = {
            "system": system,
            "receptor": receptor,
            "affinity": float(pred),
            "confidence": float(conf),
            "receptor_info": get_receptor_info(receptor),
        }

        if unc is not None:
            result["uncertainty"] = float(unc)

        return result

    def _update_mechanisms(
        self, results: Dict, category: str, info: Dict, pred: float
    ) -> None:
        """Update abuse mechanisms in results."""
        mechanisms = info.get("mechanisms", [])
        if mechanisms:
            results["mechanisms"][category] = {
                "mechanisms": mechanisms,
                "contribution": float(pred),
            }

    def _update_risk_levels(
        self, results: Dict, category: str, category_result: Dict
    ) -> None:
        """Update risk levels in results."""
        if "risk_level" in category_result:
            results["risk_levels"][category] = {
                "level": category_result["risk_level"],
                "description": category_result["risk_description"],
                "probability": category_result["probability"],
                "confidence": category_result["confidence"],
            }

    def _update_withdrawal_info(self, results: Dict, category_result: Dict) -> None:
        """Update withdrawal information in results."""
        if category_result["risk_level"] in WITHDRAWAL_SEVERITY:
            severity_info = get_withdrawal_info(category_result["risk_level"])
            results["withdrawal"] = {
                "severity": category_result["risk_level"],
                "probability": category_result["probability"],
                "confidence": category_result["confidence"],
                "duration": severity_info["duration"],
                "symptoms": severity_info["symptoms"],
                "management": severity_info["management"],
            }

    def _update_predictions_with_web_data(self, results: Dict, web_data: Dict) -> None:
        """Update predictions with web data."""
        if not web_data:
            return

        # Update abuse potential information
        if "abuse_potential" in web_data:
            for category in results["abuse_potential"]:
                web_abuse = web_data["abuse_potential"].get(category["category"], {})
                if web_abuse:
                    category["web_references"] = web_abuse.get("references", [])
                    category["reported_cases"] = web_abuse.get("cases", [])
                    category["legal_status"] = web_abuse.get("legal_status", {})

        # Update receptor information
        if "receptor_data" in web_data:
            web_receptors = web_data["receptor_data"]
            for receptor in results["receptor_affinities"]:
                web_rec = web_receptors.get(receptor["receptor"], {})
                if web_rec:
                    receptor["literature_data"] = web_rec.get("literature", [])
                    receptor["binding_data"] = web_rec.get("binding", {})

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
                model_path = f"{self.model_dir}/abuse_ensemble.pt"
            else:
                model_path = f"{self.model_dir}/abuse_gnn.pt"

            self.model.load_state_dict(torch.load(model_path, map_location=self.device))
            self.logger.info(f"Loaded model from {model_path}")

        except Exception as e:
            self.logger.error(f"Error loading models: {str(e)}")
