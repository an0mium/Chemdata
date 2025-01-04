"""Enhanced ML integration for chemical compound analysis.

This module provides:
1. Model loading and management
2. Prediction generation for multiple endpoints
3. Feature extraction and processing
4. Model validation and monitoring
5. Ensemble methods and transfer learning
6. Supporting evidence collection
7. Performance tracking
8. Cache management
"""

import logging
from typing import Dict, List, Optional, Tuple, Any, Union
import numpy as np
import pandas as pd
from pathlib import Path
import torch
import torch.nn as nn
from sklearn.base import BaseEstimator
from sklearn.ensemble import (
    RandomForestClassifier,
    RandomForestRegressor,
    GradientBoostingRegressor,
    GradientBoostingClassifier,
)
from sklearn.neural_network import MLPClassifier, MLPRegressor

from ...processors.structure.ml.activity import ActivityPredictor
from ...processors.structure.ml.ensemble import EnsemblePredictor
from ...processors.structure.ml.features import FeatureExtractor
from ...processors.structure.ml.utils import ModelUtils
from ...processors.structure.descriptors import MolecularDescriptors
from ...processors.structure.pharmacophore import PharmacophoreGenerator
from ...processors.structure.similarity import StructureSimilarity
from ...web_enrichment.data_sources.swiss import SwissClient
from ...web_enrichment.data_sources.chembl import ChEMBLClient
from ...web_enrichment.data_sources.pubchem import PubChemClient


class MLIntegrator:
    """Enhanced ML integration for compound analysis."""

    def __init__(
        self,
        activity_predictor: Optional[ActivityPredictor] = None,
        ensemble_predictor: Optional[EnsemblePredictor] = None,
        feature_extractor: Optional[FeatureExtractor] = None,
        model_utils: Optional[ModelUtils] = None,
        model_config: Optional[Dict] = None,
    ):
        """Initialize ML integration.

        Args:
            activity_predictor: Activity prediction model
            ensemble_predictor: Ensemble prediction model
            feature_extractor: Feature extraction handler
            model_utils: Model utility functions
            model_config: Optional model configuration
        """
        self.logger = logging.getLogger(__name__)

        # Initialize core components
        self.activity_predictor = activity_predictor or ActivityPredictor()
        self.ensemble_predictor = ensemble_predictor or EnsemblePredictor()
        self.feature_extractor = feature_extractor or FeatureExtractor()
        self.model_utils = model_utils or ModelUtils()

        # Initialize descriptor generators
        self.descriptors = MolecularDescriptors()
        self.pharmacophore = PharmacophoreGenerator()
        self.similarity = StructureSimilarity()

        # Initialize data source clients
        self.swiss_client = SwissClient()
        self.chembl_client = ChEMBLClient()
        self.pubchem_client = PubChemClient()

        # Initialize model cache
        self._model_cache = {}
        self._feature_cache = {}

        # Performance tracking
        self.performance_metrics = {
            "binding": {},
            "activity": {},
            "toxicity": {},
            "abuse": {},
            "interactions": {},
        }

        # Initialize models
        self._initialize_models(model_config)

    def _initialize_models(self, config: Optional[Dict] = None):
        """Initialize ML models with configuration."""
        try:
            # Default config
            default_config = {
                "binding": {
                    "model_type": "ensemble",
                    "n_estimators": 100,
                    "max_depth": 10,
                },
                "activity": {
                    "model_type": "neural",
                    "hidden_layers": [512, 256, 128],
                    "dropout": 0.3,
                },
                "toxicity": {
                    "model_type": "ensemble",
                    "n_estimators": 100,
                    "max_depth": 8,
                },
                "abuse": {
                    "model_type": "ensemble",
                    "n_estimators": 100,
                    "max_depth": 8,
                },
            }

            config = config or default_config

            # Initialize models based on config
            for pred_type, model_config in config.items():
                if model_config["model_type"] == "ensemble":
                    if pred_type in ["binding", "activity"]:
                        model = RandomForestRegressor(
                            n_estimators=model_config["n_estimators"],
                            max_depth=model_config["max_depth"],
                            random_state=42,
                        )
                    else:
                        model = RandomForestClassifier(
                            n_estimators=model_config["n_estimators"],
                            max_depth=model_config["max_depth"],
                            random_state=42,
                        )
                elif model_config["model_type"] == "neural":
                    layers = model_config["hidden_layers"]
                    dropout = model_config.get("dropout", 0.3)

                    model = nn.Sequential(
                        nn.Linear(1024, layers[0]),
                        nn.ReLU(),
                        nn.Dropout(dropout),
                        *[
                            layer
                            for i in range(len(layers) - 1)
                            for layer in [
                                nn.Linear(layers[i], layers[i + 1]),
                                nn.ReLU(),
                                nn.Dropout(dropout),
                            ]
                        ],
                        nn.Linear(layers[-1], 1),
                    )

                self._model_cache[pred_type] = model

            # Load pre-trained models
            self._load_pretrained_models()

        except Exception as e:
            self.logger.error(f"Error initializing models: {str(e)}")

    def _load_pretrained_models(self):
        """Load pre-trained model weights."""
        try:
            # Load activity predictor models
            self.activity_predictor.load_models()

            # Load ensemble predictor models
            self.ensemble_predictor.load_models()

        except Exception as e:
            self.logger.error(f"Error loading pre-trained models: {str(e)}")

    def predict_binding(
        self,
        compounds: pd.DataFrame,
        targets: Optional[List[str]] = None,
        confidence_threshold: float = 0.7,
    ) -> Dict[str, Any]:
        """Predict binding affinities with confidence scores.

        Args:
            compounds: Compound DataFrame
            targets: Optional list of target proteins
            confidence_threshold: Minimum confidence threshold

        Returns:
            Dictionary of predictions and metadata
        """
        try:
            # Extract features
            features = self.feature_extractor.get_binding_features(compounds)

            # Get base predictions
            predictions = self.activity_predictor.predict_binding(
                features, targets=targets
            )

            # Get ensemble predictions
            ensemble_preds = self.ensemble_predictor.predict_binding(
                features, targets=targets
            )

            # Get confidence scores
            confidence = self.activity_predictor.get_confidence_scores(
                features, predictions
            )

            # Filter by confidence threshold
            filtered_preds = {}
            for target, pred in predictions.items():
                if confidence[target] >= confidence_threshold:
                    filtered_preds[target] = {
                        "affinity": pred,
                        "confidence": confidence[target],
                        "ensemble_prediction": ensemble_preds.get(target),
                    }

            # Add supporting evidence
            evidence = self.get_supporting_evidence(compounds, filtered_preds)

            # Get similar known binders
            similar_binders = self.get_similar_known_binders(compounds, filtered_preds)

            return {
                "predictions": filtered_preds,
                "evidence": evidence,
                "similar_binders": similar_binders,
                "features_used": features.columns.tolist(),
            }

        except Exception as e:
            self.logger.error(f"Error predicting binding: {str(e)}")
            return {}

    def predict_activity(
        self,
        compounds: pd.DataFrame,
        activity_types: Optional[List[str]] = None,
        include_pharmacophores: bool = True,
    ) -> Dict[str, Any]:
        """Predict activity profiles with pharmacophore analysis.

        Args:
            compounds: Compound DataFrame
            activity_types: Optional list of activity types
            include_pharmacophores: Whether to include pharmacophore analysis

        Returns:
            Dictionary of predictions and metadata
        """
        try:
            # Extract features
            features = self.feature_extractor.get_activity_features(compounds)

            # Get activity predictions
            predictions = self.activity_predictor.predict_activity(
                features, activity_types=activity_types
            )

            # Get ensemble predictions
            ensemble_preds = self.ensemble_predictor.predict_activity(
                features, activity_types=activity_types
            )

            # Get probabilities
            probabilities = self.activity_predictor.predict_proba(features)

            # Add pharmacophore analysis if requested
            pharmacophores = None
            if include_pharmacophores:
                pharmacophores = self.pharmacophore.analyze_compounds(compounds)

            # Get supporting evidence
            evidence = self.get_supporting_evidence(compounds, predictions)

            return {
                "predictions": predictions,
                "ensemble_predictions": ensemble_preds,
                "probabilities": probabilities,
                "pharmacophores": pharmacophores,
                "evidence": evidence,
                "features_used": features.columns.tolist(),
            }

        except Exception as e:
            self.logger.error(f"Error predicting activity: {str(e)}")
            return {}

    def predict_toxicity(
        self,
        compounds: pd.DataFrame,
        include_alerts: bool = True,
    ) -> Dict[str, Any]:
        """Predict toxicity risks with structural alerts.

        Args:
            compounds: Compound DataFrame
            include_alerts: Whether to include structural alerts

        Returns:
            Dictionary of predictions and metadata
        """
        try:
            # Extract features
            features = self.feature_extractor.get_toxicity_features(compounds)

            # Get ensemble predictions
            predictions = self.ensemble_predictor.predict_toxicity(features)

            # Get detailed risk assessment
            risks = self.ensemble_predictor.assess_toxicity_risks(features)

            # Add structural alerts if requested
            alerts = None
            if include_alerts:
                alerts = self.get_toxicity_alerts(compounds)

            # Get supporting evidence
            evidence = self.get_supporting_evidence(compounds, predictions)

            return {
                "predictions": predictions,
                "risks": risks,
                "structural_alerts": alerts,
                "evidence": evidence,
                "features_used": features.columns.tolist(),
            }

        except Exception as e:
            self.logger.error(f"Error predicting toxicity: {str(e)}")
            return {}

    def predict_abuse_potential(
        self,
        compounds: pd.DataFrame,
        include_similarities: bool = True,
    ) -> Dict[str, Any]:
        """Predict abuse potential with similarity analysis.

        Args:
            compounds: Compound DataFrame
            include_similarities: Whether to include similarity analysis

        Returns:
            Dictionary of predictions and metadata
        """
        try:
            # Extract features
            features = self.feature_extractor.get_abuse_features(compounds)

            # Get ensemble predictions
            predictions = self.ensemble_predictor.predict_abuse_potential(features)

            # Get risk factors
            risk_factors = self.ensemble_predictor.assess_abuse_risks(features)

            # Add similarity analysis if requested
            similarities = None
            if include_similarities:
                similarities = self.get_similarity_to_known_drugs(compounds)

            # Get supporting evidence
            evidence = self.get_supporting_evidence(compounds, predictions)

            return {
                "predictions": predictions,
                "risk_factors": risk_factors,
                "structural_similarities": similarities,
                "evidence": evidence,
                "features_used": features.columns.tolist(),
            }

        except Exception as e:
            self.logger.error(f"Error predicting abuse potential: {str(e)}")
            return {}

    def predict_interactions(
        self,
        compounds: pd.DataFrame,
        known_drugs: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """Predict drug interactions.

        Args:
            compounds: Compound DataFrame
            known_drugs: Optional list of known drugs to check

        Returns:
            Dictionary of predictions and metadata
        """
        try:
            # Extract features
            features = self.feature_extractor.get_interaction_features(compounds)

            # Get interaction predictions
            predictions = self.activity_predictor.predict_interactions(
                features, known_drugs=known_drugs
            )

            # Get ensemble predictions
            ensemble_preds = self.ensemble_predictor.predict_interactions(
                features, known_drugs=known_drugs
            )

            # Get supporting evidence
            evidence = self.get_supporting_evidence(compounds, predictions)

            return {
                "predictions": predictions,
                "ensemble_predictions": ensemble_preds,
                "evidence": evidence,
                "features_used": features.columns.tolist(),
            }

        except Exception as e:
            self.logger.error(f"Error predicting interactions: {str(e)}")
            return {}

    def get_supporting_evidence(
        self, compounds: pd.DataFrame, predictions: Dict
    ) -> Dict[str, List[Dict]]:
        """Get supporting evidence for predictions.

        Args:
            compounds: Compound DataFrame
            predictions: Prediction results

        Returns:
            Dictionary mapping compounds to evidence
        """
        try:
            evidence = {}

            for idx, row in compounds.iterrows():
                compound_evidence = []

                # Check structural features
                features = self.descriptors.get_descriptors(row["smiles"])
                if features:
                    compound_evidence.append(
                        {
                            "type": "structural",
                            "features": features,
                        }
                    )

                # Check pharmacophores
                pharmacophores = self.pharmacophore.get_pharmacophores(row["smiles"])
                if pharmacophores:
                    compound_evidence.append(
                        {
                            "type": "pharmacophore",
                            "patterns": pharmacophores,
                        }
                    )

                # Check similar compounds
                similar = self.similarity.find_similar_compounds(
                    row["smiles"], threshold=0.7
                )
                if similar:
                    compound_evidence.append(
                        {
                            "type": "similarity",
                            "compounds": similar,
                        }
                    )

                evidence[row["name"]] = compound_evidence

            return evidence

        except Exception as e:
            self.logger.error(f"Error getting supporting evidence: {str(e)}")
            return {}

    def get_toxicity_alerts(self, compounds: pd.DataFrame) -> Dict[str, List[str]]:
        """Get toxicity structural alerts.

        Args:
            compounds: Compound DataFrame

        Returns:
            Dictionary mapping compounds to alerts
        """
        try:
            alerts = {}

            for idx, row in compounds.iterrows():
                compound_alerts = []

                # Check for known toxic substructures
                toxic_groups = self.descriptors.get_toxic_groups(row["smiles"])
                if toxic_groups:
                    compound_alerts.extend(toxic_groups)

                # Check for reactive groups
                reactive_groups = self.descriptors.get_reactive_groups(row["smiles"])
                if reactive_groups:
                    compound_alerts.extend(reactive_groups)

                # Check for metabolism-related alerts
                metabolic_alerts = self.descriptors.get_metabolic_alerts(row["smiles"])
                if metabolic_alerts:
                    compound_alerts.extend(metabolic_alerts)

                alerts[row["name"]] = compound_alerts

            return alerts

        except Exception as e:
            self.logger.error(f"Error getting toxicity alerts: {str(e)}")
            return {}

    def get_similarity_to_known_drugs(
        self, compounds: pd.DataFrame
    ) -> Dict[str, List[Dict]]:
        """Get structural similarity to known drugs.

        Args:
            compounds: Compound DataFrame

        Returns:
            Dictionary mapping compounds to similar known drugs
        """
        try:
            similarities = {}

            for idx, row in compounds.iterrows():
                # Find similar known drugs
                similar_drugs = self.similarity.find_similar_known_drugs(
                    row["smiles"], threshold=0.6, top_n=5
                )

                if similar_drugs:
                    similarities[row["name"]] = [
                        {
                            "name": drug["name"],
                            "similarity": drug["similarity"],
                            "mechanism": drug.get("mechanism", "unknown"),
                            "scheduling": drug.get("scheduling", "unknown"),
                        }
                        for drug in similar_drugs
                    ]
                else:
                    similarities[row["name"]] = []

            return similarities

        except Exception as e:
            self.logger.error(f"Error getting drug similarities: {str(e)}")
            return {}

    def update_performance_metrics(
        self, prediction_type: str, metrics: Dict[str, float]
    ) -> None:
        """Update model performance metrics.

        Args:
            prediction_type: Type of prediction
            metrics: Performance metrics
        """
        try:
            if prediction_type in self.performance_metrics:
                self.performance_metrics[prediction_type].update(metrics)

                # Log performance update
                self.logger.info(
                    f"Updated {prediction_type} performance metrics: {metrics}"
                )

                # Check for performance degradation
                if self._check_performance_degradation(prediction_type, metrics):
                    self.logger.warning(
                        f"Performance degradation detected for {prediction_type}"
                    )

        except Exception as e:
            self.logger.error(f"Error updating performance metrics: {str(e)}")

    def _check_performance_degradation(
        self, prediction_type: str, metrics: Dict[str, float]
    ) -> bool:
        """Check for model performance degradation.

        Args:
            prediction_type: Type of prediction
            metrics: New performance metrics

        Returns:
            True if degradation detected, False otherwise
        """
        try:
            # Get historical metrics
            historical = self.performance_metrics[prediction_type]

            if not historical:
                return False

            # Check for significant degradation
            degradation_threshold = 0.1  # 10% degradation threshold

            for metric, value in metrics.items():
                if metric in historical:
                    previous = historical[metric]
                    if (previous - value) / previous > degradation_threshold:
                        return True

            return False

        except Exception as e:
            self.logger.error(f"Error checking performance degradation: {str(e)}")
            return False

    def clear_cache(self) -> None:
        """Clear model and feature caches."""
        try:
            self._model_cache.clear()
            self._feature_cache.clear()
            self.logger.info("Cleared model and feature caches")

        except Exception as e:
            self.logger.error(f"Error clearing cache: {str(e)}")
