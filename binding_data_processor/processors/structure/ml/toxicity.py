"""Toxicity prediction module for chemical structures."""

import logging
from typing import Dict, List, Optional, Union

import numpy as np
from rdkit import Chem

from .base import MLPredictor, MLPredictorConfig
from ..base import BaseStructureProcessor
from ...psychopharm.predictors.toxicity import ToxicityPredictor

logger = logging.getLogger(__name__)


class ToxicityPredictor(MLPredictor):
    """Predictor for compound toxicity using ML models."""

    def __init__(
        self,
        config: Optional[Union[Dict, MLPredictorConfig]] = None,
        **kwargs,
    ):
        """Initialize toxicity predictor.

        Args:
            config: Predictor configuration
            **kwargs: Additional arguments passed to parent classes
        """
        super().__init__(config, **kwargs)
        self.toxicity_model = ToxicityPredictor()
        self.endpoints = [
            "acute_toxicity",
            "carcinogenicity",
            "mutagenicity",
            "reproductive_toxicity",
            "hepatotoxicity",
            "cardiotoxicity",
            "neurotoxicity",
        ]

    def preprocess(self, data: Union[str, Chem.Mol]) -> np.ndarray:
        """Preprocess input data.

        Args:
            data: Input compound as SMILES string or RDKit molecule

        Returns:
            Preprocessed features
        """
        mol = self._validate_input(data)
        if mol is None:
            raise ValueError("Invalid compound input")
        return self.featurizer.get_features(mol)

    def predict(
        self,
        data: Union[str, Chem.Mol],
        endpoints: Optional[List[str]] = None,
        **kwargs,
    ) -> Dict[str, float]:
        """Predict toxicity for compound.

        Args:
            data: Input compound as SMILES string or RDKit molecule
            endpoints: List of toxicity endpoints to predict. If None, predicts all endpoints.
            **kwargs: Additional prediction arguments

        Returns:
            Dictionary mapping toxicity endpoints to predicted probabilities
        """
        if not self.is_trained:
            logger.warning("Model not trained, using pre-trained weights")

        try:
            features = self.preprocess(data)
            if endpoints is None:
                endpoints = self.endpoints

            predictions = {}
            for endpoint in endpoints:
                if endpoint not in self.endpoints:
                    logger.warning(f"Unknown endpoint: {endpoint}, skipping")
                    continue
                pred = self.toxicity_model.predict_endpoint(features, endpoint)
                predictions[endpoint] = float(pred)

            return predictions

        except Exception as e:
            logger.error(f"Error predicting toxicity: {str(e)}")
            return {}

    def train(
        self,
        train_data: List[Union[str, Chem.Mol]],
        train_labels: Dict[str, List[float]],
        val_data: Optional[List[Union[str, Chem.Mol]]] = None,
        val_labels: Optional[Dict[str, List[float]]] = None,
        endpoints: Optional[List[str]] = None,
        **kwargs,
    ) -> Dict[str, float]:
        """Train toxicity prediction models.

        Args:
            train_data: Training compounds
            train_labels: Training labels mapping endpoints to label lists
            val_data: Validation compounds
            val_labels: Validation labels mapping endpoints to label lists
            endpoints: List of endpoints to train. If None, trains all endpoints.
            **kwargs: Additional training arguments

        Returns:
            Dictionary of training metrics
        """
        try:
            if endpoints is None:
                endpoints = self.endpoints

            metrics = {}
            for endpoint in endpoints:
                if endpoint not in train_labels:
                    logger.warning(f"No training labels for endpoint: {endpoint}, skipping")
                    continue

                # Preprocess data
                X_train = np.vstack([self.preprocess(x) for x in train_data])
                y_train = np.array(train_labels[endpoint])

                if val_data is not None and val_labels is not None and endpoint in val_labels:
                    X_val = np.vstack([self.preprocess(x) for x in val_data])
                    y_val = np.array(val_labels[endpoint])
                else:
                    X_val = None
                    y_val = None

                # Train model
                endpoint_metrics = self.toxicity_model.train_endpoint(
                    X_train,
                    y_train,
                    endpoint,
                    X_val=X_val,
                    y_val=y_val,
                    **kwargs,
                )
                metrics[endpoint] = endpoint_metrics

            self.is_trained = True
            return metrics

        except Exception as e:
            logger.error(f"Error training models: {str(e)}")
            return {}

    def evaluate(
        self,
        test_data: List[Union[str, Chem.Mol]],
        test_labels: Dict[str, List[float]],
        endpoints: Optional[List[str]] = None,
        metrics: Optional[List[str]] = None,
        **kwargs,
    ) -> Dict[str, Dict[str, float]]:
        """Evaluate toxicity prediction models.

        Args:
            test_data: Test compounds
            test_labels: Test labels mapping endpoints to label lists
            endpoints: List of endpoints to evaluate. If None, evaluates all endpoints.
            metrics: List of metrics to compute
            **kwargs: Additional evaluation arguments

        Returns:
            Dictionary mapping endpoints to metric dictionaries
        """
        try:
            if endpoints is None:
                endpoints = self.endpoints

            if metrics is None:
                metrics = ["accuracy", "precision", "recall", "f1", "auc"]

            results = {}
            for endpoint in endpoints:
                if endpoint not in test_labels:
                    logger.warning(f"No test labels for endpoint: {endpoint}, skipping")
                    continue

                # Preprocess data
                X_test = np.vstack([self.preprocess(x) for x in test_data])
                y_test = np.array(test_labels[endpoint])

                # Evaluate model
                endpoint_metrics = self.toxicity_model.evaluate_endpoint(
                    X_test,
                    y_test,
                    endpoint,
                    metrics=metrics,
                    **kwargs,
                )
                results[endpoint] = endpoint_metrics

            return results

        except Exception as e:
            logger.error(f"Error evaluating models: {str(e)}")
            return {}
