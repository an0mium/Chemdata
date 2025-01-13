"""Molecular interaction prediction module."""

import logging
from typing import Dict, List, Optional, Union, Any

import numpy as np
import torch
from rdkit import Chem

from .base import ActivityPredictor, ActivityPredictorConfig


class InteractionPredictor(ActivityPredictor):
    """Predict molecular interactions."""

    def __init__(
        self,
        config: Optional[Union[Dict, ActivityPredictorConfig]] = None,
        **kwargs,
    ):
        """Initialize interaction predictor.

        Args:
            config: Predictor configuration
            **kwargs: Additional arguments passed to parent classes
        """
        super().__init__(config=config, **kwargs)
        self.logger = logging.getLogger(self.__class__.__name__)

    def predict(
        self,
        mol: Union[str, Chem.Mol],
        return_confidence: bool = False,
        threshold: float = 0.5,
        **kwargs,
    ) -> Union[Dict[str, float], tuple[Dict[str, float], Dict[str, float]]]:
        """Predict molecular interactions.

        Args:
            mol: Input molecule or SMILES
            return_confidence: Whether to return prediction confidence
            threshold: Confidence threshold for predictions
            **kwargs: Additional arguments

        Returns:
            Dictionary of interaction predictions, optionally with confidence scores
        """
        try:
            # Convert SMILES to molecule if needed
            if isinstance(mol, str):
                mol = Chem.MolFromSmiles(mol)
            if mol is None:
                raise ValueError("Invalid molecule")

            # Preprocess molecule
            features = self.preprocess(mol)
            if features is None:
                raise ValueError("Error preprocessing molecule")

            # Make predictions
            with torch.no_grad():
                outputs = self.model(features)
                probabilities = torch.sigmoid(outputs).cpu().numpy()

            # Get predictions above threshold
            predictions = {}
            confidences = {}
            for i, (interaction, prob) in enumerate(zip(self.model.interaction_types, probabilities[0])):
                if prob >= threshold:
                    predictions[interaction] = float(prob)
                    if return_confidence:
                        confidences[interaction] = float(prob)

            if return_confidence:
                return predictions, confidences
            return predictions

        except Exception as e:
            self.logger.error(f"Error predicting interactions: {str(e)}")
            if return_confidence:
                return {}, {}
            return {}

    def analyze_interactions(
        self,
        mol: Union[str, Chem.Mol],
        include_substructures: bool = True,
        **kwargs,
    ) -> Dict[str, Any]:
        """Analyze predicted interactions and relevant molecular features.

        Args:
            mol: Input molecule or SMILES
            include_substructures: Whether to analyze interaction-specific substructures
            **kwargs: Additional arguments

        Returns:
            Dictionary containing interaction analysis results
        """
        try:
            # Get predictions
            predictions, confidences = self.predict(mol, return_confidence=True, **kwargs)

            # Sort interactions by confidence
            sorted_interactions = sorted(predictions.items(), key=lambda x: x[1], reverse=True)

            analysis = {
                "predictions": predictions,
                "confidences": confidences,
                "top_interactions": [interaction for interaction, _ in sorted_interactions[:5]],
            }

            # Analyze substructures if requested
            if include_substructures and isinstance(mol, Chem.Mol):
                from rdkit.Chem import AllChem

                analysis["substructures"] = {}
                for interaction in predictions:
                    if hasattr(self.model, f"{interaction}_smarts"):
                        smarts = getattr(self.model, f"{interaction}_smarts")
                        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
                        analysis["substructures"][interaction] = len(matches)

            return analysis

        except Exception as e:
            self.logger.error(f"Error analyzing interactions: {str(e)}")
            return {}
