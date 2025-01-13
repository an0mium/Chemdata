"""Target prediction module."""

import logging
from typing import Dict, List, Optional, Union, Any

import numpy as np
import torch
from rdkit import Chem

from .base import ActivityPredictor, ActivityPredictorConfig


class TargetPredictor(ActivityPredictor):
    """Predict molecular targets for compounds."""

    def __init__(
        self,
        config: Optional[Union[Dict, ActivityPredictorConfig]] = None,
        **kwargs,
    ):
        """Initialize target predictor.

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
        """Predict molecular targets.

        Args:
            mol: Input molecule or SMILES
            return_confidence: Whether to return prediction confidence
            threshold: Confidence threshold for predictions
            **kwargs: Additional arguments

        Returns:
            Dictionary of target predictions, optionally with confidence scores
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
            for i, (target, prob) in enumerate(zip(self.model.target_names, probabilities[0])):
                if prob >= threshold:
                    predictions[target] = float(prob)
                    if return_confidence:
                        confidences[target] = float(prob)

            if return_confidence:
                return predictions, confidences
            return predictions

        except Exception as e:
            self.logger.error(f"Error predicting targets: {str(e)}")
            if return_confidence:
                return {}, {}
            return {}

    def analyze_targets(
        self,
        mol: Union[str, Chem.Mol],
        include_substructures: bool = True,
        **kwargs,
    ) -> Dict[str, Any]:
        """Analyze predicted targets and relevant molecular features.

        Args:
            mol: Input molecule or SMILES
            include_substructures: Whether to analyze target-specific substructures
            **kwargs: Additional arguments

        Returns:
            Dictionary containing target analysis results
        """
        try:
            # Get predictions
            predictions, confidences = self.predict(mol, return_confidence=True, **kwargs)

            # Sort targets by confidence
            sorted_targets = sorted(predictions.items(), key=lambda x: x[1], reverse=True)

            analysis = {
                "predictions": predictions,
                "confidences": confidences,
                "top_targets": [target for target, _ in sorted_targets[:5]],
            }

            # Analyze substructures if requested
            if include_substructures and isinstance(mol, Chem.Mol):
                from rdkit.Chem import AllChem

                analysis["substructures"] = {}
                for target in predictions:
                    if hasattr(self.model, f"{target}_smarts"):
                        smarts = getattr(self.model, f"{target}_smarts")
                        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
                        analysis["substructures"][target] = len(matches)

            return analysis

        except Exception as e:
            self.logger.error(f"Error analyzing targets: {str(e)}")
            return {}
