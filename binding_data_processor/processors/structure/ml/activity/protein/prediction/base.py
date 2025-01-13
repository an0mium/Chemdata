"""Base classes for protein structure prediction."""

from abc import ABC, abstractmethod
from typing import Dict, Optional, Any
from Bio.PDB import Structure

from binding_data_processor.core.config import ProteinAnalysisConfig


class BasePredictor(ABC):
    """Base class for protein structure predictors."""

    def __init__(self, config: Optional[ProteinAnalysisConfig] = None):
        """Initialize predictor.

        Args:
            config: Optional configuration
        """
        self.config = config or ProteinAnalysisConfig()

    @abstractmethod
    async def predict_structure(self, sequence: str, **kwargs) -> Optional[Structure]:
        """Predict protein structure from sequence.

        Args:
            sequence: Amino acid sequence
            **kwargs: Additional arguments

        Returns:
            Predicted structure or None if prediction fails
        """
        pass

    @abstractmethod
    def predict_properties(self, structure: Structure) -> Dict[str, Any]:
        """Predict protein properties from structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary containing predicted properties
        """
        pass


class BasePropertyPredictor(ABC):
    """Base class for protein property predictors."""

    @abstractmethod
    def predict(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """Predict properties from features.

        Args:
            features: Structure features

        Returns:
            Dictionary of predictions
        """
        pass

    @abstractmethod
    def calculate_confidence(self, predictions: Dict[str, Any]) -> Dict[str, float]:
        """Calculate confidence scores for predictions.

        Args:
            predictions: Dictionary of predictions

        Returns:
            Dictionary of confidence scores
        """
        pass
