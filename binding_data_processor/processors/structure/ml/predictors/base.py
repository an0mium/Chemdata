"""Base class for ML predictors."""

import os
import json
from datetime import datetime
from typing import Any, Dict, Optional


class BasePredictor:
    """Base class for ML predictors."""

    def __init__(self, model_dir: str):
        """Initialize predictor.

        Args:
            model_dir: Directory containing model files
        """
        self.model_dir = model_dir
        self._version = self._load_version()

    @property
    def version(self) -> str:
        """Get model version."""
        return self._version

    def _load_version(self) -> str:
        """Load version information from model directory."""
        try:
            version_file = os.path.join(self.model_dir, "version.json")
            if os.path.exists(version_file):
                with open(version_file) as f:
                    version_info = json.load(f)
                    return version_info.get("version", "Unknown")
            return self._generate_version()
        except Exception:
            return "Unknown"

    def _generate_version(self) -> str:
        """Generate version information for model."""
        try:
            # Get latest model file modification time
            model_files = [
                f
                for f in os.listdir(self.model_dir)
                if f.endswith((".pkl", ".pt", ".h5", ".model"))
            ]
            if not model_files:
                return "Unknown"

            latest = max(
                os.path.getmtime(os.path.join(self.model_dir, f)) for f in model_files
            )
            timestamp = datetime.fromtimestamp(latest)
            version = timestamp.strftime("%Y%m%d_%H%M%S")

            # Save version info
            version_info = {
                "version": version,
                "files": model_files,
                "generated": datetime.now().isoformat(),
            }
            version_file = os.path.join(self.model_dir, "version.json")
            with open(version_file, "w") as f:
                json.dump(version_info, f, indent=2)

            return version
        except Exception:
            return "Unknown"

    def predict(self, smiles: str) -> Dict[str, Any]:
        """Make predictions for compound.

        Args:
            smiles: SMILES string of compound

        Returns:
            Dictionary containing predictions
        """
        raise NotImplementedError("Subclasses must implement predict()")

    def batch_predict(self, smiles_list: list[str]) -> list[Dict[str, Any]]:
        """Make predictions for multiple compounds.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            List of prediction dictionaries
        """
        return [self.predict(smiles) for smiles in smiles_list]

    def _validate_smiles(self, smiles: str) -> bool:
        """Validate SMILES string.

        Args:
            smiles: SMILES string to validate

        Returns:
            Whether SMILES is valid
        """
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        return mol is not None

    def _get_confidence(self, probability: float) -> float:
        """Calculate confidence score from probability.

        Args:
            probability: Raw probability value

        Returns:
            Confidence score between 0 and 1
        """
        # Simple sigmoid-based confidence scoring
        # Could be replaced with more sophisticated methods
        if probability < 0.5:
            probability = 1 - probability
        return 2 * (probability - 0.5)

    def _format_prediction(
        self,
        label: str,
        probability: float,
        explanation: Optional[str] = None,
        mechanisms: Optional[list[str]] = None,
    ) -> Dict[str, Any]:
        """Format prediction output.

        Args:
            label: Prediction label/class
            probability: Prediction probability
            explanation: Optional explanation of prediction
            mechanisms: Optional list of mechanisms

        Returns:
            Formatted prediction dictionary
        """
        prediction = {
            "label": label,
            "probability": probability,
            "confidence": self._get_confidence(probability),
        }
        if explanation:
            prediction["explanation"] = explanation
        if mechanisms:
            prediction["mechanisms"] = mechanisms
        return prediction
