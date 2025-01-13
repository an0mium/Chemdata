"""Protein structure prediction module."""

import logging
from typing import Optional, Dict, Any, Union
from Bio.PDB import Structure

from ...alphafold import AlphaFoldPredictor, AlphaFoldConfig, AlphaFoldResult
from binding_data_processor.core.config import ProteinAnalysisConfig
from ..analysis import ProteinStructureAnalyzer

logger = logging.getLogger(__name__)


class ProteinPredictor:
    """Predicts protein structures and their properties."""

    def __init__(self, config: Optional[ProteinAnalysisConfig] = None):
        """Initialize predictor.

        Args:
            config: Optional configuration
        """
        self.config = config or ProteinAnalysisConfig()
        self.analyzer = ProteinStructureAnalyzer(self.config)
        self.logger = logging.getLogger(__name__)

        # Initialize AlphaFold predictor if enabled
        if self.config.use_alphafold:
            self.alphafold = AlphaFoldPredictor(
                config=self.config.alphafold_config,
                cache_dir=self.config.cache_dir,
            )

    async def predict_structure(self, sequence: str, use_alphafold: Optional[bool] = None, **kwargs) -> Optional[Structure]:
        """Predict protein structure from sequence.

        Args:
            sequence: Amino acid sequence
            use_alphafold: Whether to use AlphaFold (overrides config)
            **kwargs: Additional arguments passed to the predictor

        Returns:
            Predicted structure or None if prediction fails
        """
        try:
            # Determine whether to use AlphaFold
            if use_alphafold is None:
                use_alphafold = self.config.use_alphafold

            if use_alphafold:
                if not hasattr(self, "alphafold"):
                    self.logger.warning("AlphaFold not initialized")
                    return None

                # Use AlphaFold for prediction
                result = await self.alphafold.predict(sequence)
                if isinstance(result, AlphaFoldResult):
                    return result.structure
                return None

            else:
                self.logger.warning("No alternative prediction method available")
                return None

        except Exception as e:
            self.logger.error(f"Error predicting structure: {str(e)}")
            return None

    def predict_properties(self, structure: Structure) -> Dict[str, Any]:
        """Predict protein properties from structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary containing predicted properties
        """
        try:
            # Get analysis results
            analysis_results = self.analyzer.analyze_structure(structure)

            # Extract features for prediction
            features = self._extract_features(analysis_results)

            # Make predictions based on features
            predictions = self._predict(features)

            return {
                "predictions": predictions,
                "features": features,
                "confidence": self._calculate_confidence(predictions),
                "analysis": analysis_results,
            }

        except Exception as e:
            self.logger.error(f"Error predicting properties: {str(e)}")
            return {}

    def _extract_features(self, analysis_results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract prediction features from analysis results."""
        features = {}

        # Extract structural features
        if "structure" in analysis_results:
            struct_data = analysis_results["structure"]
            features.update(
                {
                    "size": struct_data.get("size", 0),
                    "surface_area": struct_data.get("surface_area", 0.0),
                    "volume": struct_data.get("volume", 0.0),
                    "compactness": struct_data.get("compactness", 0.0),
                }
            )

        # Extract dynamics features
        if "dynamics" in analysis_results:
            dyn_data = analysis_results["dynamics"]
            features.update(
                {
                    "flexibility": dyn_data.get("mean_flexibility", 0.0),
                    "num_domains": dyn_data.get("num_domains", 1),
                    "domain_contacts": dyn_data.get("total_contacts", 0),
                }
            )

        # Extract stability features
        if "stability" in analysis_results:
            stab_data = analysis_results["stability"]
            features.update(
                {
                    "energy": stab_data.get("total_energy", 0.0),
                    "clashes": stab_data.get("num_clashes", 0),
                    "quality_score": stab_data.get("quality_score", 0.0),
                }
            )

        return features

    def _predict(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """Make predictions based on extracted features."""
        predictions = {}

        # Predict stability
        stability_score = self._predict_stability(features)
        predictions["stability"] = {
            "score": stability_score,
            "category": "stable" if stability_score > 0.5 else "unstable",
        }

        # Predict flexibility
        flexibility_score = self._predict_flexibility(features)
        predictions["flexibility"] = {
            "score": flexibility_score,
            "category": "flexible" if flexibility_score > 0.5 else "rigid",
        }

        # Predict function
        function_probs = self._predict_function(features)
        predictions["function"] = {
            "probabilities": function_probs,
            "primary": max(function_probs.items(), key=lambda x: x[1])[0],
        }

        return predictions

    def _predict_stability(self, features: Dict[str, Any]) -> float:
        """Predict protein stability score."""
        # Simple heuristic-based prediction
        quality = features.get("quality_score", 0.0)
        clashes = features.get("clashes", 0)
        energy = features.get("energy", 0.0)

        # Normalize and combine factors
        stability = 0.4 * min(1.0, quality) + 0.3 * (1.0 - min(1.0, clashes / 100)) + 0.3 * (1.0 - min(1.0, abs(energy) / 1000))

        return max(0.0, min(1.0, stability))

    def _predict_flexibility(self, features: Dict[str, Any]) -> float:
        """Predict protein flexibility score."""
        # Simple heuristic-based prediction
        base_flex = features.get("flexibility", 0.0)
        domains = features.get("num_domains", 1)
        contacts = features.get("domain_contacts", 0)

        # Normalize and combine factors
        flexibility = 0.5 * min(1.0, base_flex) + 0.3 * min(1.0, (domains - 1) / 5) + 0.2 * (1.0 - min(1.0, contacts / 100))

        return max(0.0, min(1.0, flexibility))

    def _predict_function(self, features: Dict[str, Any]) -> Dict[str, float]:
        """Predict protein function probabilities."""
        # Extract relevant features
        size = features.get("size", 0)
        surface = features.get("surface_area", 0.0)
        volume = features.get("volume", 0.0)
        flexibility = features.get("flexibility", 0.0)

        # Calculate basic probabilities
        probs = {
            "enzyme": 0.0,
            "transport": 0.0,
            "signaling": 0.0,
            "structural": 0.0,
        }

        # Size-based probabilities
        if size < 200:
            probs["signaling"] += 0.3
        elif size < 500:
            probs["enzyme"] += 0.3
        else:
            probs["transport"] += 0.3

        # Surface area to volume ratio
        if surface > 0 and volume > 0:
            ratio = surface / volume
            if ratio > 0.5:
                probs["signaling"] += 0.2
            else:
                probs["structural"] += 0.2

        # Flexibility-based probabilities
        if flexibility > 0.6:
            probs["enzyme"] += 0.2
            probs["signaling"] += 0.1
        else:
            probs["structural"] += 0.2
            probs["transport"] += 0.1

        # Normalize probabilities
        total = sum(probs.values())
        if total > 0:
            probs = {k: v / total for k, v in probs.items()}

        return probs

    def _calculate_confidence(self, predictions: Dict[str, Any]) -> Dict[str, float]:
        """Calculate confidence scores for predictions."""
        confidence = {}

        # Stability confidence
        if "stability" in predictions:
            stability = predictions["stability"]["score"]
            # Higher confidence when score is far from decision boundary
            confidence["stability"] = 0.5 + abs(stability - 0.5)

        # Flexibility confidence
        if "flexibility" in predictions:
            flexibility = predictions["flexibility"]["score"]
            # Higher confidence when score is far from decision boundary
            confidence["flexibility"] = 0.5 + abs(flexibility - 0.5)

        # Function confidence
        if "function" in predictions:
            probs = predictions["function"]["probabilities"]
            if probs:
                # Higher confidence when one category dominates
                top_prob = max(probs.values())
                confidence["function"] = top_prob

        return confidence
