"""Base classes for binding prediction and analysis."""

import abc
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Union, Any

import numpy as np
import torch
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from rdkit import Chem
from rdkit.Chem import AllChem

from ..base import ActivityPredictor, ActivityPredictorConfig
from ...base import MLProcessor
from ....base import BaseStructureProcessor
from ....ml.models.graph_utils import mol_to_graph, compute_fingerprints
from ..protein import ProteinStructureAnalyzer
from .constants import BINDING_SITE_PARAMS
from .structure import StructureAnalyzer
from .pharmacophore import PharmacophoreGenerator
from .alphafold import AlphaFoldIntegrator

logger = logging.getLogger(__name__)


@dataclass
class BindingPredictorConfig(ActivityPredictorConfig):
    """Configuration for binding predictors."""

    # Additional binding-specific settings
    use_pharmacophore: bool = True
    use_alphafold: bool = True
    min_pocket_volume: float = 100.0
    min_pocket_depth: float = 4.0
    max_pocket_exposure: float = 0.8
    residue_cutoff: float = 8.0
    surface_cutoff: float = 4.0


class BindingPredictor(ActivityPredictor):
    """Base class for predicting binding properties."""

    def __init__(
        self,
        config: Optional[Union[Dict, BindingPredictorConfig]] = None,
        **kwargs,
    ):
        """Initialize binding predictor.

        Args:
            config: Predictor configuration
            **kwargs: Additional arguments passed to parent classes
        """
        if isinstance(config, dict):
            config = BindingPredictorConfig(**config)
        elif config is None:
            config = BindingPredictorConfig()

        super().__init__(config=config, **kwargs)

        # Initialize analyzers
        self.structure_analyzer = StructureAnalyzer()
        self.protein_analyzer = ProteinStructureAnalyzer()
        self.pharmacophore_generator = PharmacophoreGenerator() if config.use_pharmacophore else None
        self.alphafold_integrator = AlphaFoldIntegrator() if config.use_alphafold else None

    def predict_binding(
        self,
        mol: Union[str, Chem.Mol],
        protein: Optional[Structure] = None,
        return_confidence: bool = False,
        **kwargs,
    ) -> Union[Dict[str, float], Tuple[Dict[str, float], Dict[str, float]]]:
        """Predict binding properties.

        Args:
            mol: Input molecule or SMILES
            protein: Optional protein structure
            return_confidence: Whether to return confidence scores
            **kwargs: Additional arguments

        Returns:
            Dictionary of binding predictions, optionally with confidence scores
        """
        try:
            # Convert SMILES if needed
            if isinstance(mol, str):
                mol = Chem.MolFromSmiles(mol)
            if mol is None:
                return {} if not return_confidence else ({}, {})

            # Get molecule features
            features = self.preprocess(mol)
            if features is None:
                return {} if not return_confidence else ({}, {})

            # Get protein features if available
            protein_features = {}
            if protein is not None:
                protein_features = self.protein_analyzer.get_structure_properties(protein)

            # Get pharmacophore features
            pharm_features = {}
            if self.config.use_pharmacophore and self.pharmacophore_generator:
                pharm_features = self.pharmacophore_generator.get_features(mol)

            # Combine all features
            all_features = {
                **features,
                **protein_features,
                **pharm_features,
            }

            # Make predictions
            predictions = self._predict(all_features)
            if not return_confidence:
                return predictions

            # Get confidence scores
            confidence = self._get_confidence(all_features, predictions)
            return predictions, confidence

        except Exception as e:
            self.logger.error(f"Error predicting binding: {str(e)}")
            return {} if not return_confidence else ({}, {})

    def _predict(self, features: Dict[str, Any]) -> Dict[str, float]:
        """Make binding predictions from features.

        Args:
            features: Combined feature dictionary

        Returns:
            Dictionary of predictions
        """
        raise NotImplementedError

    def _get_confidence(
        self,
        features: Dict[str, Any],
        predictions: Dict[str, float],
    ) -> Dict[str, float]:
        """Get confidence scores for predictions.

        Args:
            features: Input features
            predictions: Model predictions

        Returns:
            Dictionary of confidence scores
        """
        raise NotImplementedError


class BindingSitePredictor(BindingPredictor):
    """Predicts protein binding sites using ML and structural analysis."""

    def predict_sites(
        self,
        structure: Structure,
        include_ml_scores: bool = True,
        min_score: float = 0.5,
        max_sites: int = 5,
    ) -> List[Dict[str, Any]]:
        """Predict potential binding sites.

        Args:
            structure: BioPython Structure object
            include_ml_scores: Whether to include ML model predictions
            min_score: Minimum confidence score threshold
            max_sites: Maximum number of sites to return

        Returns:
            List of predicted binding sites with scores and properties
        """
        try:
            # Get structure properties
            properties = self.protein_analyzer.get_structure_properties(structure)

            # Find potential pockets
            pockets = self._detect_pockets(structure, properties)

            # Score and rank pockets
            scored_pockets = self._score_pockets(
                pockets,
                properties,
                include_ml_scores=include_ml_scores,
            )

            # Filter by score and limit number
            filtered_pockets = [p for p in scored_pockets if p["score"] >= min_score]
            filtered_pockets.sort(key=lambda x: x["score"], reverse=True)
            filtered_pockets = filtered_pockets[:max_sites]

            # Add detailed analysis
            for pocket in filtered_pockets:
                pocket["properties"] = self._analyze_pocket_properties(
                    structure,
                    pocket["residues"],
                    properties,
                )

            return filtered_pockets

        except Exception as e:
            self.logger.error(f"Error predicting binding sites: {str(e)}")
            return []

    def _detect_pockets(
        self,
        structure: Structure,
        properties: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Detect potential binding pockets.

        Args:
            structure: BioPython Structure object
            properties: Structure properties

        Returns:
            List of detected pockets
        """
        try:
            # Get surface and cavity points
            surface = self.protein_analyzer.get_molecular_surface(structure)
            cavities = self.protein_analyzer.find_surface_cavities(surface)

            pockets = []
            for cavity in cavities:
                if self._meets_pocket_criteria(cavity, properties):
                    pocket = {
                        "center": cavity["center"],
                        "volume": cavity["volume"],
                        "residues": self._get_pocket_residues(cavity, structure),
                        "surface_points": cavity["surface_points"],
                        "depth": cavity["depth"],
                    }
                    pockets.append(pocket)

            return pockets

        except Exception as e:
            self.logger.error(f"Error detecting pockets: {str(e)}")
            return []

    def _meets_pocket_criteria(
        self,
        cavity: Dict[str, Any],
        properties: Dict[str, Any],
    ) -> bool:
        """Check if cavity meets criteria for binding pocket.

        Args:
            cavity: Cavity properties
            properties: Structure properties

        Returns:
            True if cavity meets criteria
        """
        try:
            # Check volume
            if cavity["volume"] < self.config.min_pocket_volume:
                return False

            # Check depth
            if cavity["depth"] < self.config.min_pocket_depth:
                return False

            # Check exposure
            if cavity.get("exposure", 1.0) > self.config.max_pocket_exposure:
                return False

            return True

        except Exception as e:
            self.logger.error(f"Error checking pocket criteria: {str(e)}")
            return False

    def _get_pocket_residues(
        self,
        cavity: Dict[str, Any],
        structure: Structure,
    ) -> List[int]:
        """Get residues forming a binding pocket.

        Args:
            cavity: Cavity properties
            structure: BioPython Structure object

        Returns:
            List of residue numbers
        """
        try:
            residues = set()
            center = cavity["center"]
            cutoff = self.config.residue_cutoff

            for residue in structure.get_residues():
                # Check if any atom is within cutoff
                for atom in residue:
                    dist = np.linalg.norm(atom.get_coord() - center)
                    if dist < cutoff:
                        residues.add(residue.get_id()[1])
                        break

            return sorted(list(residues))

        except Exception as e:
            self.logger.error(f"Error getting pocket residues: {str(e)}")
            return []

    def _score_pockets(
        self,
        pockets: List[Dict[str, Any]],
        properties: Dict[str, Any],
        include_ml_scores: bool = True,
    ) -> List[Dict[str, Any]]:
        """Score potential binding pockets.

        Args:
            pockets: List of detected pockets
            properties: Structure properties
            include_ml_scores: Whether to include ML model predictions

        Returns:
            List of scored pockets
        """
        try:
            scored_pockets = []
            for pocket in pockets:
                # Calculate geometric score
                geom_score = self._score_geometry(pocket)

                # Calculate conservation score
                cons_score = self._score_conservation(
                    pocket["residues"],
                    properties.get("conservation", {}),
                )

                # Calculate pharmacophore score if available
                pharm_score = 0.0
                if self.pharmacophore_generator:
                    pharm_score = self._score_pharmacophore(pocket, properties)

                # Calculate ML score if requested
                ml_score = 0.0
                if include_ml_scores and self.model:
                    ml_score = self._get_ml_score(pocket, properties)

                # Combine scores
                total_score = 0.3 * geom_score + 0.2 * cons_score + 0.2 * pharm_score + 0.3 * ml_score

                pocket["score"] = float(total_score)
                pocket["score_components"] = {
                    "geometry": float(geom_score),
                    "conservation": float(cons_score),
                    "pharmacophore": float(pharm_score),
                    "ml": float(ml_score),
                }
                scored_pockets.append(pocket)

            return scored_pockets

        except Exception as e:
            self.logger.error(f"Error scoring pockets: {str(e)}")
            return []

    def _score_geometry(self, pocket: Dict[str, Any]) -> float:
        """Score pocket geometry.

        Args:
            pocket: Pocket properties

        Returns:
            Geometric score (0-1)
        """
        try:
            # Score volume (prefer 200-1000 Å³)
            volume = pocket["volume"]
            if volume < 100:
                vol_score = 0.0
            elif volume < 200:
                vol_score = volume / 200
            elif volume < 1000:
                vol_score = 1.0
            else:
                vol_score = np.exp(-(volume - 1000) / 500)

            # Score depth
            depth = pocket["depth"]
            depth_score = 1.0 - np.exp(-depth / 5.0)  # Characteristic length 5Å

            # Score shape
            shape_score = self._score_pocket_shape(pocket)

            # Combine scores
            return float(0.4 * vol_score + 0.3 * depth_score + 0.3 * shape_score)

        except Exception as e:
            self.logger.error(f"Error scoring geometry: {str(e)}")
            return 0.0

    def _score_pocket_shape(self, pocket: Dict[str, Any]) -> float:
        """Score pocket shape properties.

        Args:
            pocket: Pocket properties

        Returns:
            Shape score (0-1)
        """
        try:
            if "surface_points" not in pocket:
                return 0.5

            points = pocket["surface_points"]
            if len(points) < 4:
                return 0.5

            # Calculate principal components
            centered = points - np.mean(points, axis=0)
            cov = np.cov(centered.T)
            eigenvals = np.linalg.eigvalsh(cov)
            eigenvals.sort()

            # Calculate shape descriptors
            if eigenvals[-1] > 0:
                sphericity = eigenvals[0] / eigenvals[-1]
                elongation = eigenvals[-1] / eigenvals[-2]
                flatness = eigenvals[-2] / eigenvals[-3]

                # Score shape (prefer roughly spherical pockets)
                shape_score = 0.4 * sphericity + 0.3 * (1.0 - elongation) + 0.3 * (1.0 - flatness)
                return float(shape_score)

            return 0.5

        except Exception as e:
            self.logger.error(f"Error scoring pocket shape: {str(e)}")
            return 0.5

    def _score_conservation(
        self,
        residues: List[int],
        conservation: Dict[int, float],
    ) -> float:
        """Score pocket conservation.

        Args:
            residues: List of residue numbers
            conservation: Conservation scores by residue

        Returns:
            Conservation score (0-1)
        """
        try:
            if not residues or not conservation:
                return 0.5

            scores = [conservation.get(res, 0.5) for res in residues]
            return float(np.mean(scores))

        except Exception as e:
            self.logger.error(f"Error scoring conservation: {str(e)}")
            return 0.5

    def _score_pharmacophore(
        self,
        pocket: Dict[str, Any],
        properties: Dict[str, Any],
    ) -> float:
        """Score pocket pharmacophore features.

        Args:
            pocket: Pocket properties
            properties: Structure properties

        Returns:
            Pharmacophore score (0-1)
        """
        try:
            if not self.pharmacophore_generator:
                return 0.0

            # Get pharmacophore features
            features = self.pharmacophore_generator.get_pocket_features(
                pocket["residues"],
                properties,
            )

            # Score feature complementarity
            return self.pharmacophore_generator.score_features(features)

        except Exception as e:
            self.logger.error(f"Error scoring pharmacophore: {str(e)}")
            return 0.0

    def _get_ml_score(
        self,
        pocket: Dict[str, Any],
        properties: Dict[str, Any],
    ) -> float:
        """Get ML model prediction score.

        Args:
            pocket: Pocket properties
            properties: Structure properties

        Returns:
            ML prediction score (0-1)
        """
        try:
            if not self.model:
                return 0.0

            # Extract features
            features = self._extract_ml_features(pocket, properties)

            # Get prediction
            score = self.model.predict_proba(features)[0, 1]
            return float(score)

        except Exception as e:
            self.logger.error(f"Error getting ML score: {str(e)}")
            return 0.0

    def _extract_ml_features(
        self,
        pocket: Dict[str, Any],
        properties: Dict[str, Any],
    ) -> np.ndarray:
        """Extract features for ML model.

        Args:
            pocket: Pocket properties
            properties: Structure properties

        Returns:
            Feature array
        """
        try:
            features = []

            # Geometric features
            features.extend(
                [
                    pocket["volume"],
                    pocket["depth"],
                    len(pocket["residues"]),
                ]
            )

            # Conservation features
            if "conservation" in properties:
                cons_scores = [properties["conservation"].get(res, 0.5) for res in pocket["residues"]]
                features.extend(
                    [
                        np.mean(cons_scores),
                        np.std(cons_scores),
                        np.min(cons_scores),
                        np.max(cons_scores),
                    ]
                )
            else:
                features.extend([0.5, 0.0, 0.5, 0.5])

            # Shape features
            if "surface_points" in pocket:
                points = pocket["surface_points"]
                if len(points) >= 4:
                    centered = points - np.mean(points, axis=0)
                    cov = np.cov(centered.T)
                    eigenvals = np.linalg.eigvalsh(cov)
                    eigenvals.sort()
                    features.extend(
                        [
                            eigenvals[0] / eigenvals[-1],  # Sphericity
                            eigenvals[-1] / eigenvals[-2],  # Elongation
                            eigenvals[-2] / eigenvals[-3],  # Flatness
                        ]
                    )
                else:
                    features.extend([0.5, 1.0, 1.0])
            else:
                features.extend([0.5, 1.0, 1.0])

            # Pharmacophore features if available
            if self.pharmacophore_generator:
                pharm_features = self.pharmacophore_generator.get_pocket_features(
                    pocket["residues"],
                    properties,
                )
                features.extend(pharm_features)

            return np.array(features).reshape(1, -1)

        except Exception as e:
            self.logger.error(f"Error extracting ML features: {str(e)}")
            return np.zeros((1, 10))


class BindingModePredictor(BindingPredictor):
    """Predicts binding modes and poses."""

    def predict_modes(
        self,
        mol: Union[str, Chem.Mol],
        site: Optional[Dict[str, Any]] = None,
        return_scores: bool = False,
        **kwargs,
    ) -> Union[List[Dict[str, Any]], Tuple[List[Dict[str, Any]], List[float]]]:
        """Predict binding modes for molecule.

        Args:
            mol: Input molecule or SMILES
            site: Optional binding site to predict modes for
            return_scores: Whether to return mode scores
            **kwargs: Additional arguments

        Returns:
            List of binding mode predictions, optionally with scores
        """
        try:
            # Convert SMILES if needed
            if isinstance(mol, str):
                mol = Chem.MolFromSmiles(mol)
            if mol is None:
                return [] if not return_scores else ([], [])

            # Generate 3D conformer if needed
            if not mol.GetNumConformers():
                AllChem.EmbedMolecule(mol, randomSeed=self.config.random_seed)
                AllChem.MMFFOptimizeMolecule(mol)

            # Get pharmacophore features
            pharm_features = {}
            if self.pharmacophore_generator:
                pharm_features = self.pharmacophore_generator.get_features(mol)

            # Predict binding modes
            modes = self._predict_binding_modes(mol, site, pharm_features)

            # Score modes
            scored_modes = []
            scores = []
            for mode in modes:
                score = self._score_binding_mode(mode, site, pharm_features)
                mode["score"] = score
                scored_modes.append(mode)
                scores.append(score)

            # Sort by score
            sorted_idx = np.argsort(scores)[::-1]
            scored_modes = [scored_modes[i] for i in sorted_idx]
            scores = [scores[i] for i in sorted_idx]

            if return_scores:
                return scored_modes, scores
            return scored_modes

        except Exception as e:
            self.logger.error(f"Error predicting binding modes: {str(e)}")
            return [] if not return_scores else ([], [])

    def _predict_binding_modes(
        self,
        mol: Chem.Mol,
        site: Optional[Dict[str, Any]],
        pharm_features: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Predict possible binding modes.

        Args:
            mol: Input molecule
            site: Optional binding site
            pharm_features: Pharmacophore features

        Returns:
            List of predicted binding modes
        """
        raise NotImplementedError

    def _score_binding_mode(
        self,
        mode: Dict[str, Any],
        site: Optional[Dict[str, Any]],
        pharm_features: Dict[str, Any],
    ) -> float:
        """Score predicted binding mode.

        Args:
            mode: Binding mode
            site: Optional binding site
            pharm_features: Pharmacophore features

        Returns:
            Mode score (0-1)
        """
        raise NotImplementedError
