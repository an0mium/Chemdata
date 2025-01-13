"""Pharmacophore detection and analysis functionality."""

from typing import Dict, List, Optional, Tuple, Union
import logging

from rdkit import Chem

from .base import PharmacophoreGenerator
from .features import PharmacophoreFeature
from .conformers import ConformerGenerator
from .alignment import MolecularAligner


class PharmacophoreDetector:
    """Class for detecting and analyzing pharmacophore features in molecules."""

    def __init__(self):
        """Initialize pharmacophore detector with required components."""
        self.generator = PharmacophoreGenerator()
        self.conformer_generator = ConformerGenerator()
        self.aligner = MolecularAligner()
        self.logger = logging.getLogger(__name__)

    def detect_features(
        self,
        mol: Chem.Mol,
        include_vectors: bool = True,
        generate_conformers: bool = True,
        optimize_conformers: bool = True,
    ) -> Dict[str, List[PharmacophoreFeature]]:
        """Detect pharmacophore features in a molecule.

        Args:
            mol: Input molecule
            include_vectors: Whether to calculate feature vectors
            generate_conformers: Whether to generate 3D conformers if needed
            optimize_conformers: Whether to optimize generated conformers

        Returns:
            Dictionary mapping feature types to lists of features
        """
        try:
            # Generate conformers if needed
            if generate_conformers and not mol.GetNumConformers():
                mol = self.conformer_generator.generate(mol, optimize=optimize_conformers)
                if mol is None:
                    raise ValueError("Failed to generate conformers")

            # Detect features
            features = self.generator.generate_features(
                mol,
                include_vectors=include_vectors,
            )

            return features

        except Exception as e:
            self.logger.error(f"Error detecting pharmacophore features: {str(e)}")
            return {}

    def align_to_template(
        self,
        mol: Chem.Mol,
        template: Chem.Mol,
        mol_conf_id: int = 0,
        template_conf_id: int = 0,
        symmetry_aware: bool = True,
    ) -> Tuple[float, Chem.Mol]:
        """Align a molecule to a template molecule.

        Args:
            mol: Molecule to align
            template: Template molecule
            mol_conf_id: Conformer ID for mol
            template_conf_id: Conformer ID for template
            symmetry_aware: Whether to consider molecular symmetry

        Returns:
            Tuple of (RMSD, aligned molecule)
        """
        try:
            result = self.aligner.align_mol_to_template(
                mol,
                template,
                mol_conf_id=mol_conf_id,
                template_conf_id=template_conf_id,
                symmetry_aware=symmetry_aware,
            )
            return result.rmsd, result.aligned_mol

        except Exception as e:
            self.logger.error(f"Error aligning to template: {str(e)}")
            return float("inf"), mol

    def align_to_pharmacophore(
        self,
        mol: Chem.Mol,
        pharm_points: List[Tuple[float, float, float]],
        weights: Optional[List[float]] = None,
        feature_types: Optional[List[str]] = None,
        conf_id: int = 0,
    ) -> Tuple[float, Chem.Mol]:
        """Align a molecule to pharmacophore points.

        Args:
            mol: Molecule to align
            pharm_points: List of 3D pharmacophore point coordinates
            weights: Optional weights for pharmacophore points
            feature_types: Optional feature types for points
            conf_id: Conformer ID to use

        Returns:
            Tuple of (RMSD, aligned molecule)
        """
        try:
            result = self.aligner.align_to_pharmacophore(
                mol,
                pharm_points,
                weights=weights,
                feature_types=feature_types,
                conf_id=conf_id,
            )
            return result.rmsd, result.aligned_mol

        except Exception as e:
            self.logger.error(f"Error aligning to pharmacophore: {str(e)}")
            return float("inf"), mol

    def calculate_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        conf_id1: int = 0,
        conf_id2: int = 0,
    ) -> Dict[str, float]:
        """Calculate similarity metrics between two molecules.

        Args:
            mol1: First molecule
            mol2: Second molecule
            conf_id1: Conformer ID for first molecule
            conf_id2: Conformer ID for second molecule

        Returns:
            Dictionary of similarity scores
        """
        try:
            # Calculate shape similarity
            shape_score = self.aligner.calculate_shape_similarity(
                mol1,
                mol2,
                conf_id1=conf_id1,
                conf_id2=conf_id2,
            )

            # Calculate volume overlap
            volume_score = self.aligner.calculate_volume_overlap(
                mol1,
                mol2,
                conf_id1=conf_id1,
                conf_id2=conf_id2,
            )

            # Calculate feature overlap
            feature_score = self.aligner.calculate_feature_overlap(
                mol1,
                mol2,
                conf_id1=conf_id1,
                conf_id2=conf_id2,
            )

            return {
                "shape_similarity": shape_score,
                "volume_overlap": volume_score,
                "feature_overlap": feature_score,
                "combined_score": (shape_score + volume_score + feature_score) / 3,
            }

        except Exception as e:
            self.logger.error(f"Error calculating similarity: {str(e)}")
            return {
                "shape_similarity": 0.0,
                "volume_overlap": 0.0,
                "feature_overlap": 0.0,
                "combined_score": 0.0,
            }


__all__ = [
    "PharmacophoreDetector",
    "PharmacophoreGenerator",
    "PharmacophoreFeature",
    "ConformerGenerator",
    "MolecularAligner",
]
