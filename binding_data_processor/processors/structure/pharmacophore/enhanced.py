"""Enhanced pharmacophore detection with integrated structure analysis capabilities."""

import logging
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from ..mixins import (
    StructureValidationMixin,
    ConformerGenerationMixin,
    PharmacophoreFeatureMixin,
    StructureAlignmentMixin,
    StructureVisualizationMixin,
)
from .base import PharmacophoreGenerator
from .features import PharmacophoreFeature


class EnhancedPharmacophoreGenerator(
    PharmacophoreGenerator,
    StructureValidationMixin,
    ConformerGenerationMixin,
    PharmacophoreFeatureMixin,
    StructureAlignmentMixin,
    StructureVisualizationMixin,
):
    """Enhanced pharmacophore generator with integrated capabilities."""

    def __init__(self):
        """Initialize enhanced pharmacophore generator."""
        super().__init__()
        self.logger = logging.getLogger(__name__)

    def generate_from_smiles(self, smiles: str) -> List[PharmacophoreFeature]:
        """Generate pharmacophore features from SMILES string.

        Args:
            smiles: SMILES string of molecule

        Returns:
            List of pharmacophore features
        """
        mol = self.validate_structure(smiles)
        if mol is None:
            return []

        # Generate conformer if needed
        if not mol.GetNumConformers():
            mol = self.generate_conformers(mol, n_confs=1)
            if mol is None:
                return []

        return self.generate(mol)

    def align_pharmacophores(
        self,
        ref_mol: Chem.Mol,
        probe_mol: Chem.Mol,
    ) -> Tuple[float, List[PharmacophoreFeature], List[PharmacophoreFeature]]:
        """Align pharmacophores of two molecules.

        Args:
            ref_mol: Reference molecule
            probe_mol: Probe molecule to align

        Returns:
            Tuple of (RMSD, ref features, aligned probe features)
        """
        if ref_mol is None or probe_mol is None:
            return float("inf"), [], []

        # Generate pharmacophores
        ref_features = self.generate(ref_mol)
        probe_features = self.generate(probe_mol)

        # Align structures
        rmsd, aligned_mol = self.align_structures(
            ref_mol,
            probe_mol,
            ref_features,
            probe_features,
        )

        # Generate features for aligned structure
        aligned_features = self.generate(aligned_mol)

        return rmsd, ref_features, aligned_features

    def compare_pharmacophores(
        self,
        features1: List[PharmacophoreFeature],
        features2: List[PharmacophoreFeature],
        distance_cutoff: float = 2.0,
    ) -> Dict[str, float]:
        """Compare two pharmacophore models.

        Args:
            features1: First set of pharmacophore features
            features2: Second set of pharmacophore features
            distance_cutoff: Maximum distance for matching features

        Returns:
            Dictionary of similarity metrics
        """
        try:
            if not features1 or not features2:
                return {
                    "overlap": 0.0,
                    "coverage": 0.0,
                    "similarity": 0.0,
                }

            # Group features by type
            features1_by_type = self._group_by_type(features1)
            features2_by_type = self._group_by_type(features2)

            # Calculate matches for each feature type
            matches = 0
            total_features = 0
            for feature_type in set(features1_by_type.keys()) | set(features2_by_type.keys()):
                f1 = features1_by_type.get(feature_type, [])
                f2 = features2_by_type.get(feature_type, [])

                if not f1 or not f2:
                    total_features += len(f1) + len(f2)
                    continue

                # Find matching features based on distance
                matched = set()
                for i, feat1 in enumerate(f1):
                    for j, feat2 in enumerate(f2):
                        if j in matched:
                            continue
                        dist = np.linalg.norm(np.array(feat1.position) - np.array(feat2.position))
                        if dist <= distance_cutoff:
                            matches += 1
                            matched.add(j)
                            break

                total_features += len(f1) + len(f2)

            # Calculate similarity metrics
            overlap = matches / min(len(features1), len(features2))
            coverage = matches / max(len(features1), len(features2))
            similarity = 2 * matches / total_features

            return {
                "overlap": float(overlap),
                "coverage": float(coverage),
                "similarity": float(similarity),
            }

        except Exception as e:
            self.logger.error(f"Error comparing pharmacophores: {str(e)}")
            return {
                "overlap": 0.0,
                "coverage": 0.0,
                "similarity": 0.0,
            }

    def visualize_pharmacophore(
        self,
        mol: Chem.Mol,
        features: List[PharmacophoreFeature],
        show_vectors: bool = True,
    ) -> Optional[str]:
        """Generate pharmacophore visualization.

        Args:
            mol: Molecule to visualize
            features: Pharmacophore features to highlight
            show_vectors: Whether to show feature vectors

        Returns:
            SVG string of visualization
        """
        try:
            # Prepare atom highlights
            highlight_atoms = {}
            for feature in features:
                color = self._get_feature_color(feature.feature_type)
                for atom_idx in feature.atoms:
                    highlight_atoms[atom_idx] = color

            # Generate visualization
            return self.visualize_features(mol, features)

        except Exception as e:
            self.logger.error(f"Error visualizing pharmacophore: {str(e)}")
            return None

    def _group_by_type(
        self,
        features: List[PharmacophoreFeature],
    ) -> Dict[str, List[PharmacophoreFeature]]:
        """Group features by type."""
        grouped = {}
        for feature in features:
            if feature.feature_type not in grouped:
                grouped[feature.feature_type] = []
            grouped[feature.feature_type].append(feature)
        return grouped

    def _get_feature_color(self, feature_type: str) -> Tuple[float, float, float]:
        """Get color for feature type."""
        colors = {
            "hbd": (0, 1, 0),  # Green
            "hba": (1, 0, 0),  # Red
            "aromatic": (1, 1, 0),  # Yellow
            "hydrophobe": (0, 0, 1),  # Blue
            "positive": (1, 0, 1),  # Magenta
            "negative": (0, 1, 1),  # Cyan
            "metal": (0.7, 0.7, 0.7),  # Gray
            "xbond": (1, 0.5, 0),  # Orange
            "stereocenter": (0.5, 0, 0.5),  # Purple
            "conjugated": (0.8, 0.4, 0),  # Brown
        }
        return colors.get(feature_type, (0.5, 0.5, 0.5))  # Default gray
