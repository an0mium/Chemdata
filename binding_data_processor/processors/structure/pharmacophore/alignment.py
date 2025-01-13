"""Molecular alignment and shape-based methods for pharmacophore analysis."""

from typing import Dict, List, Optional, Tuple, Union
import logging
from dataclasses import dataclass

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import rdShapeHelpers, rdDepictor
from rdkit.Geometry import Point3D
from rdkit.Chem import rdMolTransforms

from .conformers import ConformerGenerator


@dataclass
class AlignmentResult:
    """Results from molecular alignment."""

    rmsd: float
    shape_score: float
    volume_score: float
    feature_score: float
    aligned_mol: Chem.Mol
    rotation: np.ndarray
    translation: np.ndarray


class MolecularAligner:
    """Enhanced class for aligning molecules and calculating shape similarity."""

    def __init__(self, use_crippen: bool = True, use_features: bool = True):
        """Initialize the molecular aligner.

        Args:
            use_crippen: Whether to use Crippen contributions for alignment
            use_features: Whether to use pharmacophore features for alignment
        """
        self.conformer_generator = ConformerGenerator()
        self.use_crippen = use_crippen
        self.use_features = use_features
        self.logger = logging.getLogger(__name__)

    def align_mol_to_template(
        self,
        mol: Chem.Mol,
        template: Chem.Mol,
        mol_conf_id: int = 0,
        template_conf_id: int = 0,
        symmetry_aware: bool = True,
    ) -> AlignmentResult:
        """Align a molecule to a template molecule with enhanced options.

        Args:
            mol: Molecule to align
            template: Template molecule
            mol_conf_id: Conformer ID for mol
            template_conf_id: Conformer ID for template
            symmetry_aware: Whether to consider molecular symmetry

        Returns:
            AlignmentResult object
        """
        try:
            # Generate conformers if needed
            if not mol.GetNumConformers():
                mol = self.conformer_generator.generate(mol)
            if not template.GetNumConformers():
                template = self.conformer_generator.generate(template)

            if mol is None or template is None:
                raise ValueError("Failed to generate conformers")

            # Create copy for alignment
            aligned_mol = Chem.Mol(mol)

            # Calculate Crippen contributions if enabled
            if self.use_crippen:
                crippen_mol = rdMolDescriptors._CalcCrippenContribs(mol)
                crippen_template = rdMolDescriptors._CalcCrippenContribs(template)
            else:
                crippen_mol = crippen_template = None

            # Perform alignment
            if symmetry_aware:
                rmsd, (rotation, translation) = self._align_with_symmetry(
                    aligned_mol,
                    template,
                    mol_conf_id,
                    template_conf_id,
                    crippen_mol,
                    crippen_template,
                )
            else:
                rmsd = AllChem.AlignMol(
                    aligned_mol,
                    template,
                    prbCid=mol_conf_id,
                    refCid=template_conf_id,
                )
                rotation, translation = self._get_alignment_transform(aligned_mol, template)

            # Calculate additional scores
            shape_score = self.calculate_shape_similarity(
                aligned_mol,
                template,
                mol_conf_id,
                template_conf_id,
            )
            volume_score = self.calculate_volume_overlap(
                aligned_mol,
                template,
                mol_conf_id,
                template_conf_id,
            )
            feature_score = self.calculate_feature_overlap(
                aligned_mol,
                template,
                mol_conf_id,
                template_conf_id,
            )

            return AlignmentResult(
                rmsd=rmsd,
                shape_score=shape_score,
                volume_score=volume_score,
                feature_score=feature_score,
                aligned_mol=aligned_mol,
                rotation=rotation,
                translation=translation,
            )

        except Exception as e:
            self.logger.error(f"Error in molecule alignment: {str(e)}")
            raise

    def align_to_pharmacophore(
        self,
        mol: Chem.Mol,
        pharm_points: List[Tuple[float, float, float]],
        weights: Optional[List[float]] = None,
        feature_types: Optional[List[str]] = None,
        conf_id: int = 0,
    ) -> AlignmentResult:
        """Align a molecule to pharmacophore points with enhanced options.

        Args:
            mol: Molecule to align
            pharm_points: List of 3D pharmacophore point coordinates
            weights: Optional weights for pharmacophore points
            feature_types: Optional feature types for points
            conf_id: Conformer ID to use

        Returns:
            AlignmentResult object
        """
        try:
            if not weights:
                weights = [1.0] * len(pharm_points)

            # Convert points to numpy array
            points = np.array(pharm_points)

            # Generate conformers if needed
            if not mol.GetNumConformers():
                mol = self.conformer_generator.generate(mol)

            if mol is None:
                raise ValueError("Failed to generate conformers")

            # Create copy for alignment
            aligned_mol = Chem.Mol(mol)

            # Get molecule points (use atom positions)
            conf = aligned_mol.GetConformer(conf_id)
            mol_points = []
            for i in range(aligned_mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                mol_points.append([pos.x, pos.y, pos.z])
            mol_points = np.array(mol_points)

            # Calculate optimal rotation and translation
            rotation, translation = self._kabsch_align(mol_points, points, weights)

            # Apply transformation
            conf = aligned_mol.GetConformer(conf_id)
            rdMolTransforms.TransformConformer(conf, rotation, translation)

            # Calculate scores
            rmsd = self._weighted_rmsd(mol_points, points, weights)
            shape_score = 0.0  # No reference molecule for shape comparison
            volume_score = 0.0  # No reference molecule for volume comparison
            feature_score = (
                self._calculate_feature_match(
                    aligned_mol,
                    pharm_points,
                    feature_types,
                    weights,
                )
                if feature_types
                else 0.0
            )

            return AlignmentResult(
                rmsd=rmsd,
                shape_score=shape_score,
                volume_score=volume_score,
                feature_score=feature_score,
                aligned_mol=aligned_mol,
                rotation=rotation,
                translation=translation,
            )

        except Exception as e:
            self.logger.error(f"Error in pharmacophore alignment: {str(e)}")
            raise

    def calculate_shape_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        conf_id1: int = 0,
        conf_id2: int = 0,
        grid_spacing: float = 0.5,
    ) -> float:
        """Calculate enhanced shape similarity between two molecules.

        Args:
            mol1: First molecule
            mol2: Second molecule
            conf_id1: Conformer ID for first molecule
            conf_id2: Conformer ID for second molecule
            grid_spacing: Grid spacing for shape calculation

        Returns:
            Shape Tanimoto score (0-1)
        """
        try:
            # Generate conformers if needed
            if not mol1.GetNumConformers():
                mol1 = self.conformer_generator.generate(mol1)
            if not mol2.GetNumConformers():
                mol2 = self.conformer_generator.generate(mol2)

            if mol1 is None or mol2 is None:
                raise ValueError("Failed to generate conformers")

            # Calculate shape similarity using different methods
            tanimoto = 1.0 - rdShapeHelpers.ShapeTanimotoDist(
                mol1,
                mol2,
                confId1=conf_id1,
                confId2=conf_id2,
                gridSpacing=grid_spacing,
            )
            protrude = 1.0 - rdShapeHelpers.ShapeProtrudeDist(
                mol1,
                mol2,
                confId1=conf_id1,
                confId2=conf_id2,
                gridSpacing=grid_spacing,
            )

            # Combine scores (weighted average)
            return 0.7 * tanimoto + 0.3 * protrude

        except Exception as e:
            self.logger.error(f"Error calculating shape similarity: {str(e)}")
            return 0.0

    def calculate_volume_overlap(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        conf_id1: int = 0,
        conf_id2: int = 0,
    ) -> float:
        """Calculate volume overlap between two molecules.

        Args:
            mol1: First molecule
            mol2: Second molecule
            conf_id1: Conformer ID for first molecule
            conf_id2: Conformer ID for second molecule

        Returns:
            Volume overlap score (0-1)
        """
        try:
            # Calculate molecular volumes
            vol1 = AllChem.ComputeMolVolume(mol1, confId=conf_id1)
            vol2 = AllChem.ComputeMolVolume(mol2, confId=conf_id2)

            # Calculate overlap volume using shape helper
            overlap = rdShapeHelpers.ShapeProtrudeDist(
                mol1,
                mol2,
                confId1=conf_id1,
                confId2=conf_id2,
            )

            # Calculate Tanimoto-like score
            if vol1 + vol2 - overlap > 0:
                return overlap / (vol1 + vol2 - overlap)
            return 0.0

        except Exception as e:
            self.logger.error(f"Error calculating volume overlap: {str(e)}")
            return 0.0

    def calculate_feature_overlap(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        conf_id1: int = 0,
        conf_id2: int = 0,
    ) -> float:
        """Calculate pharmacophore feature overlap between molecules.

        Args:
            mol1: First molecule
            mol2: Second molecule
            conf_id1: Conformer ID for first molecule
            conf_id2: Conformer ID for second molecule

        Returns:
            Feature overlap score (0-1)
        """
        try:
            if not self.use_features:
                return 0.0

            # Get conformers
            conf1 = mol1.GetConformer(conf_id1)
            conf2 = mol2.GetConformer(conf_id2)

            # Calculate feature vectors (using Crippen contributions as proxy)
            feat1 = rdMolDescriptors._CalcCrippenContribs(mol1)
            feat2 = rdMolDescriptors._CalcCrippenContribs(mol2)

            # Calculate spatial feature overlap
            overlap = 0.0
            total = 0.0

            for (idx1, f1), (idx2, f2) in zip(enumerate(feat1), enumerate(feat2)):
                if abs(f1[0] - f2[0]) < 0.1:  # Similar feature type
                    # Get 3D positions
                    pos1 = conf1.GetAtomPosition(idx1)
                    pos2 = conf2.GetAtomPosition(idx2)

                    # Calculate distance
                    dist = pos1.Distance(pos2)
                    if dist < 2.0:  # Distance threshold
                        overlap += np.exp(-dist)
                    total += 1.0

            return overlap / total if total > 0 else 0.0

        except Exception as e:
            self.logger.error(f"Error calculating feature overlap: {str(e)}")
            return 0.0

    def _align_with_symmetry(
        self,
        mol: Chem.Mol,
        template: Chem.Mol,
        mol_conf_id: int,
        template_conf_id: int,
        crippen_mol: Optional[List] = None,
        crippen_template: Optional[List] = None,
    ) -> Tuple[float, Tuple[np.ndarray, np.ndarray]]:
        """Perform symmetry-aware alignment.

        Args:
            mol: Molecule to align
            template: Template molecule
            mol_conf_id: Conformer ID for mol
            template_conf_id: Conformer ID for template
            crippen_mol: Optional Crippen contributions for mol
            crippen_template: Optional Crippen contributions for template

        Returns:
            Tuple of (RMSD, transform)
        """
        # Get symmetry classes
        sym_mol = rdMolDescriptors.GetSymmSSSR(mol)
        sym_template = rdMolDescriptors.GetSymmSSSR(template)

        best_rmsd = float("inf")
        best_rotation = None
        best_translation = None

        # Try different symmetry-equivalent alignments
        for sym_op in self._generate_symmetry_ops(sym_mol, sym_template):
            # Create copy for this alignment attempt
            test_mol = Chem.Mol(mol)

            # Apply symmetry operation
            self._apply_symmetry_op(test_mol, sym_op, mol_conf_id)

            # Perform alignment
            rmsd = AllChem.AlignMol(
                test_mol,
                template,
                prbCid=mol_conf_id,
                refCid=template_conf_id,
            )

            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_rotation, best_translation = self._get_alignment_transform(test_mol, template)

        if best_rotation is None or best_translation is None:
            raise ValueError("Failed to find valid alignment")

        # Apply best transformation to original molecule
        conf = mol.GetConformer(mol_conf_id)
        rdMolTransforms.TransformConformer(conf, best_rotation, best_translation)

        return best_rmsd, (best_rotation, best_translation)

    def _generate_symmetry_ops(self, sym1: List, sym2: List) -> List[Tuple[np.ndarray, np.ndarray]]:
        """Generate symmetry operations for alignment.

        Args:
            sym1: Symmetry classes for first molecule
            sym2: Symmetry classes for second molecule

        Returns:
            List of possible symmetry transforms
        """
        ops = [(np.eye(3), np.zeros(3))]  # Identity transform

        # Add rotations around symmetry axes
        for axis in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
            for angle in [90, 180, 270]:
                # Create rotation matrix for this axis and angle
                theta = np.radians(angle)
                c = np.cos(theta)
                s = np.sin(theta)
                x, y, z = axis
                rotation = np.array(
                    [
                        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
                        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
                        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
                    ]
                )
                ops.append((rotation, np.zeros(3)))

        return ops

    def _apply_symmetry_op(
        self,
        mol: Chem.Mol,
        transform: Tuple[np.ndarray, np.ndarray],
        conf_id: int = 0,
    ) -> None:
        """Apply symmetry operation to molecule.

        Args:
            mol: Molecule to transform
            transform: Symmetry transform to apply
            conf_id: Conformer ID to transform
        """
        rotation, translation = transform
        conf = mol.GetConformer(conf_id)
        rdMolTransforms.TransformConformer(conf, rotation, translation)

    def _get_alignment_transform(
        self,
        aligned: Chem.Mol,
        reference: Chem.Mol,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Extract alignment transform between molecules.

        Args:
            aligned: Aligned molecule
            reference: Reference molecule

        Returns:
            Transform that aligns the molecules
        """
        # Get conformer positions
        conf_aligned = aligned.GetConformer()
        conf_ref = reference.GetConformer()

        # Calculate transform using three non-collinear points
        points_aligned = []
        points_ref = []
        for i in range(min(3, aligned.GetNumAtoms())):
            points_aligned.append(conf_aligned.GetAtomPosition(i))
            points_ref.append(conf_ref.GetAtomPosition(i))

        # Convert points to numpy arrays
        p1 = np.array([p.x for p in points_aligned])
        p2 = np.array([p.y for p in points_aligned])
        p3 = np.array([p.z for p in points_aligned])
        points1 = np.vstack([p1, p2, p3]).T

        p1 = np.array([p.x for p in points_ref])
        p2 = np.array([p.y for p in points_ref])
        p3 = np.array([p.z for p in points_ref])
        points2 = np.vstack([p1, p2, p3]).T

        # Calculate centroids
        centroid1 = np.mean(points1, axis=0)
        centroid2 = np.mean(points2, axis=0)

        # Center the points
        points1_centered = points1 - centroid1
        points2_centered = points2 - centroid2

        # Calculate covariance matrix
        H = points1_centered.T @ points2_centered

        # SVD
        U, S, Vt = np.linalg.svd(H)

        # Calculate rotation matrix
        R = Vt.T @ U.T

        # Ensure right-handed coordinate system
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        # Calculate translation
        t = centroid2 - R @ centroid1

        return R, t

    def _kabsch_align(
        self,
        P: np.ndarray,
        Q: np.ndarray,
        weights: Optional[List[float]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Kabsch algorithm for optimal rotation and translation.

        Args:
            P: Points to align (N x 3)
            Q: Reference points (N x 3)
            weights: Optional point weights

        Returns:
            Tuple of (rotation matrix, translation vector)
        """
        if weights is None:
            weights = np.ones(len(P))
        weights = np.array(weights)

        # Center the points
        centroid_P = np.average(P, axis=0, weights=weights)
        centroid_Q = np.average(Q, axis=0, weights=weights)

        P_centered = P - centroid_P
        Q_centered = Q - centroid_Q

        # Calculate covariance matrix
        C = np.dot(P_centered.T, weights[:, None] * Q_centered)

        # SVD
        V, S, W = np.linalg.svd(C)

        # Ensure right-handed coordinate system
        d = np.linalg.det(np.dot(W, V.T))
        if d < 0:
            W[-1] *= -1

        # Calculate rotation matrix
        rotation = np.dot(W, V.T)

        # Calculate translation
        translation = centroid_Q - np.dot(centroid_P, rotation)

        return rotation, translation

    def _weighted_rmsd(
        self,
        P: np.ndarray,
        Q: np.ndarray,
        weights: Optional[List[float]] = None,
    ) -> float:
        """Calculate weighted RMSD between point sets.

        Args:
            P: First set of points
            Q: Second set of points
            weights: Optional point weights

        Returns:
            Weighted RMSD value
        """
        if weights is None:
            weights = np.ones(len(P))
        weights = np.array(weights)

        # Calculate squared differences
        diff_sq = np.sum((P - Q) ** 2, axis=1)

        # Calculate weighted average
        weighted_msd = np.average(diff_sq, weights=weights)

        return np.sqrt(weighted_msd)

    def _calculate_feature_match(
        self,
        mol: Chem.Mol,
        points: List[Tuple[float, float, float]],
        feature_types: List[str],
        weights: List[float],
    ) -> float:
        """Calculate how well molecule matches pharmacophore features.

        Args:
            mol: Molecule to check
            points: Feature point coordinates
            feature_types: Feature types for points
            weights: Feature weights

        Returns:
            Feature match score (0-1)
        """
        try:
            if not mol.GetNumConformers():
                return 0.0

            conf = mol.GetConformer()
            total_score = 0.0
            total_weight = sum(weights)

            for point, feat_type, weight in zip(points, feature_types, weights):
                # Find closest atom
                min_dist = float("inf")
                for i in range(mol.GetNumAtoms()):
                    pos = conf.GetAtomPosition(i)
                    dist = Point3D(*point).Distance(pos)
                    min_dist = min(min_dist, dist)

                # Score based on distance
                if min_dist < 2.0:  # Distance threshold
                    score = weight * np.exp(-min_dist)
                    total_score += score

            return total_score / total_weight if total_weight > 0 else 0.0

        except Exception as e:
            self.logger.error(f"Error calculating feature match: {str(e)}")
            return 0.0
