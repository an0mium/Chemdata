"""Mixins for sharing functionality between structure processing modules."""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDepictor, rdFMCS

logger = logging.getLogger(__name__)


class StructureValidationMixin:
    """Mixin for structure validation and standardization."""

    def validate_structure(self, smiles: str) -> Optional[Chem.Mol]:
        """Validate and standardize chemical structure.

        Args:
            smiles: SMILES string to validate

        Returns:
            Validated RDKit molecule or None if invalid
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            mol = Chem.AddHs(mol)
            try:
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception as e:
                logger.warning(f"3D generation failed: {str(e)}")

            Chem.SanitizeMol(mol)
            return mol

        except Exception as e:
            logger.error(f"Structure validation error: {str(e)}")
            return None


class ConformerGenerationMixin:
    """Mixin for 3D conformer generation and optimization."""

    def generate_conformers(
        self,
        mol: Chem.Mol,
        n_confs: int = 10,
        optimize: bool = True,
        random_seed: int = 42,
    ) -> Optional[Chem.Mol]:
        """Generate multiple conformers.

        Args:
            mol: Input molecule
            n_confs: Number of conformers to generate
            optimize: Whether to optimize conformers
            random_seed: Random seed for reproducibility

        Returns:
            Molecule with generated conformers
        """
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMultipleConfs(
                mol,
                numConfs=n_confs,
                randomSeed=random_seed,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True,
            )

            if optimize:
                for conf_id in range(mol.GetNumConformers()):
                    AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)

            return mol

        except Exception as e:
            logger.error(f"Conformer generation error: {str(e)}")
            return None


class PharmacophoreFeatureMixin:
    """Mixin for pharmacophore feature detection."""

    FEATURE_DEFINITIONS = {
        "donor": "[!#6;!H0]-[!#6]",
        "acceptor": "[$([!#6;+0]);!$([F,Cl,Br,I]);!$([o,s,nX3]);!$([Nv5,Pv5,Sv4,Sv6])]",
        "aromatic": "a",
        "hydrophobe": "[#6+0!$([#6](~[#7,#8,#9])~[#7,#8,#9])]",
        "positive": "[+,+2,+3,+4]",
        "negative": "[-,-2,-3,-4]",
        "hbd": "[N,O,S;H1,H2]-[!$(*=[O,N,P,S])]",
        "hba": "[$([O,S;H0;v2]),$([O,S;-])]",
    }

    def detect_features(
        self,
        mol: Chem.Mol,
        include_vectors: bool = True,
    ) -> List[Dict]:
        """Detect pharmacophore features.

        Args:
            mol: Input molecule
            include_vectors: Whether to include feature vectors

        Returns:
            List of detected features
        """
        try:
            features = []
            for name, smarts in self.FEATURE_DEFINITIONS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern is None:
                    continue
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    feature = {
                        "type": name,
                        "atoms": match,
                        "position": self._get_feature_position(mol, match),
                    }
                    if include_vectors and name in ["hbd", "hba", "donor", "acceptor"]:
                        feature["vector"] = self._get_feature_vector(mol, match, name)
                    features.append(feature)

            # Add ring-based features
            rings = mol.GetRingInfo().AtomRings()
            for ring in rings:
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    features.append(
                        {
                            "type": "aromatic",
                            "atoms": ring,
                            "position": self._get_feature_position(mol, ring),
                        }
                    )

            return features

        except Exception as e:
            logger.error(f"Feature detection error: {str(e)}")
            return []

    def _get_feature_position(
        self,
        mol: Chem.Mol,
        atoms: Tuple[int, ...],
    ) -> np.ndarray:
        """Get 3D position of feature center."""
        if not mol.GetNumConformers():
            return np.zeros(3)

        positions = []
        conf = mol.GetConformer()
        for atom_idx in atoms:
            pos = conf.GetAtomPosition(atom_idx)
            positions.append([pos.x, pos.y, pos.z])

        return np.mean(positions, axis=0)

    def _get_feature_vector(
        self,
        mol: Chem.Mol,
        atoms: Tuple[int, ...],
        feature_type: str,
    ) -> Optional[np.ndarray]:
        """Get directional vector for feature."""
        try:
            if not mol.GetNumConformers():
                return None

            conf = mol.GetConformer()

            if feature_type in ["hbd", "donor"]:
                # Vector from heavy atom to H
                heavy_idx = atoms[0]
                h_idx = atoms[1]
                start = conf.GetAtomPosition(heavy_idx)
                end = conf.GetAtomPosition(h_idx)
            elif feature_type in ["hba", "acceptor"]:
                # Vector from neighbor to acceptor
                acc_idx = atoms[0]
                nbr_idx = mol.GetAtomWithIdx(acc_idx).GetNeighbors()[0].GetIdx()
                start = conf.GetAtomPosition(nbr_idx)
                end = conf.GetAtomPosition(acc_idx)
            else:
                return None

            vector = np.array([end.x - start.x, end.y - start.y, end.z - start.z])
            norm = np.linalg.norm(vector)
            if norm > 0:
                return vector / norm
            return None

        except Exception:
            return None


class StructureAlignmentMixin:
    """Mixin for structure alignment and RMSD calculation."""

    def align_structures(
        self,
        ref_mol: Chem.Mol,
        probe_mol: Chem.Mol,
        ref_features: Optional[List[Dict]] = None,
        probe_features: Optional[List[Dict]] = None,
    ) -> Tuple[float, Chem.Mol]:
        """Align structures using features or MCS.

        Args:
            ref_mol: Reference molecule
            probe_mol: Probe molecule to align
            ref_features: Optional reference features for feature-based alignment
            probe_features: Optional probe features for feature-based alignment

        Returns:
            Tuple of (RMSD, aligned molecule)
        """
        try:
            if ref_features and probe_features:
                return self._align_by_features(ref_mol, probe_mol, ref_features, probe_features)
            return self._align_by_mcs(ref_mol, probe_mol)

        except Exception as e:
            logger.error(f"Alignment error: {str(e)}")
            return float("inf"), probe_mol

    def _align_by_features(
        self,
        ref_mol: Chem.Mol,
        probe_mol: Chem.Mol,
        ref_features: List[Dict],
        probe_features: List[Dict],
    ) -> Tuple[float, Chem.Mol]:
        """Align using pharmacophore features."""
        try:
            if not ref_mol.GetNumConformers() or not probe_mol.GetNumConformers():
                return float("inf"), probe_mol

            # Match features by type
            ref_coords = []
            probe_coords = []
            for ref_feat in ref_features:
                for probe_feat in probe_features:
                    if ref_feat["type"] == probe_feat["type"]:
                        ref_coords.append(ref_feat["position"])
                        probe_coords.append(probe_feat["position"])

            if len(ref_coords) < 3:
                return self._align_by_mcs(ref_mol, probe_mol)

            # Calculate transformation
            ref_coords = np.array(ref_coords)
            probe_coords = np.array(probe_coords)
            rotation, translation = self._get_transformation(probe_coords, ref_coords)

            # Apply transformation
            conf = probe_mol.GetConformer()
            for i in range(probe_mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                coords = np.array([pos.x, pos.y, pos.z])
                new_coords = rotation.dot(coords) + translation
                conf.SetAtomPosition(i, new_coords)

            rmsd = self._calculate_rmsd(ref_coords, probe_coords)
            return rmsd, probe_mol

        except Exception as e:
            logger.error(f"Feature alignment error: {str(e)}")
            return float("inf"), probe_mol

    def _align_by_mcs(
        self,
        ref_mol: Chem.Mol,
        probe_mol: Chem.Mol,
    ) -> Tuple[float, Chem.Mol]:
        """Align using maximum common substructure."""
        try:
            # Find MCS
            mcs = rdFMCS.FindMCS([ref_mol, probe_mol])
            if mcs is None:
                return float("inf"), probe_mol

            # Get matching atoms
            pattern = Chem.MolFromSmarts(mcs.smartsString)
            ref_match = ref_mol.GetSubstructMatch(pattern)
            probe_match = probe_mol.GetSubstructMatch(pattern)

            if not ref_match or not probe_match:
                return float("inf"), probe_mol

            # Get coordinates
            ref_conf = ref_mol.GetConformer()
            probe_conf = probe_mol.GetConformer()
            ref_coords = []
            probe_coords = []
            for ref_idx, probe_idx in zip(ref_match, probe_match):
                ref_pos = ref_conf.GetAtomPosition(ref_idx)
                probe_pos = probe_conf.GetAtomPosition(probe_idx)
                ref_coords.append([ref_pos.x, ref_pos.y, ref_pos.z])
                probe_coords.append([probe_pos.x, probe_pos.y, probe_pos.z])

            # Calculate and apply transformation
            ref_coords = np.array(ref_coords)
            probe_coords = np.array(probe_coords)
            rotation, translation = self._get_transformation(probe_coords, ref_coords)

            for i in range(probe_mol.GetNumAtoms()):
                pos = probe_conf.GetAtomPosition(i)
                coords = np.array([pos.x, pos.y, pos.z])
                new_coords = rotation.dot(coords) + translation
                probe_conf.SetAtomPosition(i, new_coords)

            rmsd = self._calculate_rmsd(ref_coords, probe_coords)
            return rmsd, probe_mol

        except Exception as e:
            logger.error(f"MCS alignment error: {str(e)}")
            return float("inf"), probe_mol

    def _get_transformation(
        self,
        coords1: np.ndarray,
        coords2: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Get rotation and translation to align coordinates."""
        # Center coordinates
        center1 = np.mean(coords1, axis=0)
        center2 = np.mean(coords2, axis=0)
        coords1_centered = coords1 - center1
        coords2_centered = coords2 - center2

        # Calculate rotation matrix
        covariance = coords1_centered.T.dot(coords2_centered)
        U, S, Vt = np.linalg.svd(covariance)
        rotation = Vt.T.dot(U.T)

        # Ensure right-handed coordinate system
        if np.linalg.det(rotation) < 0:
            Vt[-1] *= -1
            rotation = Vt.T.dot(U.T)

        # Calculate translation
        translation = center2 - rotation.dot(center1)

        return rotation, translation

    def _calculate_rmsd(
        self,
        coords1: np.ndarray,
        coords2: np.ndarray,
    ) -> float:
        """Calculate RMSD between coordinate sets."""
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff * diff, axis=1)))


class StructureVisualizationMixin:
    """Mixin for structure visualization."""

    def visualize_structure(
        self,
        mol: Chem.Mol,
        highlight_atoms: Optional[List[int]] = None,
        highlight_bonds: Optional[List[int]] = None,
        highlight_colors: Optional[Dict[int, Tuple[float, float, float]]] = None,
        size: Tuple[int, int] = (400, 400),
    ) -> Optional[str]:
        """Generate 2D depiction.

        Args:
            mol: Molecule to visualize
            highlight_atoms: Optional list of atoms to highlight
            highlight_bonds: Optional list of bonds to highlight
            highlight_colors: Optional color mapping for highlighted atoms
            size: Image size (width, height)

        Returns:
            SVG string of visualization
        """
        try:
            rdDepictor.Compute2DCoords(mol)
            drawer = Draw.rdDepictor.MolDraw2DSVG(size[0], size[1])
            drawer.drawOptions().addStereoAnnotation = True
            drawer.drawOptions().addAtomIndices = False

            drawer.DrawMolecule(
                mol,
                highlightAtoms=highlight_atoms,
                highlightBonds=highlight_bonds,
                highlightAtomColors=highlight_colors,
            )
            drawer.FinishDrawing()
            return drawer.GetDrawingText()

        except Exception as e:
            logger.error(f"Visualization error: {str(e)}")
            return None

    def visualize_features(
        self,
        mol: Chem.Mol,
        features: List[Dict],
        size: Tuple[int, int] = (400, 400),
    ) -> Optional[str]:
        """Visualize pharmacophore features.

        Args:
            mol: Molecule to visualize
            features: Features to highlight
            size: Image size (width, height)

        Returns:
            SVG string of visualization
        """
        try:
            # Color scheme for features
            colors = {
                "donor": (0, 1, 0),  # Green
                "acceptor": (1, 0, 0),  # Red
                "aromatic": (1, 1, 0),  # Yellow
                "hydrophobe": (0, 0, 1),  # Blue
                "positive": (1, 0, 1),  # Magenta
                "negative": (0, 1, 1),  # Cyan
                "hbd": (0, 1, 0),  # Green
                "hba": (1, 0, 0),  # Red
            }

            # Collect atoms to highlight
            highlight_atoms = {}
            for feature in features:
                color = colors.get(feature["type"], (0.5, 0.5, 0.5))
                for atom_idx in feature["atoms"]:
                    highlight_atoms[atom_idx] = color

            return self.visualize_structure(
                mol,
                highlight_atoms=list(highlight_atoms.keys()),
                highlight_colors=highlight_atoms,
                size=size,
            )

        except Exception as e:
            logger.error(f"Feature visualization error: {str(e)}")
            return None
