"""Pharmacophore feature generation and analysis."""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures

logger = logging.getLogger(__name__)


class PharmacophoreGenerator:
    """Generate and analyze pharmacophore features from molecules."""

    def __init__(
        self,
        feature_factory: Optional[str] = None,
        include_3d: bool = True,
        use_conformers: bool = True,
    ):
        """Initialize pharmacophore generator.

        Args:
            feature_factory: Path to feature factory definition file
            include_3d: Whether to include 3D pharmacophore features
            use_conformers: Whether to generate conformers for 3D features
        """
        self.include_3d = include_3d
        self.use_conformers = use_conformers

        # Initialize feature factory
        if feature_factory:
            self.factory = ChemicalFeatures.BuildFeatureFactoryFromFile(feature_factory)
        else:
            self.factory = ChemicalFeatures.BuildFeatureFactory()

        # Define standard feature types
        self.feature_types = {
            "Donor": "H-bond donor",
            "Acceptor": "H-bond acceptor",
            "Aromatic": "Aromatic ring",
            "Hydrophobe": "Hydrophobic region",
            "LumpedHydrophobe": "Merged hydrophobic regions",
            "PosIonizable": "Positive ionizable",
            "NegIonizable": "Negative ionizable",
        }

    def generate(self, mols: List[Chem.Mol]) -> Dict[str, np.ndarray]:
        """Generate pharmacophore features for molecules.

        Args:
            mols: List of RDKit molecules

        Returns:
            Dictionary mapping feature types to feature arrays
        """
        try:
            features = {}

            for mol in mols:
                # Generate 2D features
                mol_features = self._get_2d_features(mol)

                # Generate 3D features if requested
                if self.include_3d:
                    mol_3d_features = self._get_3d_features(mol)
                    mol_features.update(mol_3d_features)

                # Add to overall features
                for feat_type, feat_array in mol_features.items():
                    if feat_type not in features:
                        features[feat_type] = []
                    features[feat_type].append(feat_array)

            # Convert lists to arrays
            for feat_type in features:
                features[feat_type] = np.array(features[feat_type])

            return features

        except Exception as e:
            logger.error(f"Error generating pharmacophore features: {str(e)}")
            return {}

    def _get_2d_features(self, mol: Chem.Mol) -> Dict[str, np.ndarray]:
        """Get 2D pharmacophore features.

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary mapping feature types to feature arrays
        """
        try:
            features = {}

            # Get chemical features
            feats = self.factory.GetFeaturesForMol(mol)

            # Group by feature type
            for feat_type in self.feature_types:
                type_feats = []
                for feat in feats:
                    if feat.GetFamily() == feat_type:
                        # Get feature properties
                        atoms = feat.GetAtomIds()
                        pos = feat.GetPos()
                        props = {"atoms": atoms, "position": pos, "type": feat_type}
                        type_feats.append(props)

                if type_feats:
                    features[feat_type] = np.array(type_feats)

            return features

        except Exception as e:
            logger.error(f"Error getting 2D features: {str(e)}")
            return {}

    def _get_3d_features(self, mol: Chem.Mol) -> Dict[str, np.ndarray]:
        """Get 3D pharmacophore features.

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary mapping feature types to feature arrays
        """
        try:
            features = {}

            # Generate conformer if needed
            if self.use_conformers and not mol.GetNumConformers():
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)

            if mol.GetNumConformers():
                # Get 3D features
                conf = mol.GetConformer()

                for feat_type in self.feature_types:
                    type_feats = []

                    # Get features of this type
                    feats = self.factory.GetFeaturesForMol(mol)
                    for feat in feats:
                        if feat.GetFamily() == feat_type:
                            # Get 3D properties
                            atoms = feat.GetAtomIds()
                            pos = feat.GetPos()
                            normal = self._get_feature_normal(mol, atoms, conf)

                            props = {"atoms": atoms, "position": pos, "normal": normal, "type": feat_type}
                            type_feats.append(props)

                    if type_feats:
                        features[f"{feat_type}_3D"] = np.array(type_feats)

            return features

        except Exception as e:
            logger.error(f"Error getting 3D features: {str(e)}")
            return {}

    def _get_feature_normal(self, mol: Chem.Mol, atoms: List[int], conf: Chem.Conformer) -> np.ndarray:
        """Calculate normal vector for a pharmacophore feature.

        Args:
            mol: RDKit molecule
            atoms: Atom indices for feature
            conf: Molecule conformer

        Returns:
            Normal vector as numpy array
        """
        try:
            if len(atoms) < 3:
                return np.zeros(3)

            # Get coordinates
            coords = []
            for idx in atoms[:3]:
                pos = conf.GetAtomPosition(idx)
                coords.append([pos.x, pos.y, pos.z])
            coords = np.array(coords)

            # Calculate normal
            v1 = coords[1] - coords[0]
            v2 = coords[2] - coords[0]
            normal = np.cross(v1, v2)

            # Normalize
            norm = np.linalg.norm(normal)
            if norm > 0:
                normal = normal / norm

            return normal

        except Exception as e:
            logger.error(f"Error calculating feature normal: {str(e)}")
            return np.zeros(3)

    def get_info(self) -> Dict[str, Dict]:
        """Get information about available feature types."""
        info = {}

        # 2D features
        for feat_type, desc in self.feature_types.items():
            info[feat_type] = {"type": "2D", "description": desc}

        # 3D features
        if self.include_3d:
            for feat_type, desc in self.feature_types.items():
                info[f"{feat_type}_3D"] = {"type": "3D", "description": f"3D {desc.lower()}"}

        return info
