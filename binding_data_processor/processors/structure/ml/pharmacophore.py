"""Pharmacophore feature detection and analysis."""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    ChemicalFeatures,
)

from .base import MLProcessor


class PharmacophoreGenerator(MLProcessor):
    """Generate pharmacophore features."""

    # Pharmacophore feature definitions
    FEATURE_TYPES = {
        "donor": {
            "description": "Hydrogen bond donor",
            "smarts": ["[!H0;N,O,S]"],
            "family": "Donor",
        },
        "acceptor": {
            "description": "Hydrogen bond acceptor",
            "smarts": ["[N,O,S;H0]"],
            "family": "Acceptor",
        },
        "aromatic": {
            "description": "Aromatic ring center",
            "smarts": ["a1aaaaa1", "a1aaaa1"],
            "family": "Aromatic",
        },
        "hydrophobic": {
            "description": "Hydrophobic group",
            "smarts": [
                "[C;!$(C=[O,N,S]);!$(C#[N,C]);!$([C@H]([NH2])[C](=O)[OH])][C;!$(C=[O,N,S]);!$(C#[N,C]);!$([C@H]([NH2])[C](=O)[OH])]",
                "c1ccccc1",
                "[CH3][CH3]",
                "[CH2][CH2][CH2]",
            ],
            "family": "Hydrophobe",
        },
        "positive": {
            "description": "Positive ionizable",
            "smarts": [
                "[NH2;!$(NC=O)]",
                "[NH;!$(NC=O)]",
                "[NH0;+;!$(NC=O)]",
                "[NH2][CX4]",
                "C(=N)N",
            ],
            "family": "PosIonizable",
        },
        "negative": {
            "description": "Negative ionizable",
            "smarts": [
                "[OH;$(O[C,S,P]=O)]",
                "[O-;$(O[C,S,P]=O)]",
                "[O;$(O[C,S,P]=O)]",
                "[CO2H]",
                "[SO3H]",
                "[PO3H]",
            ],
            "family": "NegIonizable",
        },
    }

    def __init__(
        self,
        feature_types: Optional[List[str]] = None,
        include_counts: bool = True,
        include_distances: bool = True,
        max_distance: float = 20.0,
    ):
        """Initialize pharmacophore generator.

        Args:
            feature_types: Types of features to detect
            include_counts: Include feature counts
            include_distances: Include feature distances
            max_distance: Maximum distance between features
        """
        super().__init__()
        self.feature_types = feature_types or list(self.FEATURE_TYPES.keys())
        self.include_counts = include_counts
        self.include_distances = include_distances
        self.max_distance = max_distance

        # Initialize feature factory
        self._setup_factory()

    def _setup_factory(self):
        """Set up pharmacophore feature factory."""
        try:
            params = []
            for feat_type in self.feature_types:
                if feat_type in self.FEATURE_TYPES:
                    for smarts in self.FEATURE_TYPES[feat_type]["smarts"]:
                        params.append(
                            (
                                self.FEATURE_TYPES[feat_type]["family"],
                                smarts,
                                feat_type,
                            )
                        )

            self.factory = ChemicalFeatures.BuildFeatureFactoryFromParams(params)
            self.logger.info("Initialized pharmacophore feature factory")

        except Exception as e:
            self.logger.error(f"Error setting up feature factory: {str(e)}")
            self.factory = None

    def generate(
        self,
        mols: List[Chem.Mol],
        feature_types: Optional[List[str]] = None,
    ) -> Dict[str, np.ndarray]:
        """Generate pharmacophore features for molecules.

        Args:
            mols: List of RDKit molecules
            feature_types: Types of features to generate

        Returns:
            Dictionary mapping feature types to feature arrays
        """
        if feature_types is None:
            feature_types = self.feature_types

        if self.factory is None:
            self.logger.error("Feature factory not initialized")
            return {}

        features = {}
        try:
            # Generate feature counts
            if self.include_counts:
                counts = []
                for mol in mols:
                    try:
                        mol_features = self.factory.GetFeaturesForMol(mol)
                        feat_counts = {feat_type: 0 for feat_type in feature_types}
                        for feat in mol_features:
                            feat_type = feat.GetFamily().lower()
                            if feat_type in feat_counts:
                                feat_counts[feat_type] += 1
                        counts.append(list(feat_counts.values()))
                    except Exception as e:
                        self.logger.debug(
                            f"Error getting features for molecule: {str(e)}"
                        )
                        counts.append([0] * len(feature_types))
                features["counts"] = np.array(counts)

            # Generate feature distances
            if self.include_distances:
                distances = []
                for mol in mols:
                    try:
                        mol_features = self.factory.GetFeaturesForMol(mol)
                        feat_positions = {}
                        for feat in mol_features:
                            feat_type = feat.GetFamily().lower()
                            if feat_type in feature_types:
                                if feat_type not in feat_positions:
                                    feat_positions[feat_type] = []
                                feat_positions[feat_type].append(feat.GetPos())

                        # Calculate minimum distances between feature types
                        dist_matrix = []
                        for i, type1 in enumerate(feature_types):
                            for j, type2 in enumerate(feature_types[i:], i):
                                min_dist = self.max_distance
                                if type1 in feat_positions and type2 in feat_positions:
                                    for pos1 in feat_positions[type1]:
                                        for pos2 in feat_positions[type2]:
                                            dist = np.linalg.norm(
                                                np.array(pos1) - np.array(pos2)
                                            )
                                            min_dist = min(min_dist, dist)
                                dist_matrix.append(min_dist)
                        distances.append(dist_matrix)
                    except Exception as e:
                        self.logger.debug(f"Error calculating distances: {str(e)}")
                        n_dists = (len(feature_types) * (len(feature_types) + 1)) // 2
                        distances.append([self.max_distance] * n_dists)
                features["distances"] = np.array(distances)

            return features

        except Exception as e:
            self.logger.error(f"Error generating pharmacophore features: {str(e)}")
            return {}

    def get_info(self) -> Dict[str, Dict]:
        """Get information about available features."""
        info = {}
        for feat_type in self.feature_types:
            if feat_type in self.FEATURE_TYPES:
                info[feat_type] = {
                    "description": self.FEATURE_TYPES[feat_type]["description"],
                    "smarts": self.FEATURE_TYPES[feat_type]["smarts"],
                    "family": self.FEATURE_TYPES[feat_type]["family"],
                }
        return info

    def get_feature_positions(
        self, mol: Chem.Mol, feature_type: str
    ) -> List[Tuple[float, float, float]]:
        """Get 3D positions of features of a given type.

        Args:
            mol: RDKit molecule
            feature_type: Type of feature

        Returns:
            List of (x, y, z) coordinates
        """
        if self.factory is None or feature_type not in self.FEATURE_TYPES:
            return []

        try:
            positions = []
            mol_features = self.factory.GetFeaturesForMol(mol)
            for feat in mol_features:
                if feat.GetFamily().lower() == feature_type:
                    positions.append(feat.GetPos())
            return positions
        except Exception as e:
            self.logger.error(f"Error getting feature positions: {str(e)}")
            return []
