"""Enhanced molecular fingerprint generation with ML integration."""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem,
    MACCSkeys,
    rdFingerprintGenerator,
    rdMolDescriptors,
    rdReducedGraphs,
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

from .base import MLProcessor
from ..descriptors import MolecularDescriptors


class FingerprintGenerator(MLProcessor):
    """Generate and analyze molecular fingerprints with ML capabilities."""

    # Enhanced fingerprint type definitions
    FINGERPRINT_TYPES = {
        "morgan": {
            "description": "Morgan (ECFP) circular fingerprints",
            "generator": AllChem.GetMorganFingerprintAsBitVect,
            "params": {
                "radius": 2,
                "nBits": 2048,
                "useChirality": True,
                "useBondTypes": True,
                "useFeatures": False,
            },
            "ml_importance": 0.9,  # Relative importance for ML models
        },
        "morgan_feat": {
            "description": "Morgan feature-based fingerprints",
            "generator": AllChem.GetMorganFingerprintAsBitVect,
            "params": {
                "radius": 2,
                "nBits": 2048,
                "useChirality": True,
                "useBondTypes": True,
                "useFeatures": True,
            },
            "ml_importance": 0.85,
        },
        "rdkit": {
            "description": "RDKit topological fingerprints",
            "generator": Chem.RDKFingerprint,
            "params": {
                "minPath": 1,
                "maxPath": 7,
                "fpSize": 2048,
                "nBitsPerHash": 2,
                "useHs": True,
                "tgtDensity": 0.3,
            },
            "ml_importance": 0.7,
        },
        "maccs": {
            "description": "MACCS structural keys",
            "generator": rdMolDescriptors.GetMACCSKeysFingerprint,
            "params": {},
            "ml_importance": 0.6,
        },
        "pattern": {
            "description": "Substructure pattern fingerprints",
            "generator": rdMolDescriptors.PatternFingerprint,
            "params": {
                "fpSize": 2048,
                "tautomerFingerprints": True,
            },
            "ml_importance": 0.5,
        },
        "layered": {
            "description": "Layered fingerprint",
            "generator": Chem.LayeredFingerprint,
            "params": {"fpSize": 2048},
            "ml_importance": 0.4,
        },
        "reduced": {
            "description": "Reduced graph fingerprint",
            "generator": rdReducedGraphs.GetErGFingerprint,
            "params": {"fuzzIncrement": 0.3, "sanitize": True},
            "ml_importance": 0.3,
        },
        "topological_torsion": {
            "description": "Topological torsion fingerprints",
            "generator": rdMolDescriptors.GetTopologicalTorsionFingerprintAsIntVect,
            "params": {
                "targetSize": 4,
                "includeChirality": True,
            },
            "ml_importance": 0.6,
        },
        "atom_pair": {
            "description": "Atom pair fingerprints",
            "generator": rdMolDescriptors.GetAtomPairFingerprintAsIntVect,
            "params": {
                "minLength": 1,
                "maxLength": 30,
                "includeChirality": True,
            },
            "ml_importance": 0.55,
        },
    }

    def __init__(
        self,
        fingerprint_types: Optional[List[str]] = None,
        params: Optional[Dict[str, Dict]] = None,
        use_ml: bool = True,
        descriptor_generator: Optional[MolecularDescriptors] = None,
    ):
        """Initialize enhanced fingerprint generator.

        Args:
            fingerprint_types: Types of fingerprints to generate
            params: Custom parameters for fingerprint types
            use_ml: Whether to use ML enhancements
            descriptor_generator: Optional molecular descriptor generator
        """
        super().__init__()
        self.fingerprint_types = fingerprint_types or [
            "morgan",
            "morgan_feat",
            "rdkit",
            "maccs",
        ]
        self.params = params or {}
        self.use_ml = use_ml
        self.descriptor_generator = descriptor_generator

        # Initialize components
        self._setup_generators()
        if use_ml:
            self._setup_ml_components()

    def _setup_generators(self):
        """Set up fingerprint generators with enhanced error handling."""
        self.generators = {}
        for fp_type in self.fingerprint_types:
            if fp_type in self.FINGERPRINT_TYPES:
                # Merge default and custom parameters
                params = self.FINGERPRINT_TYPES[fp_type]["params"].copy()
                if fp_type in self.params:
                    params.update(self.params[fp_type])
                self.generators[fp_type] = (
                    self.FINGERPRINT_TYPES[fp_type]["generator"],
                    params,
                )

    def _setup_ml_components(self):
        """Initialize ML components for fingerprint analysis."""
        self.scaler = StandardScaler()
        self.classifier = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
        )
        self.feature_importance = {}

    def generate(
        self,
        mols: List[Chem.Mol],
        fingerprint_types: Optional[List[str]] = None,
        include_descriptors: bool = True,
    ) -> Dict[str, np.ndarray]:
        """Generate enhanced fingerprints for molecules.

        Args:
            mols: List of RDKit molecules
            fingerprint_types: Types of fingerprints to generate
            include_descriptors: Whether to include molecular descriptors

        Returns:
            Dictionary mapping fingerprint types to feature arrays
        """
        if fingerprint_types is None:
            fingerprint_types = self.fingerprint_types

        fingerprints = {}
        try:
            # Generate fingerprints
            for fp_type in fingerprint_types:
                if fp_type in self.generators:
                    generator, params = self.generators[fp_type]
                    fp_list = []
                    for mol in mols:
                        try:
                            fp = generator(mol, **params)
                            # Convert to numpy array
                            if isinstance(
                                fp,
                                (
                                    DataStructs.ExplicitBitVect,
                                    DataStructs.SparseBitVect,
                                ),
                            ):
                                arr = np.zeros((0,), dtype=np.int8)
                                DataStructs.ConvertToNumpyArray(fp, arr)
                            else:
                                arr = np.array(list(fp))
                            fp_list.append(arr)
                        except Exception as e:
                            self.logger.debug(
                                f"Error generating {fp_type} fingerprint: {str(e)}"
                            )
                            size = (
                                params.get("nBits", 2048) if fp_type != "maccs" else 166
                            )
                            fp_list.append(np.zeros(size, dtype=np.int8))
                    fingerprints[fp_type] = np.vstack(fp_list)

            # Add molecular descriptors if requested
            if include_descriptors and self.descriptor_generator:
                try:
                    descriptors = self.descriptor_generator.calculate(mols)
                    if descriptors is not None:
                        fingerprints["descriptors"] = descriptors
                except Exception as e:
                    self.logger.error(f"Error calculating descriptors: {str(e)}")

            return fingerprints

        except Exception as e:
            self.logger.error(f"Error generating fingerprints: {str(e)}")
            return {}

    def train_classifier(
        self,
        mols: List[Chem.Mol],
        labels: List[int],
        fingerprint_types: Optional[List[str]] = None,
    ) -> float:
        """Train ML classifier on fingerprints.

        Args:
            mols: Training molecules
            labels: Binary activity labels
            fingerprint_types: Types of fingerprints to use

        Returns:
            ROC AUC score
        """
        if not self.use_ml:
            return 0.0

        try:
            # Generate fingerprints
            fps = self.generate(mols, fingerprint_types)
            if not fps:
                return 0.0

            # Combine fingerprints
            features = []
            for fp_type, fp_array in fps.items():
                if fp_type != "descriptors":
                    importance = self.FINGERPRINT_TYPES[fp_type]["ml_importance"]
                    features.append(fp_array * importance)
            if "descriptors" in fps:
                features.append(self.scaler.fit_transform(fps["descriptors"]))

            X = np.hstack(features)
            y = np.array(labels)

            # Train classifier
            self.classifier.fit(X, y)
            y_pred = self.classifier.predict_proba(X)[:, 1]
            score = roc_auc_score(y, y_pred)

            # Store feature importance
            importance = self.classifier.feature_importances_
            start = 0
            for fp_type, fp_array in fps.items():
                size = fp_array.shape[1]
                self.feature_importance[fp_type] = np.mean(
                    importance[start : start + size]
                )
                start += size

            return score

        except Exception as e:
            self.logger.error(f"Error training classifier: {str(e)}")
            return 0.0

    def predict_activity(
        self,
        mols: List[Chem.Mol],
        fingerprint_types: Optional[List[str]] = None,
    ) -> np.ndarray:
        """Predict activity probabilities for molecules.

        Args:
            mols: Molecules to predict
            fingerprint_types: Types of fingerprints to use

        Returns:
            Array of activity probabilities
        """
        if not self.use_ml:
            return np.array([])

        try:
            # Generate fingerprints
            fps = self.generate(mols, fingerprint_types)
            if not fps:
                return np.array([])

            # Combine fingerprints
            features = []
            for fp_type, fp_array in fps.items():
                if fp_type != "descriptors":
                    importance = self.FINGERPRINT_TYPES[fp_type]["ml_importance"]
                    features.append(fp_array * importance)
            if "descriptors" in fps:
                features.append(self.scaler.transform(fps["descriptors"]))

            X = np.hstack(features)
            return self.classifier.predict_proba(X)[:, 1]

        except Exception as e:
            self.logger.error(f"Error predicting activity: {str(e)}")
            return np.array([])

    def calculate_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        fingerprint_type: str = "morgan",
        metric: str = "tanimoto",
    ) -> Optional[float]:
        """Calculate similarity between molecules.

        Args:
            mol1: First molecule
            mol2: Second molecule
            fingerprint_type: Type of fingerprint to use
            metric: Similarity metric (tanimoto, dice, cosine)

        Returns:
            Similarity score
        """
        if fingerprint_type not in self.generators:
            return None

        try:
            generator, params = self.generators[fingerprint_type]
            fp1 = generator(mol1, **params)
            fp2 = generator(mol2, **params)

            if metric == "tanimoto":
                return DataStructs.TanimotoSimilarity(fp1, fp2)
            elif metric == "dice":
                return DataStructs.DiceSimilarity(fp1, fp2)
            elif metric == "cosine":
                return DataStructs.CosineSimilarity(fp1, fp2)
            else:
                return DataStructs.TanimotoSimilarity(fp1, fp2)

        except Exception as e:
            self.logger.error(f"Error calculating similarity: {str(e)}")
            return None

    def get_info(self) -> Dict[str, Dict]:
        """Get information about available fingerprint types."""
        info = {}
        for fp_type in self.fingerprint_types:
            if fp_type in self.FINGERPRINT_TYPES:
                info[fp_type] = {
                    "description": self.FINGERPRINT_TYPES[fp_type]["description"],
                    "params": self.generators[fp_type][1].copy(),
                    "ml_importance": self.FINGERPRINT_TYPES[fp_type]["ml_importance"],
                    "feature_importance": self.feature_importance.get(fp_type, 0.0),
                }
        return info
