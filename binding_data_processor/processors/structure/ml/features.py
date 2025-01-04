"""Molecular feature extraction and processing with ML capabilities."""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    MACCSkeys,
    rdDecomposition,
    rdMolDescriptors,
    rdReducedGraphs,
)
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data

from .base import MLProcessor
from .descriptors import DescriptorGenerator
from .fingerprints import FingerprintGenerator
from .graphs import GraphFeatureGenerator
from .pharmacophore import PharmacophoreGenerator
from .utils import mol_to_graph, compute_fingerprints


class EnhancedFeatureExtractor(MLProcessor):
    """Extract and process molecular features with ML capabilities."""

    def __init__(
        self,
        feature_types: Optional[List[str]] = None,
        fingerprint_radius: int = 2,
        include_3d: bool = False,
        normalize: bool = True,
        use_gpu: bool = True,
    ):
        """Initialize feature extractor.

        Args:
            feature_types: Types of features to compute
            fingerprint_radius: Radius for Morgan fingerprints
            include_3d: Whether to include 3D descriptors
            normalize: Whether to normalize features
            use_gpu: Whether to use GPU acceleration
        """
        super().__init__()
        self.feature_types = feature_types or [
            "morgan",
            "maccs",
            "rdkit",
            "descriptors",
            "pharmacophore",
            "graph",
        ]
        self.fingerprint_radius = fingerprint_radius
        self.include_3d = include_3d
        self.normalize = normalize
        self.use_gpu = use_gpu and torch.cuda.is_available()

        # Initialize feature generators
        self.fingerprint_gen = FingerprintGenerator(
            radius=fingerprint_radius,
            use_features=True,
            use_chirality=True,
        )
        self.descriptor_gen = DescriptorGenerator()
        self.pharmacophore_gen = PharmacophoreGenerator()
        self.graph_gen = GraphFeatureGenerator()

        # Initialize normalizers
        self.scalers = {}

        # Load pre-trained models if available
        self._load_models()

    def _load_models(self):
        """Load pre-trained feature extraction models."""
        try:
            model_dir = Path(__file__).parent / "models"
            if model_dir.exists():
                # Load feature extraction models
                for model_path in model_dir.glob("*.pt"):
                    model_name = model_path.stem
                    if model_name in self.feature_types:
                        model = torch.load(model_path, map_location=self.device)
                        setattr(self, f"{model_name}_model", model)
                        self.logger.info(f"Loaded {model_name} model")
        except Exception as e:
            self.logger.error(f"Error loading models: {str(e)}")

    def extract_features(
        self,
        mols: List[Chem.Mol],
        feature_types: Optional[List[str]] = None,
    ) -> Dict[str, Union[np.ndarray, Data]]:
        """Extract features from molecules.

        Args:
            mols: List of RDKit molecules
            feature_types: Types of features to extract

        Returns:
            Dictionary mapping feature types to feature arrays
        """
        if feature_types is None:
            feature_types = self.feature_types

        features = {}
        try:
            # Extract basic features
            if "fingerprints" in feature_types:
                features["fingerprints"] = self.fingerprint_gen.generate(mols)

            if "descriptors" in feature_types:
                features["descriptors"] = self.descriptor_gen.generate(mols)

            if "pharmacophore" in feature_types:
                features["pharmacophore"] = self.pharmacophore_gen.generate(mols)

            if "graph" in feature_types:
                features["graph"] = self.graph_gen.generate(mols)

            # Extract learned features if models available
            for feat_type in feature_types:
                if hasattr(self, f"{feat_type}_model"):
                    model = getattr(self, f"{feat_type}_model")
                    features[feat_type] = self._extract_learned_features(mols, model)

            # Normalize features if requested
            if self.normalize:
                features = self._normalize_features(features)

            return features

        except Exception as e:
            self.logger.error(f"Error extracting features: {str(e)}")
            return {}

    def _extract_learned_features(
        self,
        mols: List[Chem.Mol],
        model: nn.Module,
    ) -> np.ndarray:
        """Extract features using a learned model.

        Args:
            mols: List of RDKit molecules
            model: PyTorch model for feature extraction

        Returns:
            Array of learned features
        """
        try:
            model.eval()
            features = []

            with torch.no_grad():
                for mol in mols:
                    # Convert molecule to graph
                    graph = mol_to_graph(mol)
                    if self.use_gpu:
                        graph = graph.to(self.device)

                    # Extract features
                    feat = model(graph)
                    features.append(feat.cpu().numpy())

            return np.vstack(features)

        except Exception as e:
            self.logger.error(f"Error extracting learned features: {str(e)}")
            return np.array([])

    def _normalize_features(
        self,
        features: Dict[str, Union[np.ndarray, Data]],
    ) -> Dict[str, Union[np.ndarray, Data]]:
        """Normalize feature arrays.

        Args:
            features: Dictionary of features

        Returns:
            Dictionary of normalized features
        """
        try:
            normalized = {}
            for feat_type, feat_array in features.items():
                if isinstance(feat_array, np.ndarray):
                    if feat_type not in self.scalers:
                        self.scalers[feat_type] = StandardScaler()
                        self.scalers[feat_type].fit(feat_array)
                    normalized[feat_type] = self.scalers[feat_type].transform(
                        feat_array
                    )
                else:
                    normalized[feat_type] = feat_array
            return normalized

        except Exception as e:
            self.logger.error(f"Error normalizing features: {str(e)}")
            return features

    def save_scalers(self, path: Union[str, Path]):
        """Save feature scalers."""
        try:
            import joblib

            scaler_path = Path(path)
            scaler_path.mkdir(parents=True, exist_ok=True)

            for feat_type, scaler in self.scalers.items():
                joblib.dump(scaler, scaler_path / f"{feat_type}_scaler.joblib")
            self.logger.info(f"Saved scalers to {path}")

        except Exception as e:
            self.logger.error(f"Error saving scalers: {str(e)}")

    def load_scalers(self, path: Union[str, Path]):
        """Load feature scalers."""
        try:
            import joblib

            scaler_path = Path(path)

            for scaler_file in scaler_path.glob("*_scaler.joblib"):
                feat_type = scaler_file.stem.replace("_scaler", "")
                self.scalers[feat_type] = joblib.load(scaler_file)
            self.logger.info(f"Loaded scalers from {path}")

        except Exception as e:
            self.logger.error(f"Error loading scalers: {str(e)}")

    def get_feature_info(self) -> Dict[str, Dict]:
        """Get information about available features."""
        info = {}

        # Basic features
        info["fingerprints"] = self.fingerprint_gen.get_info()
        info["descriptors"] = self.descriptor_gen.get_info()
        info["pharmacophore"] = self.pharmacophore_gen.get_info()
        info["graph"] = self.graph_gen.get_info()

        # Learned features
        for feat_type in self.feature_types:
            if hasattr(self, f"{feat_type}_model"):
                model = getattr(self, f"{feat_type}_model")
                info[feat_type] = {
                    "type": "learned",
                    "model": model.__class__.__name__,
                    "output_dim": model.output_dim,
                }

        return info
