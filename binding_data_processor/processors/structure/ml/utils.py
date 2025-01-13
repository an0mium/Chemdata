"""ML utilities for molecular data processing.

This module provides:
1. Comprehensive molecular featurization
2. Advanced data preprocessing
3. Model evaluation and metrics
4. Uncertainty estimation
5. Integration with web app and database
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union, Any

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Fragments,
    MACCSkeys,
    rdMolDescriptors,
    rdRGroupDecomposition,
)
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    precision_score,
    r2_score,
    recall_score,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Batch, Data
from torch_geometric.nn import global_add_pool, global_max_pool, global_mean_pool

from ....models.core import CompoundData
from ..descriptors import DescriptorCalculator
from ..pharmacophore import PharmacophoreGenerator
from ..similarity import SimilarityProcessor
from .core_utils import mol_to_graph, compute_fingerprints


class DataPreprocessor:
    """Data preprocessing utilities."""

    def __init__(self):
        """Initialize preprocessor."""
        self.logger = logging.getLogger(__name__)
        self.scaler = StandardScaler()

    def normalize(self, data: np.ndarray, fit: bool = True) -> np.ndarray:
        """Normalize data using StandardScaler.

        Args:
            data: Input data array
            fit: Whether to fit scaler on data

        Returns:
            Normalized data array
        """
        if fit:
            return self.scaler.fit_transform(data)
        return self.scaler.transform(data)

    def split_data(self, X: np.ndarray, y: np.ndarray, test_size: float = 0.2, random_state: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Split data into train and test sets.

        Args:
            X: Feature matrix
            y: Target values
            test_size: Fraction of data to use for testing
            random_state: Random seed

        Returns:
            X_train, X_test, y_train, y_test arrays
        """
        from sklearn.model_selection import train_test_split

        return train_test_split(X, y, test_size=test_size, random_state=random_state)


class FeatureExtractor:
    """Feature extraction utilities."""

    def __init__(self):
        """Initialize extractor."""
        self.logger = logging.getLogger(__name__)
        self.featurizer = MolecularFeaturizer()

    def extract_features(
        self,
        mols: List[Chem.Mol],
        feature_types: Optional[List[str]] = None,
    ) -> Dict[str, np.ndarray]:
        """Extract features from molecules.

        Args:
            mols: List of RDKit molecules
            feature_types: Types of features to extract

        Returns:
            Dictionary of feature arrays
        """
        features = {}
        for mol in mols:
            mol_features = self.featurizer.generate_features(mol, feature_types=feature_types)
            for k, v in mol_features.items():
                if k not in features:
                    features[k] = []
                features[k].append(v)

        # Convert lists to arrays
        for k in features:
            features[k] = np.array(features[k])

        return features


class ModelEvaluator:
    """Model evaluation utilities."""

    def __init__(self):
        """Initialize evaluator."""
        self.logger = logging.getLogger(__name__)

    def evaluate_classifier(
        self,
        y_true: np.ndarray,
        y_pred: np.ndarray,
        y_prob: Optional[np.ndarray] = None,
    ) -> Dict[str, float]:
        """Evaluate classifier performance.

        Args:
            y_true: True labels
            y_pred: Predicted labels
            y_prob: Predicted probabilities

        Returns:
            Dictionary of metrics
        """
        metrics = {
            "accuracy": accuracy_score(y_true, y_pred),
            "precision": precision_score(y_true, y_pred, average="weighted"),
            "recall": recall_score(y_true, y_pred, average="weighted"),
            "f1": f1_score(y_true, y_pred, average="weighted"),
        }

        if y_prob is not None:
            metrics["roc_auc"] = roc_auc_score(y_true, y_prob, multi_class="ovr")

        return metrics

    def evaluate_regressor(
        self,
        y_true: np.ndarray,
        y_pred: np.ndarray,
    ) -> Dict[str, float]:
        """Evaluate regressor performance.

        Args:
            y_true: True values
            y_pred: Predicted values

        Returns:
            Dictionary of metrics
        """
        return {
            "mae": mean_absolute_error(y_true, y_pred),
            "mse": mean_squared_error(y_true, y_pred),
            "rmse": mean_squared_error(y_true, y_pred, squared=False),
            "r2": r2_score(y_true, y_pred),
        }


class UncertaintyEstimator:
    """Uncertainty estimation utilities."""

    def __init__(self):
        """Initialize estimator."""
        self.logger = logging.getLogger(__name__)

    def monte_carlo_dropout(
        self,
        model: nn.Module,
        X: torch.Tensor,
        n_samples: int = 100,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Estimate uncertainty using MC dropout.

        Args:
            model: PyTorch model with dropout
            X: Input tensor
            n_samples: Number of forward passes

        Returns:
            Mean and std of predictions
        """
        model.train()  # Enable dropout

        preds = []
        for _ in range(n_samples):
            with torch.no_grad():
                pred = model(X)
                preds.append(pred.cpu().numpy())

        preds = np.stack(preds)
        return np.mean(preds, axis=0), np.std(preds, axis=0)

    def ensemble_uncertainty(
        self,
        predictions: List[np.ndarray],
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Estimate uncertainty from ensemble predictions.

        Args:
            predictions: List of predictions from ensemble models

        Returns:
            Mean and std of predictions
        """
        predictions = np.stack(predictions)
        return np.mean(predictions, axis=0), np.std(predictions, axis=0)


class MolecularFeaturizer:
    """Molecular feature extraction."""

    def __init__(self, model_dir: Optional[Union[str, Path]] = None):
        """Initialize featurizer.

        Args:
            model_dir: Optional directory for model checkpoints
        """
        self.logger = logging.getLogger(__name__)
        self.model_dir = Path(model_dir) if model_dir else None
        self.descriptor_calculator = DescriptorCalculator()
        self.pharmacophore_generator = PharmacophoreGenerator()
        self.similarity_calculator = SimilarityProcessor()

        # Load structural alerts and SMARTS patterns
        self.structural_alerts = self._load_structural_alerts()
        self.functional_groups = self._load_functional_groups()
        self.toxicophores = self._load_toxicophores()
        self.psychophores = self._load_psychophores()

    def generate_features(
        self,
        mol: Chem.Mol,
        include_3d: bool = False,
        feature_types: Optional[List[str]] = None,
        batch_size: Optional[int] = None,
    ) -> Dict[str, Union[np.ndarray, Data, CompoundData]]:
        """Generate comprehensive molecular features.

        Args:
            mol: RDKit molecule
            include_3d: Whether to include 3D features
            feature_types: Types of features to generate
            batch_size: Optional batch size for graph features

        Returns:
            Dictionary of feature arrays and objects
        """
        if mol is None:
            return {}

        features = {}
        try:
            # Generate descriptors and fragments
            features.update(self._generate_descriptors_and_fragments(mol, feature_types))

            # Generate decomposition and fingerprints
            features.update(self._generate_decomposition_and_fingerprints(mol, feature_types))

            # Generate pharmacophore and graph features
            features.update(self._generate_pharmacophore_and_graph(mol, feature_types, batch_size))

            # Generate 3D and structural features
            features.update(self._generate_3d_and_structure(mol, include_3d, feature_types))

            # Create CompoundData object
            features["compound"] = self._create_compound_data(mol)

            return features

        except Exception as e:
            self.logger.error(f"Error generating features: {str(e)}")
            return features

    def _generate_descriptors_and_fragments(self, mol: Chem.Mol, feature_types: Optional[List[str]] = None) -> Dict:
        """Generate descriptors and fragment features."""
        features = {}

        # Generate descriptors
        if feature_types is None or "descriptors" in feature_types:
            features.update(self._generate_descriptors(mol))

        # Generate fragment features
        if feature_types is None or "fragments" in feature_types:
            features.update(self._generate_fragment_features(mol))

        return features

    def _generate_decomposition_and_fingerprints(self, mol: Chem.Mol, feature_types: Optional[List[str]] = None) -> Dict:
        """Generate decomposition and fingerprint features."""
        features = {}

        # Generate decomposition features
        if feature_types is None or "decomposition" in feature_types:
            features["decomposition"] = self._analyze_decomposition(mol)

        # Generate fingerprints
        if feature_types is None or "fingerprints" in feature_types:
            features.update({"fingerprints": compute_fingerprints(mol)})

        return features

    def _generate_pharmacophore_and_graph(
        self,
        mol: Chem.Mol,
        feature_types: Optional[List[str]] = None,
        batch_size: Optional[int] = None,
    ) -> Dict:
        """Generate pharmacophore and graph features."""
        features = {}

        # Generate pharmacophore features
        if feature_types is None or "pharmacophore" in feature_types:
            features["pharmacophore"] = self.pharmacophore_generator.generate_features(mol)

        # Generate graph features
        if feature_types is None or "graph" in feature_types:
            features.update(self._generate_graph_data(mol, batch_size))

        return features

    def _generate_3d_and_structure(
        self,
        mol: Chem.Mol,
        include_3d: bool = False,
        feature_types: Optional[List[str]] = None,
    ) -> Dict:
        """Generate 3D and structural features."""
        features = {}

        # Generate 3D features
        if include_3d and (feature_types is None or "3d" in feature_types):
            features["3d"] = self._generate_3d_features(mol)

        # Generate structural analysis
        features.update(self._analyze_structure(mol))

        return features

    def _generate_descriptors(self, mol: Chem.Mol) -> Dict:
        """Generate molecular descriptors."""
        descriptors = self.descriptor_calculator.calculate(mol)
        descriptors.update(
            {
                "num_rings": rdMolDescriptors.CalcNumRings(mol),
                "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
                "num_aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
                "num_saturated_rings": rdMolDescriptors.CalcNumSaturatedRings(mol),
                "num_heterocycles": rdMolDescriptors.CalcNumHeterocycles(mol),
                "num_spiro_atoms": rdMolDescriptors.CalcNumSpiroAtoms(mol),
                "num_bridgeheads": rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
            }
        )
        return {"descriptors": descriptors}

    def _generate_fragment_features(self, mol: Chem.Mol) -> Dict:
        """Generate fragment-based features."""
        return {
            "fragments": {
                "num_acid_groups": Fragments.fr_Al_COO(mol),
                "num_alcohols": Fragments.fr_alcohol(mol),
                "num_amides": Fragments.fr_amide(mol),
                "num_amines": Fragments.fr_amine(mol),
                "num_aromatic_rings": Fragments.fr_Ar_N(mol),
                "num_esters": Fragments.fr_ester(mol),
                "num_ethers": Fragments.fr_ether(mol),
                "num_ketones": Fragments.fr_ketone(mol),
                "num_phenols": Fragments.fr_phenol(mol),
            }
        }

    def _generate_graph_data(self, mol: Chem.Mol, batch_size: Optional[int] = None) -> Dict:
        """Generate graph-based features."""
        graph_data = mol_to_graph(mol)
        features = {"graph": graph_data}

        if batch_size:
            batch = Batch.from_data_list([graph_data])
            features["graph"] = batch
            features["graph_pooled"] = {
                "sum": global_add_pool(graph_data.x, graph_data.batch),
                "max": global_max_pool(graph_data.x, graph_data.batch),
                "mean": global_mean_pool(graph_data.x, graph_data.batch),
            }

        return features

    def _analyze_structure(self, mol: Chem.Mol) -> Dict:
        """Analyze molecular structure."""
        return {
            "alerts": self._check_structural_alerts(mol),
            "functional_groups": self._identify_functional_groups(mol),
            "toxicophores": self._identify_toxicophores(mol),
            "psychophores": self._identify_psychophores(mol),
        }

    def _create_compound_data(self, mol: Chem.Mol) -> CompoundData:
        """Create CompoundData object from molecule."""
        compound = CompoundData(
            name=mol.GetProp("_Name") if mol.HasProp("_Name") else "",
            smiles=Chem.MolToSmiles(mol),
        )
        compound.molecular_weight = Descriptors.ExactMolWt(mol)
        compound.logp = Descriptors.MolLogP(mol)
        compound.tpsa = Descriptors.TPSA(mol)
        compound.hbd = rdMolDescriptors.CalcNumHBD(mol)
        compound.hba = rdMolDescriptors.CalcNumHBA(mol)
        compound.rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        return compound

    def _analyze_decomposition(self, mol: Chem.Mol) -> Dict:
        """Analyze molecule decomposition."""
        decomp = {}

        try:
            # Get Murcko scaffold
            scaffold = AllChem.MurckoScaffoldSmiles(mol=mol, includeChirality=True)
            decomp["scaffold_smiles"] = scaffold

            # Get largest ring system
            ring_info = mol.GetRingInfo()
            if ring_info.NumRings() > 0:
                # Get atoms in rings
                ring_atoms = set()
                for ring in ring_info.AtomRings():
                    ring_atoms.update(ring)

                # Create substructure from ring atoms
                ring_mol = Chem.PathToSubmol(mol, list(ring_atoms))
                if ring_mol:
                    decomp["largest_ring_system"] = Chem.MolToSmiles(ring_mol)

            # Identify R-groups using RGroupDecomposition
            params = rdRGroupDecomposition.RGroupDecompositionParameters()
            params.removeHydrogensPostMatch = True

            # Use scaffold as core pattern
            core = Chem.MolFromSmiles(scaffold)
            if core:
                decomp_obj = rdRGroupDecomposition.RGroupDecomposition(core, params)
                decomp_obj.Add(mol)
                if decomp_obj.Process():
                    rgroups = decomp_obj.GetRGroupsAsColumns(asSmiles=True)[0]
                    decomp["rgroups"] = {k: v for k, v in rgroups.items() if k != "Core"}

        except Exception as e:
            self.logger.error(f"Error in decomposition analysis: {str(e)}")

        return decomp

    def _generate_3d_features(self, mol: Chem.Mol) -> Dict[str, np.ndarray]:
        """Generate 3D molecular features."""
        features = {}

        try:
            # Generate 3D conformation if not present
            if mol.GetNumConformers() == 0:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
                AllChem.MMFFOptimizeMolecule(mol)

            # Calculate 3D descriptors
            features["radius_of_gyration"] = rdMolDescriptors.CalcRadiusOfGyration(mol)
            features["inertial_shape_factor"] = rdMolDescriptors.CalcInertialShapeFactor(mol)
            features["npr1"] = rdMolDescriptors.CalcNPR1(mol)
            features["npr2"] = rdMolDescriptors.CalcNPR2(mol)
            features["pmi1"] = rdMolDescriptors.CalcPMI1(mol)
            features["pmi2"] = rdMolDescriptors.CalcPMI2(mol)
            features["pmi3"] = rdMolDescriptors.CalcPMI3(mol)
            features["spherocity"] = rdMolDescriptors.CalcSpherocityIndex(mol)
            features["asphericity"] = rdMolDescriptors.CalcAsphericity(mol)
            features["eccentricity"] = rdMolDescriptors.CalcEccentricity(mol)

        except Exception as e:
            self.logger.error(f"Error generating 3D features: {str(e)}")

        return features

    def _check_structural_alerts(self, mol: Chem.Mol) -> List[str]:
        """Check molecule for structural alerts."""
        alerts = []
        for name, smarts in self.structural_alerts.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                alerts.append(name)
        return alerts

    def _identify_functional_groups(self, mol: Chem.Mol) -> List[str]:
        """Identify functional groups."""
        groups = []
        for name, smarts in self.functional_groups.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                groups.append(name)
        return groups

    def _identify_toxicophores(self, mol: Chem.Mol) -> List[str]:
        """Identify toxicophores."""
        toxicophores = []
        for name, smarts in self.toxicophores.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                toxicophores.append(name)
        return toxicophores

    def _identify_psychophores(self, mol: Chem.Mol) -> List[str]:
        """Identify psychophores."""
        psychophores = []
        for name, smarts in self.psychophores.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                psychophores.append(name)
        return psychophores

    def _load_structural_alerts(self) -> Dict[str, str]:
        """Load structural alert SMARTS patterns."""
        return {
            "michael_acceptor": "[C;H1,H2]=C[C;H1]=O",
            "epoxide": "[O;R1][C;R1][C;R1]",
            "acyl_halide": "[C;H1](=O)[F,Cl,Br,I]",
            "hydrazine": "[N;H1][N;H1,H2]",
            "beta_lactam": "[C;R1]1[C;R1](=O)[N;R1][C;R1]1",
            "phosphoramide": "[P](=O)([O,N])[N]",
            "aromatic_nitro": "[N+](=O)[O-]",
            "sulfonate_ester": "OS(=O)(=O)[C,c]",
            "alkyl_halide": "[C;!R][F,Cl,Br,I]",
            "peroxide": "[O;H1][O;H1]",
            "imine": "[C;H1]=[N;H1]",
            "thiol": "[S;H1]",
            "aldehyde": "[C;H1](=O)[C,H]",
            "quinone": "[O;H1]c1c([O;H1])cccc1",
            "hydroquinone": "[O;H1]c1ccc([O;H1])cc1",
            "anhydride": "[C;H1](=O)O[C;H1](=O)",
            "azide": "[N-][N+]#N",
            "aziridine": "[N;R1][C;R1][C;R1]",
            "thioester": "[S;H1][C;H1](=O)",
            "sulfonyl_halide": "S(=O)(=O)[F,Cl,Br,I]",
        }

    def _load_functional_groups(self) -> Dict[str, str]:
        """Load functional group SMARTS patterns."""
        return {
            "alcohol": "[OH]",
            "phenol": "[OH]c",
            "carboxylic_acid": "[OH]C=O",
            "amine": "[NH2]",
            "amide": "NC=O",
            "ester": "COC=O",
            "ether": "[OR]",
            "ketone": "[#6]C(=O)[#6]",
            "aldehyde": "[CH]=[O]",
            "alkene": "C=C",
            "alkyne": "C#C",
            "aromatic": "c1ccccc1",
            "heterocycle": "[!#6;R]",
            "nitro": "[N+](=O)[O-]",
            "nitrile": "C#N",
            "sulfide": "[SX2]",
            "sulfoxide": "[SX3](=O)",
            "sulfone": "[SX4](=O)(=O)",
            "phosphate": "[P](=O)([O])[O]",
        }

    def _load_toxicophores(self) -> Dict[str, str]:
        """Load toxicophore SMARTS patterns."""
        return {
            "alkylating_agent": "[C;!R][F,Cl,Br,I]",
            "michael_acceptor": "[C;H1,H2]=C[C;H1]=O",
            "epoxide": "[O;R1][C;R1][C;R1]",
            "aziridine": "[N;R1][C;R1][C;R1]",
            "quinone": "[O;H1]c1c([O;H1])cccc1",
            "aromatic_nitro": "[N+](=O)[O-]",
            "aromatic_amine": "c[NH2]",
            "hydrazine": "[N;H1][N;H1,H2]",
            "aliphatic_halide": "[C;!R][F,Cl,Br,I]",
            "peroxide": "[O;H1][O;H1]",
            "acyl_halide": "[C;H1](=O)[F,Cl,Br,I]",
            "thiocarbonyl": "[C;X3]=[S;X1]",
            "isocyanate": "[N;X2]=[C;X2]=[O;X1]",
            "isothiocyanate": "[N;X2]=[C;X2]=[S;X1]",
            "alpha_haloketone": "[C;X4][C;X3](=O)[F,Cl,Br,I]",
        }

    def _load_psychophores(self) -> Dict[str, str]:
        """Load psychophore SMARTS patterns."""
        return {
            "phenethylamine": "NCCc1ccccc1",
            "tryptamine": "NCCc1c[nH]c2ccccc12",
            "amphetamine": "CC(N)Cc1ccccc1",
            "cathinone": "CC(N)C(=O)c1ccccc1",
            "phenylpiperazine": "c1ccccc1N1CCNCC1",
            "phenylpiperidine": "c1ccccc1N1CCCCC1",
            "benzofuran": "c1ccc2c(c1)cco2",
            "indole": "c1ccc2c(c1)[nH]cc2",
            "benzodioxole": "c1ccc2c(c1)OCO2",
            "quinoline": "c1ccc2c(c1)cccn2",
            "isoquinoline": "c1ccc2c(c1)ccnc2",
            "tropane": "CN1[C@H]2CC[C@@H]1CC2",
            "morphinan": "c1ccc2c(c1)CC1C3CCCC(C2)N3C1",
            "ergoline": "CN1C[C@@H](C=C2[C@H]1Cc1c[nH]c3cccc2c13)C",
            "lysergamide": "CN1C[C@@H](C=C2[C@H]1Cc1c[nH]c3cccc2c13)C(=O)N",
        }
