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
from typing import Dict, List, Optional, Set, Tuple, Union

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
    rdDecomposition,
    rdMolDescriptors,
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
from torch_geometric.data import Batch, Data
from torch_geometric.nn import global_add_pool, global_max_pool, global_mean_pool

from ....models.core import CompoundData
from ..descriptors import DescriptorCalculator
from ..pharmacophore import PharmacophoreGenerator
from ..similarity import MolecularSimilarity
from .activity import ActivityPredictor


class MolecularFeaturizer:
    """Molecular feature extraction."""

    def __init__(self, model_dir: Optional[Union[str, Path]] = None):
        """
        Initialize featurizer.

        Args:
            model_dir: Optional directory for model checkpoints
        """
        self.logger = logging.getLogger(__name__)
        self.model_dir = Path(model_dir) if model_dir else None
        self.descriptor_calculator = DescriptorCalculator()
        self.pharmacophore_generator = PharmacophoreGenerator()
        self.similarity_calculator = MolecularSimilarity()
        self.activity_predictor = ActivityPredictor()

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
        """
        Generate comprehensive molecular features.

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
            features.update(
                self._generate_pharmacophore_and_graph(mol, feature_types, batch_size)
            )

            # Generate 3D and structural features
            features.update(
                self._generate_3d_and_structure(mol, include_3d, feature_types)
            )

            # Generate activity predictions
            features["predicted_activities"] = (
                self.activity_predictor.predict_activities(mol)
            )

            # Create CompoundData object
            features["compound"] = self._create_compound_data(mol)

            return features

        except Exception as e:
            self.logger.error(f"Error generating features: {str(e)}")
            return features

    def _generate_descriptors_and_fragments(
        self, mol: Chem.Mol, feature_types: Optional[List[str]] = None
    ) -> Dict:
        """Generate descriptors and fragment features."""
        features = {}

        # Generate descriptors
        if feature_types is None or "descriptors" in feature_types:
            features.update(self._generate_descriptors(mol))

        # Generate fragment features
        if feature_types is None or "fragments" in feature_types:
            features.update(self._generate_fragment_features(mol))

        return features

    def _generate_decomposition_and_fingerprints(
        self, mol: Chem.Mol, feature_types: Optional[List[str]] = None
    ) -> Dict:
        """Generate decomposition and fingerprint features."""
        features = {}

        # Generate decomposition features
        if feature_types is None or "decomposition" in feature_types:
            features["decomposition"] = self._analyze_decomposition(mol)

        # Generate fingerprints
        if feature_types is None or "fingerprints" in feature_types:
            features.update(self._generate_fingerprints(mol))

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
            features["pharmacophore"] = (
                self.pharmacophore_generator.generate_features(mol)
            )

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
        descriptors.update({
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "num_aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
            "num_saturated_rings": rdMolDescriptors.CalcNumSaturatedRings(mol),
            "num_heterocycles": rdMolDescriptors.CalcNumHeterocycles(mol),
            "num_spiro_atoms": rdMolDescriptors.CalcNumSpiroAtoms(mol),
            "num_bridgeheads": rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
        })
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

    def _generate_fingerprints(self, mol: Chem.Mol) -> Dict:
        """Generate molecular fingerprints."""
        return {
            "morgan": self._generate_morgan_fingerprint(mol),
            "maccs": self._generate_maccs_fingerprint(mol),
            "rdkit": self._generate_rdkit_fingerprint(mol),
            "pattern": self._generate_pattern_fingerprint(mol),
            "layered": self._generate_layered_fingerprint(mol),
        }

    def _generate_graph_data(
        self, mol: Chem.Mol, batch_size: Optional[int] = None
    ) -> Dict:
        """Generate graph-based features."""
        graph_data = self._generate_graph_features(mol)
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
        
        # Get Murcko scaffold
        scaffold = rdDecomposition.GetScaffoldForMol(mol)
        decomp["scaffold_smiles"] = Chem.MolToSmiles(scaffold) if scaffold else ""
        
        # Get framework
        framework = rdDecomposition.GetFrameworkForMol(mol)
        decomp["framework_smiles"] = Chem.MolToSmiles(framework) if framework else ""
        
        # Get ring systems
        ring_systems = rdDecomposition.GetRingSystems(mol)
        decomp["ring_systems"] = [
            Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            for smiles in ring_systems
        ]
    def _generate_morgan_fingerprint(
        self, mol: Chem.Mol, radius: int = 3, nbits: int = 2048
    ) -> np.ndarray:
        """Generate Morgan fingerprint with chirality."""
        return np.array(
            list(
                AllChem.GetMorganFingerprintAsBitVect(
                    mol,
                    radius,
                    nBits=nbits,
                    useChirality=True,
                    useFeatures=True,
                )
            )
        )

    def _generate_maccs_fingerprint(self, mol: Chem.Mol) -> np.ndarray:
        """Generate MACCS keys fingerprint."""
        return np.array(list(MACCSkeys.GenMACCSKeys(mol)))

    def _generate_rdkit_fingerprint(
        self, mol: Chem.Mol, nbits: int = 2048
    ) -> np.ndarray:
        """Generate RDKit topological fingerprint."""
        return np.array(
            list(
                Chem.RDKFingerprint(
                    mol,
                    fpSize=nbits,
                    minPath=1,
                    maxPath=7,
                    useHs=True,
                )
            )
        )

    def _generate_pattern_fingerprint(
        self, mol: Chem.Mol, nbits: int = 2048
    ) -> np.ndarray:
        """Generate pattern fingerprint."""
        return np.array(
            list(
                Chem.PatternFingerprint(
                    mol,
                    fpSize=nbits,
                    tautomerFingerprints=True,
                )
            )
        )

    def _generate_layered_fingerprint(
        self, mol: Chem.Mol, nbits: int = 2048
    ) -> np.ndarray:
        """Generate layered fingerprint."""
        return np.array(
            list(
                Chem.LayeredFingerprint(
                    mol,
                    fpSize=nbits,
                    layerFlags=0xFFFFFFFF,
                )
            )
        )

    def _generate_graph_features(self, mol: Chem.Mol) -> Data:
        """Generate molecular graph features."""
        # Node features
        atomic_nums = []
        aromatic = []
        sp = []
        sp2 = []
        sp3 = []
        num_hs = []
        formal_charge = []
        radical_electrons = []
        in_ring = []
        chirality = []

        for atom in mol.GetAtoms():
            atomic_nums.append(atom.GetAtomicNum())
            aromatic.append(1 if atom.GetIsAromatic() else 0)
            hybridization = atom.GetHybridization()
            sp.append(1 if hybridization == Chem.HybridizationType.SP else 0)
            sp2.append(1 if hybridization == Chem.HybridizationType.SP2 else 0)
            sp3.append(1 if hybridization == Chem.HybridizationType.SP3 else 0)
            num_hs.append(atom.GetTotalNumHs())
            formal_charge.append(atom.GetFormalCharge())
            radical_electrons.append(atom.GetNumRadicalElectrons())
            in_ring.append(1 if atom.IsInRing() else 0)
            if atom.HasProp("_CIPCode"):
                chirality.append(1 if atom.GetProp("_CIPCode") == "R" else -1)
            else:
                chirality.append(0)

        x = torch.tensor(
            [
                atomic_nums,
                aromatic,
                sp,
                sp2,
                sp3,
                num_hs,
                formal_charge,
                radical_electrons,
                in_ring,
                chirality,
            ],
            dtype=torch.float,
        ).t()

        # Edge features
        edge_indices = []
        edge_attrs = []

        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            edge_indices += [[i, j], [j, i]]

            # Bond features
            bond_type = bond.GetBondType()
            bond_conjugated = bond.GetIsConjugated()
            bond_in_ring = bond.IsInRing()
            bond_stereo = bond.GetStereo()

            edge_attr = [
                int(bond_type == Chem.BondType.SINGLE),
                int(bond_type == Chem.BondType.DOUBLE),
                int(bond_type == Chem.BondType.TRIPLE),
                int(bond_type == Chem.BondType.AROMATIC),
                int(bond_conjugated),
                int(bond_in_ring),
                int(bond_stereo > Chem.BondStereo.STEREONONE),
            ]

            edge_attrs += [edge_attr, edge_attr]

        edge_index = torch.tensor(edge_indices, dtype=torch.long).t()
        edge_attr = torch.tensor(edge_attrs, dtype=torch.float)

        return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)

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
            features["inertial_shape_factor"] = (
                rdMolDescriptors.CalcInertialShapeFactor(mol)
            )
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


class DataPreprocessor:
    """Data preprocessing utilities."""

    def __init__(self):
        """Initialize preprocessor."""
        self.logger = logging.getLogger(__name__)

    def normalize_features(
        self,
        features: np.ndarray,
        method: str = "standard",
        params: Optional[Dict] = None,
    ) -> Tuple[np.ndarray, Dict]:
        """
        Normalize feature values.

        Args:
            features: Feature array
            method: Normalization method (standard/minmax/robust)
            params: Optional normalization parameters

        Returns:
            Normalized features and parameters
        """
        try:
            if method == "standard":
                if params is None:
                    mean = np.mean(features, axis=0)
                    std = np.std(features, axis=0)
                    params = {"mean": mean, "std": std}
                features = (features - params["mean"]) / (params["std"] + 1e-8)

            elif method == "minmax":
                if params is None:
                    min_val = np.min(features, axis=0)
                    max_val = np.max(features, axis=0)
                    params = {"min": min_val, "max": max_val}
                features = (features - params["min"]) / (
                    params["max"] - params["min"] + 1e-8
                )

            elif method == "robust":
                if params is None:
                    q1 = np.percentile(features, 25, axis=0)
                    q3 = np.percentile(features, 75, axis=0)
                    iqr = q3 - q1
                    params = {"q1": q1, "q3": q3, "iqr": iqr}
                features = (features - params["q1"]) / (params["iqr"] + 1e-8)

            return features, params

        except Exception as e:
            self.logger.error(f"Error normalizing features: {str(e)}")
            return features, {}

    def handle_missing_values(
        self,
        features: np.ndarray,
        strategy: str = "mean",
        fill_value: Optional[float] = None,
    ) -> Tuple[np.ndarray, Dict]:
        """
        Handle missing values in features.

        Args:
            features: Feature array
            strategy: Handling strategy (mean/median/constant)
            fill_value: Value for constant strategy

        Returns:
            Processed features and parameters
        """
        try:
            mask = np.isnan(features)
            params = {}

            if strategy == "mean":
                params["fill_values"] = np.nanmean(features, axis=0)
            elif strategy == "median":
                params["fill_values"] = np.nanmedian(features, axis=0)
            elif strategy == "constant":
                params["fill_values"] = (
                    np.full(features.shape[1], fill_value)
                    if fill_value is not None
                    else np.zeros(features.shape[1])
                )

            features = features.copy()
            for i in range(features.shape[1]):
                features[mask[:, i], i] = params["fill_values"][i]

            return features, params

        except Exception as e:
            self.logger.error(f"Error handling missing values: {str(e)}")
            return features, {}


class ModelEvaluator:
    """Model evaluation utilities."""

    def __init__(self):
        """Initialize evaluator."""
        self.logger = logging.getLogger(__name__)

    def compute_metrics(
        self,
        y_true: np.ndarray,
        y_pred: np.ndarray,
        task_type: str = "regression",
    ) -> Dict[str, float]:
        """
        Compute evaluation metrics.

        Args:
            y_true: True values
            y_pred: Predicted values
            task_type: Task type (regression/classification)

        Returns:
            Dictionary of metrics
        """
        try:
            metrics = {}

            if task_type == "regression":
                metrics["mse"] = float(mean_squared_error(y_true, y_pred))
                metrics["rmse"] = float(np.sqrt(metrics["mse"]))
                metrics["mae"] = float(mean_absolute_error(y_true, y_pred))
                metrics["r2"] = float(r2_score(y_true, y_pred))
                metrics["pearson"] = float(
                    np.corrcoef(y_true.flatten(), y_pred.flatten())[0, 1]
                )
                metrics["spearman"] = float(
                    pd.Series(y_true.flatten()).corr(
                        pd.Series(y_pred.flatten()), method="spearman"
                    )
                )
            else:
                y_pred_class = (y_pred > 0.5).astype(int)
                metrics["accuracy"] = float(accuracy_score(y_true, y_pred_class))
                metrics["precision"] = float(precision_score(y_true, y_pred_class))
                metrics["recall"] = float(recall_score(y_true, y_pred_class))
                metrics["f1"] = float(f1_score(y_true, y_pred_class))
                try:
                    metrics["auc"] = float(roc_auc_score(y_true, y_pred))
                except (ValueError, TypeError) as e:
                    self.logger.warning(f"Could not compute AUC: {str(e)}")
                    metrics["auc"] = 0.0

            return metrics

        except Exception as e:
            self.logger.error(f"Error computing metrics: {str(e)}")
            return {}


class UncertaintyEstimator:
    """Model uncertainty estimation."""

    def __init__(self):
        """Initialize estimator."""
        self.logger = logging.getLogger(__name__)

    def estimate(
        self,
        model: nn.Module,
        inputs: torch.Tensor,
        num_samples: int = 10,
        method: str = "dropout",
    ) -> np.ndarray:
        """
        Estimate prediction uncertainty.

        Args:
            model: PyTorch model
            inputs: Input tensor
            num_samples: Number of Monte Carlo samples
            method: Uncertainty method (dropout/ensemble/bootstrap)

        Returns:
            Uncertainty estimates
        """
        try:
            if method == "dropout":
                return self._dropout_uncertainty(model, inputs, num_samples)
            elif method == "ensemble":
                return self._ensemble_uncertainty(model, inputs, num_samples)
            elif method == "bootstrap":
                return self._bootstrap_uncertainty(model, inputs, num_samples)
            else:
                raise ValueError(f"Unknown uncertainty method: {method}")

        except Exception as e:
            self.logger.error(f"Error estimating uncertainty: {str(e)}")
            return np.zeros(inputs.shape[0])

    def _dropout_uncertainty(
        self,
        model: nn.Module,
        inputs: torch.Tensor,
        num_samples: int,
    ) -> np.ndarray:
        """Estimate uncertainty using MC dropout."""
        model.train()  # Enable dropout
        predictions = []

        with torch.no_grad():
            for _ in range(num_samples):
                outputs = model(inputs)
                predictions.append(outputs.cpu().numpy())

        predictions = np.stack(predictions)
        return np.std(predictions, axis=0)

    def _ensemble_uncertainty(
        self,
        model: nn.Module,
        inputs: torch.Tensor,
        num_samples: int,
    ) -> np.ndarray:
        """Estimate uncertainty using deep ensembles."""
        predictions = []

        with torch.no_grad():
            for _ in range(num_samples):
                # Add random noise to model parameters
                for param in model.parameters():
                    noise = torch.randn_like(param) * 0.1
                    param.data += noise

                outputs = model(inputs)
                predictions.append(outputs.cpu().numpy())

                # Restore original parameters
                for param in model.parameters():
                    param.data -= noise

        predictions = np.stack(predictions)
        return np.std(predictions, axis=0)

    def _bootstrap_uncertainty(
        self,
        model: nn.Module,
        inputs: torch.Tensor,
        num_samples: int,
    ) -> np.ndarray:
        """Estimate uncertainty using bootstrapping."""
        predictions = []
        batch_size = inputs.shape[0]

        with torch.no_grad():
            for _ in range(num_samples):
                # Sample with replacement
                indices = np.random.choice(batch_size, size=batch_size, replace=True)
                batch = inputs[indices]
                outputs = model(batch)
                predictions.append(outputs.cpu().numpy())

        predictions = np.stack(predictions)
        return np.std(predictions, axis=0)
