"""Base structure processing functionality.

This module provides core functionality for:
1. Structure validation and standardization 
2. SMILES/InChI conversion and handling
3. Structure cleaning and normalization
4. Fingerprint generation and similarity
5. Substructure searching and highlighting
6. Descriptor calculation
7. ML-enhanced structure analysis
8. Pharmacophore detection
9. Toxicophore analysis
10. Activity prediction
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Descriptors3D,
    MACCSkeys,
    Draw,
    rdDeprotect,
    rdRGroupDecomposition,
    rdMolDescriptors,
    rdMolTransforms,
    rdFMCS,
    rdMolAlign,
)
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Descriptors import MoleculeDescriptors

# Optional ML imports
try:
    import torch
    from torch_geometric.data import Data
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.cluster import DBSCAN
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.svm import SVC, SVR

    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False

logger = logging.getLogger(__name__)


class BaseStructureProcessor:
    """Base class for molecular structure processing."""

    # Common fingerprint parameters
    FP_PARAMS = {
        "morgan": {
            "radius": 2,
            "n_bits": 2048,
            "use_features": True,
        },
        "maccs": {},
        "topological": {
            "n_bits": 2048,
            "min_path": 1,
            "max_path": 7,
        },
        "rdkit": {
            "min_path": 1,
            "max_path": 7,
            "fp_size": 2048,
        },
    }

    def __init__(self, config: Optional[Dict] = None):
        """Initialize structure processor."""
        self.logger = logging.getLogger(__name__)
        self.config = config or {}

        # Initialize ML components if available
        if ML_AVAILABLE:
            self.scaler = StandardScaler()
            self.pca = PCA(n_components=0.95)
            self.clustering = DBSCAN(eps=0.5, min_samples=5)

            # Initialize ML models
            self.activity_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
            self.affinity_regressor = RandomForestRegressor(n_estimators=100, random_state=42)
            self.toxicity_classifier = RandomForestClassifier(n_estimators=100, random_state=42)

    def validate_structure(self, smiles: str) -> Optional[Chem.Mol]:
        """Validate and standardize chemical structure."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate 3D conformation
            try:
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)
            except (ValueError, RuntimeError) as e:
                self.logger.warning(f"Could not generate 3D conformation: {str(e)}")

            return mol

        except Exception as e:
            self.logger.error(f"Error validating structure: {str(e)}")
            return None

    def clean_structure(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Clean and normalize structure."""
        try:
            if mol is None:
                return None

            # Remove hydrogens
            mol = Chem.RemoveHs(mol)

            # Sanitize molecule
            Chem.SanitizeMol(mol)

            # Canonicalize atom order
            mol = Chem.RenumberAtoms(mol, Chem.CanonicalRankAtoms(mol))

            # Add hydrogens back
            mol = Chem.AddHs(mol)

            return mol

        except Exception as e:
            self.logger.error(f"Structure cleaning error: {str(e)}")
            return None

    def standardize_smiles(self, smiles: str) -> Optional[str]:
        """Generate standardized canonical SMILES."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception as e:
            self.logger.error(f"SMILES standardization error: {str(e)}")
            return None

    def standardize_inchi(self, inchi: str) -> Optional[str]:
        """Generate standardized InChI."""
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                return None
            return Chem.MolToInchi(mol)
        except Exception as e:
            self.logger.error(f"InChI standardization error: {str(e)}")
            return None

    def generate_conformer(
        self,
        mol: Chem.Mol,
        num_confs: int = 1,
        optimize: bool = True,
    ) -> Optional[Chem.Mol]:
        """Generate 3D conformer(s)."""
        try:
            if mol is None:
                return None

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate conformers
            conf_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_confs,
                randomSeed=42,
                useRandomCoords=True,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True,
                enforceChirality=True,
            )

            if not conf_ids:
                return None

            # Optimize conformers
            if optimize:
                for conf_id in conf_ids:
                    if AllChem.MMFFHasAllMoleculeParams(mol):
                        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
                    else:
                        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)

            return mol

        except Exception as e:
            self.logger.error(f"Conformer generation error: {str(e)}")
            return None

    def generate_fingerprint(
        self,
        mol: Chem.Mol,
        fp_type: str = "morgan",
        params: Optional[Dict] = None,
    ) -> Optional[Union[DataStructs.ExplicitBitVect, np.ndarray]]:
        """Generate molecular fingerprint."""
        try:
            if mol is None:
                return None

            # Get default parameters
            fp_params = self.FP_PARAMS.get(fp_type, {}).copy()
            if params:
                fp_params.update(params)

            # Generate fingerprint
            if fp_type == "morgan":
                return AllChem.GetMorganFingerprintAsBitVect(
                    mol,
                    radius=fp_params["radius"],
                    nBits=fp_params["n_bits"],
                    useFeatures=fp_params["use_features"],
                )
            elif fp_type == "maccs":
                return MACCSkeys.GenMACCSKeys(mol)
            elif fp_type == "topological":
                return rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                    mol,
                    minLength=fp_params["min_path"],
                    maxLength=fp_params["max_path"],
                    nBits=fp_params["n_bits"],
                )
            elif fp_type == "rdkit":
                return Chem.RDKFingerprint(
                    mol,
                    minPath=fp_params["min_path"],
                    maxPath=fp_params["max_path"],
                    fpSize=fp_params["fp_size"],
                )
            else:
                self.logger.warning(f"Unknown fingerprint type: {fp_type}")
                return None

        except Exception as e:
            self.logger.error(f"Error generating fingerprint: {str(e)}")
            return None

    def calculate_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        method: str = "tanimoto",
        fp_type: str = "morgan",
    ) -> Optional[float]:
        """Calculate similarity between molecules."""
        try:
            if mol1 is None or mol2 is None:
                return None

            # Generate fingerprints
            fp1 = self.generate_fingerprint(mol1, fp_type)
            fp2 = self.generate_fingerprint(mol2, fp_type)
            if fp1 is None or fp2 is None:
                return None

            # Calculate similarity
            if method == "tanimoto":
                return DataStructs.TanimotoSimilarity(fp1, fp2)
            elif method == "dice":
                return DataStructs.DiceSimilarity(fp1, fp2)
            elif method == "cosine":
                return DataStructs.CosineSimilarity(fp1, fp2)
            else:
                self.logger.warning(f"Unknown similarity method: {method}")
                return None

        except Exception as e:
            self.logger.error(f"Error calculating similarity: {str(e)}")
            return None

    def calculate_descriptors(
        self,
        mol: Chem.Mol,
        include_3d: bool = False,
    ) -> Dict[str, float]:
        """Calculate molecular descriptors."""
        try:
            if mol is None:
                return {}

            descriptors = {}

            # Basic descriptors
            descriptors["mw"] = Descriptors.ExactMolWt(mol)
            descriptors["logp"] = Descriptors.MolLogP(mol)
            descriptors["hbd"] = Descriptors.NumHDonors(mol)
            descriptors["hba"] = Descriptors.NumHAcceptors(mol)
            descriptors["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
            descriptors["aromatic_rings"] = Descriptors.NumAromaticRings(mol)
            descriptors["heavy_atoms"] = mol.GetNumHeavyAtoms()
            descriptors["tpsa"] = Descriptors.TPSA(mol)

            # Constitutional descriptors
            descriptors["atoms"] = mol.GetNumAtoms()
            descriptors["bonds"] = mol.GetNumBonds()
            descriptors["rings"] = Descriptors.RingCount(mol)
            descriptors["aliphatic_rings"] = Descriptors.NumAliphaticRings(mol)

            # Topological descriptors
            descriptors["bertz"] = Descriptors.BertzCT(mol)
            descriptors["wiener_index"] = Descriptors.WienerIndex(mol)
            descriptors["balaban_j"] = Descriptors.BalabanJ(mol)

            # Electronic descriptors
            descriptors["max_abs_charge"] = Descriptors.MaxAbsPartialCharge(mol)
            descriptors["max_pos_charge"] = Descriptors.MaxPartialCharge(mol)
            descriptors["max_neg_charge"] = Descriptors.MinPartialCharge(mol)

            # 3D descriptors if requested and conformer exists
            if include_3d and mol.GetNumConformers() > 0:
                descriptors["asphericity"] = Descriptors3D.Asphericity(mol)
                descriptors["eccentricity"] = Descriptors3D.Eccentricity(mol)
                descriptors["inertial_shape_factor"] = Descriptors3D.InertialShapeFactor(mol)
                descriptors["npr1"] = Descriptors3D.NPR1(mol)
                descriptors["npr2"] = Descriptors3D.NPR2(mol)
                descriptors["pmi1"] = Descriptors3D.PMI1(mol)
                descriptors["pmi2"] = Descriptors3D.PMI2(mol)
                descriptors["pmi3"] = Descriptors3D.PMI3(mol)
                descriptors["radius_of_gyration"] = Descriptors3D.RadiusOfGyration(mol)
                descriptors["spherocity_index"] = Descriptors3D.SpherocityIndex(mol)

            return descriptors

        except Exception as e:
            self.logger.error(f"Error calculating descriptors: {str(e)}")
            return {}

    def get_graph_features(self, mol: Chem.Mol) -> Optional[Data]:
        """Convert molecule to graph features for GNN."""
        if not ML_AVAILABLE:
            self.logger.warning("ML functionality not available")
            return None

        try:
            # Node features
            atom_features = []
            for atom in mol.GetAtoms():
                features = [
                    atom.GetAtomicNum(),
                    atom.GetTotalDegree(),
                    atom.GetFormalCharge(),
                    atom.GetTotalNumHs(),
                    atom.GetIsAromatic(),
                    atom.GetMass(),
                    atom.GetExplicitValence(),
                    atom.GetImplicitValence(),
                    int(atom.IsInRing()),
                    int(atom.IsInRingSize(5)),
                    int(atom.IsInRingSize(6)),
                ]
                atom_features.append(features)

            # Edge features
            edge_indices = []
            edge_features = []
            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                edge_indices += [[i, j], [j, i]]
                features = [
                    bond.GetBondTypeAsDouble(),
                    bond.GetIsAromatic(),
                    bond.IsInRing(),
                    bond.IsInRingSize(5),
                    bond.IsInRingSize(6),
                    float(bond.GetIsConjugated()),
                    float(bond.GetIsRotatable()),
                ]
                edge_features += [features, features]

            return Data(
                x=torch.tensor(atom_features, dtype=torch.float),
                edge_index=torch.tensor(edge_indices, dtype=torch.long).t().contiguous(),
                edge_attr=torch.tensor(edge_features, dtype=torch.float),
            )

        except Exception as e:
            self.logger.error(f"Error generating graph features: {str(e)}")
            return None

    def predict_activity(
        self,
        mol: Chem.Mol,
        activity_type: str,
    ) -> Optional[Dict[str, float]]:
        """Predict biological activity."""
        if not ML_AVAILABLE:
            self.logger.warning("ML functionality not available")
            return None

        try:
            # Calculate descriptors
            descriptors = self.calculate_descriptors(mol)
            if not descriptors:
                return None

            # Convert to feature vector
            X = np.array(list(descriptors.values())).reshape(1, -1)

            # Scale features
            X_scaled = self.scaler.transform(X)

            # Make predictions
            if activity_type == "classification":
                prob = self.activity_classifier.predict_proba(X_scaled)[0]
                return {
                    "active_probability": float(prob[1]),
                    "inactive_probability": float(prob[0]),
                }
            elif activity_type == "regression":
                pred = self.affinity_regressor.predict(X_scaled)[0]
                return {"predicted_affinity": float(pred)}
            else:
                self.logger.warning(f"Unknown activity type: {activity_type}")
                return None

        except Exception as e:
            self.logger.error(f"Error predicting activity: {str(e)}")
            return None

    def predict_toxicity(self, mol: Chem.Mol) -> Optional[Dict[str, float]]:
        """Predict toxicity risks."""
        if not ML_AVAILABLE:
            self.logger.warning("ML functionality not available")
            return None

        try:
            # Calculate descriptors
            descriptors = self.calculate_descriptors(mol)
            if not descriptors:
                return None

            # Convert to feature vector
            X = np.array(list(descriptors.values())).reshape(1, -1)

            # Scale features
            X_scaled = self.scaler.transform(X)

            # Make predictions
            prob = self.toxicity_classifier.predict_proba(X_scaled)[0]
            return {
                "toxic_probability": float(prob[1]),
                "nontoxic_probability": float(prob[0]),
            }

        except Exception as e:
            self.logger.error(f"Error predicting toxicity: {str(e)}")
            return None

    @staticmethod
    def _get_functional_group_smarts() -> Dict[str, str]:
        """Get SMARTS patterns for functional groups."""
        return {
            "[OH]": "hydroxyl",
            "[NH2]": "primary_amine",
            "[NH1][!$(*C=[O,N,S])]": "secondary_amine",
            "[NH0][!$(*C=[O,N,S])]": "tertiary_amine",
            "[CX3](=[OX1])[OX2H1]": "carboxylic_acid",
            "[CX3](=[OX1])[OX2H0]": "carboxylate",
            "[CX3](=[OX1])[NX3H2]": "primary_amide",
            "[CX3](=[OX1])[NX3H1]": "secondary_amide",
            "[CX3](=[OX1])[NX3H0]": "tertiary_amide",
            "[$([CX3]([#6])[#6])]=O": "ketone",
            "[$([CX3H][#6])]=O": "aldehyde",
            "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]": "phosphate",
            "[SX4](=[OX1])(=[OX1])([OX2H])[OX2H]": "sulfonic_acid",
            "[#6][SX2][#6]": "thioether",
            "[#6][SX2H]": "thiol",
            "[#6][PX3][#6]": "phosphine",
            "[CX3](=[OX1])[F,Cl,Br,I]": "acid_halide",
            "[NX3][CX3](=[OX1])[#6]": "amide",
            "[CX3](=[OX1])[OX2][CX4]": "ester",
            "[CX3](=[OX1])[OX2][CX3](=[OX1])": "anhydride",
            "[NX3H2,NX4H3][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]": "amino_acid",
        }
