"""Core utilities for molecular structure processing.

This module provides core functionality used across multiple modules:
1. Graph conversion utilities
2. Fingerprint computation
3. Shared molecular operations
"""

import logging
from typing import Dict, List, Optional, Union

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors
from torch_geometric.data import Data

logger = logging.getLogger(__name__)

# Atom feature dimensions
ATOM_FEATURES = {
    "atomic_num": list(range(1, 119)),  # Atomic numbers 1-118
    "degree": [0, 1, 2, 3, 4, 5, 6],
    "formal_charge": [-3, -2, -1, 0, 1, 2, 3],
    "chiral_tag": [0, 1, 2, 3],
    "hybridization": [
        Chem.rdchem.HybridizationType.SP,
        Chem.rdchem.HybridizationType.SP2,
        Chem.rdchem.HybridizationType.SP3,
        Chem.rdchem.HybridizationType.SP3D,
        Chem.rdchem.HybridizationType.SP3D2,
    ],
    "num_h": [0, 1, 2, 3, 4],
    "implicit_valence": [0, 1, 2, 3, 4, 5, 6],
    "aromatic": [0, 1],
}

# Bond feature dimensions
BOND_FEATURES = {
    "bond_type": [
        Chem.rdchem.BondType.SINGLE,
        Chem.rdchem.BondType.DOUBLE,
        Chem.rdchem.BondType.TRIPLE,
        Chem.rdchem.BondType.AROMATIC,
    ],
    "conjugated": [0, 1],
    "in_ring": [0, 1],
    "stereo": [0, 1, 2, 3, 4, 5],
}


def _one_hot(val, choices: list) -> list:
    """Create one-hot encoding."""
    encoding = [0] * len(choices)
    try:
        idx = choices.index(val)
        encoding[idx] = 1
    except:
        pass
    return encoding


def mol_to_graph(
    mol: Chem.Mol,
    add_hs: bool = True,
    compute_distances: bool = False,
    return_pyg: bool = False,
) -> Union[Dict[str, torch.Tensor], Data]:
    """Convert molecule to graph representation.

    Args:
        mol: Input molecule
        add_hs: Whether to add explicit hydrogens
        compute_distances: Whether to compute pairwise atomic distances
        return_pyg: Whether to return PyTorch Geometric Data object

    Returns:
        Either dictionary containing graph tensors or PyG Data object
    """
    if mol is None:
        return None

    if add_hs:
        mol = Chem.AddHs(mol)

    # Get atom features
    atom_features = []
    for atom in mol.GetAtoms():
        features = []
        features.append(_one_hot(atom.GetAtomicNum(), ATOM_FEATURES["atomic_num"]))
        features.append(_one_hot(atom.GetDegree(), ATOM_FEATURES["degree"]))
        features.append(_one_hot(atom.GetFormalCharge(), ATOM_FEATURES["formal_charge"]))
        features.append(_one_hot(atom.GetChiralTag(), ATOM_FEATURES["chiral_tag"]))
        features.append(_one_hot(atom.GetHybridization(), ATOM_FEATURES["hybridization"]))
        features.append(_one_hot(atom.GetTotalNumHs(), ATOM_FEATURES["num_h"]))
        features.append(_one_hot(atom.GetImplicitValence(), ATOM_FEATURES["implicit_valence"]))
        features.append([atom.GetIsAromatic()])
        atom_features.append(torch.cat([torch.tensor(f) for f in features]))

    # Get bond features
    edge_indices = []
    edge_features = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        edge_indices += [[i, j], [j, i]]

        features = []
        features.append(_one_hot(bond.GetBondType(), BOND_FEATURES["bond_type"]))
        features.append([bond.GetIsConjugated()])
        features.append([bond.IsInRing()])
        features.append(_one_hot(bond.GetStereo(), BOND_FEATURES["stereo"]))
        edge_features += [torch.cat([torch.tensor(f) for f in features])] * 2

    # Convert to tensors
    x = torch.stack(atom_features)
    edge_index = torch.tensor(edge_indices).t().contiguous()
    edge_attr = torch.stack(edge_features)

    # Optionally compute distance matrix
    distances = None
    if compute_distances and mol.GetNumConformers() > 0:
        pos = torch.tensor(mol.GetConformer().GetPositions())
        distances = torch.cdist(pos, pos)

    if return_pyg:
        data = Data(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_attr,
        )
        if distances is not None:
            data.distances = distances
        return data
    else:
        data = {
            "node_features": x,
            "edge_index": edge_index,
            "edge_features": edge_attr,
        }
        if distances is not None:
            data["distances"] = distances
        return data


def compute_fingerprints(
    mol: Chem.Mol,
    fp_types: Optional[List[str]] = None,
    as_tensor: bool = True,
) -> Union[Dict[str, np.ndarray], Dict[str, torch.Tensor]]:
    """Compute molecular fingerprints.

    Args:
        mol: Input molecule
        fp_types: Types of fingerprints to compute
        as_tensor: Whether to return PyTorch tensors

    Returns:
        Dictionary mapping fingerprint types to feature vectors
    """
    if mol is None:
        return None

    if fp_types is None:
        fp_types = [
            "morgan",
            "maccs",
            "topological",
            "atom_pairs",
            "torsion",
            "estate",
        ]

    fps = {}
    for fp_type in fp_types:
        try:
            if fp_type == "morgan":
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
                arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps["morgan"] = arr

            elif fp_type == "maccs":
                fp = AllChem.GetMACCSKeysFingerprint(mol)
                arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps["maccs"] = arr

            elif fp_type == "topological":
                fp = Chem.RDKFingerprint(mol, fpSize=2048)
                arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps["topological"] = arr

            elif fp_type == "atom_pairs":
                fp = AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=2048)
                arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps["atom_pairs"] = arr

            elif fp_type == "torsion":
                fp = AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=2048)
                arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps["torsion"] = arr

            elif fp_type == "estate":
                fps["estate"] = np.array(rdMolDescriptors.CalcEStateIndices(mol), dtype=np.float32)

        except Exception as e:
            logger.warning(f"Error computing {fp_type} fingerprint: {str(e)}")
            continue

    if as_tensor:
        return {k: torch.from_numpy(v) for k, v in fps.items()}
    return fps


"""Graph utilities for molecular structure processing.

This module provides utilities for:
1. Converting molecules to graph representations
2. Computing molecular fingerprints and descriptors
3. Graph feature extraction and processing
4. 3D structure analysis and manipulation
5. Molecular similarity calculations
"""

import logging
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    rdMolDescriptors,
    rdMolTransforms,
    rdShapeHelpers,
)
from torch_geometric.data import Data

logger = logging.getLogger(__name__)

# Atom feature dimensions
ATOM_FEATURES = {
    "atomic_num": list(range(1, 119)),  # Atomic numbers 1-118
    "degree": [0, 1, 2, 3, 4, 5, 6],
    "formal_charge": [-3, -2, -1, 0, 1, 2, 3],
    "chiral_tag": [0, 1, 2, 3],
    "hybridization": [
        Chem.rdchem.HybridizationType.SP,
        Chem.rdchem.HybridizationType.SP2,
        Chem.rdchem.HybridizationType.SP3,
        Chem.rdchem.HybridizationType.SP3D,
        Chem.rdchem.HybridizationType.SP3D2,
    ],
    "num_h": [0, 1, 2, 3, 4],
    "implicit_valence": [0, 1, 2, 3, 4, 5, 6],
    "aromatic": [0, 1],
}

# Bond feature dimensions
BOND_FEATURES = {
    "bond_type": [
        Chem.rdchem.BondType.SINGLE,
        Chem.rdchem.BondType.DOUBLE,
        Chem.rdchem.BondType.TRIPLE,
        Chem.rdchem.BondType.AROMATIC,
    ],
    "conjugated": [0, 1],
    "in_ring": [0, 1],
    "stereo": [0, 1, 2, 3, 4, 5],
}


def mol_to_graph(
    mol: Chem.Mol,
    add_hs: bool = True,
    compute_distances: bool = False,
    return_pyg: bool = False,
) -> Union[Dict[str, torch.Tensor], Data]:
    """Convert molecule to graph representation.

    Args:
        mol: Input molecule
        add_hs: Whether to add explicit hydrogens
        compute_distances: Whether to compute pairwise atomic distances
        return_pyg: Whether to return PyTorch Geometric Data object

    Returns:
        Either dictionary containing graph tensors or PyG Data object
    """
    if mol is None:
        return None

    if add_hs:
        mol = Chem.AddHs(mol)

    # Get atom features
    atom_features = []
    for atom in mol.GetAtoms():
        features = []
        features.append(_one_hot(atom.GetAtomicNum(), ATOM_FEATURES["atomic_num"]))
        features.append(_one_hot(atom.GetDegree(), ATOM_FEATURES["degree"]))
        features.append(_one_hot(atom.GetFormalCharge(), ATOM_FEATURES["formal_charge"]))
        features.append(_one_hot(atom.GetChiralTag(), ATOM_FEATURES["chiral_tag"]))
        features.append(_one_hot(atom.GetHybridization(), ATOM_FEATURES["hybridization"]))
        features.append(_one_hot(atom.GetTotalNumHs(), ATOM_FEATURES["num_h"]))
        features.append(_one_hot(atom.GetImplicitValence(), ATOM_FEATURES["implicit_valence"]))
        features.append([atom.GetIsAromatic()])
        atom_features.append(torch.cat([torch.tensor(f) for f in features]))

    # Get bond features
    edge_indices = []
    edge_features = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        edge_indices += [[i, j], [j, i]]

        features = []
        features.append(_one_hot(bond.GetBondType(), BOND_FEATURES["bond_type"]))
        features.append([bond.GetIsConjugated()])
        features.append([bond.IsInRing()])
        features.append(_one_hot(bond.GetStereo(), BOND_FEATURES["stereo"]))
        edge_features += [torch.cat([torch.tensor(f) for f in features])] * 2

    # Convert to tensors
    x = torch.stack(atom_features)
    edge_index = torch.tensor(edge_indices).t().contiguous()
    edge_attr = torch.stack(edge_features)

    # Optionally compute distance matrix
    distances = None
    if compute_distances and mol.GetNumConformers() > 0:
        pos = torch.tensor(mol.GetConformer().GetPositions())
        distances = torch.cdist(pos, pos)

    if return_pyg:
        data = Data(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_attr,
        )
        if distances is not None:
            data.distances = distances
        return data
    else:
        data = {
            "node_features": x,
            "edge_index": edge_index,
            "edge_features": edge_attr,
        }
        if distances is not None:
            data["distances"] = distances
        return data


def compute_fingerprints(
    mol: Chem.Mol,
    fp_types: Optional[List[str]] = None,
    as_tensor: bool = True,
) -> Union[Dict[str, np.ndarray], Dict[str, torch.Tensor]]:
    """Compute molecular fingerprints.

    Args:
        mol: Input molecule
        fp_types: Types of fingerprints to compute
        as_tensor: Whether to return PyTorch tensors

    Returns:
        Dictionary mapping fingerprint types to feature vectors
    """
    if mol is None:
        return None

    if fp_types is None:
        fp_types = [
            "morgan",
            "maccs",
            "topological",
            "atom_pairs",
            "torsion",
            "estate",
        ]

    fps = {}
    for fp_type in fp_types:
        try:
            if fp_type == "morgan":
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
                arr = np.zeros((1,))
                AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
                fps["morgan"] = arr

            elif fp_type == "maccs":
                fp = AllChem.GetMACCSKeysFingerprint(mol)
                arr = np.zeros((1,))
                AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
                fps["maccs"] = arr

            elif fp_type == "topological":
                fp = Chem.RDKFingerprint(mol, fpSize=2048)
                arr = np.zeros((1,))
                AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
                fps["topological"] = arr

            elif fp_type == "atom_pairs":
                fp = AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=2048)
                arr = np.zeros((1,))
                AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
                fps["atom_pairs"] = arr

            elif fp_type == "torsion":
                fp = AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=2048)
                arr = np.zeros((1,))
                AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
                fps["torsion"] = arr

            elif fp_type == "estate":
                fps["estate"] = np.array(rdMolDescriptors.CalcEStateIndices(mol), dtype=np.float32)

        except Exception as e:
            logger.warning(f"Error computing {fp_type} fingerprint: {str(e)}")
            continue

    if as_tensor:
        return {k: torch.from_numpy(v) for k, v in fps.items()}
    return fps


def compute_3d_descriptors(mol: Chem.Mol, conformer_id: int = -1) -> Dict[str, float]:
    """Compute 3D molecular descriptors.

    Args:
        mol: Input molecule with 3D coordinates
        conformer_id: Conformer ID to use

    Returns:
        Dictionary of 3D descriptors
    """
    if mol is None or mol.GetNumConformers() == 0:
        return {}

    try:
        # Shape descriptors
        descriptors = {
            "asphericity": rdShapeHelpers.ComputeASPH(mol, confId=conformer_id),
            "eccentricity": rdShapeHelpers.ComputeECCEN(mol, confId=conformer_id),
            "inertial_shape_factor": rdShapeHelpers.ComputeISF(mol, confId=conformer_id),
            "npr1": rdShapeHelpers.ComputeNPR1(mol, confId=conformer_id),
            "npr2": rdShapeHelpers.ComputeNPR2(mol, confId=conformer_id),
            "pmi1": rdShapeHelpers.ComputePMI1(mol, confId=conformer_id),
            "pmi2": rdShapeHelpers.ComputePMI2(mol, confId=conformer_id),
            "pmi3": rdShapeHelpers.ComputePMI3(mol, confId=conformer_id),
            "radius_of_gyration": rdShapeHelpers.ComputeRadiusOfGyration(mol, confId=conformer_id),
            "spherocity": rdShapeHelpers.ComputeSPH(mol, confId=conformer_id),
        }

        # Add volume and surface area if available
        if hasattr(rdShapeHelpers, "ComputeMolVolume"):
            descriptors["volume"] = rdShapeHelpers.ComputeMolVolume(mol, confId=conformer_id)
        if hasattr(rdShapeHelpers, "ComputeMolSurfaceArea"):
            descriptors["surface_area"] = rdShapeHelpers.ComputeMolSurfaceArea(mol, confId=conformer_id)

        return descriptors

    except Exception as e:
        logger.error(f"Error computing 3D descriptors: {str(e)}")
        return {}


def compute_shape_similarity(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    conf_id1: int = -1,
    conf_id2: int = -1,
) -> float:
    """Compute 3D shape similarity between molecules.

    Args:
        mol1: First molecule
        mol2: Second molecule
        conf_id1: Conformer ID for first molecule
        conf_id2: Conformer ID for second molecule

    Returns:
        Shape Tanimoto score
    """
    try:
        score = rdShapeHelpers.ShapeTanimotoDist(mol1, mol2, confId1=conf_id1, confId2=conf_id2)
        return float(score)
    except Exception as e:
        logger.error(f"Error computing shape similarity: {str(e)}")
        return 0.0


def apply_transformation(
    mol: Chem.Mol,
    matrix: Union[np.ndarray, torch.Tensor],
    conf_id: int = -1,
) -> bool:
    """Apply transformation matrix to molecule coordinates.

    Args:
        mol: Input molecule
        matrix: 4x4 transformation matrix
        conf_id: Conformer ID to transform

    Returns:
        Success status
    """
    try:
        if isinstance(matrix, torch.Tensor):
            matrix = matrix.numpy()

        conf = mol.GetConformer(conf_id)
        transform = rdMolTransforms.Transform3D()
        for i in range(3):
            for j in range(4):
                transform.SetElement(i, j, float(matrix[i, j]))
        rdMolTransforms.TransformConformer(conf, transform)
        return True

    except Exception as e:
        logger.error(f"Error applying transformation: {str(e)}")
        return False


def _one_hot(val, choices: list) -> list:
    """Create one-hot encoding."""
    encoding = [0] * len(choices)
    try:
        idx = choices.index(val)
        encoding[idx] = 1
    except:
        pass
    return encoding
