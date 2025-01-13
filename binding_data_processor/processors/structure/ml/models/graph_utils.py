"""Graph utilities for molecular structure processing.

This module provides utilities for:
1. Converting molecules to graph representations
2. Computing molecular fingerprints and descriptors
3. Graph feature extraction and processing
4. 3D structure analysis and manipulation
5. Molecular similarity calculations
"""

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D


def mol_to_graph(mol: Chem.Mol) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Convert RDKit molecule to graph representation.

    Args:
        mol: RDKit molecule

    Returns:
        Tuple of (node features, edge index, edge features)
    """
    # Get atom features
    atoms = mol.GetAtoms()
    num_atoms = len(atoms)

    # Node features
    atomic_nums = []
    aromatic = []
    hybridization = []
    num_hs = []

    for atom in atoms:
        atomic_nums.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        hybridization.append(atom.GetHybridization())
        num_hs.append(atom.GetTotalNumHs())

    # Convert to tensor
    x = torch.tensor([atomic_nums, aromatic, hybridization, num_hs], dtype=torch.float).t()

    # Get edge features
    src = []
    dst = []
    edge_features = []

    for bond in mol.GetBonds():
        # Get atoms in bond
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()

        # Add edges in both directions
        src += [start, end]
        dst += [end, start]

        # Get bond features
        bond_type = bond.GetBondType()
        is_conj = bond.GetIsConjugated()
        is_ring = bond.IsInRing()

        # Add same features for both directions
        edge_features.extend([[bond_type, is_conj, is_ring]] * 2)

    edge_index = torch.tensor([src, dst], dtype=torch.long)
    edge_attr = torch.tensor(edge_features, dtype=torch.float)

    return x, edge_index, edge_attr


def compute_fingerprints(mol: Chem.Mol, radius: int = 2, nBits: int = 2048) -> torch.Tensor:
    """Compute Morgan fingerprints for molecule.

    Args:
        mol: RDKit molecule
        radius: Fingerprint radius
        nBits: Number of bits in fingerprint

    Returns:
        Binary fingerprint tensor
    """
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    arr = np.zeros((1,))
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    return torch.from_numpy(arr).float()


def compute_3d_descriptors(mol: Chem.Mol) -> Dict[str, float]:
    """Compute 3D molecular descriptors.

    Args:
        mol: RDKit molecule with 3D coordinates

    Returns:
        Dictionary of descriptor names and values
    """
    # Ensure molecule has 3D coordinates
    if not mol.GetNumConformers():
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

    descriptors = {
        "asphericity": Descriptors3D.Asphericity(mol),
        "eccentricity": Descriptors3D.Eccentricity(mol),
        "inertial_shape_factor": Descriptors3D.InertialShapeFactor(mol),
        "npr1": Descriptors3D.NPR1(mol),
        "npr2": Descriptors3D.NPR2(mol),
        "pmi1": Descriptors3D.PMI1(mol),
        "pmi2": Descriptors3D.PMI2(mol),
        "pmi3": Descriptors3D.PMI3(mol),
        "radius_of_gyration": Descriptors3D.RadiusOfGyration(mol),
        "spherocity_index": Descriptors3D.SpherocityIndex(mol),
    }

    # Add surface area and volume descriptors
    try:
        descriptors.update(
            {
                "surface_area": AllChem.ComputeMolSurf(mol),
                "volume": AllChem.ComputeMolVolume(mol),
            }
        )
    except:
        # Some molecules may fail surface/volume calculations
        descriptors.update(
            {
                "surface_area": 0.0,
                "volume": 0.0,
            }
        )

    return descriptors
