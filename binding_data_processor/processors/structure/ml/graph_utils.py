"""Graph utilities for molecular data processing."""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import torch
from rdkit import Chem
from torch_geometric.data import Data

logger = logging.getLogger(__name__)


def mol_to_graph(mol: Chem.Mol) -> Data:
    """Convert RDKit molecule to PyTorch Geometric graph.

    Args:
        mol: RDKit molecule

    Returns:
        PyTorch Geometric Data object
    """
    try:
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

    except Exception as e:
        logger.error(f"Error converting molecule to graph: {str(e)}")
        return None


def compute_fingerprints(mol: Chem.Mol) -> Dict[str, np.ndarray]:
    """Compute molecular fingerprints.

    Args:
        mol: RDKit molecule

    Returns:
        Dictionary of fingerprint arrays
    """
    try:
        from rdkit.Chem import AllChem, MACCSkeys

        fingerprints = {}

        # Morgan fingerprint
        morgan = AllChem.GetMorganFingerprintAsBitVect(
            mol,
            radius=3,
            nBits=2048,
            useChirality=True,
            useFeatures=True,
        )
        fingerprints["morgan"] = np.array(list(morgan))

        # MACCS keys
        maccs = MACCSkeys.GenMACCSKeys(mol)
        fingerprints["maccs"] = np.array(list(maccs))

        # RDKit fingerprint
        rdkit = Chem.RDKFingerprint(
            mol,
            fpSize=2048,
            minPath=1,
            maxPath=7,
            useHs=True,
        )
        fingerprints["rdkit"] = np.array(list(rdkit))

        # Pattern fingerprint
        pattern = Chem.PatternFingerprint(
            mol,
            fpSize=2048,
            tautomerFingerprints=True,
        )
        fingerprints["pattern"] = np.array(list(pattern))

        # Layered fingerprint
        layered = Chem.LayeredFingerprint(
            mol,
            fpSize=2048,
            layerFlags=0xFFFFFFFF,
        )
        fingerprints["layered"] = np.array(list(layered))

        return fingerprints

    except Exception as e:
        logger.error(f"Error computing fingerprints: {str(e)}")
        return {}
