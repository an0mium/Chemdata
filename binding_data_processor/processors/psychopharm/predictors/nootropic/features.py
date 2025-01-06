"""Feature extraction utilities for nootropic prediction."""

from typing import Dict, List, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from .....models.core import CompoundData


def extract_fingerprints(mol: Chem.Mol) -> np.ndarray:
    """Extract molecular fingerprints."""
    # Morgan fingerprints (ECFP4)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return np.array(fp)


def extract_descriptors(mol: Chem.Mol) -> np.ndarray:
    """Extract molecular descriptors."""
    descriptors = []

    # Physical properties
    descriptors.extend(
        [
            Descriptors.ExactMolWt(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.TPSA(mol),
            Descriptors.MolLogP(mol),
        ]
    )

    # Topological descriptors
    descriptors.extend(
        [
            Descriptors.BertzCT(mol),
            Descriptors.Chi0n(mol),
            Descriptors.Chi1n(mol),
            Descriptors.Chi2n(mol),
            Descriptors.Chi3n(mol),
            Descriptors.Chi4n(mol),
        ]
    )

    # Constitutional descriptors
    descriptors.extend(
        [
            rdMolDescriptors.CalcNumRings(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
            rdMolDescriptors.CalcNumAliphaticRings(mol),
            rdMolDescriptors.CalcNumSaturatedRings(mol),
        ]
    )

    return np.array(descriptors)


def extract_enhanced_features(
    mol: Chem.Mol,
    fingerprints: Optional[np.ndarray] = None,
    descriptors: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Extract enhanced features combining fingerprints and descriptors."""
    if fingerprints is None:
        fingerprints = extract_fingerprints(mol)
    if descriptors is None:
        descriptors = extract_descriptors(mol)

    # Combine base features
    features = np.concatenate([fingerprints, descriptors])

    # Add interaction terms
    interactions = []
    for i in range(len(descriptors)):
        for j in range(i + 1, len(descriptors)):
            interactions.append(descriptors[i] * descriptors[j])

    # Add polynomial terms
    poly_terms = []
    for d in descriptors:
        poly_terms.extend([d**2, d**3])

    return np.concatenate(
        [
            features,
            np.array(interactions),
            np.array(poly_terms),
        ]
    )


def extract_all_features(
    compound: CompoundData,
    feature_types: List[str],
) -> Dict[str, np.ndarray]:
    """Extract all requested feature types for a compound."""
    mol = Chem.MolFromSmiles(compound.smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {compound.smiles}")

    features = {}
    fingerprints = None
    descriptors = None

    for feature_type in feature_types:
        if feature_type == "fingerprints":
            fingerprints = extract_fingerprints(mol)
            features[feature_type] = fingerprints

        elif feature_type == "descriptors":
            descriptors = extract_descriptors(mol)
            features[feature_type] = descriptors

        elif feature_type == "enhanced":
            features[feature_type] = extract_enhanced_features(mol, fingerprints, descriptors)

        else:
            raise ValueError(f"Unknown feature type: {feature_type}")

    return features
