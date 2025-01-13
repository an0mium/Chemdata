"""Features module for toxicity prediction.

This module provides functionality for extracting and processing features
related to toxicity prediction, including:
- Molecular descriptors
- Structural features
- Chemical properties
- Toxicophores
"""

import logging
from typing import Dict, List, Optional, Tuple, Union
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

logger = logging.getLogger(__name__)

# Toxicophore SMARTS patterns
TOXICOPHORES = {
    "michael_acceptor": "[C;H0,H1]=[C;H2]C(=O)[!N]",
    "epoxide": "C1OC1",
    "aziridine": "C1NC1",
    "halogen": "[Cl,Br,I,F]",
    "aldehyde": "[CH;D2]=O",
    "hydrazine": "[NH][NH2]",
    "nitro": "[N+](=O)[O-]",
    "nitroso": "N=O",
    "peroxide": "OO",
    "thiocyanate": "SC#N",
    "isocyanate": "N=C=O",
    "isothiocyanate": "N=C=S",
    "acyl_halide": "C(=O)[Cl,Br,I,F]",
    "sulfonyl_halide": "S(=O)(=O)[Cl,Br,I,F]",
    "phosphoryl_halide": "P(=O)[Cl,Br,I,F]",
}


def calculate_molecular_descriptors(mol: Chem.Mol) -> Dict[str, float]:
    """Calculate molecular descriptors for toxicity prediction.

    Args:
        mol: RDKit molecule object

    Returns:
        Dictionary of descriptor names and values
    """
    try:
        descriptors = {
            "MW": Descriptors.ExactMolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "TPSA": Descriptors.TPSA(mol),
            "HBA": rdMolDescriptors.CalcNumHBA(mol),
            "HBD": rdMolDescriptors.CalcNumHBD(mol),
            "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "Rings": rdMolDescriptors.CalcNumRings(mol),
            "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "HeteroAtoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
            "SaturatedRings": rdMolDescriptors.CalcNumSaturatedRings(mol),
            "AliphaticRings": rdMolDescriptors.CalcNumAliphaticRings(mol),
            "HeavyAtoms": mol.GetNumHeavyAtoms(),
        }

        # Add topological descriptors
        descriptors.update(
            {
                "BertzCT": Descriptors.BertzCT(mol),
                "Chi0v": Descriptors.Chi0v(mol),
                "Chi1v": Descriptors.Chi1v(mol),
                "Chi2v": Descriptors.Chi2v(mol),
                "Chi3v": Descriptors.Chi3v(mol),
                "Chi4v": Descriptors.Chi4v(mol),
                "HallKierAlpha": Descriptors.HallKierAlpha(mol),
                "Kappa1": Descriptors.Kappa1(mol),
                "Kappa2": Descriptors.Kappa2(mol),
                "Kappa3": Descriptors.Kappa3(mol),
            }
        )

        return descriptors

    except Exception as e:
        logger.error(f"Error calculating molecular descriptors: {str(e)}")
        return {}


def detect_toxicophores(mol: Chem.Mol) -> Dict[str, int]:
    """Detect toxicophore patterns in molecule.

    Args:
        mol: RDKit molecule object

    Returns:
        Dictionary mapping toxicophore names to counts
    """
    try:
        results = {}
        for name, smarts in TOXICOPHORES.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None:
                matches = mol.GetSubstructMatches(pattern)
                results[name] = len(matches)
            else:
                logger.warning(f"Invalid SMARTS pattern for {name}")
                results[name] = 0
        return results

    except Exception as e:
        logger.error(f"Error detecting toxicophores: {str(e)}")
        return {name: 0 for name in TOXICOPHORES}


def calculate_fingerprints(
    mol: Chem.Mol,
    fp_type: str = "morgan",
    radius: int = 2,
    n_bits: int = 2048,
) -> np.ndarray:
    """Calculate molecular fingerprints.

    Args:
        mol: RDKit molecule object
        fp_type: Type of fingerprint ("morgan", "maccs", "topological")
        radius: Radius for Morgan fingerprints
        n_bits: Number of bits in fingerprint

    Returns:
        Numpy array of fingerprint bits
    """
    try:
        if fp_type == "morgan":
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        elif fp_type == "maccs":
            fp = AllChem.GetMACCSKeysFingerprint(mol)
        elif fp_type == "topological":
            fp = Chem.RDKFingerprint(mol, fpSize=n_bits)
        else:
            raise ValueError(f"Unknown fingerprint type: {fp_type}")

        return np.array(fp)

    except Exception as e:
        logger.error(f"Error calculating fingerprints: {str(e)}")
        return np.zeros(n_bits)


def extract_descriptors(mol: Chem.Mol, enhanced: bool = False) -> Dict[str, float]:
    """Extract molecular descriptors for toxicity prediction.

    This is a wrapper around calculate_molecular_descriptors for consistency
    with other feature extraction functions.

    Args:
        mol: RDKit molecule object
        enhanced: Whether to include enhanced descriptors

    Returns:
        Dictionary of descriptor names and values
    """
    descriptors = calculate_molecular_descriptors(mol)

    if enhanced:
        # Add enhanced descriptors
        descriptors.update(
            {
                "FractionCSP3": Descriptors.FractionCSP3(mol),
                "NumSpiroAtoms": rdMolDescriptors.CalcNumSpiroAtoms(mol),
                "NumBridgeheadAtoms": rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
                "NumAtomStereoCenters": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
                "NumUnspecifiedAtomStereoCenters": rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol),
                "RingCount": rdMolDescriptors.CalcNumRings(mol),
                "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
                "NumSaturatedRings": rdMolDescriptors.CalcNumSaturatedRings(mol),
                "NumAliphaticRings": rdMolDescriptors.CalcNumAliphaticRings(mol),
                "NumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                "NumHAcceptors": rdMolDescriptors.CalcNumHBA(mol),
                "NumHDonors": rdMolDescriptors.CalcNumHBD(mol),
                "MolLogP": Descriptors.MolLogP(mol),
                "MolMR": Descriptors.MolMR(mol),
                "TPSA": Descriptors.TPSA(mol),
                "LabuteASA": Descriptors.LabuteASA(mol),
                "PEOE_VSA1": Descriptors.PEOE_VSA1(mol),
                "PEOE_VSA2": Descriptors.PEOE_VSA2(mol),
                "PEOE_VSA3": Descriptors.PEOE_VSA3(mol),
                "SMR_VSA1": Descriptors.SMR_VSA1(mol),
                "SMR_VSA2": Descriptors.SMR_VSA2(mol),
                "SMR_VSA3": Descriptors.SMR_VSA3(mol),
                "SlogP_VSA1": Descriptors.SlogP_VSA1(mol),
                "SlogP_VSA2": Descriptors.SlogP_VSA2(mol),
                "SlogP_VSA3": Descriptors.SlogP_VSA3(mol),
            }
        )

    return descriptors


def extract_fingerprints(
    mol: Chem.Mol,
    fp_type: str = "morgan",
    radius: int = 2,
    n_bits: int = 2048,
) -> np.ndarray:
    """Extract molecular fingerprints (wrapper around calculate_fingerprints).

    Args:
        mol: RDKit molecule object
        fp_type: Type of fingerprint ("morgan", "maccs", "topological")
        radius: Radius for Morgan fingerprints
        n_bits: Number of bits in fingerprint

    Returns:
        Numpy array of fingerprint bits
    """
    return calculate_fingerprints(mol, fp_type, radius, n_bits)


def extract_features(
    mol: Chem.Mol,
    include_descriptors: bool = True,
    include_toxicophores: bool = True,
    include_fingerprints: bool = True,
    fp_type: str = "morgan",
    fp_radius: int = 2,
    fp_bits: int = 2048,
) -> Dict[str, Union[float, int, np.ndarray]]:
    """Extract all features for toxicity prediction.

    Args:
        mol: RDKit molecule object
        include_descriptors: Whether to include molecular descriptors
        include_toxicophores: Whether to include toxicophore counts
        include_fingerprints: Whether to include fingerprints
        fp_type: Type of fingerprint
        fp_radius: Radius for Morgan fingerprints
        fp_bits: Number of fingerprint bits

    Returns:
        Dictionary containing all requested features
    """
    features = {}

    if include_descriptors:
        features.update(calculate_molecular_descriptors(mol))

    if include_toxicophores:
        features.update(detect_toxicophores(mol))

    if include_fingerprints:
        features["fingerprint"] = calculate_fingerprints(mol, fp_type, fp_radius, fp_bits)

    return features


def extract_all_features(mol: Chem.Mol, enhanced: bool = False) -> Dict[str, Union[float, int, np.ndarray]]:
    """Extract all possible features for toxicity prediction.

    This is a convenience wrapper around extract_features that enables all feature types.

    Args:
        mol: RDKit molecule object
        enhanced: Whether to include enhanced features

    Returns:
        Dictionary containing all available features
    """
    features = extract_features(
        mol,
        include_descriptors=True,
        include_toxicophores=True,
        include_fingerprints=True,
    )

    if enhanced:
        # Add enhanced descriptors
        enhanced_descriptors = extract_descriptors(mol, enhanced=True)
        features.update(enhanced_descriptors)

        # Add enhanced fingerprints
        features["ecfp6"] = calculate_fingerprints(mol, fp_type="morgan", radius=3, n_bits=2048)
        features["fcfp4"] = calculate_fingerprints(mol, fp_type="morgan", radius=2, n_bits=2048, useFeatures=True)
        features["maccs"] = calculate_fingerprints(mol, fp_type="maccs")

        # Add enhanced toxicophores
        features.update(detect_toxicophores(mol))

    return features


def extract_enhanced_features(mol: Chem.Mol) -> Dict[str, Union[float, int, np.ndarray]]:
    """Extract enhanced features for toxicity prediction.

    This is a convenience wrapper around extract_all_features with enhanced=True.

    Args:
        mol: RDKit molecule object

    Returns:
        Dictionary containing all enhanced features
    """
    return extract_all_features(mol, enhanced=True)
