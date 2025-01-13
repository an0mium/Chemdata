"""Feature extraction utilities for nootropic prediction.

This module provides comprehensive feature extraction for nootropic effect prediction,
including:
1. Molecular fingerprints and descriptors
2. Pharmacophore features
3. Binding site features
4. Literature-derived features
5. Community data features
6. Enhanced interaction and polynomial features
"""

import logging
from typing import Dict, List, Optional, Any, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

from .....models.core import CompoundData


def extract_nootropic_features(
    compound: CompoundData,
    feature_types: List[str],
    models: Optional[Dict[str, Dict[str, Any]]] = None,
    device: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """Extract all requested feature types for a compound.

    Args:
        compound: Compound to extract features for
        feature_types: Types of features to extract
        models: Optional dictionary of loaded models for model-based features
        device: Optional device specification for model-based features

    Returns:
        Dictionary of feature names to feature arrays
    """
    logger = logging.getLogger(__name__)
    features = {}

    try:
        mol = Chem.MolFromSmiles(compound.smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {compound.smiles}")

        # Track intermediate results to avoid recomputation
        fingerprints = None
        descriptors = None
        pharmacophore = None

        for feature_type in feature_types:
            if feature_type == "molecular":
                fingerprints = extract_fingerprints(mol)
                descriptors = extract_descriptors(mol)
                features["molecular"] = np.concatenate([fingerprints, descriptors])

            elif feature_type == "pharmacophore":
                pharmacophore = extract_pharmacophore_features(mol)
                features["pharmacophore"] = pharmacophore

            elif feature_type == "binding":
                binding = extract_binding_features(mol, models, device)
                features["binding"] = binding

            elif feature_type == "literature":
                literature = extract_literature_features(compound, models)
                features["literature"] = literature

            elif feature_type == "community":
                community = extract_community_features(compound, models)
                features["community"] = community

            elif feature_type == "enhanced":
                enhanced = extract_features(
                    mol,
                    fingerprints=fingerprints,
                    descriptors=descriptors,
                    pharmacophore=pharmacophore,
                )
                features["enhanced"] = enhanced

            else:
                logger.warning(f"Unknown feature type: {feature_type}")

    except Exception as e:
        logger.error(f"Error extracting features: {str(e)}")
        raise

    return features


def extract_fingerprints(mol: Chem.Mol) -> np.ndarray:
    """Extract molecular fingerprints.

    Extracts multiple fingerprint types and combines them:
    1. Morgan fingerprints (ECFP4)
    2. MACCS keys
    3. Topological torsion fingerprints
    4. Atom pair fingerprints
    """
    # Morgan fingerprints (ECFP4)
    ecfp4 = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))

    # MACCS keys
    maccs = np.array(AllChem.GetMACCSKeysFingerprint(mol))

    # Topological torsion
    tt = np.array(AllChem.GetTopologicalTorsionFingerprint(mol))

    # Atom pairs
    ap = np.array(AllChem.GetAtomPairFingerprint(mol))

    return np.concatenate([ecfp4, maccs, tt[:100], ap[:100]])  # Truncate to control dimensionality


def extract_descriptors(mol: Chem.Mol) -> np.ndarray:
    """Extract molecular descriptors.

    Extracts comprehensive molecular descriptors including:
    1. Physical properties
    2. Topological descriptors
    3. Constitutional descriptors
    4. Electronic descriptors
    5. Geometric descriptors (if 3D conformer available)
    """
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
            Descriptors.MolMR(mol),
            Descriptors.FractionCSP3(mol),
            Descriptors.NumAliphaticRings(mol),
            Descriptors.NumAromaticRings(mol),
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
            Descriptors.HallKierAlpha(mol),
            Descriptors.Kappa1(mol),
            Descriptors.Kappa2(mol),
            Descriptors.Kappa3(mol),
        ]
    )

    # Constitutional descriptors
    descriptors.extend(
        [
            rdMolDescriptors.CalcNumRings(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
            rdMolDescriptors.CalcNumAliphaticRings(mol),
            rdMolDescriptors.CalcNumSaturatedRings(mol),
            rdMolDescriptors.CalcNumHeterocycles(mol),
            rdMolDescriptors.CalcNumSpiroAtoms(mol),
            rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
        ]
    )

    # Electronic descriptors
    descriptors.extend(
        [
            Descriptors.MaxPartialCharge(mol),
            Descriptors.MinPartialCharge(mol),
            Descriptors.MaxAbsPartialCharge(mol),
            Descriptors.MinAbsPartialCharge(mol),
        ]
    )

    # Geometric descriptors if 3D conformer available
    if mol.GetNumConformers() > 0:
        descriptors.extend(
            [
                Descriptors.NPR1(mol),
                Descriptors.NPR2(mol),
                Descriptors.RadiusOfGyration(mol),
                Descriptors.InertialShapeFactor(mol),
                Descriptors.Asphericity(mol),
                Descriptors.Eccentricity(mol),
            ]
        )

    return np.array(descriptors)


def extract_pharmacophore_features(mol: Chem.Mol) -> np.ndarray:
    """Extract pharmacophore features.

    Uses Gobbi 2D pharmacophore fingerprints that encode:
    1. Hydrogen bond donors/acceptors
    2. Basic/acidic groups
    3. Aromatic rings
    4. Hydrophobic regions
    5. Charged groups
    """
    # Generate 3D conformer if needed
    if not mol.GetNumConformers():
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

    # Generate Gobbi 2D pharmacophore fingerprint
    pharm2d = Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)
    return np.array(pharm2d)


def extract_binding_features(
    mol: Chem.Mol,
    models: Optional[Dict[str, Dict[str, Any]]] = None,
    device: Optional[str] = None,
) -> np.ndarray:
    """Extract binding site features.

    Uses loaded models to predict:
    1. Binding site interactions
    2. Protein-ligand contacts
    3. Binding pocket complementarity
    4. Key interaction points
    """
    if not models or "binding" not in models:
        return np.array([])

    try:
        binding_features = []
        for model_name, model in models["binding"].items():
            if device:
                model = model.to(device)
            pred = model.predict_binding(mol)
            binding_features.extend(pred)
        return np.array(binding_features)
    except Exception as e:
        logger.error(f"Error extracting binding features: {str(e)}")
        return np.array([])


def extract_literature_features(
    compound: CompoundData,
    models: Optional[Dict[str, Dict[str, Any]]] = None,
) -> np.ndarray:
    """Extract literature-derived features.

    Analyzes literature data for:
    1. Reported mechanisms
    2. Clinical effects
    3. Safety profiles
    4. Structure-activity relationships
    """
    if not models or "literature" not in models:
        return np.array([])

    try:
        literature_features = []
        for model_name, model in models["literature"].items():
            pred = model.analyze_literature(compound)
            literature_features.extend(pred)
        return np.array(literature_features)
    except Exception as e:
        logger.error(f"Error extracting literature features: {str(e)}")
        return np.array([])


def extract_community_features(
    compound: CompoundData,
    models: Optional[Dict[str, Dict[str, Any]]] = None,
) -> np.ndarray:
    """Extract community data features.

    Analyzes community reports for:
    1. Subjective effects
    2. Usage patterns
    3. Safety reports
    4. Interaction reports
    """
    if not models or "community" not in models:
        return np.array([])

    try:
        community_features = []
        for model_name, model in models["community"].items():
            pred = model.analyze_community_data(compound)
            community_features.extend(pred)
        return np.array(community_features)
    except Exception as e:
        logger.error(f"Error extracting community features: {str(e)}")
        return np.array([])


def extract_features(
    mol: Chem.Mol,
    fingerprints: Optional[np.ndarray] = None,
    descriptors: Optional[np.ndarray] = None,
    pharmacophore: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Extract enhanced features combining multiple feature types.

    Generates:
    1. Interaction terms between descriptors
    2. Polynomial terms of descriptors
    3. Combined fingerprint and pharmacophore features
    4. Derived features from multiple sources
    """
    features = []

    # Get base features if not provided
    if fingerprints is None:
        fingerprints = extract_fingerprints(mol)
    if descriptors is None:
        descriptors = extract_descriptors(mol)
    if pharmacophore is None:
        pharmacophore = extract_pharmacophore_features(mol)

    # Combine base features
    features.extend(fingerprints)
    features.extend(descriptors)
    features.extend(pharmacophore)

    # Add interaction terms between descriptors
    for i in range(len(descriptors)):
        for j in range(i + 1, len(descriptors)):
            features.append(descriptors[i] * descriptors[j])

    # Add polynomial terms of descriptors
    for d in descriptors:
        features.extend([d**2, d**3])

    # Add derived features
    features.extend(
        [
            np.mean(fingerprints),
            np.std(fingerprints),
            np.mean(descriptors),
            np.std(descriptors),
            np.mean(pharmacophore),
            np.std(pharmacophore),
        ]
    )

    return np.array(features)
