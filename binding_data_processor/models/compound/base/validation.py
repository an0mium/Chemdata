"""Validation utilities for compound data.

This module provides validation functions for:
1. Chemical identifiers (SMILES, InChI, CAS)
2. Database IDs (PubChem, ChEMBL, etc.)
3. Chemical properties (MW, LogP, etc.)
4. Structural features
"""

import re
from typing import Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from .types import CompoundType, LegalStatus, PsychoactiveClass, NootropicMechanism, BBBPermeability


class ValidationError(Exception):
    """Raised when compound data validation fails."""

    pass


def validate_smiles(smiles: str) -> Tuple[bool, Optional[str]]:
    """
    Validate SMILES string using RDKit.

    Args:
        smiles: SMILES string to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not smiles or not isinstance(smiles, str):
        return False, "Empty or invalid SMILES string"

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES format"

        # Validate atom types
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                return False, f"Invalid atom type at position {atom.GetIdx()}"

        # Validate valence
        try:
            Chem.SanitizeMol(mol)
        except ValueError as e:
            return False, f"Structure validation failed: {str(e)}"

        return True, None

    except Exception as e:
        return False, f"SMILES parsing failed: {str(e)}"


def validate_inchi(inchi: str) -> Tuple[bool, Optional[str]]:
    """
    Validate InChI string using RDKit.

    Args:
        inchi: InChI string to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not inchi or not isinstance(inchi, str):
        return False, "Empty or invalid InChI string"

    if not inchi.startswith("InChI=1S/"):
        return False, "Invalid InChI format - must start with 'InChI=1S/'"

    try:
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            return False, "Invalid InChI format"

        # Validate structure
        try:
            Chem.SanitizeMol(mol)
        except ValueError as e:
            return False, f"Structure validation failed: {str(e)}"

        return True, None

    except Exception as e:
        return False, f"InChI parsing failed: {str(e)}"


def validate_cas_number(cas: str) -> Tuple[bool, Optional[str]]:
    """
    Validate CAS Registry Number format and checksum.

    Args:
        cas: CAS number to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not cas or not isinstance(cas, str):
        return False, "Empty or invalid CAS number"

    # Validate format
    pattern = r"^\d{1,7}-\d{2}-\d$"
    if not re.match(pattern, cas):
        return False, "Invalid CAS number format"

    # Validate checksum
    try:
        numbers = cas.replace("-", "")
        check_digit = int(numbers[-1])
        numbers = numbers[:-1]

        total = sum(int(num) * (i + 1) for i, num in enumerate(reversed(numbers)))

        if (total % 10) != check_digit:
            return False, "Invalid CAS number checksum"

        return True, None

    except Exception as e:
        return False, f"CAS number validation failed: {str(e)}"


def validate_pubchem_id(cid: str) -> Tuple[bool, Optional[str]]:
    """
    Validate PubChem Compound ID format.

    Args:
        cid: PubChem CID to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not cid or not isinstance(cid, str):
        return False, "Empty or invalid PubChem CID"

    try:
        cid_int = int(cid)
        if cid_int <= 0:
            return False, "PubChem CID must be positive"
        return True, None
    except ValueError:
        return False, "PubChem CID must be numeric"


def validate_chembl_id(chembl_id: str) -> Tuple[bool, Optional[str]]:
    """
    Validate ChEMBL ID format.

    Args:
        chembl_id: ChEMBL ID to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not chembl_id or not isinstance(chembl_id, str):
        return False, "Empty or invalid ChEMBL ID"

    pattern = r"^CHEMBL\d+$"
    if not re.match(pattern, chembl_id):
        return False, "Invalid ChEMBL ID format"

    return True, None


def validate_molecular_weight(mw: float) -> Tuple[bool, Optional[str]]:
    """
    Validate molecular weight value.

    Args:
        mw: Molecular weight to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(mw, (int, float)):
        return False, "Molecular weight must be numeric"

    if mw <= 0:
        return False, "Molecular weight must be positive"

    if mw > 5000:  # Arbitrary upper limit
        return False, "Molecular weight seems too large"

    return True, None


def validate_logp(logp: float) -> Tuple[bool, Optional[str]]:
    """
    Validate LogP value.

    Args:
        logp: LogP value to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(logp, (int, float)):
        return False, "LogP must be numeric"

    if abs(logp) > 20:  # Arbitrary limit
        return False, "LogP value seems unrealistic"

    return True, None


def validate_structure(mol: Chem.Mol) -> Tuple[bool, Optional[str]]:
    """
    Validate chemical structure using RDKit.

    Args:
        mol: RDKit molecule to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if mol is None:
        return False, "Invalid molecule"

    try:
        # Check for valid atoms
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                return False, f"Invalid atom type at position {atom.GetIdx()}"

        # Validate valence
        try:
            Chem.SanitizeMol(mol)
        except ValueError as e:
            return False, f"Structure validation failed: {str(e)}"

        # Check for disconnected fragments
        if len(Chem.GetMolFrags(mol)) > 1:
            return False, "Structure contains disconnected fragments"

        # Check for unusual bond types
        for bond in mol.GetBonds():
            if bond.GetBondType() not in [
                Chem.BondType.SINGLE,
                Chem.BondType.DOUBLE,
                Chem.BondType.TRIPLE,
                Chem.BondType.AROMATIC,
            ]:
                return False, f"Unusual bond type at position {bond.GetIdx()}"

        return True, None

    except Exception as e:
        return False, f"Structure validation failed: {str(e)}"


def validate_properties(mol: Chem.Mol) -> Tuple[bool, Optional[str]]:
    """
    Validate chemical properties using RDKit.

    Args:
        mol: RDKit molecule to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if mol is None:
        return False, "Invalid molecule"

    try:
        # Calculate basic properties
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)

        # Validate ranges
        if not 0 < mw < 5000:
            return False, f"Suspicious molecular weight: {mw}"

        if abs(logp) > 20:
            return False, f"Suspicious LogP value: {logp}"

        if tpsa < 0:
            return False, f"Invalid TPSA value: {tpsa}"

        if hbd < 0 or hba < 0:
            return False, "Invalid H-bond donor/acceptor count"

        if rotatable < 0:
            return False, "Invalid rotatable bond count"

        return True, None

    except Exception as e:
        return False, f"Property calculation failed: {str(e)}"


def validate_compound_type(compound_type: str) -> Tuple[bool, Optional[str]]:
    """
    Validate compound type enum value.

    Args:
        compound_type: Compound type to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        CompoundType(compound_type)
        return True, None
    except ValueError:
        return False, f"Invalid compound type: {compound_type}"


def validate_legal_status(status: str) -> Tuple[bool, Optional[str]]:
    """
    Validate legal status enum value.

    Args:
        status: Legal status to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        LegalStatus(status)
        return True, None
    except ValueError:
        return False, f"Invalid legal status: {status}"


def validate_psychoactive_class(pclass: str) -> Tuple[bool, Optional[str]]:
    """
    Validate psychoactive class enum value.

    Args:
        pclass: Psychoactive class to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        PsychoactiveClass(pclass)
        return True, None
    except ValueError:
        return False, f"Invalid psychoactive class: {pclass}"


def validate_nootropic_mechanism(mechanism: str) -> Tuple[bool, Optional[str]]:
    """
    Validate nootropic mechanism enum value.

    Args:
        mechanism: Nootropic mechanism to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        NootropicMechanism(mechanism)
        return True, None
    except ValueError:
        return False, f"Invalid nootropic mechanism: {mechanism}"


def validate_bbb_permeability(permeability: str) -> Tuple[bool, Optional[str]]:
    """
    Validate BBB permeability enum value.

    Args:
        permeability: BBB permeability to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        BBBPermeability(permeability)
        return True, None
    except ValueError:
        return False, f"Invalid BBB permeability: {permeability}"


def validate_confidence(confidence: float) -> Tuple[bool, Optional[str]]:
    """
    Validate confidence score.

    Args:
        confidence: Confidence score to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(confidence, (int, float)):
        return False, "Confidence must be numeric"

    if not 0 <= confidence <= 1:
        return False, "Confidence must be between 0 and 1"

    return True, None


def validate_probability(prob: float) -> Tuple[bool, Optional[str]]:
    """
    Validate probability value.

    Args:
        prob: Probability value to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(prob, (int, float)):
        return False, "Probability must be numeric"

    if not 0 <= prob <= 1:
        return False, "Probability must be between 0 and 1"

    return True, None


def validate_affinity_value(value: float) -> Tuple[bool, Optional[str]]:
    """
    Validate binding affinity value.

    Args:
        value: Affinity value to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(value, (int, float)):
        return False, "Affinity value must be numeric"

    if value < 0:
        return False, "Affinity value must be non-negative"

    if value > 1e6:  # Arbitrary upper limit
        return False, "Affinity value seems too large"

    return True, None


def validate_affinity_type(atype: str) -> Tuple[bool, Optional[str]]:
    """
    Validate binding affinity type.

    Args:
        atype: Affinity type to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    valid_types = {"Ki", "IC50", "EC50", "Kd"}
    if atype not in valid_types and atype != "N/A":
        return False, f"Invalid affinity type: {atype}"
    return True, None


def validate_affinity_unit(unit: str) -> Tuple[bool, Optional[str]]:
    """
    Validate binding affinity unit.

    Args:
        unit: Affinity unit to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    valid_units = {"nM", "uM", "mM", "pM"}
    if unit not in valid_units and unit != "N/A":
        return False, f"Invalid affinity unit: {unit}"
    return True, None
