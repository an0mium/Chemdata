"""Utilities for extracting and cleaning compound names from web sources."""

import re
from typing import Dict, List, Optional, Set, Tuple


def clean_name(name: str) -> str:
    """Clean a compound name by removing unwanted characters and standardizing format.

    Args:
        name: Raw compound name string

    Returns:
        Cleaned compound name string
    """
    if not name:
        return ""

    # Convert to lowercase
    name = name.lower()

    # Remove special characters but keep hyphens, numbers and important chemical notation
    name = re.sub(r"[^a-z0-9\-\+\(\)\[\]\{\}]", " ", name)

    # Replace multiple spaces with single space
    name = re.sub(r"\s+", " ", name)

    # Remove leading/trailing whitespace
    name = name.strip()

    return name


def standardize_name(name: str) -> str:
    """Standardize a compound name to a canonical form.

    Handles common variations in chemical nomenclature.

    Args:
        name: Compound name to standardize

    Returns:
        Standardized compound name
    """
    # Clean the name first
    std_name = clean_name(name)

    # Remove common prefixes
    prefixes = ["the ", "a ", "an "]
    for prefix in prefixes:
        if std_name.startswith(prefix):
            std_name = std_name[len(prefix) :]

    # Standardize salt forms
    salt_forms = {r"hydrochloride\b": "hcl", r"hydrogen chloride\b": "hcl", r"sulfate\b": "so4", r"sulphate\b": "so4", r"acetate\b": "oac", r"phosphate\b": "po4"}

    for pattern, replacement in salt_forms.items():
        std_name = re.sub(pattern, replacement, std_name)

    return std_name.strip()


def extract_identifiers(text: str) -> Tuple[Set[str], Set[str], Set[str], Set[str]]:
    """Extract chemical identifiers from text.

    Looks for:
    - Chemical names (e.g. "serotonin", "dopamine")
    - CAS numbers (e.g. "50-67-9")
    - InChI Keys (e.g. "QZAYGJVTTNCVMB-UHFFFAOYSA-N")
    - SMILES strings

    Args:
        text: Text to extract identifiers from

    Returns:
        Tuple of (chemical_names, cas_numbers, inchi_keys, smiles)
    """
    if not text:
        return set(), set(), set(), set()

    chemical_names: Set[str] = set()
    cas_numbers: Set[str] = set()
    inchi_keys: Set[str] = set()
    smiles: Set[str] = set()

    # Extract CAS numbers (e.g. 50-67-9)
    cas_pattern = r"\b\d{1,7}-\d{2}-\d\b"
    cas_numbers.update(re.findall(cas_pattern, text))

    # Extract InChI Keys (27 characters, all uppercase letters and numbers)
    inchi_pattern = r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b"
    inchi_keys.update(re.findall(inchi_pattern, text))

    # Extract potential SMILES strings
    smiles_pattern = r"\b[A-Za-z0-9@\+\-\[\]\(\)\{\}/\\=#$\.]+\b"
    for match in re.finditer(smiles_pattern, text):
        if looks_like_smiles(match.group()):
            smiles.add(match.group())

    # Extract chemical names
    # Split on common delimiters and clean each potential name
    potential_names = re.split(r"[;,\n\t]", text)
    for name in potential_names:
        cleaned = clean_name(name)
        if cleaned and len(cleaned) > 2:  # Avoid very short strings
            chemical_names.add(cleaned)

    return chemical_names, cas_numbers, inchi_keys, smiles


def looks_like_smiles(text: str) -> bool:
    """Basic check if a string looks like it could be a SMILES string.

    This is a basic heuristic check. For proper SMILES validation,
    use a chemistry toolkit like RDKit.

    Args:
        text: String to check

    Returns:
        Boolean indicating if the string looks like SMILES
    """
    if not text or len(text) < 2:
        return False

    # Should contain at least one atom symbol
    if not re.search(r"[CNOPS]", text):
        return False

    # Check for balanced parentheses/brackets
    if text.count("(") != text.count(")"):
        return False
    if text.count("[") != text.count("]"):
        return False

    # Should not contain invalid characters
    if re.search(r"[^A-Za-z0-9@\+\-\[\]\(\)\{\}/\\=#$\.]", text):
        return False

    # Should have valid atom symbols
    atoms = re.findall(r"[A-Z][a-z]?", text)
    valid_atoms = {"C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "B", "Si", "Se", "H"}
    return all(atom in valid_atoms for atom in atoms)


def get_name_variants(name: str) -> List[str]:
    """Generate common variants of a compound name.

    Args:
        name: Base compound name

    Returns:
        List of name variants
    """
    variants = set()

    # Add original name
    variants.add(name)

    # Add cleaned name
    cleaned = clean_name(name)
    variants.add(cleaned)

    # Add standardized name
    standardized = standardize_name(name)
    variants.add(standardized)

    # Handle common abbreviations and variations
    abbrev_map = {"hydrochloride": "hcl", "hydrogen chloride": "hcl", "sulfate": "so4", "sulphate": "so4", "acetate": "oac", "phosphate": "po4"}

    for full, abbrev in abbrev_map.items():
        if full in standardized:
            variants.add(standardized.replace(full, abbrev))

    # Handle common stereochemistry notations
    if "r-" in standardized:
        variants.add(standardized.replace("r-", "(r)-"))
    if "s-" in standardized:
        variants.add(standardized.replace("s-", "(s)-"))

    return sorted(list(variants))


def extract_chemical_properties(text: str) -> Dict[str, str]:
    """Extract chemical property information from text.

    Looks for common chemical properties like:
    - Molecular weight
    - Melting point
    - Boiling point
    - Density
    - LogP
    etc.

    Args:
        text: Text to extract properties from

    Returns:
        Dictionary of property names and values
    """
    properties = {}

    # Common property patterns
    patterns = {
        "molecular_weight": r"(?:molecular weight|mw|mol wt)[:\s]+(\d+\.?\d*)",
        "melting_point": r"(?:melting point|mp)[:\s]+(\d+\.?\d*)",
        "boiling_point": r"(?:boiling point|bp)[:\s]+(\d+\.?\d*)",
        "density": r"density[:\s]+(\d+\.?\d*)",
        "logp": r"(?:logp|log p)[:\s]+([-]?\d+\.?\d*)",
        "pka": r"pka[:\s]+([-]?\d+\.?\d*)",
        "ph": r"ph[:\s]+([-]?\d+\.?\d*)",
    }

    for prop, pattern in patterns.items():
        match = re.search(pattern, text.lower())
        if match:
            properties[prop] = match.group(1)

    return properties
