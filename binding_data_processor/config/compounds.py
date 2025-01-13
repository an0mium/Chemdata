"""Main module for compound configuration and data access."""

from typing import Dict, List, Union, Any
from .compounds.proteins import PROTEINS
from .compounds.biomolecules import BIOMOLECULES
from .compounds.custom_compounds import CUSTOM_COMPOUNDS


def get_protein_data(protein_name: str) -> Union[str, List[str], None]:
    """Get UniProt ID(s) for a protein.

    Args:
        protein_name: Name of the protein to look up

    Returns:
        UniProt ID(s) for the protein, or None if not found
    """
    return PROTEINS.get(protein_name)


def get_biomolecule_data(name: str) -> Union[str, None]:
    """Get CAS number for a biomolecule.

    Args:
        name: Name of the biomolecule

    Returns:
        CAS number if found, None otherwise
    """
    return BIOMOLECULES.get(name)


def get_custom_compound_data(category: str, compound_id: str = None) -> Union[Dict[str, Any], List[Dict[str, Any]], None]:
    """Get data for custom compounds.

    Args:
        category: Category of compounds (e.g. "EMERGING_THREATS")
        compound_id: Optional specific compound ID within the category

    Returns:
        If compound_id provided: Dictionary of compound data if found, None otherwise
        If no compound_id: List of all compounds in category if found, None otherwise
    """
    if category not in CUSTOM_COMPOUNDS:
        return None

    if compound_id is None:
        return CUSTOM_COMPOUNDS[category]

    for compound in CUSTOM_COMPOUNDS[category]:
        if compound.get("name") == compound_id:
            return compound
    return None


def list_proteins() -> List[str]:
    """Get list of all protein names."""
    return list(PROTEINS.keys())


def list_biomolecules() -> List[str]:
    """Get list of all biomolecule names."""
    return list(BIOMOLECULES.keys())


def list_custom_compound_categories() -> List[str]:
    """Get list of all custom compound categories."""
    return list(CUSTOM_COMPOUNDS.keys())


def list_custom_compounds(category: str) -> List[str]:
    """Get list of all compound names in a category.

    Args:
        category: Category to list compounds from

    Returns:
        List of compound names in the category
    """
    if category not in CUSTOM_COMPOUNDS:
        return []
    return [c["name"] for c in CUSTOM_COMPOUNDS[category]]


def get_all_compounds() -> Dict[str, Any]:
    """Get complete compound data dictionary.

    Returns:
        Dictionary containing all protein, biomolecule and custom compound data
    """
    return {"proteins": PROTEINS, "biomolecules": BIOMOLECULES, "custom_compounds": CUSTOM_COMPOUNDS}
