"""Protein interface analysis functionality."""

from typing import Dict, List, Any, Optional, Tuple
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue

from ...binding.structure.geometry import calculate_surface_area
from .utils import get_ca_atoms


def analyze_interfaces(
    structure: Structure,
    chains: Optional[List[str]] = None,
    cutoff: float = 5.0,
    min_contacts: int = 5,
) -> Dict[str, Any]:
    """Analyze protein chain interfaces.

    Args:
        structure: BioPython Structure object
        chains: Optional list of chain IDs to analyze. If None, analyze all chains.
        cutoff: Distance cutoff for interface contacts in Angstroms
        min_contacts: Minimum number of contacts to consider an interface

    Returns:
        Dictionary containing:
        - interfaces: List of interface descriptions
        - contact_maps: Contact maps for each interface
        - areas: Interface surface areas
        - statistics: Overall interface statistics
    """
    try:
        # Get chains to analyze
        if chains is None:
            chains = [chain.id for chain in structure.get_chains()]

        interfaces = []
        contact_maps = {}
        areas = {}

        # Analyze each chain pair
        for i, chain1 in enumerate(chains[:-1]):
            for chain2 in chains[i + 1 :]:
                interface = _analyze_chain_interface(
                    structure,
                    chain1,
                    chain2,
                    cutoff=cutoff,
                    min_contacts=min_contacts,
                )

                if interface:
                    interface_id = f"{chain1}_{chain2}"
                    interfaces.append(interface)
                    contact_maps[interface_id] = interface["contacts"]
                    areas[interface_id] = interface["area"]

        # Calculate overall statistics
        statistics = _calculate_interface_statistics(interfaces)

        return {
            "interfaces": interfaces,
            "contact_maps": contact_maps,
            "areas": areas,
            "statistics": statistics,
        }

    except Exception as e:
        logger.error(f"Error analyzing interfaces: {str(e)}")
        return {}


def _analyze_chain_interface(
    structure: Structure,
    chain1: str,
    chain2: str,
    cutoff: float = 5.0,
    min_contacts: int = 5,
) -> Optional[Dict[str, Any]]:
    """Analyze interface between two chains.

    Args:
        structure: BioPython Structure object
        chain1: First chain ID
        chain2: Second chain ID
        cutoff: Distance cutoff for contacts
        min_contacts: Minimum contacts for interface

    Returns:
        Interface description dictionary or None if no significant interface
    """
    try:
        # Get chain atoms
        atoms1 = [atom for atom in structure[0][chain1].get_atoms()]
        atoms2 = [atom for atom in structure[0][chain2].get_atoms()]

        # Find contacts
        contacts = []
        for atom1 in atoms1:
            res1 = atom1.get_parent()
            for atom2 in atoms2:
                res2 = atom2.get_parent()
                diff = atom1.get_coord() - atom2.get_coord()
                dist = np.sqrt(np.sum(diff * diff))
                if dist <= cutoff:
                    contacts.append((res1.get_id()[1], res2.get_id()[1]))

        if len(contacts) < min_contacts:
            return None

        # Calculate interface area
        interface_area = _calculate_interface_area(
            structure[0][chain1],
            structure[0][chain2],
        )

        # Get interface residues
        interface_residues1 = sorted(list(set(c[0] for c in contacts)))
        interface_residues2 = sorted(list(set(c[1] for c in contacts)))

        return {
            "chain1": chain1,
            "chain2": chain2,
            "contacts": contacts,
            "area": float(interface_area),
            "residues1": interface_residues1,
            "residues2": interface_residues2,
            "n_contacts": len(contacts),
            "n_residues": len(interface_residues1) + len(interface_residues2),
        }

    except Exception as e:
        logger.error(f"Error analyzing chain interface: {str(e)}")
        return None


def _calculate_interface_area(chain1, chain2) -> float:
    """Calculate interface surface area between two chains.

    Args:
        chain1: First BioPython Chain object
        chain2: Second BioPython Chain object

    Returns:
        Interface area in Å²
    """
    try:
        # Calculate surface areas
        area_1 = calculate_surface_area([atom.get_coord() for atom in chain1.get_atoms()])
        area_2 = calculate_surface_area([atom.get_coord() for atom in chain2.get_atoms()])
        area_complex = calculate_surface_area([atom.get_coord() for atom in chain1.get_atoms()] + [atom.get_coord() for atom in chain2.get_atoms()])

        # Interface area is the sum of individual areas minus complex area
        interface_area = area_1 + area_2 - area_complex

        return float(interface_area)

    except Exception as e:
        logger.error(f"Error calculating interface area: {str(e)}")
        return 0.0


def _calculate_interface_statistics(
    interfaces: List[Dict[str, Any]],
) -> Dict[str, float]:
    """Calculate overall interface statistics.

    Args:
        interfaces: List of interface descriptions

    Returns:
        Dictionary of interface statistics
    """
    if not interfaces:
        return {}

    areas = [interface["area"] for interface in interfaces]
    n_contacts = [interface["n_contacts"] for interface in interfaces]
    n_residues = [interface["n_residues"] for interface in interfaces]

    return {
        "mean_area": float(np.mean(areas)),
        "total_area": float(np.sum(areas)),
        "min_area": float(np.min(areas)),
        "max_area": float(np.max(areas)),
        "mean_contacts": float(np.mean(n_contacts)),
        "mean_residues": float(np.mean(n_residues)),
    }
