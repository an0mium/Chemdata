"""Utility functions for protein structure analysis."""

import logging
import numpy as np
from Bio.PDB import Atom, Residue, Structure, Chain
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from scipy.spatial import ConvexHull
from typing import Dict, List, Optional, Tuple, Union

logger = logging.getLogger(__name__)

# Standard vdW radii in Angstroms
VDW_RADII = {
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
    "P": 1.8,
    "H": 1.2,
    "F": 1.47,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
}

# Property scales
RESIDUE_PROPERTIES = {
    "hydrophobicity": {  # Kyte-Doolittle
        "ILE": 4.5,
        "VAL": 4.2,
        "LEU": 3.8,
        "PHE": 2.8,
        "CYS": 2.5,
        "MET": 1.9,
        "ALA": 1.8,
        "GLY": -0.4,
        "THR": -0.7,
        "SER": -0.8,
        "TRP": -0.9,
        "TYR": -1.3,
        "PRO": -1.6,
        "HIS": -3.2,
        "GLU": -3.5,
        "GLN": -3.5,
        "ASP": -3.5,
        "ASN": -3.5,
        "LYS": -3.9,
        "ARG": -4.5,
    },
    "charge": {  # Net charge at pH 7
        "ARG": 1,
        "LYS": 1,
        "ASP": -1,
        "GLU": -1,
        "HIS": 0.5,  # Can be protonated
    },
    "volume": {  # Residue volumes in Å³
        "GLY": 60.1,
        "ALA": 88.6,
        "SER": 89.0,
        "CYS": 108.5,
        "ASP": 111.1,
        "PRO": 112.7,
        "ASN": 114.1,
        "THR": 116.1,
        "GLU": 138.4,
        "VAL": 140.0,
        "GLN": 143.8,
        "HIS": 153.2,
        "MET": 162.9,
        "ILE": 166.7,
        "LEU": 166.7,
        "LYS": 168.6,
        "ARG": 173.4,
        "PHE": 189.9,
        "TYR": 193.6,
        "TRP": 227.8,
    },
    "reference_sasa": {  # Reference SASA values for fully exposed residues
        "ALA": 113.0,
        "ARG": 241.0,
        "ASN": 158.0,
        "ASP": 151.0,
        "CYS": 140.0,
        "GLN": 189.0,
        "GLU": 183.0,
        "GLY": 85.0,
        "HIS": 194.0,
        "ILE": 182.0,
        "LEU": 180.0,
        "LYS": 211.0,
        "MET": 204.0,
        "PHE": 218.0,
        "PRO": 143.0,
        "SER": 122.0,
        "THR": 146.0,
        "TRP": 259.0,
        "TYR": 229.0,
        "VAL": 160.0,
    },
}


def get_atom_radius(atom: Atom.Atom) -> float:
    """Get van der Waals radius for atom."""
    return VDW_RADII.get(atom.element, 1.5)  # Default radius


def calculate_center_of_mass(atoms: List[Atom.Atom], weighted: bool = True) -> np.ndarray:
    """Calculate center of mass for a set of atoms.

    Args:
        atoms: List of BioPython Atom objects
        weighted: Whether to weight by atomic mass

    Returns:
        Center of mass coordinates as numpy array
    """
    try:
        if not atoms:
            return np.zeros(3)

        coords = []
        weights = []

        for atom in atoms:
            coords.append(atom.get_coord())
            weights.append(atom.mass if weighted and hasattr(atom, "mass") else 1.0)

        coords = np.array(coords)
        weights = np.array(weights)

        com = np.average(coords, weights=weights, axis=0)
        return com

    except Exception as e:
        logger.error(f"Error calculating center of mass: {str(e)}")
        return np.zeros(3)


def calculate_radius_of_gyration(atoms: List[Atom.Atom]) -> float:
    """Calculate radius of gyration for a set of atoms."""
    try:
        if not atoms:
            return 0.0

        coords = []
        masses = []

        for atom in atoms:
            coords.append(atom.get_coord())
            masses.append(atom.mass if hasattr(atom, "mass") else 1.0)

        coords = np.array(coords)
        masses = np.array(masses)

        # Calculate center of mass
        com = calculate_center_of_mass(atoms)

        # Calculate radius of gyration
        rg2 = np.sum(masses * np.sum((coords - com) ** 2, axis=1))
        total_mass = np.sum(masses)

        return np.sqrt(rg2 / total_mass)

    except Exception as e:
        logger.error(f"Error calculating radius of gyration: {str(e)}")
        return 0.0


def get_residue_property(residue: Residue.Residue, property_name: str, default: float = 0.0) -> float:
    """Get residue property from property scales."""
    if property_name not in RESIDUE_PROPERTIES:
        logger.warning(f"Unknown property scale: {property_name}")
        return default
    return RESIDUE_PROPERTIES[property_name].get(residue.get_resname(), default)


def calculate_torsion(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """Calculate torsion angle between 4 points in degrees."""
    try:
        v1 = p2 - p1
        v2 = p3 - p2
        v3 = p4 - p3

        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)

        n1_norm = normalize_vector(n1)
        n2_norm = normalize_vector(n2)

        x = np.dot(n1_norm, n2_norm)
        y = np.dot(np.cross(n1_norm, normalize_vector(v2)), n2_norm)

        return np.degrees(np.arctan2(y, x))

    except Exception as e:
        logger.error(f"Error calculating torsion angle: {str(e)}")
        return 0.0


def get_residue_neighbors(residue: Residue.Residue, structure: Structure.Structure, cutoff: float = 10.0) -> List[Residue.Residue]:
    """Get neighboring residues within cutoff distance."""
    try:
        neighbors = []
        if not residue.has_id("CA"):
            return neighbors

        res_ca = residue["CA"].get_coord()
        for model in structure:
            for chain in model:
                for other_res in chain:
                    if other_res != residue and other_res.has_id("CA"):
                        other_ca = other_res["CA"].get_coord()
                        dist = np.linalg.norm(res_ca - other_ca)
                        if dist <= cutoff:
                            neighbors.append(other_res)
        return neighbors
    except Exception as e:
        logger.error(f"Error getting residue neighbors: {str(e)}")
        return []


def calculate_surface_area(coords: np.ndarray, radii: np.ndarray, probe_radius: float = 1.4) -> float:
    """Calculate solvent accessible surface area using Shrake-Rupley algorithm."""
    try:
        from Bio.PDB.SASA import ShrakeRupley

        sr = ShrakeRupley()
        return sr.compute(coords, radii + probe_radius)
    except Exception as e:
        logger.error(f"Error calculating surface area: {str(e)}")
        return 0.0


def calculate_surface_exposure(residue: Residue.Residue, structure: Structure.Structure) -> float:
    """Calculate relative surface exposure of residue."""
    try:
        # Calculate SASA for residue
        coords = []
        radii = []
        for atom in residue:
            coords.append(atom.get_coord())
            radii.append(get_atom_radius(atom))

        if not coords:
            return 0.0

        residue_sasa = calculate_surface_area(np.array(coords), np.array(radii))

        # Get reference SASA for residue type
        ref_sasa = get_residue_property(residue, "reference_sasa")

        return residue_sasa / ref_sasa if ref_sasa > 0 else 0.0

    except Exception as e:
        logger.error(f"Error calculating surface exposure: {str(e)}")
        return 0.0


def calculate_residue_depth(residue: Residue.Residue, surface_points: np.ndarray) -> float:
    """Calculate residue depth from protein surface."""
    try:
        if not residue.has_id("CA"):
            return 0.0

        ca_coord = residue["CA"].get_coord()
        distances = np.linalg.norm(surface_points - ca_coord, axis=1)
        return float(np.min(distances))

    except Exception as e:
        logger.error(f"Error calculating residue depth: {str(e)}")
        return 0.0


def get_secondary_structure(structure: Structure.Structure) -> Dict[int, str]:
    """Get secondary structure assignments using DSSP."""
    try:
        dssp_dict = dssp_dict_from_pdb_file(structure.id)[0]
        ss_map = {}

        for residue in structure.get_residues():
            res_id = (residue.get_parent().id, residue.id[1])
            if res_id in dssp_dict:
                ss_map[residue.id[1]] = dssp_dict[res_id][2]
            else:
                ss_map[residue.id[1]] = "-"

        return ss_map

    except Exception as e:
        logger.error(f"Error getting secondary structure: {str(e)}")
        return {}


def calculate_interface_area(chain1: Chain.Chain, chain2: Chain.Chain) -> float:
    """Calculate interface surface area between chains."""
    try:
        # Calculate SASA for individual chains
        coords1 = [atom.get_coord() for atom in chain1.get_atoms()]
        radii1 = [get_atom_radius(atom) for atom in chain1.get_atoms()]
        sasa1 = calculate_surface_area(np.array(coords1), np.array(radii1))

        coords2 = [atom.get_coord() for atom in chain2.get_atoms()]
        radii2 = [get_atom_radius(atom) for atom in chain2.get_atoms()]
        sasa2 = calculate_surface_area(np.array(coords2), np.array(radii2))

        # Calculate SASA for complex
        coords = coords1 + coords2
        radii = radii1 + radii2
        complex_sasa = calculate_surface_area(np.array(coords), np.array(radii))

        # Interface area is the difference
        return (sasa1 + sasa2 - complex_sasa) / 2

    except Exception as e:
        logger.error(f"Error calculating interface area: {str(e)}")
        return 0.0


def calculate_cavity_volume(points: np.ndarray) -> float:
    """Calculate volume of cavity from surface points."""
    try:
        hull = ConvexHull(points)
        return hull.volume
    except Exception as e:
        logger.error(f"Error calculating cavity volume: {str(e)}")
        return 0.0


def normalize_vector(v: np.ndarray) -> np.ndarray:
    """Normalize vector to unit length."""
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def get_atom_coords(structure: Structure.Structure) -> Tuple[np.ndarray, np.ndarray]:
    """Get atom coordinates and masses from structure."""
    try:
        coords = []
        masses = []

        for atom in structure.get_atoms():
            coords.append(atom.get_coord())
            masses.append(atom.mass if hasattr(atom, "mass") else 1.0)

        return np.array(coords), np.array(masses)

    except Exception as e:
        logger.error(f"Error getting atom coordinates: {str(e)}")
        return np.array([]), np.array([])


def get_backbone_atoms(residue: Residue.Residue) -> Tuple[Optional[Atom.Atom], Optional[Atom.Atom], Optional[Atom.Atom]]:
    """Get backbone N, CA, C atoms from residue."""
    n = residue["N"] if "N" in residue else None
    ca = residue["CA"] if "CA" in residue else None
    c = residue["C"] if "C" in residue else None
    return n, ca, c


def get_sequence_from_structure(structure: Structure.Structure, chain_id: Optional[str] = None) -> str:
    """Extract amino acid sequence from structure."""
    sequence = []
    for model in structure:
        for chain in model:
            if chain_id is None or chain.id == chain_id:
                for residue in chain:
                    if residue.get_resname() in RESIDUE_PROPERTIES["hydrophobicity"]:
                        sequence.append(residue.get_resname())
    return "".join(sequence)


def calculate_distance_matrix(coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
    """Calculate pairwise distance matrix between two sets of coordinates."""
    try:
        diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
        return np.sqrt(np.sum(diff * diff, axis=2))
    except Exception as e:
        logger.error(f"Error calculating distance matrix: {str(e)}")
        return np.array([[]])


def find_contacts(atoms1: List[Atom.Atom], atoms2: List[Atom.Atom], cutoff: float = 5.0) -> List[Tuple[Atom.Atom, Atom.Atom, float]]:
    """Find contacting atom pairs within cutoff distance."""
    try:
        contacts = []
        coords1 = np.array([a.get_coord() for a in atoms1])
        coords2 = np.array([a.get_coord() for a in atoms2])

        dist_matrix = calculate_distance_matrix(coords1, coords2)
        contact_indices = np.where(dist_matrix <= cutoff)

        for i, j in zip(*contact_indices):
            contacts.append((atoms1[i], atoms2[j], dist_matrix[i, j]))

        return contacts

    except Exception as e:
        logger.error(f"Error finding contacts: {str(e)}")
        return []
