"""Utility functions for protein structure analysis."""

import logging
from typing import Dict, List, Tuple, Optional, Union, Any
import numpy as np
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB.SASA import ShrakeRupley
from scipy.spatial import ConvexHull, Delaunay

logger = logging.getLogger(__name__)


def get_residue_property(residue: Residue, property_name: str) -> Any:
    """Get property value for a residue.

    Args:
        residue: BioPython Residue object
        property_name: Name of property to get. Options:
            - hydrophobicity: Hydrophobicity score
            - charge: Residue charge
            - volume: Residue volume in Å³
            - surface_area: Solvent accessible surface area in Å²
            - secondary_structure: Secondary structure type (H, E, C)
            - coordinates: Array of atomic coordinates
            - center: Center of mass coordinates
            - backbone_angles: (phi, psi) angles in degrees
            - packing_density: Local packing density
            - conservation: Conservation score if available
            - depth: Depth from protein surface in Å
            - contacts: List of contacting residue IDs
            - interface: Whether residue is in an interface

    Returns:
        Property value, type depends on property
    """
    try:
        if property_name == "hydrophobicity":
            return _get_hydrophobicity(residue)
        elif property_name == "charge":
            return _get_charge(residue)
        elif property_name == "volume":
            return _calculate_residue_volume(residue)
        elif property_name == "surface_area":
            return _calculate_residue_surface_area(residue)
        elif property_name == "secondary_structure":
            return _get_secondary_structure(residue)
        elif property_name == "coordinates":
            return _get_residue_coordinates(residue)
        elif property_name == "center":
            return _calculate_residue_center(residue)
        elif property_name == "backbone_angles":
            return _calculate_backbone_angles(residue)
        elif property_name == "packing_density":
            return _calculate_packing_density(residue)
        elif property_name == "conservation":
            return _get_conservation_score(residue)
        elif property_name == "depth":
            # Get surface atoms first
            structure = residue.get_parent().get_parent().get_parent()
            surface_atoms = get_surface_atoms(structure)
            if not surface_atoms:
                return 0.0
            surface_points = np.array([atom.get_coord() for atom in surface_atoms])

            # Calculate depth directly here to avoid circular imports
            try:
                # Calculate residue center
                residue_coords = []
                for atom in residue:
                    residue_coords.append(atom.get_coord())
                if not residue_coords:
                    return 0.0

                residue_center = np.mean(residue_coords, axis=0)

                # Calculate minimum distance to any surface point
                min_distance = float("inf")
                for point in surface_points:
                    dist = np.linalg.norm(point - residue_center)
                    min_distance = min(min_distance, dist)

                return min_distance

            except Exception as e:
                logger.error(f"Error calculating residue depth: {str(e)}")
                return 0.0

        elif property_name == "contacts":
            return _get_residue_contacts(residue)
        elif property_name == "interface":
            return _is_interface_residue(residue)
        else:
            logger.warning(f"Unknown residue property: {property_name}")
            return None
    except Exception as e:
        logger.error(f"Error getting residue property {property_name}: {str(e)}")
        return None


def get_ca_atoms(structure: Structure) -> List[Atom]:
    """Get alpha carbon atoms from protein structure.

    Args:
        structure: BioPython Structure object

    Returns:
        List of CA atoms
    """
    try:
        ca_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if "CA" in residue:
                        ca_atoms.append(residue["CA"])
        return ca_atoms
    except Exception as e:
        logger.error(f"Error getting CA atoms: {str(e)}")
        return []


def get_atom_radius(atom: Atom) -> float:
    """Get van der Waals radius for an atom.

    Args:
        atom: BioPython Atom object

    Returns:
        Van der Waals radius in Angstroms
    """
    try:
        # Standard vdW radii in Angstroms
        radii = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8, "H": 1.2, "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98}
        return radii.get(atom.element, 1.5)  # Default radius
    except Exception as e:
        logger.error(f"Error getting atom radius: {str(e)}")
        return 1.5


def get_atom_coords(structure: Structure) -> np.ndarray:
    """Get coordinates of all atoms in structure.

    Args:
        structure: BioPython Structure object

    Returns:
        Array of shape (n_atoms, 3) containing xyz coordinates
    """
    try:
        coords = []
        for atom in structure.get_atoms():
            coords.append(atom.get_coord())
        return np.array(coords)
    except Exception as e:
        logger.error(f"Error getting atom coordinates: {str(e)}")
        return np.array([])


def get_backbone_atoms(residue: Residue) -> Dict[str, Atom]:
    """Get backbone atoms (N, CA, C, O) for a residue.

    Args:
        residue: BioPython Residue object

    Returns:
        Dictionary mapping atom names to Atom objects
    """
    try:
        backbone = {}
        for atom in residue:
            if atom.get_name() in ["N", "CA", "C", "O"]:
                backbone[atom.get_name()] = atom
        return backbone
    except Exception as e:
        logger.error(f"Error getting backbone atoms: {str(e)}")
        return {}


def calculate_center_of_mass(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> np.ndarray:
    """Calculate center of mass for a set of coordinates.

    Args:
        coords: Array of shape (n_atoms, 3) containing xyz coordinates
        masses: Optional array of shape (n_atoms,) containing atomic masses

    Returns:
        Array of shape (3,) containing center of mass coordinates
    """
    try:
        if masses is None:
            masses = np.ones(len(coords))
        return np.average(coords, weights=masses, axis=0)
    except Exception as e:
        logger.error(f"Error calculating center of mass: {str(e)}")
        return np.zeros(3)


def calculate_radius_of_gyration(coords: np.ndarray, masses: Optional[np.ndarray] = None) -> float:
    """Calculate radius of gyration.

    Args:
        coords: Array of shape (n_atoms, 3) containing xyz coordinates
        masses: Optional array of shape (n_atoms,) containing atomic masses

    Returns:
        Radius of gyration in Angstroms
    """
    try:
        if masses is None:
            masses = np.ones(len(coords))

        com = calculate_center_of_mass(coords, masses)
        r2 = np.sum(masses * np.sum((coords - com) ** 2, axis=1))
        total_mass = np.sum(masses)
        return np.sqrt(r2 / total_mass)
    except Exception as e:
        logger.error(f"Error calculating radius of gyration: {str(e)}")
        return 0.0


def calculate_sasa(coords: np.ndarray, radii: np.ndarray, probe_radius: float = 1.4) -> float:
    """Calculate solvent accessible surface area.

    Args:
        coords: Array of atomic coordinates (N x 3)
        radii: Array of atomic radii (N)
        probe_radius: Probe radius (water = 1.4 Å)

    Returns:
        SASA in Å²
    """
    try:
        sr = ShrakeRupley()
        return sr.compute(coords, radii + probe_radius)
    except Exception as e:
        logger.error(f"Error calculating SASA: {str(e)}")
        return 0.0


def calculate_distance_matrix(coords: np.ndarray) -> np.ndarray:
    """Calculate pairwise distance matrix between coordinates.

    Args:
        coords: Array of shape (n_points, 3) containing xyz coordinates

    Returns:
        Array of shape (n_points, n_points) containing pairwise distances
    """
    try:
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        return np.sqrt(np.sum(diff * diff, axis=-1))
    except Exception as e:
        logger.error(f"Error calculating distance matrix: {str(e)}")
        return np.array([])


def find_contacts(coords1: np.ndarray, coords2: np.ndarray, cutoff: float = 5.0) -> List[Tuple[int, int]]:
    """Find contacts between two sets of coordinates within cutoff distance.

    Args:
        coords1: Array of shape (n1, 3) containing xyz coordinates
        coords2: Array of shape (n2, 3) containing xyz coordinates
        cutoff: Distance cutoff in Angstroms

    Returns:
        List of tuples containing indices of contacting points
    """
    try:
        diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
        dist = np.sqrt(np.sum(diff * diff, axis=-1))
        contacts = np.where(dist <= cutoff)
        return list(zip(contacts[0], contacts[1]))
    except Exception as e:
        logger.error(f"Error finding contacts: {str(e)}")
        return []


def find_cavities(surface_points: np.ndarray, min_radius: float = 2.0) -> List[Dict[str, Any]]:
    """Find cavities in molecular surface using alpha shapes.

    Args:
        surface_points: Array of surface point coordinates (N x 3)
        min_radius: Minimum cavity radius in Angstroms

    Returns:
        List of cavity properties (center, radius, volume, etc)
    """
    try:
        cavities = []
        tri = Delaunay(surface_points)

        for simplex in tri.simplices:
            center = surface_points[simplex].mean(axis=0)
            radius = np.linalg.norm(surface_points[simplex[0]] - center)

            if radius > min_radius:
                cavity = {
                    "center": center,
                    "radius": radius,
                    "volume": (4 / 3) * np.pi * radius**3,
                    "surface_points": surface_points[simplex],
                    "depth": calculate_cavity_depth(center, surface_points),
                }
                cavities.append(cavity)

        return cavities
    except Exception as e:
        logger.error(f"Error finding cavities: {str(e)}")
        return []


def calculate_cavity_depth(center: np.ndarray, surface_points: np.ndarray) -> float:
    """Calculate depth of cavity from surface.

    Args:
        center: Cavity center coordinates
        surface_points: Surface point coordinates

    Returns:
        Cavity depth in Angstroms
    """
    try:
        distances = np.linalg.norm(surface_points - center, axis=1)
        return float(np.min(distances))
    except Exception as e:
        logger.error(f"Error calculating cavity depth: {str(e)}")
        return 0.0


def calculate_volume(coords: np.ndarray) -> float:
    """Calculate volume of point cloud using convex hull.

    Args:
        coords: Array of point coordinates (N x 3)

    Returns:
        Volume in Å³
    """
    try:
        hull = ConvexHull(coords)
        return hull.volume
    except Exception as e:
        logger.error(f"Error calculating volume: {str(e)}")
        return 0.0


def get_sequence_from_structure(structure: Structure) -> str:
    """Extract amino acid sequence from structure.

    Args:
        structure: BioPython Structure object

    Returns:
        Amino acid sequence as string
    """
    try:
        sequence = ""
        for model in structure:
            for chain in model:
                for residue in chain:
                    if "CA" in residue:  # Only standard amino acids
                        sequence += residue.get_resname()
        return sequence
    except Exception as e:
        logger.error(f"Error getting sequence: {str(e)}")
        return ""


def get_secondary_structure(structure: Structure) -> Dict[str, float]:
    """Get secondary structure composition.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary mapping SS types to fractions
    """
    try:
        dssp = dssp_dict_from_pdb_file(structure.id)[0]
        ss_counts = {"H": 0, "B": 0, "E": 0, "G": 0, "I": 0, "T": 0, "S": 0}

        for residue in dssp:
            ss = dssp[residue][2]
            if ss in ss_counts:
                ss_counts[ss] += 1

        total = sum(ss_counts.values())
        if total > 0:
            return {k: v / total for k, v in ss_counts.items()}
        return ss_counts
    except Exception as e:
        logger.error(f"Error getting secondary structure: {str(e)}")
        return {}


def get_surface_atoms(structure: Structure, probe_radius: float = 1.4) -> List[Atom]:
    """Get surface-exposed atoms using SASA calculation.

    Args:
        structure: BioPython Structure object
        probe_radius: Probe radius in Angstroms

    Returns:
        List of surface-exposed atoms
    """
    try:
        coords = []
        radii = []
        atoms = []
        for atom in structure.get_atoms():
            coords.append(atom.get_coord())
            radii.append(get_atom_radius(atom))
            atoms.append(atom)

        coords = np.array(coords)
        radii = np.array(radii)

        sr = ShrakeRupley()
        sasa = sr.compute_atoms(coords, radii + probe_radius)

        return [atom for atom, area in zip(atoms, sasa) if area > 0.0]
    except Exception as e:
        logger.error(f"Error getting surface atoms: {str(e)}")
        return []


def get_interface_residues(chain1: Chain, chain2: Chain, cutoff: float = 5.0) -> List[Tuple[Residue, Residue]]:
    """Get pairs of residues in contact between two chains.

    Args:
        chain1: First BioPython Chain object
        chain2: Second BioPython Chain object
        cutoff: Distance cutoff in Angstroms

    Returns:
        List of residue pairs in contact
    """
    try:
        contacts = []
        for res1 in chain1:
            for res2 in chain2:
                for atom1 in res1:
                    for atom2 in res2:
                        if np.linalg.norm(atom1.get_coord() - atom2.get_coord()) <= cutoff:
                            contacts.append((res1, res2))
                            break
                    else:
                        continue
                    break
        return contacts
    except Exception as e:
        logger.error(f"Error getting interface residues: {str(e)}")
        return []


def calculate_torsion_angle(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """Calculate torsion angle between 4 points.

    Args:
        p1, p2, p3, p4: Arrays of shape (3,) containing xyz coordinates

    Returns:
        Torsion angle in degrees
    """
    try:
        v1 = p2 - p1
        v2 = p3 - p2
        v3 = p4 - p3

        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)

        n1 = n1 / np.linalg.norm(n1)
        n2 = n2 / np.linalg.norm(n2)

        x = np.dot(n1, n2)
        y = np.dot(np.cross(n1, v2 / np.linalg.norm(v2)), n2)
        angle = np.arctan2(y, x)

        return np.degrees(angle)
    except Exception as e:
        logger.error(f"Error calculating torsion angle: {str(e)}")
        return 0.0


def normalize_vector(v: np.ndarray) -> np.ndarray:
    """Normalize vector to unit length.

    Args:
        v: Array of shape (n,) to normalize

    Returns:
        Normalized vector
    """
    try:
        norm = np.linalg.norm(v)
        if norm == 0:
            return v
        return v / norm
    except Exception as e:
        logger.error(f"Error normalizing vector: {str(e)}")
        return v


def calculate_interface_area(chain1: Chain, chain2: Chain, probe_radius: float = 1.4) -> float:
    """Calculate interface surface area between two chains.

    Args:
        chain1: First BioPython Chain object
        chain2: Second BioPython Chain object
        probe_radius: Probe radius in Angstroms (default 1.4 for water)

    Returns:
        Interface area in Å²
    """
    try:
        # Calculate SASA for each chain separately and combined
        sr = ShrakeRupley()

        # Get coordinates and radii for chain1
        coords1 = []
        radii1 = []
        for atom in chain1.get_atoms():
            coords1.append(atom.get_coord())
            radii1.append(get_atom_radius(atom))
        coords1 = np.array(coords1)
        radii1 = np.array(radii1)

        # Get coordinates and radii for chain2
        coords2 = []
        radii2 = []
        for atom in chain2.get_atoms():
            coords2.append(atom.get_coord())
            radii2.append(get_atom_radius(atom))
        coords2 = np.array(coords2)
        radii2 = np.array(radii2)

        # Calculate SASA for each chain separately
        sasa1 = sr.compute(coords1, radii1 + probe_radius)
        sasa2 = sr.compute(coords2, radii2 + probe_radius)

        # Calculate SASA for combined chains
        coords_combined = np.vstack([coords1, coords2])
        radii_combined = np.concatenate([radii1, radii2])
        sasa_combined = sr.compute(coords_combined, radii_combined + probe_radius)

        # Interface area is the sum of individual SASAs minus combined SASA
        interface_area = float(np.sum(sasa1) + np.sum(sasa2) - np.sum(sasa_combined))
        return max(0.0, interface_area)  # Ensure non-negative

    except Exception as e:
        logger.error(f"Error calculating interface area: {str(e)}")
        return 0.0


def calculate_cavity_volume(surface_points: np.ndarray, probe_radius: float = 1.4) -> float:
    """Calculate volume of a cavity defined by surface points.

    Args:
        surface_points: Array of surface point coordinates (N x 3)
        probe_radius: Probe radius in Angstroms (default 1.4 for water)

    Returns:
        Cavity volume in Å³
    """
    try:
        if len(surface_points) < 4:
            return 0.0

        # Calculate convex hull
        hull = ConvexHull(surface_points)
        total_volume = hull.volume

        # Calculate the void volume by subtracting the volume of surface atoms
        n_points = len(surface_points)
        point_volume = (4 / 3) * np.pi * probe_radius**3
        void_volume = total_volume - (n_points * point_volume)

        return max(0.0, void_volume)  # Ensure non-negative

    except Exception as e:
        logger.error(f"Error calculating cavity volume: {str(e)}")
        return 0.0


def calculate_surface_exposure(residue: Residue, probe_radius: float = 1.4) -> float:
    """Calculate solvent accessible surface area for a residue.

    Args:
        residue: BioPython Residue object
        probe_radius: Probe radius in Angstroms (default 1.4 for water)

    Returns:
        SASA in Å²
    """
    try:
        # Get coordinates and radii for residue atoms
        coords = []
        radii = []
        for atom in residue:
            coords.append(atom.get_coord())
            radii.append(get_atom_radius(atom))

        coords = np.array(coords)
        radii = np.array(radii)

        if len(coords) < 1:
            return 0.0

        # Calculate SASA using Shrake-Rupley algorithm
        sr = ShrakeRupley()
        sasa = sr.compute(coords, radii + probe_radius)
        return float(np.sum(sasa))

    except Exception as e:
        logger.error(f"Error calculating surface exposure: {str(e)}")
        return 0.0


# Helper functions for get_residue_property
def _get_hydrophobicity(residue: Residue) -> float:
    """Get hydrophobicity value for residue."""
    hydrophobicity = {
        "ALA": 1.8,
        "ARG": -4.5,
        "ASN": -3.5,
        "ASP": -3.5,
        "CYS": 2.5,
        "GLN": -3.5,
        "GLU": -3.5,
        "GLY": -0.4,
        "HIS": -3.2,
        "ILE": 4.5,
        "LEU": 3.8,
        "LYS": -3.9,
        "MET": 1.9,
        "PHE": 2.8,
        "PRO": -1.6,
        "SER": -0.8,
        "THR": -0.7,
        "TRP": -0.9,
        "TYR": -1.3,
        "VAL": 4.2,
    }
    return hydrophobicity.get(residue.get_resname(), 0.0)


def _get_charge(residue: Residue) -> float:
    """Get charge value for residue."""
    charge = {"ARG": 1.0, "LYS": 1.0, "ASP": -1.0, "GLU": -1.0, "HIS": 0.1}
    return charge.get(residue.get_resname(), 0.0)


def _calculate_residue_volume(residue: Residue) -> float:
    """Calculate volume of residue using convex hull."""
    try:
        coords = [atom.get_coord() for atom in residue]
        if len(coords) < 4:
            return 0.0
        hull = ConvexHull(coords)
        return float(hull.volume)
    except Exception as e:
        logger.error(f"Error calculating residue volume: {str(e)}")
        return 0.0


def _calculate_residue_surface_area(residue: Residue) -> float:
    """Calculate surface area of residue."""
    try:
        coords = [atom.get_coord() for atom in residue]
        if len(coords) < 4:
            return 0.0
        hull = ConvexHull(coords)
        return float(hull.area)
    except Exception as e:
        logger.error(f"Error calculating residue surface area: {str(e)}")
        return 0.0


def _get_secondary_structure(residue: Residue) -> str:
    """Get secondary structure type for residue."""
    try:
        structure = residue.get_parent().get_parent().get_parent()
        dssp = dssp_dict_from_pdb_file(structure.id)[0]
        key = (residue.get_parent().id, residue.id[1])
        if key in dssp:
            return dssp[key][2]
        return "C"  # Coil as default
    except Exception as e:
        logger.error(f"Error getting secondary structure: {str(e)}")
        return "C"


def _get_residue_coordinates(residue: Residue) -> np.ndarray:
    """Get coordinates of all atoms in residue."""
    try:
        return np.array([atom.get_coord() for atom in residue])
    except Exception as e:
        logger.error(f"Error getting residue coordinates: {str(e)}")
        return np.array([])


def _calculate_residue_center(residue: Residue) -> np.ndarray:
    """Calculate center of mass of residue."""
    try:
        coords = [atom.get_coord() for atom in residue]
        return np.mean(coords, axis=0)
    except Exception as e:
        logger.error(f"Error calculating residue center: {str(e)}")
        return np.zeros(3)


def _calculate_backbone_angles(residue: Residue) -> Tuple[Optional[float], Optional[float]]:
    """Calculate phi/psi angles for residue."""
    try:
        prev_res = residue.get_previous_residue()
        next_res = residue.get_next_residue()

        phi = None
        if prev_res and "C" in prev_res and all(atom in residue for atom in ["N", "CA", "C"]):
            phi = calculate_torsion_angle(prev_res["C"].get_coord(), residue["N"].get_coord(), residue["CA"].get_coord(), residue["C"].get_coord())

        psi = None
        if next_res and "N" in next_res and all(atom in residue for atom in ["N", "CA", "C"]):
            psi = calculate_torsion_angle(residue["N"].get_coord(), residue["CA"].get_coord(), residue["C"].get_coord(), next_res["N"].get_coord())

        return phi, psi
    except Exception as e:
        logger.error(f"Error calculating backbone angles: {str(e)}")
        return None, None


def _calculate_packing_density(residue: Residue) -> float:
    """Calculate local packing density around residue."""
    try:
        structure = residue.get_parent().get_parent().get_parent()
        center = _calculate_residue_center(residue)
        radius = 10.0  # Å

        n_neighbors = 0
        for other_res in structure.get_residues():
            if other_res != residue:
                other_center = _calculate_residue_center(other_res)
                dist = np.linalg.norm(center - other_center)
                if dist <= radius:
                    n_neighbors += 1

        volume = (4 / 3) * np.pi * radius**3
        return float(n_neighbors / volume)
    except Exception as e:
        logger.error(f"Error calculating packing density: {str(e)}")
        return 0.0


def _get_conservation_score(residue: Residue) -> Optional[float]:
    """Get conservation score for residue if available."""
    try:
        # This would need to be implemented based on available conservation data
        return None
    except Exception as e:
        logger.error(f"Error getting conservation score: {str(e)}")
        return None


def _get_residue_contacts(residue: Residue) -> List[int]:
    """Get list of residues in contact with this residue."""
    try:
        structure = residue.get_parent().get_parent().get_parent()
        cutoff = 5.0  # Å
        contacts = []

        center = _calculate_residue_center(residue)
        for other_res in structure.get_residues():
            if other_res != residue:
                other_center = _calculate_residue_center(other_res)
                if np.linalg.norm(center - other_center) <= cutoff:
                    contacts.append(other_res.id[1])

        return contacts
    except Exception as e:
        logger.error(f"Error getting residue contacts: {str(e)}")
        return []


def _is_interface_residue(residue: Residue) -> bool:
    """Check if residue is in an interface between chains."""
    try:
        chain = residue.get_parent()
        structure = chain.get_parent().get_parent()
        cutoff = 5.0  # Å

        for other_chain in structure[0]:
            if other_chain != chain:
                contacts = get_interface_residues(chain, other_chain, cutoff)
                if any(res1 == residue or res2 == residue for res1, res2 in contacts):
                    return True

        return False
    except Exception as e:
        logger.error(f"Error checking interface residue: {str(e)}")
        return False
