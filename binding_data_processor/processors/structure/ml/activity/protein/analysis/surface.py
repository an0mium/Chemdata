"""Surface analysis functionality for protein structures."""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
from Bio.PDB import Structure, Model, Chain, Residue, Atom

from .utils import (
    get_atom_coords,
    calculate_sasa,
    normalize_vector,
)
from binding_data_processor.core.config import ProteinAnalysisConfig

logger = logging.getLogger(__name__)


def analyze_surface(structure: Structure, config: Optional[ProteinAnalysisConfig] = None) -> Dict[str, Any]:
    """Analyze protein surface properties.

    Args:
        structure: BioPython Structure object
        config: Optional configuration object

    Returns:
        Dictionary containing surface analysis results
    """
    if config is None:
        config = ProteinAnalysisConfig()

    try:
        # Get atomic coordinates and radii
        coords = []
        radii = []
        atoms = []
        for atom in structure.get_atoms():
            coords.append(atom.get_coord())
            radii.append(config.atom_radii.get(atom.element, 1.5))
            atoms.append(atom)

        coords = np.array(coords)
        radii = np.array(radii)

        # Calculate total SASA
        total_area = calculate_sasa(coords, radii)

        # Calculate component areas
        polar_area = calculate_polar_area(coords, radii, atoms)
        hydrophobic_area = calculate_hydrophobic_area(coords, radii, atoms)

        # Get exposed residues
        exposed_residues = identify_exposed_residues(structure, threshold=config.residue_exposure_threshold)

        # Calculate surface properties
        properties = {
            "total_area": total_area,
            "polar_area": polar_area,
            "hydrophobic_area": hydrophobic_area,
            "polar_ratio": polar_area / total_area if total_area > 0 else 0.0,
            "hydrophobic_ratio": hydrophobic_area / total_area if total_area > 0 else 0.0,
            "exposed_residues": exposed_residues,
            "charge": calculate_surface_charge(structure, exposed_residues, config),
            "hydrophobicity": calculate_surface_hydrophobicity(structure, exposed_residues, config),
        }

        return properties

    except Exception as e:
        logger.error(f"Error analyzing surface: {str(e)}")
        return {}


def calculate_polar_area(coords: np.ndarray, radii: np.ndarray, atoms: List[Atom]) -> float:
    """Calculate polar surface area.

    Args:
        coords: Array of shape (n_atoms, 3) containing xyz coordinates
        radii: Array of shape (n_atoms,) containing atomic radii
        atoms: List of BioPython Atom objects

    Returns:
        Polar surface area in square Angstroms
    """
    try:
        # Get polar atoms (N, O)
        polar_mask = np.array([atom.element in ["N", "O"] for atom in atoms])

        if not np.any(polar_mask):
            return 0.0

        polar_coords = coords[polar_mask]
        polar_radii = radii[polar_mask]

        return calculate_sasa(polar_coords, polar_radii)

    except Exception as e:
        logger.error(f"Error calculating polar area: {str(e)}")
        return 0.0


def calculate_hydrophobic_area(coords: np.ndarray, radii: np.ndarray, atoms: List[Atom]) -> float:
    """Calculate hydrophobic surface area.

    Args:
        coords: Array of shape (n_atoms, 3) containing xyz coordinates
        radii: Array of shape (n_atoms,) containing atomic radii
        atoms: List of BioPython Atom objects

    Returns:
        Hydrophobic surface area in square Angstroms
    """
    try:
        # Get carbon atoms not bonded to N/O
        hydrophobic_mask = []
        for atom in atoms:
            if atom.element == "C":
                is_hydrophobic = True
                # Check bonded atoms in residue
                for other in atom.get_parent():
                    if other.element in ["N", "O"]:
                        is_hydrophobic = False
                        break
                hydrophobic_mask.append(is_hydrophobic)
            else:
                hydrophobic_mask.append(False)

        hydrophobic_mask = np.array(hydrophobic_mask)
        if not np.any(hydrophobic_mask):
            return 0.0

        hydrophobic_coords = coords[hydrophobic_mask]
        hydrophobic_radii = radii[hydrophobic_mask]

        return calculate_sasa(hydrophobic_coords, hydrophobic_radii)

    except Exception as e:
        logger.error(f"Error calculating hydrophobic area: {str(e)}")
        return 0.0


def identify_exposed_residues(structure: Structure, threshold: float = 2.8) -> List[int]:
    """Identify solvent-exposed residues.

    Args:
        structure: BioPython Structure object
        threshold: Minimum SASA for exposed residue (Å²)

    Returns:
        List of exposed residue numbers
    """
    try:
        exposed = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    # Calculate SASA for residue
                    coords = []
                    radii = []
                    for atom in residue:
                        coords.append(atom.get_coord())
                        radii.append(1.4 + 1.5)  # Probe radius + vdW radius

                    if coords:
                        sasa = calculate_sasa(np.array(coords), np.array(radii))
                        if sasa > threshold:
                            exposed.append(residue.get_id()[1])

        return exposed

    except Exception as e:
        logger.error(f"Error identifying exposed residues: {str(e)}")
        return []


def calculate_surface_charge(structure: Structure, exposed_residues: List[int], config: ProteinAnalysisConfig) -> float:
    """Calculate net surface charge.

    Args:
        structure: BioPython Structure object
        exposed_residues: List of exposed residue numbers
        config: Configuration object with residue properties

    Returns:
        Net surface charge
    """
    try:
        charge = 0.0
        charge_scale = config.residue_properties["charge"]

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[1] in exposed_residues:
                        charge += charge_scale.get(residue.get_resname(), 0.0)

        return charge

    except Exception as e:
        logger.error(f"Error calculating surface charge: {str(e)}")
        return 0.0


def calculate_surface_hydrophobicity(structure: Structure, exposed_residues: List[int], config: ProteinAnalysisConfig) -> float:
    """Calculate average surface hydrophobicity.

    Args:
        structure: BioPython Structure object
        exposed_residues: List of exposed residue numbers
        config: Configuration object with residue properties

    Returns:
        Average hydrophobicity of exposed residues
    """
    try:
        hydrophobicity = []
        hydrophobicity_scale = config.residue_properties["hydrophobicity"]

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[1] in exposed_residues:
                        h = hydrophobicity_scale.get(residue.get_resname(), 0.0)
                        hydrophobicity.append(h)

        if hydrophobicity:
            return float(np.mean(hydrophobicity))
        return 0.0

    except Exception as e:
        logger.error(f"Error calculating surface hydrophobicity: {str(e)}")
        return 0.0
