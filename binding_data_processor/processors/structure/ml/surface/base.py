"""Base surface analysis functionality for protein structures."""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Any, Set
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from scipy.spatial import ConvexHull, Voronoi
from ...activity.protein.analysis.utils import calculate_sasa

logger = logging.getLogger(__name__)


class BaseSurfaceAnalyzer:
    """Base class for analyzing molecular surfaces."""

    # Standard van der Waals radii in Angstroms
    ATOM_RADII = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8, "H": 1.2, "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98}

    # Kyte-Doolittle hydrophobicity scale
    HYDROPHOBICITY = {
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
    }

    # Residue charges at physiological pH
    CHARGE = {
        "ARG": 1,  # Positive
        "LYS": 1,  # Positive
        "ASP": -1,  # Negative
        "GLU": -1,  # Negative
        "HIS": 0.1,  # Slightly positive
    }

    def __init__(self):
        """Initialize surface analyzer."""
        self.logger = logging.getLogger(self.__class__.__name__)

    def get_surface_atoms(
        self,
        structure: Structure,
        probe_radius: float = 1.4,
    ) -> Set[Atom]:
        """Get surface atoms using rolling ball algorithm.

        Args:
            structure: BioPython Structure object
            probe_radius: Probe radius in Å

        Returns:
            Set of surface atoms
        """
        try:
            surface_atoms = set()
            all_atoms = list(structure.get_atoms())
            coords = np.array([atom.get_coord() for atom in all_atoms])

            # For each atom, check if it's accessible to probe
            for i, atom in enumerate(all_atoms):
                is_surface = False
                atom_coord = coords[i]

                # Generate points on sphere around atom
                theta = np.linspace(0, np.pi, 10)
                phi = np.linspace(0, 2 * np.pi, 20)
                probe_points = []
                for t in theta:
                    for p in phi:
                        x = probe_radius * np.sin(t) * np.cos(p)
                        y = probe_radius * np.sin(t) * np.sin(p)
                        z = probe_radius * np.cos(t)
                        probe_points.append(atom_coord + np.array([x, y, z]))

                # Check if any probe position is accessible
                probe_coords = np.array(probe_points)
                distances = np.linalg.norm(probe_coords[:, None] - coords, axis=2)
                if np.any(np.all(distances > probe_radius, axis=1)):
                    surface_atoms.add(atom)

            return surface_atoms

        except Exception as e:
            self.logger.error(f"Error getting surface atoms: {str(e)}")
            return set()

    def calculate_surface_area(self, atoms: Set[Atom]) -> float:
        """Calculate solvent accessible surface area.

        Args:
            atoms: Set of atoms

        Returns:
            Surface area in Å²
        """
        try:
            if not atoms:
                return 0.0

            coords = []
            radii = []
            for atom in atoms:
                coords.append(atom.get_coord())
                radii.append(self.ATOM_RADII.get(atom.element, 1.5))

            coords = np.array(coords)
            radii = np.array(radii)

            return calculate_sasa(coords, radii)

        except Exception as e:
            self.logger.error(f"Error calculating surface area: {str(e)}")
            return 0.0

    def calculate_volume(self, coords: np.ndarray) -> float:
        """Calculate volume of point cloud.

        Args:
            coords: Array of shape (n_atoms, 3) containing xyz coordinates

        Returns:
            Volume in Å³
        """
        try:
            if len(coords) < 4:
                return 0.0
            hull = ConvexHull(coords)
            return hull.volume
        except Exception as e:
            self.logger.error(f"Error calculating volume: {str(e)}")
            return 0.0

    def get_surface_residues(self, surface_atoms: Set[Atom]) -> List[Residue]:
        """Get residues containing surface atoms.

        Args:
            surface_atoms: Set of surface atoms

        Returns:
            List of surface residues
        """
        try:
            surface_residues = set()
            for atom in surface_atoms:
                surface_residues.add(atom.get_parent())
            return list(surface_residues)
        except Exception as e:
            self.logger.error(f"Error getting surface residues: {str(e)}")
            return []

    def calculate_hydrophobicity(self, residues: List[Residue]) -> float:
        """Calculate average hydrophobicity of residues.

        Args:
            residues: List of residues

        Returns:
            Average hydrophobicity score
        """
        try:
            scores = [self.HYDROPHOBICITY.get(res.get_resname(), 0.0) for res in residues]
            return float(np.mean(scores)) if scores else 0.0
        except Exception as e:
            self.logger.error(f"Error calculating hydrophobicity: {str(e)}")
            return 0.0

    def calculate_charge(self, residues: List[Residue]) -> float:
        """Calculate net charge of residues.

        Args:
            residues: List of residues

        Returns:
            Net charge
        """
        try:
            charges = [self.CHARGE.get(res.get_resname(), 0.0) for res in residues]
            return float(sum(charges))
        except Exception as e:
            self.logger.error(f"Error calculating charge: {str(e)}")
            return 0.0

    def find_surface_cavities(
        self,
        surface_atoms: Set[Atom],
        probe_radius: float = 1.4,
    ) -> List[Dict[str, Any]]:
        """Find cavities in molecular surface.

        Args:
            surface_atoms: Set of surface atoms
            probe_radius: Probe radius in Å

        Returns:
            List of cavity properties
        """
        try:
            if not surface_atoms:
                return []

            # Get surface coordinates
            coords = np.array([atom.get_coord() for atom in surface_atoms])

            # Calculate Voronoi diagram
            vor = Voronoi(coords)
            cavities = []

            # Check each Voronoi vertex
            for vertex in vor.vertices:
                # Calculate distances to surface atoms
                distances = np.linalg.norm(coords - vertex, axis=1)
                min_dist = np.min(distances)

                # Check if vertex represents a cavity
                if min_dist > probe_radius:
                    cavity = {
                        "center": vertex,
                        "radius": float(min_dist),
                        "volume": float(4 / 3 * np.pi * min_dist**3),
                        "surrounding_atoms": [atom for i, atom in enumerate(surface_atoms) if distances[i] < min_dist + probe_radius],
                    }
                    cavities.append(cavity)

            return cavities

        except Exception as e:
            self.logger.error(f"Error finding surface cavities: {str(e)}")
            return []

    def calculate_cavity_depth(self, center: np.ndarray, surface_points: np.ndarray) -> float:
        """Calculate depth of cavity from surface.

        Args:
            center: Cavity center coordinates
            surface_points: Surface point coordinates

        Returns:
            Cavity depth in Å
        """
        try:
            distances = np.linalg.norm(surface_points - center, axis=1)
            return float(np.min(distances))
        except Exception as e:
            self.logger.error(f"Error calculating cavity depth: {str(e)}")
            return 0.0

    def identify_exposed_residues(self, structure: Structure, threshold: float = 2.8) -> List[int]:
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
                            radii.append(1.4 + self.ATOM_RADII.get(atom.element, 1.5))

                        if coords:
                            sasa = calculate_sasa(np.array(coords), np.array(radii))
                            if sasa > threshold:
                                exposed.append(residue.get_id()[1])

            return exposed

        except Exception as e:
            self.logger.error(f"Error identifying exposed residues: {str(e)}")
            return []

    def _get_site_atoms(self, structure: Structure, site_residues: List[int]) -> Set[Atom]:
        """Get atoms from specified residues.

        Args:
            structure: BioPython Structure object
            site_residues: List of residue numbers

        Returns:
            Set of atoms from specified residues
        """
        try:
            site_atoms = set()
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[1] in site_residues:
                            for atom in residue:
                                site_atoms.add(atom)
            return site_atoms

        except Exception as e:
            self.logger.error(f"Error getting site atoms: {str(e)}")
            return set()
