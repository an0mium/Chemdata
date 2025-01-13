"""Geometric analysis functionality for protein structures."""

import logging
from typing import Dict, List, Optional, Tuple, Any, Union
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_angle, calc_dihedral
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)


def calculate_surface_area(coords: Union[List[List[float]], np.ndarray]) -> float:
    """Calculate surface area using convex hull.

    Args:
        coords: List of coordinates or numpy array of shape (N, 3)

    Returns:
        Surface area in Å²
    """
    try:
        coords = np.array(coords)
        if len(coords) < 4:
            return 0.0

        hull = ConvexHull(coords)
        return float(hull.area)

    except Exception as e:
        logger.error(f"Error calculating surface area: {str(e)}")
        return 0.0


class GeometryCalculator:
    """Handles geometric calculations and analysis for protein structures."""

    def __init__(self):
        """Initialize geometry calculator."""
        self.logger = logging.getLogger(__name__)

    def get_properties(self, structure: Structure) -> Dict[str, Any]:
        """Calculate comprehensive geometric properties of structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of geometric properties
        """
        try:
            # Get coordinates and basic properties
            coords = self.get_coordinates(structure)

            properties = {
                # Basic geometric properties
                "coords": coords,
                "center_of_mass": self.calculate_center_of_mass(structure),
                "radius_of_gyration": self.calculate_radius_of_gyration(structure),
                "principal_axes": self.calculate_principal_axes(structure),
                "inertia_tensor": self.calculate_inertia_tensor(structure),
                "dimensions": self.calculate_dimensions(structure),
                # Volume and surface properties
                "volume": self.calculate_volume(structure),
                "surface_area": self.calculate_surface_area(structure),
                "packing_density": self.calculate_packing_density(structure),
                # Shape analysis
                "shape_descriptors": self.calculate_shape_descriptors(structure),
                # Structure analysis
                "secondary_structure": self.analyze_secondary_structure(structure),
                "backbone_angles": [],
            }

            # Calculate backbone angles
            for residue in structure.get_residues():
                phi, psi = self._calculate_backbone_angles(residue)
                properties["backbone_angles"].append((phi, psi))

            return properties

        except Exception as e:
            self.logger.error(f"Error calculating geometric properties: {str(e)}")
            return {}

    def analyze_site(
        self,
        structure: Structure,
        site_residues: List[int],
    ) -> Dict[str, Any]:
        """Analyze geometric properties of binding site.

        Args:
            structure: BioPython Structure object
            site_residues: List of residue numbers in binding site

        Returns:
            Dictionary of geometric properties
        """
        try:
            # Get residue objects
            residues = [structure[0]["A"][res_num] for res_num in site_residues]

            properties = {
                # Basic properties
                "volume": self.calculate_volume(residues),
                "surface_area": self.calculate_surface_area(residues),
                "depth": self._calculate_depth(residues),
                # Shape properties
                "shape": self.calculate_shape_descriptors(residues),
                "center": self.calculate_center_of_mass(residues),
                "dimensions": self.calculate_dimensions(residues),
                # Packing properties
                "packing": self._calculate_packing(residues),
                "packing_density": self.calculate_packing_density(residues),
                # Structure properties
                "secondary_structure": self._analyze_site_structure(residues),
            }

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing site geometry: {str(e)}")
            return {}

    def get_coordinates(self, structure: Structure) -> np.ndarray:
        """Get all atom coordinates from structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Array of coordinates
        """
        try:
            coords = []
            for atom in structure.get_atoms():
                coords.append(atom.get_coord())
            return np.array(coords)

        except Exception as e:
            self.logger.error(f"Error getting coordinates: {str(e)}")
            return np.array([])

    def analyze_residue(self, residue: Residue) -> Dict[str, Any]:
        """Analyze geometric properties of residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Dictionary of residue properties
        """
        try:
            properties = {
                # Position and volume
                "center": self._calculate_residue_center(residue),
                "volume": self._calculate_residue_volume(residue),
                "surface_area": self._calculate_residue_surface_area(residue),
                # Chemical properties
                "hydrophobicity": self.get_hydrophobicity(residue),
                "charge": self.get_charge(residue),
                # Packing properties
                "packing_density": self._calculate_residue_packing_density(residue),
                "electrostatics": self._calculate_residue_electrostatics(residue),
                # Shape properties
                "shape": self._calculate_residue_shape(residue),
                # Structure
                "secondary_structure": self.get_secondary_structure(residue),
            }
            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing residue: {str(e)}")
            return {}

    def calculate_center_of_mass(self, structure: Structure) -> np.ndarray:
        """Calculate center of mass of structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Array of center coordinates
        """
        try:
            coords = []
            masses = []
            for atom in structure.get_atoms():
                coords.append(atom.get_coord())
                masses.append(self._get_atomic_mass(atom))

            coords = np.array(coords)
            masses = np.array(masses)

            return np.average(coords, weights=masses, axis=0)

        except Exception as e:
            self.logger.error(f"Error calculating center of mass: {str(e)}")
            return np.zeros(3)

    def calculate_radius_of_gyration(self, structure: Structure) -> float:
        """Calculate radius of gyration of structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Radius of gyration in Å
        """
        try:
            com = self.calculate_center_of_mass(structure)

            # Calculate distances to center of mass
            distances = []
            masses = []
            for atom in structure.get_atoms():
                dist = np.linalg.norm(atom.get_coord() - com)
                distances.append(dist)
                masses.append(self._get_atomic_mass(atom))

            distances = np.array(distances)
            masses = np.array(masses)

            # Calculate Rg
            rg = np.sqrt(np.sum(masses * distances**2) / np.sum(masses))
            return float(rg)

        except Exception as e:
            self.logger.error(f"Error calculating radius of gyration: {str(e)}")
            return 0.0

    def calculate_principal_axes(self, structure: Structure) -> np.ndarray:
        """Calculate principal axes of structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Array of principal axes vectors
        """
        try:
            coords = self.get_coordinates(structure)
            if len(coords) < 3:
                return np.eye(3)

            # Center coordinates
            com = np.mean(coords, axis=0)
            centered = coords - com

            # Calculate covariance matrix
            cov = np.cov(centered.T)

            # Get eigenvectors (principal axes)
            eigenvals, eigenvecs = np.linalg.eigh(cov)

            # Sort by eigenvalue in descending order
            idx = np.argsort(eigenvals)[::-1]
            return eigenvecs[:, idx]

        except Exception as e:
            self.logger.error(f"Error calculating principal axes: {str(e)}")
            return np.eye(3)

    def calculate_inertia_tensor(self, structure: Structure) -> np.ndarray:
        """Calculate inertia tensor of structure.

        Args:
            structure: BioPython Structure object

        Returns:
            3x3 inertia tensor matrix
        """
        try:
            coords = self.get_coordinates(structure)
            if len(coords) < 3:
                return np.zeros((3, 3))

            # Center coordinates
            com = np.mean(coords, axis=0)
            centered = coords - com

            # Calculate inertia tensor
            tensor = np.zeros((3, 3))
            for coord in centered:
                r2 = np.sum(coord * coord)
                tensor += np.diag([r2, r2, r2]) - np.outer(coord, coord)

            return tensor

        except Exception as e:
            self.logger.error(f"Error calculating inertia tensor: {str(e)}")
            return np.zeros((3, 3))

    def calculate_dimensions(self, structure: Structure) -> Dict[str, float]:
        """Calculate dimensions of structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of dimensions
        """
        try:
            coords = self.get_coordinates(structure)
            if len(coords) < 2:
                return {"x": 0.0, "y": 0.0, "z": 0.0}

            # Get min/max along each axis
            min_coords = np.min(coords, axis=0)
            max_coords = np.max(coords, axis=0)
            dimensions = max_coords - min_coords

            return {
                "x": float(dimensions[0]),
                "y": float(dimensions[1]),
                "z": float(dimensions[2]),
            }

        except Exception as e:
            self.logger.error(f"Error calculating dimensions: {str(e)}")
            return {"x": 0.0, "y": 0.0, "z": 0.0}

    def calculate_volume(self, structure: Structure) -> float:
        """Calculate volume of structure using convex hull.

        Args:
            structure: BioPython Structure object

        Returns:
            Volume in Å³
        """
        try:
            coords = self.get_coordinates(structure)
            if len(coords) < 4:
                return 0.0

            hull = ConvexHull(coords)
            return float(hull.volume)

        except Exception as e:
            self.logger.error(f"Error calculating volume: {str(e)}")
            return 0.0

    def calculate_surface_area(self, structure: Structure) -> float:
        """Calculate surface area of structure using convex hull.

        Args:
            structure: BioPython Structure object

        Returns:
            Surface area in Å²
        """
        try:
            coords = self.get_coordinates(structure)
            if len(coords) < 4:
                return 0.0

            hull = ConvexHull(coords)
            return float(hull.area)

        except Exception as e:
            self.logger.error(f"Error calculating surface area: {str(e)}")
            return 0.0

    def calculate_packing_density(self, structure: Structure) -> float:
        """Calculate packing density of structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Packing density (0-1)
        """
        try:
            volume = self.calculate_volume(structure)
            if volume == 0:
                return 0.0

            # Calculate van der Waals volume
            vdw_volume = 0.0
            for atom in structure.get_atoms():
                radius = self._get_vdw_radius(atom)
                vdw_volume += 4 / 3 * np.pi * radius**3

            return float(vdw_volume / volume)

        except Exception as e:
            self.logger.error(f"Error calculating packing density: {str(e)}")
            return 0.0

    def calculate_shape_descriptors(self, structure: Structure) -> Dict[str, float]:
        """Calculate shape descriptors of structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of shape descriptors
        """
        try:
            coords = self.get_coordinates(structure)
            if len(coords) < 4:
                return {
                    "sphericity": 0.0,
                    "asphericity": 0.0,
                    "eccentricity": 0.0,
                    "elongation": 0.0,
                    "flatness": 0.0,
                }

            # Calculate principal components
            centered = coords - np.mean(coords, axis=0)
            cov = np.cov(centered.T)
            eigenvals = np.linalg.eigvals(cov)
            eigenvals.sort()

            # Shape descriptors
            a, b, c = eigenvals
            descriptors = {
                "sphericity": float(c / a),  # 1 for perfect sphere
                "asphericity": float((a - (b + c) / 2) / a),  # 0 for perfect sphere
                "eccentricity": float(np.sqrt(1 - c / a)),  # 0 for perfect sphere
                "elongation": float(a / b),  # 1 for sphere
                "flatness": float(b / c),  # 1 for sphere
            }

            return descriptors

        except Exception as e:
            self.logger.error(f"Error calculating shape descriptors: {str(e)}")
            return {
                "sphericity": 0.0,
                "asphericity": 0.0,
                "eccentricity": 0.0,
                "elongation": 0.0,
                "flatness": 0.0,
            }

    def analyze_secondary_structure(self, structure: Structure) -> Dict[str, float]:
        """Analyze secondary structure composition.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of secondary structure percentages
        """
        try:
            ss_counts = {"H": 0, "E": 0, "C": 0}  # Helix, Sheet, Coil
            total = 0

            for residue in structure.get_residues():
                ss = self.get_secondary_structure(residue)
                ss_counts[ss] += 1
                total += 1

            if total == 0:
                return {"H": 0.0, "E": 0.0, "C": 0.0}

            return {k: float(v) / total * 100 for k, v in ss_counts.items()}

        except Exception as e:
            self.logger.error(f"Error analyzing secondary structure: {str(e)}")
            return {"H": 0.0, "E": 0.0, "C": 0.0}

    def get_secondary_structure(self, residue: Residue) -> str:
        """Get secondary structure assignment for residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Secondary structure type (H, E, C)
        """
        try:
            # Calculate phi/psi angles
            phi, psi = self._calculate_backbone_angles(residue)

            if phi is None or psi is None:
                return "C"  # Coil

            # Ramachandran plot regions
            if -140 < phi < -60 and -70 < psi < -15:
                return "H"  # Alpha helix
            elif -150 < phi < -50 and 100 < psi < 180:
                return "E"  # Beta sheet
            else:
                return "C"  # Coil

        except Exception as e:
            self.logger.error(f"Error getting secondary structure: {str(e)}")
            return "C"

    def _calculate_backbone_angles(self, residue: Residue) -> Tuple[Optional[float], Optional[float]]:
        """Calculate backbone phi/psi angles.

        Args:
            residue: BioPython Residue object

        Returns:
            Tuple of (phi angle, psi angle) in degrees
        """
        try:
            # Get required atoms
            if not all(atom in residue for atom in ["N", "CA", "C"]):
                return None, None

            # Get previous and next residues
            prev_res = residue.get_previous_residue()
            next_res = residue.get_next_residue()

            # Calculate phi angle (requires previous residue)
            phi = None
            if prev_res and "C" in prev_res:
                phi = calc_dihedral(
                    prev_res["C"].get_vector(),
                    residue["N"].get_vector(),
                    residue["CA"].get_vector(),
                    residue["C"].get_vector(),
                )

            # Calculate psi angle (requires next residue)
            psi = None
            if next_res and "N" in next_res:
                psi = calc_dihedral(
                    residue["N"].get_vector(),
                    residue["CA"].get_vector(),
                    residue["C"].get_vector(),
                    next_res["N"].get_vector(),
                )

            return phi, psi

        except Exception as e:
            self.logger.error(f"Error calculating backbone angles: {str(e)}")
            return None, None

    def _calculate_depth(self, residues: List[Residue]) -> float:
        """Calculate depth of residue selection from protein surface.

        Args:
            residues: List of residues

        Returns:
            Depth in Å
        """
        try:
            # Get residue centers
            centers = []
            for res in residues:
                center = np.mean([atom.get_coord() for atom in res], axis=0)
                centers.append(center)
            centers = np.array(centers)

            # Get structure surface
            structure = residues[0].get_parent().get_parent().get_parent()
            surface_atoms = []
            for atom in structure.get_atoms():
                if not atom.is_disordered():
                    surface_atoms.append(atom.get_coord())
            surface_atoms = np.array(surface_atoms)

            if len(surface_atoms) == 0:
                return 0.0

            # Calculate minimum distance to surface
            min_dist = float("inf")
            for center in centers:
                distances = np.linalg.norm(surface_atoms - center, axis=1)
                min_dist = min(min_dist, np.min(distances))

            return min_dist

        except Exception as e:
            self.logger.error(f"Error calculating depth: {str(e)}")
            return 0.0

    def _calculate_packing(self, residues: List[Residue]) -> Dict[str, float]:
        """Calculate packing density and efficiency.

        Args:
            residues: List of residues

        Returns:
            Dictionary of packing properties
        """
        try:
            # Get coordinates and volumes
            coords = []
            volumes = []
            for res in residues:
                for atom in res:
                    coords.append(atom.get_coord())
                    # Approximate atom volume from radius
                    radius = self._get_vdw_radius(atom)
                    volumes.append(4 / 3 * np.pi * radius**3)
            coords = np.array(coords)
            volumes = np.array(volumes)

            if len(coords) < 4:
                return {}

            # Calculate total volume and sum of atomic volumes
            hull = ConvexHull(coords)
            total_volume = hull.volume
            atomic_volume = np.sum(volumes)

            return {
                "packing_density": atomic_volume / total_volume,
                "void_volume": total_volume - atomic_volume,
                "efficiency": atomic_volume / (total_volume * len(coords)),
            }

        except Exception as e:
            self.logger.error(f"Error calculating packing: {str(e)}")
            return {}

    def _analyze_site_structure(self, residues: List[Residue]) -> Dict[str, float]:
        """Analyze secondary structure composition of site.

        Args:
            residues: List of residues

        Returns:
            Dictionary of structure properties
        """
        try:
            ss_counts = {"H": 0, "E": 0, "C": 0}  # Helix, Sheet, Coil
            total = 0

            for res in residues:
                ss = self.get_secondary_structure(res)
                ss_counts[ss] += 1
                total += 1

            if total == 0:
                return {"H": 0.0, "E": 0.0, "C": 0.0}

            return {k: v / total for k, v in ss_counts.items()}

        except Exception as e:
            self.logger.error(f"Error analyzing site structure: {str(e)}")
            return {"H": 0.0, "E": 0.0, "C": 0.0}

    def _get_atomic_mass(self, atom) -> float:
        """Get atomic mass for an atom.

        Args:
            atom: BioPython Atom object

        Returns:
            Atomic mass in Da
        """
        masses = {
            "H": 1.008,
            "C": 12.011,
            "N": 14.007,
            "O": 15.999,
            "S": 32.06,
            "P": 30.974,
        }
        return masses.get(atom.element, 0.0)

    def _get_vdw_radius(self, atom) -> float:
        """Get van der Waals radius for atom.

        Args:
            atom: BioPython Atom object

        Returns:
            Radius in Å
        """
        radii = {
            "H": 1.2,
            "C": 1.7,
            "N": 1.55,
            "O": 1.52,
            "S": 1.8,
            "P": 1.8,
            "F": 1.47,
            "Cl": 1.75,
            "Br": 1.85,
            "I": 1.98,
        }
        return radii.get(atom.element, 1.5)

    def get_hydrophobicity(self, residue: Residue) -> float:
        """Get hydrophobicity value for residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Hydrophobicity score
        """
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

    def get_charge(self, residue: Residue) -> float:
        """Get charge value for residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Charge value
        """
        charge = {
            "ARG": 1.0,
            "LYS": 1.0,
            "ASP": -1.0,
            "GLU": -1.0,
            "HIS": 0.1,
        }
        return charge.get(residue.get_resname(), 0.0)

    def _calculate_residue_center(self, residue: Residue) -> np.ndarray:
        """Calculate geometric center of residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Array of center coordinates
        """
        try:
            coords = []
            for atom in residue:
                coords.append(atom.get_coord())
            return np.mean(coords, axis=0)
        except Exception as e:
            self.logger.error(f"Error calculating residue center: {str(e)}")
            return np.zeros(3)

    def _calculate_residue_volume(self, residue: Residue) -> float:
        """Calculate volume of residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Volume in Å³
        """
        try:
            coords = []
            for atom in residue:
                coords.append(atom.get_coord())
            coords = np.array(coords)

            if len(coords) < 4:
                return 0.0

            hull = ConvexHull(coords)
            return float(hull.volume)

        except Exception as e:
            self.logger.error(f"Error calculating residue volume: {str(e)}")
            return 0.0

    def _calculate_residue_surface_area(self, residue: Residue) -> float:
        """Calculate surface area of residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Surface area in Å²
        """
        try:
            coords = []
            for atom in residue:
                coords.append(atom.get_coord())
            coords = np.array(coords)

            if len(coords) < 4:
                return 0.0

            hull = ConvexHull(coords)
            return float(hull.area)

        except Exception as e:
            self.logger.error(f"Error calculating residue surface area: {str(e)}")
            return 0.0

    def _calculate_residue_packing_density(self, residue: Residue) -> float:
        """Calculate packing density of residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Packing density (0-1)
        """
        try:
            volume = self._calculate_residue_volume(residue)
            if volume == 0:
                return 0.0

            # Calculate van der Waals volume
            vdw_volume = 0.0
            for atom in residue:
                radius = self._get_vdw_radius(atom)
                vdw_volume += 4 / 3 * np.pi * radius**3

            return float(vdw_volume / volume)

        except Exception as e:
            self.logger.error(f"Error calculating residue packing density: {str(e)}")
            return 0.0

    def _calculate_residue_electrostatics(self, residue: Residue) -> Dict[str, float]:
        """Calculate electrostatic properties of residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Dictionary of electrostatic properties
        """
        try:
            charge = self.get_charge(residue)
            return {
                "charge": charge,
                "charge_density": charge / self._calculate_residue_volume(residue) if charge != 0 else 0.0,
            }

        except Exception as e:
            self.logger.error(f"Error calculating residue electrostatics: {str(e)}")
            return {"charge": 0.0, "charge_density": 0.0}

    def _calculate_residue_shape(self, residue: Residue) -> Dict[str, float]:
        """Calculate shape properties of residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Dictionary of shape properties
        """
        try:
            coords = []
            for atom in residue:
                coords.append(atom.get_coord())
            coords = np.array(coords)

            if len(coords) < 4:
                return {}

            return self.calculate_shape_descriptors(coords)

        except Exception as e:
            self.logger.error(f"Error calculating residue shape: {str(e)}")
            return {}
