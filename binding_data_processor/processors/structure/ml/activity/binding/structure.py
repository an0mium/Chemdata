"""Structure analysis functionality for binding site prediction."""

import logging
from typing import Dict, List, Optional, Tuple, Union, Any
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_angle, calc_dihedral
from Bio.PDB.SASA import calculate_sasa
from Bio.PDB.Polypeptide import three_to_one, is_aa
from Bio.Align import substitution_matrices
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import cdist

from .constants import (
    HYDROPHOBICITY,
    VOLUME,
    CHARGE,
    POCKET_PARAMS,
    GEOMETRY_PARAMS,
    ENERGY_PARAMS,
)

logger = logging.getLogger(__name__)


class StructureAnalyzer:
    """Analyzes protein structures for binding site prediction."""

    def __init__(self):
        """Initialize structure analyzer."""
        self.logger = logging.getLogger(__name__)
        self.blosum = substitution_matrices.load("BLOSUM62")

    def get_structure_properties(
        self,
        structure: Structure,
        include_conservation: bool = True,
        include_dynamics: bool = True,
    ) -> Dict[str, Any]:
        """Extract comprehensive properties from protein structure.

        Args:
            structure: BioPython Structure object
            include_conservation: Whether to include conservation analysis
            include_dynamics: Whether to include dynamic properties

        Returns:
            Dictionary of structure properties
        """
        try:
            properties = {
                "residues": [],
                "coords": [],
                "hydrophobicity": [],
                "charge": [],
                "volume": [],
                "secondary_structure": [],
                "surface_accessibility": [],
                "b_factors": [],
                "residue_depth": [],
            }

            # Get sequence for conservation analysis
            sequence = ""
            residue_ids = []

            # Analyze each residue
            for residue in structure.get_residues():
                if not is_aa(residue):
                    continue

                res_id = residue.get_id()
                res_name = residue.get_resname()

                # Store residue and coordinates
                properties["residues"].append(residue)
                for atom in residue:
                    properties["coords"].append(atom.get_coord())

                # Basic properties
                properties["hydrophobicity"].append(HYDROPHOBICITY.get(res_name, 0.0))
                properties["charge"].append(CHARGE.get(res_name, 0.0))
                properties["volume"].append(VOLUME.get(res_name, 0.0))

                # Secondary structure
                ss = self._get_secondary_structure(residue)
                properties["secondary_structure"].append(ss)

                # Surface accessibility
                sasa = calculate_sasa(residue)
                properties["surface_accessibility"].append(float(sasa))

                # B-factors and depth
                b_factors = [atom.get_bfactor() for atom in residue]
                properties["b_factors"].append(float(np.mean(b_factors)))
                properties["residue_depth"].append(self._calculate_residue_depth(residue))

                # Add to sequence
                sequence += three_to_one(res_name)
                residue_ids.append(res_id[1])

            # Convert lists to arrays
            for key in ["coords", "hydrophobicity", "charge", "volume", "b_factors"]:
                properties[key] = np.array(properties[key])

            # Optional analyses
            if include_conservation:
                properties["conservation"] = self._analyze_conservation(sequence, residue_ids)

            if include_dynamics:
                properties["dynamics"] = self._analyze_dynamics(structure, properties["b_factors"])

            return properties

        except Exception as e:
            self.logger.error(f"Error getting structure properties: {str(e)}")
            return {}

    def find_pockets(
        self,
        coords: np.ndarray,
        min_volume: float = POCKET_PARAMS["min_volume"],
        probe_radius: float = POCKET_PARAMS["probe_radius"],
    ) -> List[Dict[str, Any]]:
        """Find potential binding pockets using geometric analysis.

        Args:
            coords: Atomic coordinates
            min_volume: Minimum pocket volume in Å³
            probe_radius: Probe radius for surface calculation

        Returns:
            List of detected pockets with properties
        """
        try:
            # Calculate alpha shape
            alpha_shape = self._calculate_alpha_shape(coords, probe_radius)

            # Get Voronoi diagram
            vor = Voronoi(coords)

            # Find cavities using alpha shape and Voronoi vertices
            cavities = []
            for v in vor.vertices:
                if self._is_cavity(v, alpha_shape, coords, probe_radius):
                    cavity = {
                        "center": v,
                        "volume": self._estimate_cavity_volume(v, vor, coords),
                        "residues": self._get_cavity_residues(v, coords),
                        "depth": self._calculate_cavity_depth(v, coords),
                        "exposure": self._calculate_exposure(v, coords),
                    }
                    if cavity["volume"] >= min_volume:
                        cavities.append(cavity)

            return cavities

        except Exception as e:
            self.logger.error(f"Error detecting pockets: {str(e)}")
            return []

    def _get_secondary_structure(self, residue: Residue) -> str:
        """Get secondary structure assignment for residue.

        Args:
            residue: BioPython residue object

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
            residue: BioPython residue object

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

    def _calculate_residue_depth(self, residue: Residue) -> float:
        """Calculate residue depth from protein surface.

        Args:
            residue: BioPython residue object

        Returns:
            Depth in Å
        """
        try:
            # Get residue center
            coords = []
            for atom in residue:
                coords.append(atom.get_coord())
            center = np.mean(coords, axis=0)

            # Get surface atoms
            structure = residue.get_parent().get_parent().get_parent()
            surface_atoms = self._get_surface_atoms(structure)

            if len(surface_atoms) > 0:
                # Calculate minimum distance to surface
                distances = np.linalg.norm(surface_atoms - center, axis=1)
                return float(np.min(distances))
            return 0.0

        except Exception as e:
            self.logger.error(f"Error calculating residue depth: {str(e)}")
            return 0.0

    def _get_surface_atoms(self, structure: Structure) -> np.ndarray:
        """Get surface atoms using SASA calculation.

        Args:
            structure: BioPython Structure object

        Returns:
            Array of surface atom coordinates
        """
        try:
            surface_atoms = []
            for residue in structure.get_residues():
                if not is_aa(residue):
                    continue

                # Calculate SASA for each atom
                for atom in residue:
                    sasa = calculate_sasa(residue, atom=atom.get_id())
                    if sasa > 1.0:  # Exposed atom
                        surface_atoms.append(atom.get_coord())

            return np.array(surface_atoms)

        except Exception as e:
            self.logger.error(f"Error getting surface atoms: {str(e)}")
            return np.array([])

    def _analyze_conservation(
        self,
        sequence: str,
        residue_ids: List[int],
        window_size: int = 5,
    ) -> Dict[int, float]:
        """Analyze sequence conservation using sliding window.

        Args:
            sequence: Protein sequence
            residue_ids: List of residue numbers
            window_size: Window size for conservation calculation

        Returns:
            Dictionary mapping residue IDs to conservation scores
        """
        try:
            conservation = {}

            # Calculate conservation scores
            for i, res_id in enumerate(residue_ids):
                # Get window
                start = max(0, i - window_size // 2)
                end = min(len(sequence), i + window_size // 2 + 1)
                window = sequence[start:end]

                # Calculate conservation score
                score = self._calculate_window_conservation(window)
                conservation[res_id] = score

            return conservation

        except Exception as e:
            self.logger.error(f"Error analyzing conservation: {str(e)}")
            return {}

    def _calculate_window_conservation(self, window: str) -> float:
        """Calculate conservation score for sequence window.

        Args:
            window: Sequence window

        Returns:
            Conservation score (0-1)
        """
        try:
            # Calculate average substitution score
            scores = []
            for i in range(len(window)):
                for j in range(i + 1, len(window)):
                    score = self.blosum[window[i]][window[j]]
                    scores.append(score)

            if not scores:
                return 0.0

            # Normalize score
            avg_score = np.mean(scores)
            max_score = max(self.blosum[aa][aa] for aa in self.blosum)
            return float(avg_score / max_score)

        except Exception as e:
            self.logger.error(f"Error calculating conservation: {str(e)}")
            return 0.0

    def _analyze_dynamics(
        self,
        structure: Structure,
        b_factors: np.ndarray,
    ) -> Dict[str, float]:
        """Analyze dynamic properties using B-factors.

        Args:
            structure: BioPython Structure object
            b_factors: Array of B-factor values

        Returns:
            Dictionary of dynamic properties
        """
        try:
            # Calculate overall statistics
            avg_bfactor = float(np.mean(b_factors))
            std_bfactor = float(np.std(b_factors))

            # Calculate relative flexibility
            all_bfactors = []
            for atom in structure.get_atoms():
                all_bfactors.append(atom.get_bfactor())
            relative_flex = avg_bfactor / np.mean(all_bfactors)

            return {
                "average_bfactor": avg_bfactor,
                "bfactor_std": std_bfactor,
                "relative_flexibility": float(relative_flex),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing dynamics: {str(e)}")
            return {
                "average_bfactor": 0.0,
                "bfactor_std": 0.0,
                "relative_flexibility": 0.0,
            }

    def _calculate_alpha_shape(
        self,
        coords: np.ndarray,
        probe_radius: float,
    ) -> Any:
        """Calculate alpha shape of atomic coordinates.

        Args:
            coords: Atomic coordinates
            probe_radius: Probe radius for surface calculation

        Returns:
            Alpha shape object
        """
        try:
            from alphashape import alphashape

            return alphashape(coords, alpha=1.0 / probe_radius)
        except Exception as e:
            self.logger.error(f"Error calculating alpha shape: {str(e)}")
            return None

    def _is_cavity(
        self,
        point: np.ndarray,
        alpha_shape: Any,
        coords: np.ndarray,
        probe_radius: float,
    ) -> bool:
        """Check if a point represents a cavity.

        Args:
            point: Point to check
            alpha_shape: Alpha shape of protein
            coords: Atomic coordinates
            probe_radius: Probe radius

        Returns:
            True if point is in a cavity
        """
        try:
            # Check if point is inside alpha shape
            if not alpha_shape.contains(point):
                return False

            # Check distances to atoms
            distances = np.linalg.norm(coords - point, axis=1)
            return np.all(distances > probe_radius)

        except Exception as e:
            self.logger.error(f"Error checking cavity: {str(e)}")
            return False

    def _estimate_cavity_volume(
        self,
        center: np.ndarray,
        vor: Voronoi,
        coords: np.ndarray,
        max_radius: float = 10.0,
    ) -> float:
        """Estimate cavity volume using Voronoi cells.

        Args:
            center: Cavity center point
            vor: Voronoi diagram
            coords: Atomic coordinates
            max_radius: Maximum radius to consider

        Returns:
            Estimated cavity volume in Å³
        """
        try:
            # Find nearby Voronoi vertices
            distances = np.linalg.norm(vor.vertices - center, axis=1)
            nearby = vor.vertices[distances < max_radius]

            if len(nearby) < 4:
                return 0.0

            # Calculate convex hull volume
            hull = ConvexHull(nearby)
            return hull.volume

        except Exception as e:
            self.logger.error(f"Error estimating cavity volume: {str(e)}")
            return 0.0

    def _get_cavity_residues(
        self,
        center: np.ndarray,
        coords: np.ndarray,
        cutoff: float = GEOMETRY_PARAMS["pocket_overlap_dist"],
    ) -> List[int]:
        """Get residues forming a cavity.

        Args:
            center: Cavity center point
            coords: Atomic coordinates
            cutoff: Distance cutoff for residue inclusion

        Returns:
            List of residue indices
        """
        try:
            # Calculate distances to cavity center
            distances = np.linalg.norm(coords - center, axis=1)

            # Get indices of nearby atoms
            nearby = np.where(distances < cutoff)[0]

            # Convert to residue indices (approximate)
            residue_indices = list(set([i // 10 for i in nearby]))
            return residue_indices

        except Exception as e:
            self.logger.error(f"Error getting cavity residues: {str(e)}")
            return []

    def _calculate_cavity_depth(
        self,
        center: np.ndarray,
        coords: np.ndarray,
    ) -> float:
        """Calculate depth of cavity from protein surface.

        Args:
            center: Cavity center point
            coords: Atomic coordinates

        Returns:
            Cavity depth in Å
        """
        try:
            # Get surface atoms
            surface_atoms = self._get_surface_atoms(coords)

            if len(surface_atoms) == 0:
                return 0.0

            # Calculate minimum distance to surface
            distances = np.linalg.norm(surface_atoms - center, axis=1)
            return float(np.min(distances))

        except Exception as e:
            self.logger.error(f"Error calculating cavity depth: {str(e)}")
            return 0.0

    def _calculate_exposure(
        self,
        center: np.ndarray,
        coords: np.ndarray,
        num_directions: int = 100,
    ) -> float:
        """Calculate solvent exposure of cavity.

        Args:
            center: Cavity center point
            coords: Atomic coordinates
            num_directions: Number of directions to check

        Returns:
            Exposure score (0-1)
        """
        try:
            # Generate random directions
            directions = np.random.randn(num_directions, 3)
            directions /= np.linalg.norm(directions, axis=1)[:, np.newaxis]

            # Ray casting
            exposed = 0
            for direction in directions:
                ray = center + direction * np.arange(0, 20, 0.5)[:, np.newaxis]
                distances = cdist(ray, coords)
                if np.all(distances.min(axis=1) > 2.0):  # No collisions
                    exposed += 1

            return exposed / num_directions

        except Exception as e:
            self.logger.error(f"Error calculating exposure: {str(e)}")
            return 1.0
