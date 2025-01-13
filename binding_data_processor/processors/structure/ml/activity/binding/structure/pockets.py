"""Binding pocket detection and analysis functionality."""

import logging
from typing import Dict, List, Optional, Any, Set, Tuple
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_angle, calc_dihedral
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import cdist
from scipy.cluster.hierarchy import linkage, fcluster
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers

logger = logging.getLogger(__name__)


class PocketDetector:
    """Detects and analyzes protein binding pockets."""

    def __init__(self):
        """Initialize pocket detector."""
        self.logger = logging.getLogger(__name__)

    def find_pockets(
        self,
        structure: Structure,
        min_volume: float = 100.0,
        probe_radius: float = 1.4,
        properties: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Find potential binding pockets using multiple methods.

        Args:
            structure: BioPython Structure object
            min_volume: Minimum pocket volume in Å³
            probe_radius: Probe radius for surface calculation
            properties: Pre-calculated structure properties

        Returns:
            List of dictionaries containing pocket properties
        """
        try:
            # Get structure coordinates
            coords = np.array([atom.get_coord() for atom in structure.get_atoms()])

            # Find surface pockets using geometric criteria
            surface_pockets = self._find_surface_pockets(
                coords,
                probe_radius=probe_radius,
            )

            # Find cavities using Voronoi tessellation
            voronoi_pockets = self._find_voronoi_pockets(
                coords,
                min_volume=min_volume,
            )

            # Merge overlapping pockets
            merged_pockets = self._merge_pockets(
                surface_pockets + voronoi_pockets,
                coords,
            )

            # Filter by size and add properties
            filtered_pockets = []
            for pocket in merged_pockets:
                volume = self._calculate_pocket_volume(pocket["coords"])
                if volume >= min_volume:
                    pocket["volume"] = volume
                    pocket["residues"] = self._get_pocket_residues(
                        structure,
                        pocket["coords"],
                    )
                    pocket.update(
                        self._calculate_pocket_properties(
                            pocket["coords"],
                            pocket["residues"],
                            properties,
                        )
                    )
                    filtered_pockets.append(pocket)

            return filtered_pockets

        except Exception as e:
            self.logger.error(f"Error finding pockets: {str(e)}")
            return []

    def score_pockets(
        self,
        pockets: List[Dict[str, Any]],
        structure: Structure,
        properties: Dict[str, Any],
        include_scores: bool = True,
    ) -> List[Dict[str, Any]]:
        """Score potential binding pockets.

        Args:
            pockets: List of pocket dictionaries
            structure: BioPython Structure object
            properties: Structure properties
            include_scores: Whether to include detailed scoring

        Returns:
            List of scored pocket dictionaries
        """
        try:
            scored_pockets = []

            for pocket in pockets:
                # Calculate basic scores
                volume_score = self._score_volume(pocket["volume"])
                shape_score = self._score_shape(pocket["shape_descriptors"])
                conservation_score = self._score_conservation(
                    pocket["residues"],
                    properties.get("conservation", {}),
                )
                druggability_score = self._score_druggability(pocket)

                # Calculate total score
                total_score = 0.3 * volume_score + 0.3 * shape_score + 0.2 * conservation_score + 0.2 * druggability_score

                # Add scores to pocket
                scored_pocket = pocket.copy()
                scored_pocket["score"] = float(total_score)

                if include_scores:
                    scored_pocket["component_scores"] = {
                        "volume": float(volume_score),
                        "shape": float(shape_score),
                        "conservation": float(conservation_score),
                        "druggability": float(druggability_score),
                    }

                scored_pockets.append(scored_pocket)

            return scored_pockets

        except Exception as e:
            self.logger.error(f"Error scoring pockets: {str(e)}")
            return []

    def _find_surface_pockets(
        self,
        coords: np.ndarray,
        probe_radius: float = 1.4,
    ) -> List[Dict[str, Any]]:
        """Find pockets on protein surface using probe rolling.

        Args:
            coords: Atomic coordinates
            probe_radius: Probe radius in Å

        Returns:
            List of surface pocket properties
        """
        try:
            pockets = []

            # Generate probe positions on sphere
            theta = np.linspace(0, np.pi, 20)
            phi = np.linspace(0, 2 * np.pi, 40)
            probe_points = []
            for t in theta:
                for p in phi:
                    x = probe_radius * np.sin(t) * np.cos(p)
                    y = probe_radius * np.sin(t) * np.sin(p)
                    z = probe_radius * np.cos(t)
                    probe_points.append([x, y, z])
            probe_points = np.array(probe_points)

            # Roll probe over surface
            for i, coord in enumerate(coords):
                # Check each probe position
                positions = probe_points + coord
                distances = cdist(positions, coords)

                # Find accessible positions
                accessible = np.all(distances > probe_radius, axis=1)
                if np.any(accessible):
                    # Found potential pocket
                    pocket_coords = positions[accessible]
                    pocket = {
                        "coords": pocket_coords,
                        "center": np.mean(pocket_coords, axis=0),
                        "type": "surface",
                    }
                    pockets.append(pocket)

            return pockets

        except Exception as e:
            self.logger.error(f"Error finding surface pockets: {str(e)}")
            return []

    def _find_voronoi_pockets(
        self,
        coords: np.ndarray,
        min_volume: float = 100.0,
    ) -> List[Dict[str, Any]]:
        """Find pockets using Voronoi tessellation.

        Args:
            coords: Atomic coordinates
            min_volume: Minimum pocket volume in Å³

        Returns:
            List of Voronoi pocket properties
        """
        try:
            # Calculate Voronoi diagram
            vor = Voronoi(coords)
            pockets = []

            # Check each Voronoi vertex
            for vertex in vor.vertices:
                # Calculate distances to atoms
                distances = np.linalg.norm(coords - vertex, axis=1)
                min_dist = np.min(distances)

                # Check if vertex represents a cavity
                if min_dist > 2.8:  # Minimum cavity size
                    # Get surrounding atoms
                    nearby = coords[distances < min_dist + 2.8]
                    if len(nearby) >= 4:  # Minimum surrounding atoms
                        volume = 4 / 3 * np.pi * min_dist**3
                        if volume >= min_volume:
                            pocket = {
                                "coords": nearby,
                                "center": vertex,
                                "type": "voronoi",
                            }
                            pockets.append(pocket)

            return pockets

        except Exception as e:
            self.logger.error(f"Error finding Voronoi pockets: {str(e)}")
            return []

    def _merge_pockets(
        self,
        pockets: List[Dict[str, Any]],
        coords: np.ndarray,
    ) -> List[Dict[str, Any]]:
        """Merge overlapping pockets.

        Args:
            pockets: List of pocket dictionaries
            coords: Atomic coordinates

        Returns:
            List of merged pocket dictionaries
        """
        try:
            if not pockets:
                return []

            # Calculate distances between pocket centers
            centers = np.array([p["center"] for p in pockets])
            distances = cdist(centers, centers)

            # Cluster pockets
            linkage_matrix = linkage(distances, method="single")
            clusters = fcluster(linkage_matrix, t=5.0, criterion="distance")

            # Merge pockets in each cluster
            merged = []
            for cluster_id in np.unique(clusters):
                cluster_indices = np.where(clusters == cluster_id)[0]

                # Combine coordinates
                cluster_coords = np.concatenate([pockets[i]["coords"] for i in cluster_indices])

                # Remove duplicates
                cluster_coords = np.unique(cluster_coords, axis=0)

                merged.append(
                    {
                        "coords": cluster_coords,
                        "center": np.mean(cluster_coords, axis=0),
                        "type": "merged",
                    }
                )

            return merged

        except Exception as e:
            self.logger.error(f"Error merging pockets: {str(e)}")
            return []

    def _calculate_pocket_volume(self, coords: np.ndarray) -> float:
        """Calculate pocket volume.

        Args:
            coords: Pocket coordinates

        Returns:
            Volume in Å³
        """
        try:
            if len(coords) < 4:
                return 0.0

            hull = ConvexHull(coords)
            return hull.volume

        except Exception as e:
            self.logger.error(f"Error calculating pocket volume: {str(e)}")
            return 0.0

    def _get_pocket_residues(
        self,
        structure: Structure,
        coords: np.ndarray,
        cutoff: float = 4.0,
    ) -> List[int]:
        """Get residues forming pocket.

        Args:
            structure: BioPython Structure object
            coords: Pocket coordinates
            cutoff: Distance cutoff in Å

        Returns:
            List of residue numbers
        """
        try:
            pocket_residues = set()

            for residue in structure.get_residues():
                # Check if any atom is near pocket
                for atom in residue:
                    distances = np.linalg.norm(
                        coords - atom.get_coord(),
                        axis=1,
                    )
                    if np.any(distances < cutoff):
                        pocket_residues.add(residue.get_id()[1])
                        break

            return sorted(list(pocket_residues))

        except Exception as e:
            self.logger.error(f"Error getting pocket residues: {str(e)}")
            return []

    def _calculate_pocket_properties(
        self,
        coords: np.ndarray,
        residues: List[int],
        properties: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Calculate pocket properties.

        Args:
            coords: Pocket coordinates
            residues: Pocket residue numbers
            properties: Structure properties

        Returns:
            Dictionary of pocket properties
        """
        try:
            pocket_props = {
                "shape_descriptors": self._calculate_shape_descriptors(coords),
                "depth": self._calculate_pocket_depth(coords),
                "exposure": self._calculate_surface_exposure(coords),
            }

            # Add residue-based properties if available
            if properties:
                if "hydrophobicity" in properties:
                    hydrophobicity = [properties["hydrophobicity"][i] for i in residues if i < len(properties["hydrophobicity"])]
                    if hydrophobicity:
                        pocket_props["hydrophobicity"] = float(np.mean(hydrophobicity))

                if "conservation" in properties:
                    conservation = [properties["conservation"].get(res, 0.0) for res in residues]
                    if conservation:
                        pocket_props["conservation"] = float(np.mean(conservation))

            return pocket_props

        except Exception as e:
            self.logger.error(f"Error calculating pocket properties: {str(e)}")
            return {}

    def _calculate_shape_descriptors(self, coords: np.ndarray) -> Dict[str, float]:
        """Calculate shape descriptors for pocket.

        Args:
            coords: Pocket coordinates

        Returns:
            Dictionary of shape descriptors
        """
        try:
            if len(coords) < 4:
                return {}

            # Calculate principal components
            centered = coords - np.mean(coords, axis=0)
            cov = np.cov(centered.T)
            eigenvals, eigenvecs = np.linalg.eigh(cov)

            # Sort by eigenvalue
            idx = np.argsort(eigenvals)[::-1]
            eigenvals = eigenvals[idx]
            eigenvecs = eigenvecs[:, idx]

            return {
                "sphericity": float(np.min(eigenvals) / np.max(eigenvals)),
                "elongation": float(eigenvals[0] / eigenvals[1]),
                "flatness": float(eigenvals[1] / eigenvals[2]),
                "volume_over_surface": float(self._calculate_pocket_volume(coords) / ConvexHull(coords).area),
            }

        except Exception as e:
            self.logger.error(f"Error calculating shape descriptors: {str(e)}")
            return {}

    def _calculate_pocket_depth(self, coords: np.ndarray) -> float:
        """Calculate pocket depth.

        Args:
            coords: Pocket coordinates

        Returns:
            Depth in Å
        """
        try:
            if len(coords) < 2:
                return 0.0

            # Get pocket center and entrance
            center = np.mean(coords, axis=0)
            entrance = coords[np.argmax(coords[:, 2])]  # Highest z-coordinate

            return float(np.linalg.norm(center - entrance))

        except Exception as e:
            self.logger.error(f"Error calculating pocket depth: {str(e)}")
            return 0.0

    def _calculate_surface_exposure(self, coords: np.ndarray) -> float:
        """Calculate pocket surface exposure.

        Args:
            coords: Pocket coordinates

        Returns:
            Surface exposure score (0-1)
        """
        try:
            if len(coords) < 4:
                return 1.0

            # Calculate convex hull surface area
            hull = ConvexHull(coords)
            surface_area = hull.area

            # Calculate theoretical sphere surface area
            radius = np.mean(
                np.linalg.norm(
                    coords - np.mean(coords, axis=0),
                    axis=1,
                )
            )
            sphere_area = 4 * np.pi * radius * radius

            return float(surface_area / sphere_area)

        except Exception as e:
            self.logger.error(f"Error calculating surface exposure: {str(e)}")
            return 1.0

    def _score_volume(self, volume: float) -> float:
        """Score pocket volume.

        Args:
            volume: Pocket volume in Å³

        Returns:
            Volume score (0-1)
        """
        try:
            # Prefer pockets between 200-1000 Å³
            if volume < 100:
                return 0.0
            elif volume < 200:
                return volume / 200
            elif volume < 1000:
                return 1.0
            else:
                return np.exp(-(volume - 1000) / 500)

        except Exception as e:
            self.logger.error(f"Error scoring volume: {str(e)}")
            return 0.0

    def _score_shape(self, descriptors: Dict[str, float]) -> float:
        """Score pocket shape.

        Args:
            descriptors: Shape descriptor dictionary

        Returns:
            Shape score (0-1)
        """
        try:
            if not descriptors:
                return 0.0

            # Prefer roughly spherical pockets
            sphericity = descriptors.get("sphericity", 0.0)
            elongation = descriptors.get("elongation", 1.0)
            flatness = descriptors.get("flatness", 1.0)
            vol_surf = descriptors.get("volume_over_surface", 0.0)

            shape_score = 0.4 * sphericity + 0.2 * (1 - elongation) + 0.2 * (1 - flatness) + 0.2 * min(vol_surf / 2.0, 1.0)

            return float(shape_score)

        except Exception as e:
            self.logger.error(f"Error scoring shape: {str(e)}")
            return 0.0

    def _score_conservation(
        self,
        residues: List[int],
        conservation: Dict[int, float],
    ) -> float:
        """Score pocket conservation.

        Args:
            residues: List of residue numbers
            conservation: Conservation scores

        Returns:
            Conservation score (0-1)
        """
        try:
            if not residues or not conservation:
                return 0.0

            scores = [conservation.get(res, 0.0) for res in residues]

            return float(np.mean(scores))

        except Exception as e:
            self.logger.error(f"Error scoring conservation: {str(e)}")
            return 0.0

    def _score_druggability(self, pocket: Dict[str, Any]) -> float:
        """Score pocket druggability.

        Args:
            pocket: Pocket dictionary

        Returns:
            Druggability score (0-1)
        """
        try:
            # Volume score
            volume_score = self._score_volume(pocket["volume"])

            # Shape score
            shape_score = self._score_shape(pocket["shape_descriptors"])

            # Depth score
            depth = pocket.get("depth", 0.0)
            depth_score = min(depth / 10.0, 1.0)

            # Exposure score
            exposure = pocket.get("exposure", 1.0)
            exposure_score = 1.0 - exposure

            # Hydrophobicity score
            hydrophobicity = pocket.get("hydrophobicity", 0.0)
            hydro_score = (hydrophobicity + 4.5) / 9.0  # Normalize to 0-1

            # Combine scores
            druggability = 0.3 * volume_score + 0.2 * shape_score + 0.2 * depth_score + 0.15 * exposure_score + 0.15 * hydro_score

            return float(druggability)

        except Exception as e:
            self.logger.error(f"Error scoring druggability: {str(e)}")
            return 0.0
