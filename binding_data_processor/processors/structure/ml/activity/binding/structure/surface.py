"""Surface analysis functionality for protein structures."""

import logging
from typing import Dict, List, Optional, Any, Set, Tuple
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_angle, calc_dihedral
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)


class SurfaceAnalyzer:
    """Analyzes protein surface properties."""

    # Residue properties
    HYDROPHOBICITY = {
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

    CHARGE = {
        "ARG": 1,
        "LYS": 1,
        "ASP": -1,
        "GLU": -1,
        "HIS": 0.1,
    }

    def __init__(self):
        """Initialize surface analyzer."""
        self.logger = logging.getLogger(__name__)

    def analyze_surface(
        self,
        structure: Structure,
        probe_radius: float = 1.4,
        include_electrostatics: bool = True,
        include_hydrophobicity: bool = True,
    ) -> Dict[str, Any]:
        """Analyze protein surface properties.

        Args:
            structure: BioPython Structure object
            probe_radius: Probe radius for surface calculation (Å)
            include_electrostatics: Whether to calculate electrostatic properties
            include_hydrophobicity: Whether to calculate hydrophobic properties

        Returns:
            Dictionary of surface properties
        """
        try:
            # Get surface atoms and calculate basic properties
            surface_atoms = self._get_surface_atoms(structure, probe_radius)
            surface_area = self._calculate_surface_area(surface_atoms)
            surface_residues = self._get_surface_residues(surface_atoms)

            properties = {
                "surface_atoms": surface_atoms,
                "surface_area": surface_area,
                "surface_residues": surface_residues,
                "relative_surface_area": self._calculate_relative_surface_area(structure, surface_area),
                "surface_roughness": self._calculate_surface_roughness(surface_atoms),
                "cavities": self._find_surface_cavities(surface_atoms, probe_radius),
            }

            # Calculate electrostatic properties
            if include_electrostatics:
                electrostatics = self._analyze_electrostatics(surface_residues)
                properties.update(electrostatics)

            # Calculate hydrophobic properties
            if include_hydrophobicity:
                hydrophobicity = self._analyze_hydrophobicity(surface_residues)
                properties.update(hydrophobicity)

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing surface: {str(e)}")
            return {}

    def analyze_site(
        self,
        structure: Structure,
        site_residues: List[int],
    ) -> Dict[str, Any]:
        """Analyze surface properties of a specific site.

        Args:
            structure: BioPython Structure object
            site_residues: List of residue numbers in site

        Returns:
            Dictionary of site surface properties
        """
        try:
            # Get site atoms and surface properties
            site_atoms = self._get_site_atoms(structure, site_residues)
            surface_atoms = self._get_surface_atoms(structure)
            site_surface_atoms = set(site_atoms) & set(surface_atoms)

            # Calculate basic properties
            surface_area = self._calculate_surface_area(site_surface_atoms)
            total_area = self._calculate_surface_area(site_atoms)

            properties = {
                "surface_area": surface_area,
                "total_area": total_area,
                "exposure": float(surface_area / total_area if total_area > 0 else 0.0),
                "roughness": self._calculate_surface_roughness(site_surface_atoms),
            }

            # Calculate electrostatic properties
            site_residues = self._get_surface_residues(site_surface_atoms)
            electrostatics = self._analyze_electrostatics(site_residues)
            properties.update({f"site_{k}": v for k, v in electrostatics.items()})

            # Calculate hydrophobic properties
            hydrophobicity = self._analyze_hydrophobicity(site_residues)
            properties.update({f"site_{k}": v for k, v in hydrophobicity.items()})

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing site surface: {str(e)}")
            return {}

    def _get_surface_atoms(
        self,
        structure: Structure,
        probe_radius: float = 1.4,
    ) -> Set["Atom"]:
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
                distances = cdist(probe_coords, coords)
                if np.any(np.all(distances > probe_radius, axis=1)):
                    surface_atoms.add(atom)

            return surface_atoms

        except Exception as e:
            self.logger.error(f"Error getting surface atoms: {str(e)}")
            return set()

    def _calculate_surface_area(self, atoms: Set["Atom"]) -> float:
        """Calculate surface area from surface atoms.

        Args:
            atoms: Set of atoms

        Returns:
            Surface area in Å²
        """
        try:
            if not atoms:
                return 0.0

            # Get coordinates
            coords = np.array([atom.get_coord() for atom in atoms])

            if len(coords) < 4:
                return 0.0

            # Calculate surface area using convex hull
            hull = ConvexHull(coords)
            return hull.area

        except Exception as e:
            self.logger.error(f"Error calculating surface area: {str(e)}")
            return 0.0

    def _calculate_relative_surface_area(
        self,
        structure: Structure,
        surface_area: float,
    ) -> float:
        """Calculate relative surface area compared to sphere.

        Args:
            structure: BioPython Structure object
            surface_area: Calculated surface area

        Returns:
            Relative surface area (0-1)
        """
        try:
            # Get radius of equivalent sphere
            coords = np.array([atom.get_coord() for atom in structure.get_atoms()])
            center = np.mean(coords, axis=0)
            distances = np.linalg.norm(coords - center, axis=1)
            radius = np.mean(distances)

            # Calculate sphere surface area
            sphere_area = 4 * np.pi * radius * radius

            return float(surface_area / sphere_area if sphere_area > 0 else 0.0)

        except Exception as e:
            self.logger.error(f"Error calculating relative surface area: {str(e)}")
            return 0.0

    def _calculate_surface_roughness(self, atoms: Set["Atom"]) -> float:
        """Calculate surface roughness using fractal dimension.

        Args:
            atoms: Set of surface atoms

        Returns:
            Surface roughness score (higher means rougher)
        """
        try:
            if not atoms:
                return 0.0

            coords = np.array([atom.get_coord() for atom in atoms])
            if len(coords) < 4:
                return 0.0

            # Calculate fractal dimension using box counting
            scales = np.logspace(0, 2, 20)
            counts = []

            for scale in scales:
                # Create grid
                mins = np.min(coords, axis=0)
                maxs = np.max(coords, axis=0)
                bins = np.ceil((maxs - mins) / scale).astype(int)

                # Count occupied boxes
                H, _ = np.histogramdd(coords, bins=bins)
                counts.append(np.sum(H > 0))

            # Calculate fractal dimension from slope
            coeffs = np.polyfit(np.log(scales), np.log(counts), 1)
            fractal_dim = -coeffs[0]

            return float(fractal_dim)

        except Exception as e:
            self.logger.error(f"Error calculating surface roughness: {str(e)}")
            return 0.0

    def _find_surface_cavities(
        self,
        surface_atoms: Set["Atom"],
        probe_radius: float = 1.4,
    ) -> List[Dict[str, Any]]:
        """Find cavities in protein surface.

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

    def _analyze_electrostatics(
        self,
        residues: List[Residue],
    ) -> Dict[str, Any]:
        """Analyze electrostatic properties of surface residues.

        Args:
            residues: List of residues

        Returns:
            Dictionary of electrostatic properties
        """
        try:
            if not residues:
                return {}

            # Count charged residues
            charges = [self.CHARGE.get(res.get_resname(), 0.0) for res in residues]
            positive = sum(1 for c in charges if c > 0)
            negative = sum(1 for c in charges if c < 0)

            return {
                "positive_charges": positive,
                "negative_charges": negative,
                "net_charge": positive - negative,
                "charge_density": float((positive + negative) / len(residues) if residues else 0.0),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing electrostatics: {str(e)}")
            return {}

    def _analyze_hydrophobicity(
        self,
        residues: List[Residue],
    ) -> Dict[str, Any]:
        """Analyze hydrophobic properties of surface residues.

        Args:
            residues: List of residues

        Returns:
            Dictionary of hydrophobic properties
        """
        try:
            if not residues:
                return {}

            # Calculate hydrophobicity scores
            scores = [self.HYDROPHOBICITY.get(res.get_resname(), 0.0) for res in residues]

            # Identify hydrophobic patches
            hydrophobic_residues = [res for i, res in enumerate(residues) if scores[i] > 0]
            patches = self._find_hydrophobic_patches(hydrophobic_residues)

            return {
                "average_hydrophobicity": float(np.mean(scores)),
                "hydrophobic_residues": len(hydrophobic_residues),
                "hydrophobic_ratio": float(len(hydrophobic_residues) / len(residues) if residues else 0.0),
                "hydrophobic_patches": len(patches),
                "largest_patch_size": max((len(patch) for patch in patches), default=0),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing hydrophobicity: {str(e)}")
            return {}

    def _find_hydrophobic_patches(
        self,
        hydrophobic_residues: List[Residue],
        distance_cutoff: float = 8.0,
    ) -> List[List[Residue]]:
        """Find connected hydrophobic patches on surface.

        Args:
            hydrophobic_residues: List of hydrophobic residues
            distance_cutoff: Maximum distance between residues in patch

        Returns:
            List of residue lists representing patches
        """
        try:
            if not hydrophobic_residues:
                return []

            # Build contact network
            contacts = {}
            for i, res1 in enumerate(hydrophobic_residues):
                contacts[i] = set()
                for j, res2 in enumerate(hydrophobic_residues[i + 1 :], i + 1):
                    if self._check_residue_contact(res1, res2, distance_cutoff):
                        contacts[i].add(j)
                        if j not in contacts:
                            contacts[j] = set()
                        contacts[j].add(i)

            # Find connected components (patches)
            patches = []
            unvisited = set(range(len(hydrophobic_residues)))

            while unvisited:
                # Start new patch
                current = unvisited.pop()
                patch = {current}
                stack = [current]

                # Grow patch
                while stack:
                    node = stack.pop()
                    for neighbor in contacts[node]:
                        if neighbor in unvisited:
                            patch.add(neighbor)
                            stack.append(neighbor)
                            unvisited.remove(neighbor)

                patches.append([hydrophobic_residues[i] for i in patch])

            return patches

        except Exception as e:
            self.logger.error(f"Error finding hydrophobic patches: {str(e)}")
            return []

    def _check_residue_contact(
        self,
        res1: Residue,
        res2: Residue,
        cutoff: float = 8.0,
    ) -> bool:
        """Check if two residues are in contact.

        Args:
            res1: First residue
            res2: Second residue
            cutoff: Distance cutoff in Å

        Returns:
            True if residues are in contact
        """
        try:
            # Check CA distance first
            if "CA" in res1 and "CA" in res2:
                ca_dist = res1["CA"] - res2["CA"]
                if ca_dist > cutoff:
                    return False

            # Check all atom pairs
            for atom1 in res1:
                for atom2 in res2:
                    if atom1 - atom2 < cutoff:
                        return True

            return False

        except Exception as e:
            self.logger.error(f"Error checking residue contact: {str(e)}")
            return False

    def _get_site_atoms(
        self,
        structure: Structure,
        site_residues: List[int],
    ) -> Set["Atom"]:
        """Get atoms belonging to site residues.

        Args:
            structure: BioPython Structure object
            site_residues: List of residue numbers

        Returns:
            Set of atoms in site
        """
        try:
            site_atoms = set()
            for residue in structure.get_residues():
                if residue.get_id()[1] in site_residues:
                    site_atoms.update(residue.get_atoms())
            return site_atoms

        except Exception as e:
            self.logger.error(f"Error getting site atoms: {str(e)}")
            return set()

    def _get_surface_residues(self, surface_atoms: Set["Atom"]) -> List[Residue]:
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
