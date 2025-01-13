"""Binding site surface analysis functionality."""

import logging
import numpy as np
from typing import Dict, List, Optional, Any, Set
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from Bio.PDB.SASA import ShrakeRupley

from .base import BaseSurfaceAnalyzer

logger = logging.getLogger(__name__)


class BindingSurfaceAnalyzer(BaseSurfaceAnalyzer):
    """Analyzes binding site surface properties."""

    def analyze_binding_site(
        self,
        structure: Structure,
        site_residues: List[int],
        probe_radius: float = 1.4,
    ) -> Dict[str, Any]:
        """Analyze binding site surface properties.

        Args:
            structure: BioPython Structure object
            site_residues: List of residue numbers in binding site
            probe_radius: Probe radius for surface calculation (Å)

        Returns:
            Dictionary of binding site properties
        """
        try:
            # Get binding site atoms
            site_atoms = self._get_site_atoms(structure, site_residues)
            if not site_atoms:
                return {}

            # Calculate surface properties
            surface_atoms = self.get_surface_atoms(structure, probe_radius)
            site_surface_atoms = site_atoms & surface_atoms

            # Basic properties
            properties = {
                "total_atoms": len(site_atoms),
                "surface_atoms": len(site_surface_atoms),
                "total_area": self.calculate_surface_area(site_atoms),
                "surface_area": self.calculate_surface_area(site_surface_atoms),
                "volume": self._calculate_pocket_volume(site_atoms),
                "depth": self._calculate_pocket_depth(site_atoms, surface_atoms),
            }

            # Get residue properties
            site_residues = self.get_surface_residues(site_surface_atoms)
            properties.update(
                {
                    "hydrophobicity": self.calculate_hydrophobicity(site_residues),
                    "charge": self.calculate_charge(site_residues),
                    "residue_composition": self._get_residue_composition(site_residues),
                    "secondary_structure": self._analyze_pocket_secondary_structure(site_residues),
                    "conservation": self._analyze_pocket_conservation(site_residues),
                    "flexibility": self._analyze_pocket_flexibility(site_residues),
                }
            )

            # Analyze pocket shape
            properties.update(self._analyze_pocket_shape(site_atoms))

            # Find sub-pockets
            properties["sub_pockets"] = self._find_sub_pockets(site_atoms, probe_radius)

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing binding site: {str(e)}")
            return {}

    def _calculate_pocket_volume(self, atoms: Set[Atom]) -> float:
        """Calculate binding pocket volume.

        Args:
            atoms: Set of pocket atoms

        Returns:
            Pocket volume in Å³
        """
        try:
            if not atoms:
                return 0.0

            coords = np.array([atom.get_coord() for atom in atoms])
            hull = ConvexHull(coords)
            return hull.volume

        except Exception as e:
            self.logger.error(f"Error calculating pocket volume: {str(e)}")
            return 0.0

    def _calculate_pocket_depth(self, pocket_atoms: Set[Atom], surface_atoms: Set[Atom]) -> float:
        """Calculate binding pocket depth.

        Args:
            pocket_atoms: Set of pocket atoms
            surface_atoms: Set of surface atoms

        Returns:
            Pocket depth in Å
        """
        try:
            if not pocket_atoms or not surface_atoms:
                return 0.0

            # Get pocket center
            pocket_coords = np.array([atom.get_coord() for atom in pocket_atoms])
            pocket_center = np.mean(pocket_coords, axis=0)

            # Get surface points
            surface_coords = np.array([atom.get_coord() for atom in surface_atoms])

            # Calculate minimum distance from pocket center to surface
            distances = np.linalg.norm(surface_coords - pocket_center, axis=1)
            return float(np.min(distances))

        except Exception as e:
            self.logger.error(f"Error calculating pocket depth: {str(e)}")
            return 0.0

    def _get_residue_composition(self, residues: List[Residue]) -> Dict[str, int]:
        """Get amino acid composition of residues.

        Args:
            residues: List of residues

        Returns:
            Dictionary mapping residue names to counts
        """
        try:
            composition = {}
            for residue in residues:
                name = residue.get_resname()
                composition[name] = composition.get(name, 0) + 1
            return composition

        except Exception as e:
            self.logger.error(f"Error getting residue composition: {str(e)}")
            return {}

    def _analyze_pocket_secondary_structure(self, residues: List[Residue]) -> Dict[str, float]:
        """Analyze secondary structure composition of pocket.

        Args:
            residues: List of pocket residues

        Returns:
            Dictionary of secondary structure proportions
        """
        try:
            ss_counts = {"H": 0, "B": 0, "E": 0, "G": 0, "I": 0, "T": 0, "S": 0}
            total = 0

            for residue in residues:
                # Get DSSP assignment
                dssp = residue.xtra.get("SS_DSSP")
                if dssp and dssp in ss_counts:
                    ss_counts[dssp] += 1
                    total += 1

            if total > 0:
                return {k: v / total for k, v in ss_counts.items()}
            return ss_counts

        except Exception as e:
            self.logger.error(f"Error analyzing pocket secondary structure: {str(e)}")
            return {}

    def _analyze_pocket_conservation(self, residues: List[Residue]) -> float:
        """Calculate average conservation score of pocket residues.

        Args:
            residues: List of pocket residues

        Returns:
            Average conservation score (0-1)
        """
        try:
            scores = []
            for residue in residues:
                score = residue.xtra.get("CONSERVATION")
                if score is not None:
                    scores.append(score)

            return float(np.mean(scores)) if scores else 0.0

        except Exception as e:
            self.logger.error(f"Error analyzing pocket conservation: {str(e)}")
            return 0.0

    def _analyze_pocket_flexibility(self, residues: List[Residue]) -> float:
        """Calculate average flexibility of pocket residues from B-factors.

        Args:
            residues: List of pocket residues

        Returns:
            Average B-factor
        """
        try:
            b_factors = []
            for residue in residues:
                for atom in residue:
                    if atom.bfactor is not None:
                        b_factors.append(atom.bfactor)

            return float(np.mean(b_factors)) if b_factors else 0.0

        except Exception as e:
            self.logger.error(f"Error analyzing pocket flexibility: {str(e)}")
            return 0.0

    def _analyze_pocket_shape(self, atoms: Set[Atom]) -> Dict[str, float]:
        """Analyze geometric properties of binding pocket.

        Args:
            atoms: Set of pocket atoms

        Returns:
            Dictionary of shape properties
        """
        try:
            if not atoms:
                return {}

            coords = np.array([atom.get_coord() for atom in atoms])

            # Calculate principal components
            centered_coords = coords - np.mean(coords, axis=0)
            cov = np.cov(centered_coords.T)
            eigenvals, eigenvecs = np.linalg.eigh(cov)
            eigenvals = eigenvals[::-1]
            eigenvecs = eigenvecs[:, ::-1]

            # Calculate shape descriptors
            properties = {
                "sphericity": float(min(eigenvals) / max(eigenvals)),
                "elongation": float(max(eigenvals) / np.median(eigenvals)),
                "flatness": float(min(eigenvals) / np.median(eigenvals)),
            }

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing pocket shape: {str(e)}")
            return {}

    def _find_sub_pockets(
        self,
        atoms: Set[Atom],
        probe_radius: float = 1.4,
    ) -> List[Dict[str, Any]]:
        """Find sub-pockets within binding site.

        Args:
            atoms: Set of binding site atoms
            probe_radius: Probe radius for cavity detection

        Returns:
            List of sub-pocket properties
        """
        try:
            if not atoms:
                return []

            # Get atomic coordinates
            coords = np.array([atom.get_coord() for atom in atoms])

            # Use clustering to identify sub-pockets
            from sklearn.cluster import DBSCAN

            clustering = DBSCAN(eps=probe_radius * 2, min_samples=4).fit(coords)
            labels = clustering.labels_

            # Analyze each sub-pocket
            sub_pockets = []
            for label in set(labels):
                if label == -1:  # Skip noise points
                    continue

                # Get sub-pocket atoms
                mask = labels == label
                sub_pocket_coords = coords[mask]
                sub_pocket_atoms = set(list(atoms)[i] for i in np.where(mask)[0])

                properties = {
                    "center": np.mean(sub_pocket_coords, axis=0).tolist(),
                    "volume": self._calculate_pocket_volume(sub_pocket_atoms),
                    "n_atoms": len(sub_pocket_atoms),
                    "residues": self.get_surface_residues(sub_pocket_atoms),
                }
                sub_pockets.append(properties)

            return sub_pockets

        except Exception as e:
            self.logger.error(f"Error finding sub-pockets: {str(e)}")
            return []
