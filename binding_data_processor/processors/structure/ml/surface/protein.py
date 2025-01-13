"""Protein structure analysis functionality."""

import logging
import numpy as np
from typing import Dict, List, Optional, Set, Any
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from scipy.spatial import ConvexHull, Voronoi

from .base import BaseSurfaceAnalyzer

logger = logging.getLogger(__name__)


class ProteinStructureAnalyzer(BaseSurfaceAnalyzer):
    """Analyzes protein structure properties."""

    # Secondary structure types
    SS_TYPES = {
        "H": "alpha helix",
        "B": "beta bridge",
        "E": "beta sheet",
        "G": "3-10 helix",
        "I": "pi helix",
        "T": "turn",
        "S": "bend",
    }

    # Reference SASA values from Miller et al. (1987)
    REFERENCE_SASA = {
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
    }

    def __init__(self, **kwargs):
        """Initialize protein structure analyzer."""
        super().__init__()
        self.logger = logging.getLogger(self.__class__.__name__)

    def _get_site_atoms(self, structure: Structure, site_residues: List[int]) -> Set[Atom]:
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

    def identify_exposed_residues(self, structure: Structure, cutoff: float = 2.8) -> List[int]:
        """Identify solvent-exposed residues using SASA.

        Args:
            structure: BioPython Structure object
            cutoff: SASA cutoff for exposure (Å²)

        Returns:
            List of exposed residue numbers
        """
        try:
            exposed = []

            # Calculate SASA per residue
            for residue in structure.get_residues():
                coords = []
                radii = []
                for atom in residue:
                    coords.append(atom.get_coord())
                    radii.append(self.ATOM_RADII.get(atom.element, 1.5))

                if coords:
                    coords = np.array(coords)
                    radii = np.array(radii)

                    from Bio.PDB.SASA import ShrakeRupley

                    sr = ShrakeRupley()
                    sasa = sr.compute(coords, radii)

                    # Check if exposed
                    if sasa > cutoff:
                        exposed.append(residue.get_id()[1])

            return exposed

        except Exception as e:
            self.logger.error(f"Error identifying exposed residues: {str(e)}")
            return []

    def _calculate_surface_hydrophobicity(self, structure: Structure) -> float:
        """Calculate surface hydrophobicity using exposed residues.

        Args:
            structure: BioPython Structure object

        Returns:
            Average hydrophobicity of exposed residues
        """
        try:
            exposed = self.identify_exposed_residues(structure)
            if not exposed:
                return 0.0

            scores = []
            for residue in structure.get_residues():
                if residue.get_id()[1] in exposed:
                    scores.append(self.HYDROPHOBICITY.get(residue.get_resname(), 0.0))

            return float(np.mean(scores)) if scores else 0.0

        except Exception as e:
            self.logger.error(f"Error calculating surface hydrophobicity: {str(e)}")
            return 0.0

    def _calculate_surface_charge(self, structure: Structure) -> float:
        """Calculate net surface charge using exposed residues.

        Args:
            structure: BioPython Structure object

        Returns:
            Net charge of exposed residues
        """
        try:
            exposed = self.identify_exposed_residues(structure)
            if not exposed:
                return 0.0

            total_charge = 0.0
            for residue in structure.get_residues():
                if residue.get_id()[1] in exposed:
                    total_charge += self.CHARGE.get(residue.get_resname(), 0.0)

            return float(total_charge)

        except Exception as e:
            self.logger.error(f"Error calculating surface charge: {str(e)}")
            return 0.0

    def _calculate_surface_polarity(self, structure: Structure) -> float:
        """Calculate surface polarity using exposed residues.

        Args:
            structure: BioPython Structure object

        Returns:
            Ratio of polar to non-polar exposed surface area
        """
        try:
            exposed = self.identify_exposed_residues(structure)
            if not exposed:
                return 0.0

            polar_area = 0.0
            nonpolar_area = 0.0

            for residue in structure.get_residues():
                if residue.get_id()[1] in exposed:
                    # Calculate SASA for residue
                    coords = []
                    radii = []
                    for atom in residue:
                        coords.append(atom.get_coord())
                        radii.append(self.ATOM_RADII.get(atom.element, 1.5))

                    if coords:
                        coords = np.array(coords)
                        radii = np.array(radii)

                        from Bio.PDB.SASA import ShrakeRupley

                        sr = ShrakeRupley()
                        sasa = sr.compute(coords, radii)

                        # Classify residue polarity
                        if residue.get_resname() in ["ARG", "LYS", "ASP", "GLU", "ASN", "GLN", "HIS", "SER", "THR", "TYR"]:
                            polar_area += sasa
                        else:
                            nonpolar_area += sasa

            return float(polar_area / nonpolar_area if nonpolar_area > 0 else 0.0)

        except Exception as e:
            self.logger.error(f"Error calculating surface polarity: {str(e)}")
            return 0.0

    def analyze_surface(
        self, structure: Structure, probe_radius: float = 1.4, include_electrostatics: bool = True, include_hydrophobicity: bool = True, include_secondary_structure: bool = True
    ) -> Dict[str, Any]:
        """Analyze protein surface comprehensively.

        Args:
            structure: BioPython Structure object
            probe_radius: Probe radius for surface calculation (Å)
            include_electrostatics: Whether to calculate electrostatic properties
            include_hydrophobicity: Whether to calculate hydrophobic properties
            include_secondary_structure: Whether to analyze secondary structure

        Returns:
            Dictionary containing surface analysis results
        """
        try:
            # Get surface atoms and calculate basic properties
            surface_atoms = self.get_surface_atoms(structure, probe_radius)
            surface_residues = self.get_surface_residues(surface_atoms)
            surface_area = self.calculate_surface_area(surface_atoms)

            # Calculate basic properties
            properties = {
                "surface_area": surface_area,
                "surface_residues": len(surface_residues),
                "relative_surface_area": self._calculate_relative_surface_area(structure, surface_area),
                "surface_roughness": self._calculate_surface_roughness(surface_atoms),
                "surface_composition": self._analyze_surface_composition(surface_residues),
                "exposed_residues": self.identify_exposed_residues(structure),
            }

            # Calculate electrostatic properties
            if include_electrostatics:
                electrostatics = self._analyze_electrostatics(surface_residues)
                properties.update(electrostatics)
                properties["electrostatic_patches"] = self._find_charged_patches(structure, surface_residues)
                properties["surface_charge"] = self._calculate_surface_charge(structure)

            # Calculate hydrophobic properties
            if include_hydrophobicity:
                hydrophobicity = self._analyze_hydrophobicity(surface_residues)
                properties.update(hydrophobicity)
                properties["hydrophobic_patches"] = self._find_hydrophobic_patches(structure, surface_residues)
                properties["surface_hydrophobicity"] = self._calculate_surface_hydrophobicity(structure)
                properties["surface_polarity"] = self._calculate_surface_polarity(structure)

            # Analyze secondary structure
            if include_secondary_structure:
                properties["secondary_structure"] = self._analyze_surface_secondary_structure(structure, surface_residues)

            # Find and analyze cavities
            cavities = self.find_surface_cavities(surface_atoms, probe_radius)
            if cavities:
                properties["cavities"] = [
                    {
                        "volume": cavity["volume"],
                        "depth": self.calculate_cavity_depth(cavity["center"], np.array([a.get_coord() for a in surface_atoms])),
                        "surrounding_residues": len(self.get_surface_residues(set(cavity["surrounding_atoms"]))),
                    }
                    for cavity in cavities
                ]

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing surface: {str(e)}")
            return {}

    def analyze_site(self, structure: Structure, site_residues: List[int], probe_radius: float = 1.4) -> Dict[str, Any]:
        """Analyze surface properties of a specific site.

        Args:
            structure: BioPython Structure object
            site_residues: List of residue numbers in site
            probe_radius: Probe radius for surface calculation (Å)

        Returns:
            Dictionary containing site analysis results
        """
        try:
            # Get site atoms and surface properties
            site_atoms = self._get_site_atoms(structure, site_residues)
            if not site_atoms:
                return {}

            surface_atoms = self.get_surface_atoms(structure, probe_radius)
            site_surface_atoms = site_atoms.intersection(surface_atoms)
            site_residues = self.get_surface_residues(site_surface_atoms)

            # Calculate basic properties
            surface_area = self.calculate_surface_area(site_surface_atoms)
            total_area = self.calculate_surface_area(site_atoms)
            volume = self.calculate_volume(np.array([a.get_coord() for a in site_atoms]))

            properties = {
                "surface_area": surface_area,
                "total_area": total_area,
                "volume": volume,
                "exposure": float(surface_area / total_area if total_area > 0 else 0.0),
                "roughness": self._calculate_surface_roughness(site_surface_atoms),
                "residue_composition": self._analyze_surface_composition(site_residues),
                "secondary_structure": self._analyze_surface_secondary_structure(structure, site_residues),
                "relative_exposure": self._calculate_site_exposure(structure, site_residues),
                "conservation": self._analyze_site_conservation(structure, site_residues),
            }

            # Calculate electrostatic properties
            electrostatics = self._analyze_electrostatics(site_residues)
            properties.update({f"site_{k}": v for k, v in electrostatics.items()})

            # Calculate hydrophobic properties
            hydrophobicity = self._analyze_hydrophobicity(site_residues)
            properties.update({f"site_{k}": v for k, v in hydrophobicity.items()})

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing site surface: {str(e)}")
            return {}

    def _calculate_relative_surface_area(self, structure: Structure, surface_area: float) -> float:
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

    def _calculate_surface_roughness(self, atoms: Set[Atom]) -> float:
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

    def _analyze_surface_composition(self, residues: List[Residue]) -> Dict[str, float]:
        """Analyze amino acid composition of surface residues.

        Args:
            residues: List of residues

        Returns:
            Dictionary mapping residue types to their frequencies
        """
        try:
            composition = {}
            total = len(residues)
            if total == 0:
                return {}

            for res in residues:
                resname = res.get_resname()
                composition[resname] = composition.get(resname, 0) + 1

            # Convert to frequencies
            return {k: v / total for k, v in composition.items()}

        except Exception as e:
            self.logger.error(f"Error analyzing surface composition: {str(e)}")
            return {}

    def _analyze_electrostatics(self, residues: List[Residue]) -> Dict[str, Any]:
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

    def _analyze_hydrophobicity(self, residues: List[Residue]) -> Dict[str, Any]:
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

    def _find_hydrophobic_patches(self, structure: Structure, surface_residues: List[Residue], min_size: int = 3, max_distance: float = 5.0) -> List[Dict[str, Any]]:
        """Find continuous hydrophobic patches on surface.

        Args:
            structure: BioPython Structure object
            surface_residues: List of surface residues
            min_size: Minimum number of residues in patch
            max_distance: Maximum distance between residues in patch

        Returns:
            List of hydrophobic patch properties
        """
        try:
            patches = []
            visited = set()

            # Get hydrophobic residues
            hydrophobic = [r for r in surface_residues if self.HYDROPHOBICITY.get(r.get_resname(), 0) > 0]

            for res in hydrophobic:
                if res.get_id()[1] in visited:
                    continue

                # Find connected hydrophobic residues
                patch = self._grow_patch(res, hydrophobic, max_distance)
                if len(patch) >= min_size:
                    patch_coords = np.array([atom.get_coord() for r in patch for atom in r if atom.element != "H"])

                    patches.append(
                        {
                            "residues": [r.get_id()[1] for r in patch],
                            "size": len(patch),
                            "area": self.calculate_surface_area(set(atom for r in patch for atom in r)),
                            "center": np.mean(patch_coords, axis=0),
                            "hydrophobicity": np.mean([self.HYDROPHOBICITY.get(r.get_resname(), 0) for r in patch]),
                        }
                    )

                visited.update(r.get_id()[1] for r in patch)

            return patches

        except Exception as e:
            self.logger.error(f"Error finding hydrophobic patches: {str(e)}")
            return []

    def _find_charged_patches(self, structure: Structure, surface_residues: List[Residue], min_size: int = 3, max_distance: float = 5.0) -> List[Dict[str, Any]]:
        """Find continuous charged patches on surface.

        Args:
            structure: BioPython Structure object
            surface_residues: List of surface residues
            min_size: Minimum number of residues in patch
            max_distance: Maximum distance between residues in patch

        Returns:
            List of charged patch properties
        """
        try:
            patches = []
            visited = set()

            # Get charged residues
            charged = [r for r in surface_residues if abs(self.CHARGE.get(r.get_resname(), 0)) > 0]

            for res in charged:
                if res.get_id()[1] in visited:
                    continue

                # Find connected charged residues
                patch = self._grow_patch(res, charged, max_distance)
                if len(patch) >= min_size:
                    patch_coords = np.array([atom.get_coord() for r in patch for atom in r if atom.element != "H"])

                    patches.append(
                        {
                            "residues": [r.get_id()[1] for r in patch],
                            "size": len(patch),
                            "area": self.calculate_surface_area(set(atom for r in patch for atom in r)),
                            "center": np.mean(patch_coords, axis=0),
                            "net_charge": sum(self.CHARGE.get(r.get_resname(), 0) for r in patch),
                        }
                    )

                visited.update(r.get_id()[1] for r in patch)

            return patches

        except Exception as e:
            self.logger.error(f"Error finding charged patches: {str(e)}")
            return []

    def _grow_patch(self, start_residue: Residue, candidate_residues: List[Residue], max_distance: float) -> List[Residue]:
        """Grow a continuous patch of residues from a starting residue.

        Args:
            start_residue: Starting residue
            candidate_residues: List of residues to consider
            max_distance: Maximum distance between residues

        Returns:
            List of residues in patch
        """
        try:
            patch = [start_residue]
            to_check = [start_residue]

            while to_check:
                current = to_check.pop(0)
                current_ca = current["CA"].get_coord()

                for res in candidate_residues:
                    if res not in patch and "CA" in res:
                        dist = np.linalg.norm(current_ca - res["CA"].get_coord())
                        if dist <= max_distance:
                            patch.append(res)
                            to_check.append(res)

            return patch

        except Exception as e:
            self.logger.error(f"Error growing patch: {str(e)}")
            return [start_residue]

    def _analyze_surface_secondary_structure(self, structure: Structure, surface_residues: List[Residue]) -> Dict[str, float]:
        """Analyze secondary structure composition of surface residues.

        Args:
            structure: BioPython Structure object
            surface_residues: List of surface residues

        Returns:
            Dictionary mapping SS types to their frequencies
        """
        try:
            ss_counts = {ss: 0 for ss in self.SS_TYPES}
            surface_ids = {(res.get_parent().id, res.get_id()[1]) for res in surface_residues}

            dssp = dssp_dict_from_pdb_file(structure.id)[0]
            total = 0

            for key in dssp:
                if key[:2] in surface_ids:  # Chain ID and residue number
                    ss = dssp[key][2]  # Secondary structure
                    if ss in ss_counts:
                        ss_counts[ss] += 1
                        total += 1

            return {k: v / total for k, v in ss_counts.items()} if total > 0 else ss_counts

        except Exception as e:
            self.logger.error(f"Error analyzing surface secondary structure: {str(e)}")
            return {ss: 0.0 for ss in self.SS_TYPES}

    def _calculate_site_exposure(self, structure: Structure, site_residues: List[Residue]) -> float:
        """Calculate relative solvent exposure of site residues.

        Args:
            structure: BioPython Structure object
            site_residues: List of site residues

        Returns:
            Average relative exposure (0-1)
        """
        try:
            exposures = []
            for res in site_residues:
                # Calculate SASA for residue
                coords = []
                radii = []
                for atom in res:
                    coords.append(atom.get_coord())
                    radii.append(1.4 + self.ATOM_RADII.get(atom.element, 1.5))

                if coords:
                    coords = np.array(coords)
                    radii = np.array(radii)

                    from Bio.PDB.SASA import ShrakeRupley

                    sr = ShrakeRupley()
                    sasa = sr.compute(coords, radii)

                    # Get reference SASA for fully exposed residue
                    ref_sasa = self.REFERENCE_SASA.get(res.get_resname(), 0.0)
                    if ref_sasa > 0:
                        exposures.append(sasa / ref_sasa)

            return float(np.mean(exposures)) if exposures else 0.0

        except Exception as e:
            self.logger.error(f"Error calculating site exposure: {str(e)}")
            return 0.0

    def _analyze_site_conservation(self, structure: Structure, site_residues: List[Residue]) -> Dict[str, Any]:
        """Analyze evolutionary conservation of site residues.

        Args:
            structure: BioPython Structure object
            site_residues: List of site residues

        Returns:
            Dictionary of conservation metrics
        """
        try:
            # Placeholder - would need sequence alignment data
            # This could be implemented by:
            # 1. Loading pre-computed conservation scores
            # 2. Performing MSA and calculating conservation
            # 3. Using external conservation prediction services
            return {
                "average_conservation": 0.0,
                "highly_conserved": [],
                "variable": [],
            }

        except Exception as e:
            self.logger.error(f"Error analyzing site conservation: {str(e)}")
            return {}
