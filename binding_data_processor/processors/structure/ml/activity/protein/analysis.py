"""Protein structure analysis module."""

import logging
import numpy as np
from Bio.PDB import Structure, Chain, Residue, Atom, Vector, calc_dihedral
from typing import Dict, List, Optional, Tuple, Union, Any

from binding_data_processor.core.config import ProteinAnalysisConfig
from .utils import (
    get_atom_radius,
    calculate_center_of_mass,
    calculate_radius_of_gyration,
    get_residue_property,
    calculate_surface_area,
    calculate_surface_exposure,
    calculate_residue_depth,
    get_secondary_structure,
    calculate_interface_area,
    calculate_cavity_volume,
    get_atom_coords,
    get_backbone_atoms,
    get_sequence_from_structure,
    find_contacts,
)

logger = logging.getLogger(__name__)

# Property scales
HYDROPHOBICITY_SCALE = {
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

CHARGE_SCALE = {
    "ARG": 1,
    "LYS": 1,
    "ASP": -1,
    "GLU": -1,
    "HIS": 0.5,  # Can be protonated
}


class ProteinStructureAnalyzer:
    """Analyzes protein structures and their properties."""

    def __init__(self, config: Optional[ProteinAnalysisConfig] = None, structure: Optional[Structure] = None):
        """Initialize protein structure analyzer.

        Args:
            config: Configuration for analysis
            structure: Optional BioPython Structure object to analyze
        """
        self.config = config or ProteinAnalysisConfig()
        self.logger = logging.getLogger(self.__class__.__name__)

        if structure is not None:
            self.set_structure(structure)

    def set_structure(self, structure: Structure) -> None:
        """Set structure to analyze.

        Args:
            structure: BioPython Structure object
        """
        self.structure = structure
        self.coords, self.masses = get_atom_coords(structure)
        self.sequence = get_sequence_from_structure(structure)

    def get_structure_properties(self, structure: Optional[Structure] = None) -> Dict[str, Any]:
        """Get comprehensive properties of protein structure.

        Args:
            structure: Optional BioPython Structure object (uses stored structure if None)

        Returns:
            Dictionary of structure properties
        """
        try:
            if structure is not None:
                self.set_structure(structure)
            elif not hasattr(self, "structure"):
                raise ValueError("No structure provided or stored")

            properties = {
                "basic": self.get_basic_properties(),
                "surface": self.analyze_surface(),
                "secondary_structure": self.analyze_secondary_structure(),
                "pockets": self.find_binding_pockets() if self.config.analyze_pockets else [],
                "interfaces": self.analyze_interfaces() if self.config.analyze_interfaces else [],
                "dynamics": self.analyze_dynamics() if self.config.analyze_dynamics else {},
                "conservation": self._analyze_conservation() if self.config.analyze_conservation else {},
                "quality": self._calculate_quality_metrics() if self.config.analyze_quality else {},
            }
            return properties

        except Exception as e:
            self.logger.error(f"Error getting structure properties: {str(e)}")
            return {}

    def get_basic_properties(self) -> Dict[str, Any]:
        """Get basic structural properties."""
        try:
            return {
                "num_residues": len(list(self.structure.get_residues())),
                "num_atoms": len(list(self.structure.get_atoms())),
                "radius_of_gyration": calculate_radius_of_gyration(list(self.structure.get_atoms())),
                "center_of_mass": calculate_center_of_mass(list(self.structure.get_atoms())),
                "sequence_length": len(self.sequence),
            }
        except Exception as e:
            self.logger.error(f"Error getting basic properties: {str(e)}")
            return {}

    def analyze_surface(self) -> Dict[str, Any]:
        """Analyze surface properties."""
        try:
            # Calculate surface area
            surface_area = calculate_surface_area(self.coords, [get_atom_radius(a) for a in self.structure.get_atoms()])

            # Get exposed residues
            exposed_residues = []
            for residue in self.structure.get_residues():
                exposure = calculate_surface_exposure(residue, self.structure)
                if exposure > 0.2:  # Exposure threshold
                    exposed_residues.append(residue.get_id()[1])

            # Calculate surface properties
            surface_props = {
                "total_area": surface_area,
                "exposed_residues": exposed_residues,
                "num_exposed": len(exposed_residues),
                "hydrophobicity": np.mean([get_residue_property(r, "hydrophobicity") for r in self.structure.get_residues() if r.get_id()[1] in exposed_residues]),
                "charge": sum([get_residue_property(r, "charge") for r in self.structure.get_residues() if r.get_id()[1] in exposed_residues]),
                "polar_area": self._calculate_polar_area(),
                "hydrophobic_area": self._calculate_hydrophobic_area(),
            }

            return surface_props

        except Exception as e:
            self.logger.error(f"Error analyzing surface: {str(e)}")
            return {}

    def analyze_secondary_structure(self) -> Dict[str, float]:
        """Analyze secondary structure composition."""
        try:
            ss_map = get_secondary_structure(self.structure)

            # Count secondary structure elements
            ss_counts = {"H": 0, "B": 0, "E": 0, "G": 0, "I": 0, "T": 0, "S": 0, "-": 0}
            for ss in ss_map.values():
                if ss in ss_counts:
                    ss_counts[ss] += 1

            # Calculate percentages
            total = sum(ss_counts.values())
            if total > 0:
                ss_fractions = {k: v / total for k, v in ss_counts.items()}
            else:
                ss_fractions = ss_counts

            return {
                "helix": ss_fractions["H"] + ss_fractions["G"] + ss_fractions["I"],
                "sheet": ss_fractions["B"] + ss_fractions["E"],
                "turn": ss_fractions["T"],
                "coil": ss_fractions["S"] + ss_fractions["-"],
                "details": ss_fractions,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing secondary structure: {str(e)}")
            return {}

    def find_binding_pockets(self, min_volume: float = 100.0) -> List[Dict[str, Any]]:
        """Find potential binding pockets.

        Args:
            min_volume: Minimum pocket volume in Å³

        Returns:
            List of pocket properties
        """
        try:
            pockets = []

            # Get surface points
            surface_points = []
            for atom in self.structure.get_atoms():
                surface_points.append(atom.get_coord())
            surface_points = np.array(surface_points)

            # Find cavities
            from scipy.spatial import ConvexHull

            hull = ConvexHull(surface_points)

            # Analyze each cavity
            for simplex in hull.simplices:
                points = surface_points[simplex]
                volume = calculate_cavity_volume(points)

                if volume >= min_volume:
                    center = points.mean(axis=0)

                    # Get residues lining the pocket
                    pocket_residues = []
                    for residue in self.structure.get_residues():
                        if residue.has_id("CA"):
                            dist = np.linalg.norm(residue["CA"].get_coord() - center)
                            if dist < 8.0:  # Distance cutoff
                                pocket_residues.append(residue.get_id()[1])

                    pocket = {
                        "center": center,
                        "volume": volume,
                        "residues": pocket_residues,
                        "depth": calculate_residue_depth(list(self.structure.get_residues())[0], surface_points),
                        "hydrophobicity": np.mean([get_residue_property(r, "hydrophobicity") for r in self.structure.get_residues() if r.get_id()[1] in pocket_residues]),
                    }
                    pockets.append(pocket)

            return sorted(pockets, key=lambda x: x["volume"], reverse=True)

        except Exception as e:
            self.logger.error(f"Error finding binding pockets: {str(e)}")
            return []

    def analyze_interfaces(self) -> List[Dict[str, Any]]:
        """Analyze chain interfaces."""
        try:
            interfaces = []
            chains = list(self.structure.get_chains())

            for i in range(len(chains)):
                for j in range(i + 1, len(chains)):
                    chain1, chain2 = chains[i], chains[j]

                    # Find contacting residues
                    contacts = find_contacts(list(chain1.get_atoms()), list(chain2.get_atoms()))

                    if contacts:
                        # Get unique residue pairs
                        residue_pairs = set()
                        for a1, a2, _ in contacts:
                            res1 = a1.get_parent().get_id()[1]
                            res2 = a2.get_parent().get_id()[1]
                            residue_pairs.add((res1, res2))

                        # Calculate interface area
                        area = calculate_interface_area(chain1, chain2)

                        interface = {
                            "chains": (chain1.id, chain2.id),
                            "num_contacts": len(contacts),
                            "residue_pairs": list(residue_pairs),
                            "area": area,
                            "hydrophobicity": np.mean(
                                [get_residue_property(r, "hydrophobicity") for r in chain1.get_residues() if r.get_id()[1] in [p[0] for p in residue_pairs]]
                                + [get_residue_property(r, "hydrophobicity") for r in chain2.get_residues() if r.get_id()[1] in [p[1] for p in residue_pairs]]
                            ),
                        }
                        interfaces.append(interface)

            return interfaces

        except Exception as e:
            self.logger.error(f"Error analyzing interfaces: {str(e)}")
            return []

    def analyze_dynamics(self) -> Dict[str, Any]:
        """Analyze protein dynamics using B-factors."""
        try:
            b_factors = []
            residue_factors = {}

            for residue in self.structure.get_residues():
                res_factors = []
                for atom in residue:
                    if atom.bfactor is not None:
                        b_factors.append(atom.bfactor)
                        res_factors.append(atom.bfactor)
                if res_factors:
                    residue_factors[residue.get_id()[1]] = np.mean(res_factors)

            if not b_factors:
                return {}

            # Calculate statistics
            dynamics = {
                "mean_bfactor": float(np.mean(b_factors)),
                "std_bfactor": float(np.std(b_factors)),
                "min_bfactor": float(np.min(b_factors)),
                "max_bfactor": float(np.max(b_factors)),
                "residue_bfactors": residue_factors,
            }

            # Identify flexible regions (high B-factors)
            mean = np.mean(b_factors)
            std = np.std(b_factors)
            flexible_residues = [res_id for res_id, b_factor in residue_factors.items() if b_factor > mean + std]
            dynamics["flexible_residues"] = flexible_residues

            return dynamics

        except Exception as e:
            self.logger.error(f"Error analyzing dynamics: {str(e)}")
            return {}

    def _calculate_polar_area(self) -> float:
        """Calculate polar surface area."""
        try:
            polar_atoms = []
            radii = []

            for atom in self.structure.get_atoms():
                if atom.element in ["N", "O"]:
                    polar_atoms.append(atom.get_coord())
                    radii.append(get_atom_radius(atom))

            if polar_atoms:
                return calculate_surface_area(np.array(polar_atoms), np.array(radii))
            return 0.0

        except Exception as e:
            self.logger.error(f"Error calculating polar area: {str(e)}")
            return 0.0

    def _calculate_hydrophobic_area(self) -> float:
        """Calculate hydrophobic surface area."""
        try:
            hydrophobic_atoms = []
            radii = []

            for atom in self.structure.get_atoms():
                if atom.element == "C":
                    # Check if carbon is bonded to polar atoms
                    is_hydrophobic = True
                    for neighbor in atom.get_parent():
                        if neighbor.element in ["N", "O"]:
                            is_hydrophobic = False
                            break

                    if is_hydrophobic:
                        hydrophobic_atoms.append(atom.get_coord())
                        radii.append(get_atom_radius(atom))

            if hydrophobic_atoms:
                return calculate_surface_area(np.array(hydrophobic_atoms), np.array(radii))
            return 0.0

        except Exception as e:
            self.logger.error(f"Error calculating hydrophobic area: {str(e)}")
            return 0.0

    def _analyze_conservation(self) -> Dict[str, float]:
        """Analyze sequence conservation."""
        try:
            # This is a placeholder for conservation analysis
            # Would need sequence alignment data
            return {}
        except Exception as e:
            self.logger.error(f"Error analyzing conservation: {str(e)}")
            return {}

    def _calculate_quality_metrics(self) -> Dict[str, float]:
        """Calculate structure quality metrics."""
        try:
            metrics = {
                "clashes": self._count_clashes(),
                "rama_outliers": self._count_ramachandran_outliers(),
                "rotamer_outliers": self._count_rotamer_outliers(),
            }
            return metrics

        except Exception as e:
            self.logger.error(f"Error calculating quality metrics: {str(e)}")
            return {}

    def _count_clashes(self) -> int:
        """Count atomic clashes."""
        try:
            clashes = 0
            atoms = list(self.structure.get_atoms())

            for i in range(len(atoms)):
                for j in range(i + 1, len(atoms)):
                    # Skip bonded atoms
                    if atoms[i].get_parent() is atoms[j].get_parent():
                        continue

                    # Get atomic radii
                    r1 = get_atom_radius(atoms[i])
                    r2 = get_atom_radius(atoms[j])

                    # Calculate overlap
                    dist = np.linalg.norm(atoms[i].get_coord() - atoms[j].get_coord())
                    overlap = r1 + r2 - dist

                    # Count significant clashes
                    if overlap > self.config.clash_overlap:
                        clashes += 1

            return clashes

        except Exception as e:
            self.logger.error(f"Error counting clashes: {str(e)}")
            return 0

    def _count_ramachandran_outliers(self) -> int:
        """Count Ramachandran plot outliers."""
        try:
            outliers = 0

            # Ramachandran plot regions (phi, psi) in degrees
            allowed_regions = [
                # Alpha helix
                (-90, -45, -70, -15),
                # Beta sheet
                (-150, -90, 100, 170),
                # Left-handed helix
                (30, 100, -15, 50),
            ]

            for model in self.structure:
                for chain in model:
                    for i, residue in enumerate(chain):
                        if i == 0 or i == len(chain) - 1:
                            continue

                        # Get backbone atoms
                        try:
                            prev_c = chain[i - 1]["C"]
                            n = residue["N"]
                            ca = residue["CA"]
                            c = residue["C"]
                            next_n = chain[i + 1]["N"]
                        except KeyError:
                            continue

                        # Calculate phi angle
                        phi = calc_dihedral(prev_c.get_vector(), n.get_vector(), ca.get_vector(), c.get_vector())

                        # Calculate psi angle
                        psi = calc_dihedral(n.get_vector(), ca.get_vector(), c.get_vector(), next_n.get_vector())

                        # Convert to degrees
                        phi = np.degrees(phi)
                        psi = np.degrees(psi)

                        # Check if in allowed regions
                        is_allowed = False
                        for region in allowed_regions:
                            if region[0] <= phi <= region[1] and region[2] <= psi <= region[3]:
                                is_allowed = True
                                break

                        if not is_allowed:
                            outliers += 1

            return outliers

        except Exception as e:
            self.logger.error(f"Error counting Ramachandran outliers: {str(e)}")
            return 0

    def _count_rotamer_outliers(self) -> int:
        """Count rotamer outliers."""
        try:
            outliers = 0

            # Common rotamer angles (chi1, chi2) in degrees
            rotamers = {
                "LEU": [(62, 175), (-177, 65), (-65, -60)],
                "ILE": [(62, 175), (-177, 65), (-65, -60)],
                "VAL": [(175,), (65,), (-60,)],
                "PHE": [(62, 90), (-177, 80), (-65, -85)],
                "TYR": [(62, 90), (-177, 80), (-65, -85)],
                "TRP": [(62, 90), (-177, 80), (-65, -85)],
            }

            for residue in self.structure.get_residues():
                if residue.get_resname() in rotamers:
                    try:
                        # Calculate chi angles
                        chi1 = self._calculate_chi1(residue)
                        chi2 = self._calculate_chi2(residue)

                        # Check if matches any rotamer
                        is_rotamer = False
                        for rot in rotamers[residue.get_resname()]:
                            if len(rot) == 1:
                                if abs(chi1 - rot[0]) < 30:  # 30 degree tolerance
                                    is_rotamer = True
                                    break
                            else:
                                if abs(chi1 - rot[0]) < 30 and abs(chi2 - rot[1]) < 30:
                                    is_rotamer = True
                                    break

                        if not is_rotamer:
                            outliers += 1

                    except:
                        continue

            return outliers

        except Exception as e:
            self.logger.error(f"Error counting rotamer outliers: {str(e)}")
            return 0

    def _calculate_chi1(self, residue: Residue) -> float:
        """Calculate chi1 torsion angle."""
        try:
            n = residue["N"].get_coord()
            ca = residue["CA"].get_coord()
            cb = residue["CB"].get_coord()
            cg = residue["CG"].get_coord()
            return calc_dihedral(Vector(n), Vector(ca), Vector(cb), Vector(cg))
        except:
            raise ValueError("Could not calculate chi1")

    def _calculate_chi2(self, residue: Residue) -> float:
        """Calculate chi2 torsion angle."""
        try:
            ca = residue["CA"].get_coord()
            cb = residue["CB"].get_coord()
            cg = residue["CG"].get_coord()
            cd = residue["CD"].get_coord()
            return calc_dihedral(Vector(ca), Vector(cb), Vector(cg), Vector(cd))
        except:
            raise ValueError("Could not calculate chi2")
