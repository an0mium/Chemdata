"""Protein dynamics analysis functionality."""

import logging
from typing import Dict, List, Optional, Any, Tuple
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_angle, calc_dihedral
from Bio.PDB import NMA
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
import networkx as nx

logger = logging.getLogger(__name__)


class DynamicsAnalyzer:
    """Analyzes protein dynamics and flexibility."""

    def __init__(self):
        """Initialize dynamics analyzer."""
        self.logger = logging.getLogger(__name__)
        # Add domain analysis parameters
        self.domain_params = {
            "min_domain_size": 20,  # Minimum residues for a domain
            "interface_cutoff": 8.0,  # Å for interface detection
            "packing_cutoff": 4.0,  # Å for domain packing
        }

    def analyze_dynamics(
        self,
        structure: Structure,
        use_bfactors: bool = True,
        use_normal_modes: bool = True,
        use_contacts: bool = True,
        use_flexibility: bool = True,
        use_domain_motions: bool = True,
        use_correlations: bool = True,
    ) -> Dict[str, Any]:
        """Analyze protein dynamics using multiple methods.

        Args:
            structure: BioPython Structure object
            use_bfactors: Whether to use B-factors
            use_normal_modes: Whether to calculate normal modes
            use_contacts: Whether to analyze contact networks
            use_flexibility: Whether to analyze flexibility
            use_domain_motions: Whether to analyze domain motions
            use_correlations: Whether to analyze correlations

        Returns:
            Dictionary of dynamics properties
        """
        try:
            dynamics = {}

            # B-factor analysis
            if use_bfactors:
                bfactor_props = self._analyze_bfactors(structure)
                dynamics.update(bfactor_props)

            # Normal mode analysis
            if use_normal_modes:
                mode_props = self._analyze_normal_modes(structure)
                dynamics.update(mode_props)

            # Contact network analysis
            if use_contacts:
                contact_props = self._analyze_contacts(structure)
                dynamics.update(contact_props)

            # Flexibility analysis
            if use_flexibility:
                flex_props = self._analyze_flexibility(structure)
                dynamics.update(flex_props)

            # Domain motion analysis
            if use_domain_motions:
                motion_props = self._analyze_domain_motions(structure)
                dynamics.update(motion_props)

            # Correlation analysis
            if use_correlations:
                corr_props = self._analyze_correlations(dynamics)
                dynamics.update(corr_props)

            # Calculate overall dynamics metrics
            dynamics["metrics"] = self._calculate_dynamics_metrics(dynamics)

            return dynamics

        except Exception as e:
            self.logger.error(f"Error analyzing dynamics: {str(e)}")
            return {}

    def analyze_site_dynamics(
        self,
        residues: List[Residue],
        structure: Structure,
        include_domain_context: bool = True,
    ) -> Dict[str, float]:
        """Analyze dynamics of binding site residues.

        Args:
            residues: List of residues in binding site
            structure: Full structure for context
            include_domain_context: Whether to analyze domain context

        Returns:
            Dictionary of dynamics metrics
        """
        try:
            # Get B-factors for site residues
            bfactors = []
            for res in residues:
                for atom in res:
                    bfactors.append(atom.get_bfactor())

            # Calculate flexibility scores
            flexibility = self._calculate_residue_flexibility(residues, structure)

            # Calculate contact network properties
            contact_props = self._analyze_site_contacts(
                [res.get_id()[1] for res in residues],
                self._analyze_contacts(structure).get("contact_network", {}),
            )

            # Get domain context if requested
            domain_context = {}
            if include_domain_context:
                domain_context = self._analyze_site_domain_context(residues, structure)

            # Calculate site metrics
            metrics = {
                "average_bfactor": float(np.mean(bfactors)),
                "bfactor_std": float(np.std(bfactors)),
                "relative_bfactor": self._calculate_relative_bfactor(bfactors, structure),
                "flexibility_score": float(np.mean(flexibility)),
                "rigidity_score": 1.0 - float(np.mean(flexibility)),
                "variability": float(np.std(flexibility)),
                "contact_density": contact_props.get("contact_density", 0.0),
                "surface_exposure": contact_props.get("surface_exposure", 0.0),
            }

            # Add domain context metrics
            metrics.update(domain_context)

            return metrics

        except Exception as e:
            self.logger.error(f"Error analyzing site dynamics: {str(e)}")
            return {}

    def get_residue_dynamics(
        self,
        residue: Residue,
        structure: Structure,
    ) -> Dict[str, float]:
        """Get dynamics properties for single residue.

        Args:
            residue: BioPython Residue object
            structure: Full structure for context

        Returns:
            Dictionary of dynamics properties
        """
        try:
            # Get B-factors
            bfactors = [atom.get_bfactor() for atom in residue]
            avg_bfactor = float(np.mean(bfactors))

            # Calculate relative B-factor
            all_bfactors = []
            for res in structure.get_residues():
                for atom in res:
                    all_bfactors.append(atom.get_bfactor())
            relative_bfactor = (avg_bfactor - np.mean(all_bfactors)) / np.std(all_bfactors)

            # Calculate flexibility
            flexibility = self._calculate_residue_flexibility([residue], structure)[0]

            # Get contact network properties
            contacts = self._analyze_contacts(structure).get("contact_network", {})
            degree = len(contacts.get(residue.get_id()[1], set()))

            # Get domain membership
            domains = self._identify_domains(structure)
            domain_id = None
            for i, domain in enumerate(domains):
                if residue.get_id()[1] in domain:
                    domain_id = i + 1
                    break

            return {
                "bfactor": avg_bfactor,
                "relative_bfactor": float(relative_bfactor),
                "flexibility": float(flexibility),
                "rigidity": 1.0 - float(flexibility),
                "contact_degree": degree,
                "domain": domain_id,
            }

        except Exception as e:
            self.logger.error(f"Error getting residue dynamics: {str(e)}")
            return {}

    def _analyze_bfactors(self, structure: Structure) -> Dict[str, Any]:
        """Analyze B-factors across structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of B-factor analysis
        """
        try:
            # Collect B-factors
            bfactors = {
                "backbone": [],
                "sidechain": [],
                "all": [],
            }

            residue_bfactors = {}

            for residue in structure.get_residues():
                res_bfactors = []
                for atom in residue:
                    bfactor = atom.get_bfactor()
                    bfactors["all"].append(bfactor)
                    res_bfactors.append(bfactor)
                    if atom.get_name() in ["N", "CA", "C", "O"]:
                        bfactors["backbone"].append(bfactor)
                    else:
                        bfactors["sidechain"].append(bfactor)

                # Store average B-factor for residue
                if res_bfactors:
                    residue_bfactors[residue.get_id()[1]] = float(np.mean(res_bfactors))

            # Calculate statistics
            stats = {}
            for region, values in bfactors.items():
                if values:
                    stats[region] = {
                        "mean": float(np.mean(values)),
                        "std": float(np.std(values)),
                        "min": float(np.min(values)),
                        "max": float(np.max(values)),
                    }

            return {
                "bfactors": residue_bfactors,
                "statistics": stats,
                "mobility_profile": self._calculate_mobility_profile(residue_bfactors),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing B-factors: {str(e)}")
            return {}

    def _analyze_normal_modes(self, structure: Structure) -> Dict[str, Any]:
        """Calculate and analyze normal modes.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of normal mode properties
        """
        try:
            # Get CA coordinates
            coords = []
            residue_ids = []
            for residue in structure.get_residues():
                if "CA" in residue:
                    coords.append(residue["CA"].get_coord())
                    residue_ids.append(residue.get_id()[1])
            coords = np.array(coords)

            if len(coords) < 3:
                return {}

            # Build Hessian matrix (simplified elastic network model)
            cutoff = 12.0  # Å
            distances = cdist(coords, coords)
            hessian = np.zeros((3 * len(coords), 3 * len(coords)))

            # Fill Hessian using simple spring model
            for i in range(len(coords)):
                for j in range(len(coords)):
                    if i != j and distances[i, j] < cutoff:
                        diff = coords[i] - coords[j]
                        magnitude = np.linalg.norm(diff)
                        direction = diff / magnitude
                        force_constant = 1.0  # Simplified uniform force constant

                        # 3x3 block for residue pair
                        block = force_constant * np.outer(direction, direction)

                        # Add to Hessian
                        i3, j3 = i * 3, j * 3
                        hessian[i3 : i3 + 3, j3 : j3 + 3] = -block
                        hessian[i3 : i3 + 3, i3 : i3 + 3] += block

            # Calculate normal modes (lowest frequency modes)
            try:
                eigenvals, eigenvecs = np.linalg.eigh(hessian)
            except np.linalg.LinAlgError:
                return {}

            # Remove translational and rotational modes
            n_rigid = 6
            eigenvals = eigenvals[n_rigid:]
            eigenvecs = eigenvecs[:, n_rigid:]

            # Calculate residue fluctuations
            fluctuations = {}
            for i, res_id in enumerate(residue_ids):
                # Sum contributions from all modes
                total_fluct = 0.0
                for j in range(len(eigenvals)):
                    if eigenvals[j] > 1e-10:  # Avoid division by zero
                        mode = eigenvecs[i * 3 : (i + 1) * 3, j]
                        total_fluct += np.sum(mode * mode) / eigenvals[j]
                fluctuations[res_id] = float(total_fluct)

            # Calculate mode properties
            modes = []
            for i in range(min(10, len(eigenvals))):
                modes.append(
                    {
                        "frequency": float(np.sqrt(abs(eigenvals[i]))),
                        "collectivity": float(self._calculate_collectivity(eigenvecs[:, i])),
                        "amplitude": float(1.0 / np.sqrt(abs(eigenvals[i]))),
                    }
                )

            return {
                "mode_fluctuations": fluctuations,
                "mean_fluctuation": float(np.mean(list(fluctuations.values()))),
                "modes": modes,
                "total_variance": float(np.sum(1.0 / np.abs(eigenvals))),
                "collective_modes": self._analyze_collective_modes(eigenvecs),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing normal modes: {str(e)}")
            return {}

    def _analyze_contacts(self, structure: Structure) -> Dict[str, Any]:
        """Analyze residue contact network.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of contact network properties
        """
        try:
            # Build contact map
            contacts = {}
            residue_list = list(structure.get_residues())

            for i, res1 in enumerate(residue_list):
                res1_id = res1.get_id()[1]
                contacts[res1_id] = set()

                for j, res2 in enumerate(residue_list[i + 1 :], i + 1):
                    res2_id = res2.get_id()[1]

                    # Check contact between residues
                    if self._check_residue_contact(res1, res2):
                        contacts[res1_id].add(res2_id)
                        if res2_id not in contacts:
                            contacts[res2_id] = set()
                        contacts[res2_id].add(res1_id)

            # Calculate network properties
            degrees = {res_id: len(neighbors) for res_id, neighbors in contacts.items()}
            mean_degree = float(np.mean(list(degrees.values())))

            # Identify hubs (highly connected residues)
            hubs = [res_id for res_id, degree in degrees.items() if degree > mean_degree + np.std(list(degrees.values()))]

            return {
                "contact_network": contacts,
                "degrees": degrees,
                "mean_degree": mean_degree,
                "hubs": hubs,
                "clustering": self._calculate_clustering(contacts),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing contacts: {str(e)}")
            return {}

    def _analyze_flexibility(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein flexibility.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of flexibility analysis
        """
        try:
            # Calculate per-residue flexibility
            flexibility = []
            residue_ids = []

            for residue in structure.get_residues():
                score = self._calculate_residue_flexibility([residue], structure)[0]
                flexibility.append(score)
                residue_ids.append(residue.get_id()[1])

            # Identify flexible and rigid regions
            mean_flex = np.mean(flexibility)
            std_flex = np.std(flexibility)
            flexible_regions = []
            rigid_regions = []

            current_flexible = []
            current_rigid = []

            for i, (res_id, flex) in enumerate(zip(residue_ids, flexibility)):
                # Flexible region
                if flex > mean_flex + std_flex:
                    current_flexible.append(res_id)
                    if current_rigid:
                        rigid_regions.append(current_rigid)
                        current_rigid = []
                # Rigid region
                elif flex < mean_flex - std_flex:
                    current_rigid.append(res_id)
                    if current_flexible:
                        flexible_regions.append(current_flexible)
                        current_flexible = []
                # Neither
                else:
                    if current_flexible:
                        flexible_regions.append(current_flexible)
                        current_flexible = []
                    if current_rigid:
                        rigid_regions.append(current_rigid)
                        current_rigid = []

            # Add any remaining regions
            if current_flexible:
                flexible_regions.append(current_flexible)
            if current_rigid:
                rigid_regions.append(current_rigid)

            return {
                "residue_flexibility": dict(zip(residue_ids, flexibility)),
                "flexible_regions": flexible_regions,
                "rigid_regions": rigid_regions,
                "average_flexibility": float(mean_flex),
                "flexibility_std": float(std_flex),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing flexibility: {str(e)}")
            return {}

    def _calculate_dynamics_metrics(self, dynamics: Dict[str, Any]) -> Dict[str, float]:
        """Calculate overall dynamics metrics.

        Args:
            dynamics: Dictionary of dynamics data

        Returns:
            Dictionary of dynamics metrics
        """
        try:
            metrics = {}

            # B-factor based metrics
            if "statistics" in dynamics:
                stats = dynamics["statistics"]
                if "all" in stats:
                    metrics.update(
                        {
                            "average_bfactor": stats["all"]["mean"],
                            "bfactor_std": stats["all"]["std"],
                            "backbone_flexibility": stats.get("backbone", {}).get("mean", 0.0),
                            "sidechain_flexibility": stats.get("sidechain", {}).get("mean", 0.0),
                        }
                    )

            # Normal mode based metrics
            if "modes" in dynamics:
                modes = dynamics["modes"]
                if modes:
                    metrics.update(
                        {
                            "lowest_frequency": modes[0]["frequency"],
                            "average_collectivity": float(np.mean([m["collectivity"] for m in modes])),
                            "total_variance": dynamics.get("total_variance", 0.0),
                        }
                    )

            # Contact network metrics
            if "contact_network" in dynamics:
                metrics.update(
                    {
                        "average_contacts": dynamics.get("mean_degree", 0.0),
                        "num_hubs": len(dynamics.get("hubs", [])),
                    }
                )

            # Flexibility based metrics
            if "residue_flexibility" in dynamics:
                metrics.update(
                    {
                        "average_flexibility": dynamics.get("average_flexibility", 0.0),
                        "flexibility_std": dynamics.get("flexibility_std", 0.0),
                        "num_flexible_regions": len(dynamics.get("flexible_regions", [])),
                        "num_rigid_regions": len(dynamics.get("rigid_regions", [])),
                    }
                )

            return metrics

        except Exception as e:
            self.logger.error(f"Error calculating dynamics metrics: {str(e)}")
            return {}

    def _calculate_mobility_profile(
        self,
        bfactors: Dict[int, float],
    ) -> Dict[str, List[int]]:
        """Calculate mobility profile from B-factors.

        Args:
            bfactors: Dictionary mapping residue IDs to B-factors

        Returns:
            Dictionary categorizing residues by mobility
        """
        try:
            values = list(bfactors.values())
            if not values:
                return {"rigid": [], "moderate": [], "flexible": []}

            mean = np.mean(values)
            std = np.std(values)

            profile = {
                "rigid": [],  # < mean - std
                "moderate": [],  # between mean ± std
                "flexible": [],  # > mean + std
            }

            for res_id, bfactor in bfactors.items():
                if bfactor < mean - std:
                    profile["rigid"].append(res_id)
                elif bfactor > mean + std:
                    profile["flexible"].append(res_id)
                else:
                    profile["moderate"].append(res_id)

            return profile

        except Exception as e:
            self.logger.error(f"Error calculating mobility profile: {str(e)}")
            return {"rigid": [], "moderate": [], "flexible": []}

    def _analyze_collective_modes(self, eigenvecs: np.ndarray) -> Dict[str, Any]:
        """Analyze collective motions from normal modes.

        Args:
            eigenvecs: Normal mode eigenvectors

        Returns:
            Dictionary of collective motion properties
        """
        try:
            n_modes = min(10, eigenvecs.shape[1])  # Analyze top 10 modes
            collectivity = []

            for i in range(n_modes):
                mode = eigenvecs[:, i]
                # Calculate collectivity (participation ratio)
                squared = mode * mode
                sum_squared = np.sum(squared)
                if sum_squared > 0:
                    participation = np.exp(-np.sum((squared / sum_squared) * np.log(squared / sum_squared))) / len(mode)
                    collectivity.append(float(participation))

            return {
                "collectivity": collectivity,
                "mean_collectivity": float(np.mean(collectivity)) if collectivity else 0.0,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing collective modes: {str(e)}")
            return {}

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

    def _calculate_clustering(self, contacts: Dict[int, set]) -> Dict[int, float]:
        """Calculate clustering coefficients for contact network.

        Args:
            contacts: Dictionary mapping residue IDs to sets of contacting residues

        Returns:
            Dictionary mapping residue IDs to clustering coefficients
        """
        try:
            clustering = {}

            for res_id, neighbors in contacts.items():
                if len(neighbors) < 2:
                    clustering[res_id] = 0.0
                    continue

                # Count connections between neighbors
                connections = 0
                for n1 in neighbors:
                    for n2 in neighbors:
                        if n1 < n2 and n2 in contacts.get(n1, set()):
                            connections += 1

                # Calculate clustering coefficient
                max_connections = len(neighbors) * (len(neighbors) - 1) / 2
                if max_connections > 0:
                    clustering[res_id] = float(connections / max_connections)
                else:
                    clustering[res_id] = 0.0

            return clustering

        except Exception as e:
            self.logger.error(f"Error calculating clustering: {str(e)}")
            return {}

    def _analyze_site_contacts(
        self,
        site_residues: List[int],
        contact_network: Dict[int, set],
    ) -> Dict[str, Any]:
        """Analyze contact network properties for a specific site.

        Args:
            site_residues: List of residue numbers in site
            contact_network: Pre-calculated contact network

        Returns:
            Dictionary of site contact properties
        """
        try:
            # Count internal and external contacts
            internal_contacts = 0
            external_contacts = 0
            site_set = set(site_residues)

            for res_id in site_residues:
                neighbors = contact_network.get(res_id, set())
                internal_contacts += len(neighbors & site_set)
                external_contacts += len(neighbors - site_set)

            # Adjust for double counting of internal contacts
            internal_contacts //= 2

            return {
                "internal_contacts": internal_contacts,
                "external_contacts": external_contacts,
                "contact_density": float(internal_contacts / (len(site_residues) * (len(site_residues) - 1) / 2) if len(site_residues) > 1 else 0.0),
                "surface_exposure": float(external_contacts / (internal_contacts + external_contacts) if (internal_contacts + external_contacts) > 0 else 0.0),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing site contacts: {str(e)}")
            return {}

    def _calculate_residue_flexibility(
        self,
        residues: List[Residue],
        structure: Structure,
    ) -> List[float]:
        """Calculate flexibility scores for residues.

        Args:
            residues: List of residues to analyze
            structure: Full structure for context

        Returns:
            List of flexibility scores
        """
        try:
            scores = []

            for residue in residues:
                # Get local contacts
                local_contacts = self._count_local_contacts(residue, structure)

                # Get backbone angles
                phi, psi = self._calculate_backbone_angles(residue)

                # Calculate flexibility score components
                contact_score = 1.0 - (local_contacts / 12.0)  # Normalize by typical max contacts
                angle_score = 0.0
                if phi is not None and psi is not None:
                    # Higher score for non-regular secondary structure
                    if not (-140 < phi < -60 and -70 < psi < -15):  # Not helix
                        if not (-150 < phi < -50 and 100 < psi < 180):  # Not sheet
                            angle_score = 1.0

                # Combine scores (weighted average)
                flexibility = 0.6 * contact_score + 0.4 * angle_score
                scores.append(float(flexibility))

            return scores

        except Exception as e:
            self.logger.error(f"Error calculating residue flexibility: {str(e)}")
            return [0.0] * len(residues)

    def _count_local_contacts(
        self,
        residue: Residue,
        structure: Structure,
        cutoff: float = 8.0,
    ) -> int:
        """Count contacts within cutoff distance.

        Args:
            residue: Residue to analyze
            structure: Structure containing residue
            cutoff: Distance cutoff in Angstroms

        Returns:
            Number of contacts
        """
        try:
            contacts = 0
            res_atoms = [atom.get_coord() for atom in residue]

            for other_res in structure.get_residues():
                if other_res != residue:
                    other_atoms = [atom.get_coord() for atom in other_res]
                    for coord1 in res_atoms:
                        for coord2 in other_atoms:
                            dist = np.linalg.norm(coord1 - coord2)
                            if dist < cutoff:
                                contacts += 1
                                break

            return contacts

        except Exception as e:
            self.logger.error(f"Error counting contacts: {str(e)}")
            return 0

    def _calculate_backbone_angles(
        self,
        residue: Residue,
    ) -> Tuple[Optional[float], Optional[float]]:
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

    def _calculate_relative_bfactor(
        self,
        site_bfactors: List[float],
        structure: Structure,
    ) -> float:
        """Calculate relative B-factor compared to whole structure.

        Args:
            site_bfactors: List of B-factors for site
            structure: Full structure

        Returns:
            Relative B-factor score
        """
        try:
            # Get all B-factors
            all_bfactors = []
            for residue in structure.get_residues():
                for atom in residue:
                    all_bfactors.append(atom.get_bfactor())

            # Calculate Z-score
            site_mean = np.mean(site_bfactors)
            all_mean = np.mean(all_bfactors)
            all_std = np.std(all_bfactors)

            if all_std == 0:
                return 0.0

            return float((site_mean - all_mean) / all_std)

        except Exception as e:
            self.logger.error(f"Error calculating relative B-factor: {str(e)}")
            return 0.0

    def _calculate_collectivity(self, eigenvec: np.ndarray) -> float:
        """Calculate mode collectivity.

        Args:
            eigenvec: Mode eigenvector

        Returns:
            Collectivity (0-1)
        """
        try:
            # Normalize vector
            vec = eigenvec / np.sqrt(np.sum(eigenvec * eigenvec))
            vec = vec * vec

            # Calculate collectivity
            n = len(vec)
            entropy = -np.sum(vec * np.log(vec + 1e-6))
            return float(np.exp(entropy) / n)

        except Exception as e:
            self.logger.error(f"Error calculating collectivity: {str(e)}")
            return 0.0

    def _analyze_domain_motions(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein domain motions.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of domain motion properties
        """
        try:
            # Identify domains using contact topology
            domains = self._identify_domains(structure)

            # Calculate domain-level properties
            domain_properties = {}
            for i, domain in enumerate(domains):
                # Calculate domain center
                domain_coords = []
                for res_id in domain:
                    residue = structure[0]["A"][res_id]
                    for atom in residue:
                        domain_coords.append(atom.get_coord())
                center = np.mean(domain_coords, axis=0)

                # Calculate domain radius of gyration
                rg = np.sqrt(np.mean(np.sum((domain_coords - center) ** 2, axis=1)))

                # Calculate domain flexibility
                domain_residues = [structure[0]["A"][res_id] for res_id in domain]
                flexibility = np.mean(self._calculate_residue_flexibility(domain_residues, structure))

                domain_properties[f"domain_{i+1}"] = {
                    "residues": domain,
                    "center": center,
                    "radius_gyration": float(rg),
                    "flexibility": float(flexibility),
                    "size": len(domain),
                    "compactness": float(self._calculate_domain_compactness(domain_coords)),
                }

            # Calculate inter-domain properties
            inter_domain = {}
            for i, (name1, domain1) in enumerate(domain_properties.items()):
                for name2, domain2 in list(domain_properties.items())[i + 1 :]:
                    # Calculate distance between domain centers
                    center_dist = np.linalg.norm(domain1["center"] - domain2["center"])

                    # Calculate interface residues
                    interface = self._get_domain_interface(
                        domain1["residues"],
                        domain2["residues"],
                        structure,
                    )

                    # Calculate interface properties
                    interface_props = self._analyze_interface_properties(interface, structure)

                    inter_domain[f"{name1}_{name2}"] = {
                        "center_distance": float(center_dist),
                        "interface_residues": interface,
                        "interface_size": len(interface),
                        "interface_properties": interface_props,
                        "relative_orientation": self._calculate_domain_orientation(
                            domain1["center"],
                            domain2["center"],
                            structure,
                        ),
                    }

            # Calculate global domain organization
            global_properties = self._analyze_global_domain_organization(
                domain_properties,
                inter_domain,
            )

            return {
                "domains": domain_properties,
                "inter_domain": inter_domain,
                "global_properties": global_properties,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing domain motions: {str(e)}")
            return {}

    def _identify_domains(self, structure: Structure) -> List[List[int]]:
        """Identify protein domains using contact topology.

        Args:
            structure: BioPython Structure object

        Returns:
            List of domain residue lists
        """
        try:
            # Build contact network
            contacts = self._analyze_contacts(structure).get("contact_network", {})

            # Create graph
            G = nx.Graph()
            for res_id, neighbors in contacts.items():
                for neighbor in neighbors:
                    G.add_edge(res_id, neighbor)

            # Use community detection to identify domains
            communities = nx.community.louvain_communities(G)

            # Filter and merge small domains
            min_domain_size = self.domain_params["min_domain_size"]
            merged_domains = []
            current_domain = []

            for community in sorted(communities, key=len, reverse=True):
                if len(community) >= min_domain_size:
                    merged_domains.append(sorted(community))
                else:
                    current_domain.extend(community)
                    if len(current_domain) >= min_domain_size:
                        merged_domains.append(sorted(current_domain))
                        current_domain = []

            if current_domain:  # Add any remaining residues to closest domain
                if merged_domains:
                    merged_domains[0].extend(current_domain)
                else:
                    merged_domains.append(sorted(current_domain))

            return merged_domains

        except Exception as e:
            self.logger.error(f"Error identifying domains: {str(e)}")
            return []

    def _calculate_domain_compactness(self, coords: np.ndarray) -> float:
        """Calculate domain compactness using radius of gyration.

        Args:
            coords: Domain atomic coordinates

        Returns:
            Compactness score (0-1)
        """
        try:
            if len(coords) < 3:
                return 0.0

            # Calculate radius of gyration
            center = np.mean(coords, axis=0)
            rg = np.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1)))

            # Calculate theoretical minimum Rg for this number of atoms
            n_atoms = len(coords)
            min_rg = 3.0 * np.power(n_atoms, 1 / 3)  # Approximate for globular protein

            # Calculate compactness score (1 = maximally compact)
            return float(min_rg / rg if rg > 0 else 0.0)

        except Exception as e:
            self.logger.error(f"Error calculating domain compactness: {str(e)}")
            return 0.0

    def _analyze_interface_properties(
        self,
        interface_residues: List[int],
        structure: Structure,
    ) -> Dict[str, Any]:
        """Analyze properties of domain interface.

        Args:
            interface_residues: List of interface residue numbers
            structure: BioPython Structure object

        Returns:
            Dictionary of interface properties
        """
        try:
            if not interface_residues:
                return {}

            # Get interface residue objects
            residues = [structure[0]["A"][res_id] for res_id in interface_residues]

            # Calculate properties
            properties = {
                "hydrophobicity": float(np.mean([self.HYDROPHOBICITY.get(res.get_resname(), 0.0) for res in residues])),
                "charged_residues": len([res for res in residues if res.get_resname() in ["ARG", "LYS", "ASP", "GLU"]]),
                "aromatic_residues": len([res for res in residues if res.get_resname() in ["PHE", "TYR", "TRP"]]),
                "polar_residues": len([res for res in residues if res.get_resname() in ["SER", "THR", "ASN", "GLN"]]),
                "contact_density": self._analyze_site_contacts(
                    interface_residues,
                    self._analyze_contacts(structure).get("contact_network", {}),
                ).get("contact_density", 0.0),
            }

            # Calculate conservation if available
            conservation = self._calculate_conservation(residues)
            if conservation:
                properties["conservation"] = float(np.mean(list(conservation.values())))

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing interface properties: {str(e)}")
            return {}

    def _calculate_domain_orientation(
        self,
        center1: np.ndarray,
        center2: np.ndarray,
        structure: Structure,
    ) -> Dict[str, float]:
        """Calculate relative orientation between domains.

        Args:
            center1: First domain center
            center2: Second domain center
            structure: BioPython Structure object

        Returns:
            Dictionary of orientation metrics
        """
        try:
            # Calculate domain axis vector
            axis = center2 - center1
            axis_length = np.linalg.norm(axis)
            if axis_length > 0:
                axis = axis / axis_length

            # Calculate angle with principal axes
            orientation = {}

            # Get structure principal axes
            coords = np.array([atom.get_coord() for atom in structure.get_atoms()])
            centered = coords - np.mean(coords, axis=0)
            cov = np.cov(centered.T)
            eigenvals, eigenvecs = np.linalg.eigh(cov)

            # Sort by eigenvalue
            idx = eigenvals.argsort()[::-1]
            eigenvals = eigenvals[idx]
            eigenvecs = eigenvecs[:, idx]

            # Calculate angles with principal axes
            for i, vec in enumerate(eigenvecs.T):
                angle = np.arccos(np.abs(np.dot(axis, vec)))
                orientation[f"angle_axis_{i+1}"] = float(np.degrees(angle))

            # Add overall orientation metrics
            orientation.update(
                {
                    "separation": float(axis_length),
                    "alignment_score": float(1.0 - np.min(orientation.values()) / 90.0),
                }
            )

            return orientation

        except Exception as e:
            self.logger.error(f"Error calculating domain orientation: {str(e)}")
            return {}

    def _analyze_global_domain_organization(
        self,
        domain_properties: Dict[str, Dict[str, Any]],
        inter_domain: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Analyze global domain organization.

        Args:
            domain_properties: Dictionary of domain properties
            inter_domain: Dictionary of inter-domain properties

        Returns:
            Dictionary of global organization properties
        """
        try:
            if not domain_properties or not inter_domain:
                return {}

            # Calculate global properties
            properties = {
                "num_domains": len(domain_properties),
                "avg_domain_size": float(np.mean([d["size"] for d in domain_properties.values()])),
                "size_variation": float(np.std([d["size"] for d in domain_properties.values()])),
                "avg_interface_size": float(np.mean([d["interface_size"] for d in inter_domain.values()])),
                "total_interface_area": sum(d["interface_size"] for d in inter_domain.values()),
                "domain_packing": self._calculate_domain_packing(
                    domain_properties,
                    inter_domain,
                ),
            }

            # Add domain graph properties
            G = nx.Graph()
            for name1, domain1 in domain_properties.items():
                G.add_node(name1, size=domain1["size"])

            for name, props in inter_domain.items():
                d1, d2 = name.split("_")
                G.add_edge(
                    d1,
                    d2,
                    weight=props["interface_size"],
                    distance=props["center_distance"],
                )

            properties.update(
                {
                    "graph_density": float(nx.density(G)),
                    "avg_path_length": float(nx.average_shortest_path_length(G)),
                    "modularity": float(self._calculate_modularity(G)),
                }
            )

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing global domain organization: {str(e)}")
            return {}

    def _calculate_domain_packing(
        self,
        domain_properties: Dict[str, Dict[str, Any]],
        inter_domain: Dict[str, Dict[str, Any]],
    ) -> float:
        """Calculate domain packing density.

        Args:
            domain_properties: Dictionary of domain properties
            inter_domain: Dictionary of inter-domain properties

        Returns:
            Packing density score (0-1)
        """
        try:
            if not domain_properties or not inter_domain:
                return 0.0

            # Calculate total domain volume (approximated by spheres)
            total_volume = sum(4 / 3 * np.pi * (d["radius_gyration"] ** 3) for d in domain_properties.values())

            # Calculate convex hull volume of domain centers
            centers = np.array([d["center"] for d in domain_properties.values()])
            if len(centers) >= 4:
                hull = ConvexHull(centers)
                hull_volume = hull.volume
            else:
                # Approximate for fewer points
                hull_volume = total_volume * 1.5

            # Calculate packing density
            if hull_volume > 0:
                return float(total_volume / hull_volume)
            return 0.0

        except Exception as e:
            self.logger.error(f"Error calculating domain packing: {str(e)}")
            return 0.0

    def _calculate_modularity(self, G: nx.Graph) -> float:
        """Calculate graph modularity using community structure.

        Args:
            G: NetworkX graph of domain organization

        Returns:
            Modularity score (0-1)
        """
        try:
            if not G or G.number_of_nodes() < 2:
                return 0.0

            # Detect communities
            communities = nx.community.louvain_communities(G)

            # Calculate modularity
            modularity = nx.community.modularity(G, communities)

            return max(0.0, float(modularity))

        except Exception as e:
            self.logger.error(f"Error calculating modularity: {str(e)}")
            return 0.0
