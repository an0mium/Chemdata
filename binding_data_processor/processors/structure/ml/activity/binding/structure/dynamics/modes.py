"""Normal mode analysis and collective motions."""

import logging
from typing import Dict, List, Optional, Any
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from scipy.sparse.linalg import eigsh
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)


class NormalModeAnalyzer:
    """Analyzes protein normal modes and collective motions."""

    def __init__(self):
        """Initialize normal mode analyzer."""
        self.logger = logging.getLogger(__name__)
        self.settings = {
            "cutoff": 12.0,  # Å, distance cutoff for elastic network
            "force_constant": 1.0,  # Spring constant for elastic network
            "num_modes": 10,  # Number of lowest frequency modes to calculate
            "mass_weighted": True,  # Whether to use mass-weighted Hessian
            "temperature": 300,  # K, temperature for thermal fluctuations
            "k_b": 0.0019872041,  # kcal/mol/K, Boltzmann constant
            "distance_power": 0,  # Power for distance-dependent force constant
            "min_freq_cutoff": 1e-6,  # Minimum frequency cutoff
            "rigid_body_modes": 6,  # Number of rigid body modes to exclude
        }

    def analyze_normal_modes(
        self,
        structure: Structure,
        num_modes: Optional[int] = None,
        mass_weighted: Optional[bool] = None,
    ) -> Dict[str, Any]:
        """Calculate and analyze normal modes.

        Args:
            structure: BioPython Structure object
            num_modes: Number of modes to calculate (default: from settings)
            mass_weighted: Whether to use mass weighting (default: from settings)

        Returns:
            Dictionary containing normal mode analysis results
        """
        try:
            # Get CA coordinates and masses
            coords, masses = self._get_ca_coords_and_masses(structure)
            if len(coords) < 3:
                return {}

            # Build Hessian matrix
            hessian = self._build_hessian(
                coords,
                masses if mass_weighted or self.settings["mass_weighted"] else None,
            )

            # Calculate normal modes
            n_modes = num_modes or self.settings["num_modes"]
            eigenvals, eigenvecs = self._calculate_modes(hessian, n_modes)

            # Calculate mode properties
            modes = []
            for i in range(len(eigenvals)):
                mode = {
                    "frequency": float(np.sqrt(abs(eigenvals[i]))),
                    "eigenvalue": float(eigenvals[i]),
                    "collectivity": float(self._calculate_collectivity(eigenvecs[:, i])),
                    "amplitude": float(1.0 / np.sqrt(abs(eigenvals[i]))),
                    "fluctuations": self._calculate_mode_fluctuations(
                        eigenvecs[:, i],
                        eigenvals[i],
                        coords,
                    ),
                    "energy": float(0.5 * eigenvals[i]),
                }

                # Calculate thermal fluctuations
                if eigenvals[i] > self.settings["min_freq_cutoff"]:
                    mode["thermal_fluctuation"] = float(np.sqrt(self.settings["k_b"] * self.settings["temperature"] / eigenvals[i]))
                else:
                    mode["thermal_fluctuation"] = 0.0

                modes.append(mode)

            # Calculate residue fluctuations
            fluctuations = self._calculate_residue_fluctuations(
                eigenvals,
                eigenvecs,
                coords,
            )

            # Calculate correlation matrix
            correlations = self._calculate_correlation_matrix(
                eigenvals,
                eigenvecs,
                coords,
            )

            # Calculate domain motions
            domains = self._analyze_domain_motions(
                modes,
                coords,
                correlations,
            )

            # Calculate additional properties
            properties = {
                "modes": modes,
                "fluctuations": fluctuations,
                "correlations": correlations,
                "domains": domains,
                "total_variance": float(np.sum(1.0 / np.abs(eigenvals))),
                "effective_modes": self._calculate_effective_modes(eigenvals),
                "collective_modes": self._analyze_collective_modes(eigenvecs),
            }

            # Add energy distribution
            properties["energy_distribution"] = self._calculate_energy_distribution(eigenvals)

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing normal modes: {str(e)}")
            return {}

    def _get_ca_coords_and_masses(
        self,
        structure: Structure,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get CA coordinates and residue masses.

        Args:
            structure: BioPython Structure object

        Returns:
            Tuple of (coordinates array, masses array)
        """
        try:
            coords = []
            masses = []

            for residue in structure.get_residues():
                if "CA" in residue:
                    coords.append(residue["CA"].get_coord())
                    masses.append(self._get_residue_mass(residue))

            return np.array(coords), np.array(masses)

        except Exception as e:
            self.logger.error(f"Error getting coordinates and masses: {str(e)}")
            return np.array([]), np.array([])

    def _get_residue_mass(self, residue: Residue) -> float:
        """Calculate residue mass.

        Args:
            residue: BioPython Residue object

        Returns:
            Mass in atomic mass units
        """
        try:
            # Accurate masses for standard amino acids
            masses = {
                "ALA": 71.08,
                "ARG": 156.19,
                "ASN": 114.10,
                "ASP": 115.09,
                "CYS": 103.14,
                "GLN": 128.13,
                "GLU": 129.12,
                "GLY": 57.05,
                "HIS": 137.14,
                "ILE": 113.16,
                "LEU": 113.16,
                "LYS": 128.17,
                "MET": 131.19,
                "PHE": 147.18,
                "PRO": 97.12,
                "SER": 87.08,
                "THR": 101.11,
                "TRP": 186.21,
                "TYR": 163.18,
                "VAL": 99.13,
            }
            return masses.get(residue.get_resname(), 110.0)  # Average mass if unknown

        except Exception as e:
            self.logger.error(f"Error getting residue mass: {str(e)}")
            return 110.0

    def _build_hessian(
        self,
        coords: np.ndarray,
        masses: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Build Hessian matrix for elastic network model.

        Args:
            coords: Array of atomic coordinates
            masses: Optional array of atomic masses for mass weighting

        Returns:
            Hessian matrix
        """
        try:
            n_atoms = len(coords)
            hessian = np.zeros((3 * n_atoms, 3 * n_atoms))

            # Calculate pairwise distances
            distances = cdist(coords, coords)

            # Build Hessian using enhanced spring model
            for i in range(n_atoms):
                for j in range(n_atoms):
                    if i != j and distances[i, j] < self.settings["cutoff"]:
                        # Calculate direction cosines
                        diff = coords[i] - coords[j]
                        dist = distances[i, j]
                        direction = diff / dist

                        # Calculate distance-dependent force constant
                        k = self.settings["force_constant"]
                        if self.settings["distance_power"] != 0:
                            k *= (self.settings["cutoff"] / dist) ** self.settings["distance_power"]

                        # Apply mass weighting if requested
                        if masses is not None:
                            k /= np.sqrt(masses[i] * masses[j])

                        # 3x3 super-element
                        super_element = k * np.outer(direction, direction)

                        # Add to Hessian
                        i3, j3 = i * 3, j * 3
                        hessian[i3 : i3 + 3, j3 : j3 + 3] = -super_element
                        hessian[i3 : i3 + 3, i3 : i3 + 3] += super_element

            return hessian

        except Exception as e:
            self.logger.error(f"Error building Hessian: {str(e)}")
            return np.array([])

    def _calculate_modes(
        self,
        hessian: np.ndarray,
        num_modes: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate normal modes using sparse eigenvalue decomposition.

        Args:
            hessian: Hessian matrix
            num_modes: Number of modes to calculate

        Returns:
            Tuple of (eigenvalues array, eigenvectors array)
        """
        try:
            # Remove rigid body modes
            n_rigid = self.settings["rigid_body_modes"]
            n_modes = min(num_modes + n_rigid, len(hessian) - 1)

            # Calculate lowest frequency modes
            eigenvals, eigenvecs = eigsh(hessian, k=n_modes, which="SA")

            # Remove rigid body modes
            eigenvals = eigenvals[n_rigid:]
            eigenvecs = eigenvecs[:, n_rigid:]

            # Sort by frequency
            idx = np.argsort(np.abs(eigenvals))
            eigenvals = eigenvals[idx]
            eigenvecs = eigenvecs[:, idx]

            return eigenvals, eigenvecs

        except Exception as e:
            self.logger.error(f"Error calculating modes: {str(e)}")
            return np.array([]), np.array([])

    def _calculate_collectivity(self, eigenvec: np.ndarray) -> float:
        """Calculate mode collectivity (participation ratio).

        Args:
            eigenvec: Mode eigenvector

        Returns:
            Collectivity between 0 and 1
        """
        try:
            # Reshape to get per-residue components
            reshaped = eigenvec.reshape(-1, 3)
            squared = np.sum(reshaped * reshaped, axis=1)
            total = np.sum(squared)

            if total > 0:
                # Calculate participation ratio
                squared_total = np.sum((squared / total) ** 2)
                return 1.0 / (len(reshaped) * squared_total)
            return 0.0

        except Exception as e:
            self.logger.error(f"Error calculating collectivity: {str(e)}")
            return 0.0

    def _calculate_mode_fluctuations(
        self,
        eigenvec: np.ndarray,
        eigenval: float,
        coords: np.ndarray,
    ) -> Dict[int, float]:
        """Calculate residue fluctuations for a mode.

        Args:
            eigenvec: Mode eigenvector
            eigenval: Mode eigenvalue
            coords: Atomic coordinates

        Returns:
            Dictionary mapping residue indices to fluctuations
        """
        try:
            # Calculate thermal factor
            if eigenval > self.settings["min_freq_cutoff"]:
                factor = self.settings["k_b"] * self.settings["temperature"] / eigenval
            else:
                factor = 1.0

            # Calculate per-residue fluctuations
            reshaped = eigenvec.reshape(-1, 3)
            fluctuations = {}

            for i in range(len(coords)):
                fluct = factor * np.sum(reshaped[i] * reshaped[i])
                fluctuations[i] = float(fluct)

            return fluctuations

        except Exception as e:
            self.logger.error(f"Error calculating mode fluctuations: {str(e)}")
            return {}

    def _calculate_residue_fluctuations(
        self,
        eigenvals: np.ndarray,
        eigenvecs: np.ndarray,
        coords: np.ndarray,
    ) -> Dict[int, float]:
        """Calculate total residue fluctuations from all modes.

        Args:
            eigenvals: Array of eigenvalues
            eigenvecs: Array of eigenvectors
            coords: Atomic coordinates

        Returns:
            Dictionary mapping residue indices to fluctuations
        """
        try:
            fluctuations = {}

            # Calculate thermal factors
            factors = np.where(
                eigenvals > self.settings["min_freq_cutoff"],
                self.settings["k_b"] * self.settings["temperature"] / eigenvals,
                1.0,
            )

            # Sum contributions from all modes
            for i in range(len(coords)):
                total_fluct = 0.0
                for j in range(len(eigenvals)):
                    mode = eigenvecs[i * 3 : i * 3 + 3, j]
                    total_fluct += factors[j] * np.sum(mode * mode)
                fluctuations[i] = float(total_fluct)

            return fluctuations

        except Exception as e:
            self.logger.error(f"Error calculating residue fluctuations: {str(e)}")
            return {}

    def _calculate_correlation_matrix(
        self,
        eigenvals: np.ndarray,
        eigenvecs: np.ndarray,
        coords: np.ndarray,
    ) -> np.ndarray:
        """Calculate residue correlation matrix.

        Args:
            eigenvals: Array of eigenvalues
            eigenvecs: Array of eigenvectors
            coords: Atomic coordinates

        Returns:
            Correlation matrix
        """
        try:
            n_res = len(coords)
            correlations = np.zeros((n_res, n_res))

            # Calculate thermal factors
            factors = np.where(
                eigenvals > self.settings["min_freq_cutoff"],
                self.settings["k_b"] * self.settings["temperature"] / eigenvals,
                1.0,
            )

            # Calculate correlations
            for i in range(n_res):
                for j in range(i, n_res):
                    corr = 0.0
                    for k in range(len(eigenvals)):
                        mode_i = eigenvecs[i * 3 : i * 3 + 3, k]
                        mode_j = eigenvecs[j * 3 : j * 3 + 3, k]
                        corr += factors[k] * np.sum(mode_i * mode_j)
                    correlations[i, j] = corr
                    correlations[j, i] = corr

            # Normalize
            norms = np.sqrt(np.diag(correlations))
            correlations /= np.outer(norms, norms)

            return correlations

        except Exception as e:
            self.logger.error(f"Error calculating correlation matrix: {str(e)}")
            return np.array([])

    def _analyze_domain_motions(
        self,
        modes: List[Dict[str, Any]],
        coords: np.ndarray,
        correlations: np.ndarray,
    ) -> Dict[str, Any]:
        """Analyze domain motions from normal modes.

        Args:
            modes: List of mode dictionaries
            coords: Atomic coordinates
            correlations: Residue correlation matrix

        Returns:
            Dictionary of domain motion analysis
        """
        try:
            # Identify domains from correlation matrix
            domains = self._identify_domains(correlations)

            # Analyze domain motions for each mode
            domain_motions = []
            for i, mode in enumerate(modes):
                motion = self._analyze_mode_domain_motion(
                    mode,
                    domains,
                    coords,
                )
                if motion:
                    domain_motions.append(motion)

            # Identify hinge residues
            hinges = self._identify_hinge_residues(domains, correlations)

            # Calculate domain interface properties
            interfaces = self._analyze_domain_interfaces(domains, coords)

            return {
                "domains": domains,
                "motions": domain_motions,
                "hinge_residues": hinges,
                "interfaces": interfaces,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing domain motions: {str(e)}")
            return {}

    def _identify_domains(self, correlations: np.ndarray) -> List[List[int]]:
        """Identify protein domains from correlation matrix.

        Args:
            correlations: Residue correlation matrix

        Returns:
            List of domain residue lists
        """
        try:
            import networkx as nx
            from sklearn.cluster import SpectralClustering

            # Create graph from correlation matrix
            G = nx.Graph()
            n_res = len(correlations)
            for i in range(n_res):
                for j in range(i + 1, n_res):
                    if abs(correlations[i, j]) > 0.5:  # Correlation threshold
                        G.add_edge(i, j, weight=abs(correlations[i, j]))

            # Estimate optimal number of domains
            n_domains = max(2, int(np.sqrt(n_res / 50)))  # Heuristic

            # Cluster residues
            clustering = SpectralClustering(
                n_clusters=n_domains,
                affinity="precomputed",
                random_state=42,
            )
            adj_matrix = nx.adjacency_matrix(G).todense()
            labels = clustering.fit_predict(adj_matrix)

            # Group residues into domains
            domains = []
            for i in range(n_domains):
                domain = list(np.where(labels == i)[0])
                if len(domain) >= 10:  # Minimum domain size
                    domains.append(sorted(domain))

            return domains

        except Exception as e:
            self.logger.error(f"Error identifying domains: {str(e)}")
            return []

    def _analyze_mode_domain_motion(
        self,
        mode: Dict[str, Any],
        domains: List[List[int]],
        coords: np.ndarray,
    ) -> Optional[Dict[str, Any]]:
        """Analyze domain motion for a specific mode.

        Args:
            mode: Mode dictionary
            domains: List of domain residue lists
            coords: Atomic coordinates

        Returns:
            Dictionary describing domain motion or None
        """
        try:
            if not domains:
                return None

            motion = {
                "frequency": mode["frequency"],
                "collectivity": mode["collectivity"],
                "domain_displacements": [],
            }

            # Calculate domain centers
            centers = []
            for domain in domains:
                center = np.mean(coords[domain], axis=0)
                centers.append(center)

            # Analyze relative motions between domains
            for i, domain1 in enumerate(domains):
                for j, domain2 in enumerate(domains[i + 1 :], i + 1):
                    # Calculate displacement vectors
                    vec1 = np.mean(mode["fluctuations"][k] for k in domain1)
                    vec2 = np.mean(mode["fluctuations"][k] for k in domain2)
                    relative_motion = vec2 - vec1

                    # Calculate motion parameters
                    displacement = {
                        "domains": (i + 1, j + 1),
                        "magnitude": float(np.linalg.norm(relative_motion)),
                        "direction": relative_motion / np.linalg.norm(relative_motion),
                        "distance": float(np.linalg.norm(centers[j] - centers[i])),
                    }
                    motion["domain_displacements"].append(displacement)

            return motion

        except Exception as e:
            self.logger.error(f"Error analyzing mode domain motion: {str(e)}")
            return None

    def _identify_hinge_residues(
        self,
        domains: List[List[int]],
        correlations: np.ndarray,
    ) -> List[int]:
        """Identify hinge residues between domains.

        Args:
            domains: List of domain residue lists
            correlations: Residue correlation matrix

        Returns:
            List of hinge residue indices
        """
        try:
            hinges = []

            # Calculate domain boundaries
            boundaries = set()
            for domain in domains:
                boundaries.update([min(domain) - 1, max(domain) + 1])
            boundaries = sorted(list(boundaries))

            # Check residues near boundaries
            for res_id in boundaries:
                if 0 <= res_id < len(correlations):
                    # Calculate correlation difference across residue
                    left_corr = np.mean(correlations[res_id, max(0, res_id - 3) : res_id])
                    right_corr = np.mean(correlations[res_id, res_id + 1 : res_id + 4])
                    diff = abs(right_corr - left_corr)

                    if diff > 0.5:  # Correlation difference threshold
                        hinges.append(res_id)

            return sorted(hinges)

        except Exception as e:
            self.logger.error(f"Error identifying hinge residues: {str(e)}")
            return []

    def _analyze_domain_interfaces(
        self,
        domains: List[List[int]],
        coords: np.ndarray,
    ) -> List[Dict[str, Any]]:
        """Analyze interfaces between domains.

        Args:
            domains: List of domain residue lists
            coords: Atomic coordinates

        Returns:
            List of interface properties
        """
        try:
            interfaces = []

            for i, domain1 in enumerate(domains):
                for j, domain2 in enumerate(domains[i + 1 :], i + 1):
                    # Calculate interface residues
                    interface_residues = self._get_interface_residues(
                        domain1,
                        domain2,
                        coords,
                    )

                    if interface_residues:
                        interface = {
                            "domains": (i + 1, j + 1),
                            "residues": interface_residues,
                            "size": len(interface_residues),
                            "contacts": self._count_interface_contacts(
                                interface_residues,
                                coords,
                            ),
                        }
                        interfaces.append(interface)

            return interfaces

        except Exception as e:
            self.logger.error(f"Error analyzing domain interfaces: {str(e)}")
            return []

    def _get_interface_residues(
        self,
        domain1: List[int],
        domain2: List[int],
        coords: np.ndarray,
    ) -> List[int]:
        """Get residues at domain interface.

        Args:
            domain1: First domain residue list
            domain2: Second domain residue list
            coords: Atomic coordinates

        Returns:
            List of interface residue indices
        """
        try:
            interface = []
            cutoff = 8.0  # Distance cutoff for interface residues

            # Calculate distances between domain residues
            for res1 in domain1:
                for res2 in domain2:
                    dist = np.linalg.norm(coords[res1] - coords[res2])
                    if dist < cutoff:
                        interface.extend([res1, res2])

            return sorted(list(set(interface)))

        except Exception as e:
            self.logger.error(f"Error getting interface residues: {str(e)}")
            return []

    def _count_interface_contacts(
        self,
        interface_residues: List[int],
        coords: np.ndarray,
    ) -> int:
        """Count contacts between interface residues.

        Args:
            interface_residues: List of interface residue indices
            coords: Atomic coordinates

        Returns:
            Number of contacts
        """
        try:
            contacts = 0
            cutoff = 8.0  # Distance cutoff for contacts

            # Calculate pairwise distances
            for i, res1 in enumerate(interface_residues):
                for res2 in interface_residues[i + 1 :]:
                    if np.linalg.norm(coords[res1] - coords[res2]) < cutoff:
                        contacts += 1

            return contacts

        except Exception as e:
            self.logger.error(f"Error counting interface contacts: {str(e)}")
            return 0

    def _calculate_effective_modes(self, eigenvals: np.ndarray) -> float:
        """Calculate number of effective modes.

        Args:
            eigenvals: Array of eigenvalues

        Returns:
            Number of effective modes
        """
        try:
            if len(eigenvals) == 0:
                return 0.0

            # Calculate participation ratio of eigenvalues
            total = np.sum(1.0 / np.abs(eigenvals))
            squared_total = np.sum((1.0 / np.abs(eigenvals)) ** 2)

            return float(total**2 / squared_total)

        except Exception as e:
            self.logger.error(f"Error calculating effective modes: {str(e)}")
            return 0.0

    def _analyze_collective_modes(self, eigenvecs: np.ndarray) -> Dict[str, Any]:
        """Analyze collective motions from normal modes.

        Args:
            eigenvecs: Mode eigenvectors

        Returns:
            Dictionary of collective motion properties
        """
        try:
            n_modes = min(10, eigenvecs.shape[1] - self.settings["rigid_body_modes"])
            collectivity = []

            for i in range(self.settings["rigid_body_modes"], self.settings["rigid_body_modes"] + n_modes):
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
                "max_collectivity": float(np.max(collectivity)) if collectivity else 0.0,
                "min_collectivity": float(np.min(collectivity)) if collectivity else 0.0,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing collective modes: {str(e)}")
            return {}

    def _calculate_energy_distribution(self, eigenvals: np.ndarray) -> Dict[str, float]:
        """Calculate energy distribution across modes.

        Args:
            eigenvals: Array of eigenvalues

        Returns:
            Dictionary of energy distribution properties
        """
        try:
            if len(eigenvals) == 0:
                return {}

            # Calculate mode energies
            energies = 0.5 * eigenvals
            total_energy = np.sum(energies)

            # Calculate energy fractions
            energy_fractions = energies / total_energy

            return {
                "total_energy": float(total_energy),
                "mean_energy": float(np.mean(energies)),
                "energy_std": float(np.std(energies)),
                "max_energy_fraction": float(np.max(energy_fractions)),
                "min_energy_fraction": float(np.min(energy_fractions)),
            }

        except Exception as e:
            self.logger.error(f"Error calculating energy distribution: {str(e)}")
            return {}
