"""Protein dynamics analyzer module."""

import logging
from typing import Dict, List, Any, Optional, Union
import numpy as np
from Bio.PDB import Structure, Residue, Chain, Model

from binding_data_processor.core.config import ProteinAnalysisConfig
from ..utils import get_ca_atoms

logger = logging.getLogger(__name__)


class DynamicsAnalyzer:
    """Analyzer for protein dynamics including flexibility and normal modes."""

    def __init__(self, structure: Structure, config: Optional[ProteinAnalysisConfig] = None):
        """Initialize dynamics analyzer.

        Args:
            structure: BioPython Structure object to analyze
            config: Optional configuration object
        """
        self.structure = structure
        self.config = config or ProteinAnalysisConfig()
        self._cache = {}

    def get_coordinates(self, structure: Optional[Structure] = None) -> np.ndarray:
        """Get coordinates of all CA atoms in structure.

        Args:
            structure: Optional BioPython Structure object (uses self.structure if not provided)

        Returns:
            Array of shape (n_atoms, 3) containing xyz coordinates
        """
        try:
            if structure is None:
                structure = self.structure

            # Get CA atoms
            ca_atoms = get_ca_atoms(structure)
            if not ca_atoms:
                return np.array([])

            # Extract coordinates
            coords = []
            for atom in ca_atoms:
                coords.append(atom.get_coord())

            return np.array(coords)

        except Exception as e:
            logger.error(f"Error getting coordinates: {str(e)}")
            return np.array([])

    def analyze_flexibility(self) -> Dict[str, Any]:
        """Analyze protein flexibility using B-factors and other metrics.

        Returns:
            Dictionary containing flexibility analysis results including:
            - Mean, std, min, max B-factors
            - Per-residue B-factors
            - Flexible and rigid regions
            - Flexibility statistics
        """
        try:
            # Get B-factors for all atoms
            b_factors = []
            residue_b_factors = {}

            for residue in self.structure.get_residues():
                res_b_factors = []
                for atom in residue:
                    b_factor = atom.get_bfactor()
                    b_factors.append(b_factor)
                    res_b_factors.append(b_factor)

                # Calculate average B-factor per residue
                if res_b_factors:
                    res_id = residue.get_id()[1]
                    residue_b_factors[res_id] = float(np.mean(res_b_factors))

            # Calculate overall statistics
            mean_b = float(np.mean(b_factors))
            std_b = float(np.std(b_factors))

            # Identify highly flexible regions (B-factor > mean + std)
            flexible_residues = [res_id for res_id, b_factor in residue_b_factors.items() if b_factor > (mean_b + std_b)]

            # Identify rigid regions (B-factor < mean - std)
            rigid_residues = [res_id for res_id, b_factor in residue_b_factors.items() if b_factor < (mean_b - std_b)]

            return {
                "mean_b_factor": mean_b,
                "std_b_factor": std_b,
                "min_b_factor": float(min(b_factors)),
                "max_b_factor": float(max(b_factors)),
                "residue_b_factors": residue_b_factors,
                "flexible_residues": flexible_residues,
                "rigid_residues": rigid_residues,
                "n_flexible": len(flexible_residues),
                "n_rigid": len(rigid_residues),
            }

        except Exception as e:
            logger.error(f"Error analyzing flexibility: {str(e)}")
            return {}

    def analyze_normal_modes(self, n_modes: int = 10) -> Dict[str, Any]:
        """Analyze protein normal modes using elastic network model.

        Args:
            n_modes: Number of normal modes to calculate

        Returns:
            Dictionary containing normal modes analysis results including:
            - Mode frequencies
            - Mode vectors
            - Residue fluctuations
            - Mode statistics
        """
        try:
            # Get CA atoms for elastic network model
            ca_atoms = []
            ca_coords = []

            for residue in self.structure.get_residues():
                if "CA" in residue:
                    ca_atoms.append(residue["CA"])
                    ca_coords.append(residue["CA"].get_coord())

            ca_coords = np.array(ca_coords)

            # Build contact matrix (simple distance-based)
            n_res = len(ca_atoms)
            contacts = np.zeros((n_res, n_res))

            for i in range(n_res):
                for j in range(i + 1, n_res):
                    dist = np.linalg.norm(ca_coords[i] - ca_coords[j])
                    if dist < 12.0:  # Cutoff for contacts
                        contacts[i, j] = contacts[j, i] = 1

            # Calculate normal modes (simplified)
            # In practice, you would use a proper ENM implementation
            try:
                eigenvals, eigenvecs = np.linalg.eigh(contacts)
                # Sort by eigenvalue ascending (lowest frequency modes first)
                idx = eigenvals.argsort()
                eigenvals = eigenvals[idx]
                eigenvecs = eigenvecs[:, idx]

                # Get lowest frequency modes (excluding first 6 rigid body motions)
                modes = eigenvecs[:, 6 : 6 + n_modes]
                freqs = eigenvals[6 : 6 + n_modes]

                # Calculate residue fluctuations from modes
                flucts = np.sum(modes * modes, axis=1)

                # Cache correlation matrix for later use
                self._cache["correlation_matrix"] = np.dot(modes, modes.T)

                return {
                    "frequencies": [float(f) for f in freqs],
                    "modes": modes.tolist(),
                    "residue_fluctuations": [float(f) for f in flucts],
                    "n_modes": n_modes,
                }

            except np.linalg.LinAlgError:
                logger.warning("Could not calculate normal modes - using fallback")
                return {"frequencies": [], "modes": [], "residue_fluctuations": [], "n_modes": 0}

        except Exception as e:
            logger.error(f"Error analyzing normal modes: {str(e)}")
            return {}

    def analyze_site_dynamics(
        self,
        site_residues: List[Residue],
        structure: Structure,
        include_domain_context: bool = True,
        include_correlations: bool = True,
    ) -> Dict[str, Any]:
        """Analyze dynamics of specific binding site residues.

        Args:
            site_residues: List of residues in binding site
            structure: Full structure for context
            include_domain_context: Whether to include domain motion analysis
            include_correlations: Whether to include correlation analysis

        Returns:
            Dictionary containing site dynamics analysis including:
            - Site-specific dynamics
            - Local flexibility
            - Contact network
            - Domain context (optional)
            - Correlation analysis (optional)
        """
        try:
            results = {}

            # Basic site dynamics
            site_ids = [res.get_id()[1] for res in site_residues]

            # Get flexibility for site residues
            flexibility = self.analyze_flexibility()
            site_flexibility = {res_id: flexibility["residue_b_factors"].get(res_id, 0.0) for res_id in site_ids}
            results["flexibility"] = site_flexibility

            # Get contacts for site residues
            contacts = self.analyze_contacts()
            site_contacts = {res_id: contacts["residue_contacts"].get(res_id, []) for res_id in site_ids}
            results["contacts"] = site_contacts

            # Add enhanced site-specific analysis
            if include_domain_context:
                results["domain_context"] = self._analyze_domain_context(site_residues)

            if include_correlations:
                results["correlations"] = self._analyze_site_correlations(site_residues)

            return results

        except Exception as e:
            logger.error(f"Error analyzing site dynamics: {str(e)}")
            return {}

    def analyze_contacts(self, cutoff: float = 8.0) -> Dict[str, Any]:
        """Analyze residue-residue contacts.

        Args:
            cutoff: Distance cutoff for contacts in Angstroms

        Returns:
            Dictionary containing contact analysis results including:
            - Contact list
            - Per-residue contacts
            - Contact statistics
            - Network metrics
        """
        try:
            contacts = []
            residue_contacts = {}

            residues = list(self.structure.get_residues())

            for i, res1 in enumerate(residues):
                res1_id = res1.get_id()[1]
                residue_contacts[res1_id] = []

                for res2 in residues[i + 1 :]:
                    res2_id = res2.get_id()[1]

                    # Calculate minimum distance between residues
                    min_dist = float("inf")
                    for atom1 in res1:
                        for atom2 in res2:
                            dist = atom1 - atom2
                            min_dist = min(min_dist, dist)

                    if min_dist < cutoff:
                        contact = {"residue1": res1_id, "residue2": res2_id, "distance": float(min_dist)}
                        contacts.append(contact)
                        residue_contacts[res1_id].append(res2_id)
                        if res2_id not in residue_contacts:
                            residue_contacts[res2_id] = []
                        residue_contacts[res2_id].append(res1_id)

            # Calculate contact statistics
            contacts_per_res = [len(c) for c in residue_contacts.values()]

            return {
                "n_contacts": len(contacts),
                "contacts": contacts,
                "residue_contacts": residue_contacts,
                "mean_contacts": float(np.mean(contacts_per_res)),
                "max_contacts": float(max(contacts_per_res)),
                "min_contacts": float(min(contacts_per_res)),
            }

        except Exception as e:
            logger.error(f"Error analyzing contacts: {str(e)}")
            return {}

    def _analyze_domain_context(self, site_residues: List[Residue]) -> Dict[str, Any]:
        """Analyze domain context of binding site residues.

        Args:
            site_residues: List of residues in binding site

        Returns:
            Dictionary containing domain context analysis including:
            - Domain assignments
            - Domain motions
            - Interface contacts
        """
        try:
            # Get domain assignments
            domains = self._get_domain_assignments()

            # Map residues to domains
            residue_domains = {}
            for res in site_residues:
                domain = self._get_residue_domain(res, domains)
                if domain:
                    residue_domains[res.get_id()[1]] = domain

            # Analyze domain motions
            domain_motions = self._analyze_domain_motions(domains)

            return {
                "residue_domains": residue_domains,
                "domain_motions": domain_motions,
            }

        except Exception as e:
            logger.error(f"Error analyzing domain context: {str(e)}")
            return {}

    def _analyze_site_correlations(self, site_residues: List[Residue]) -> Dict[str, Any]:
        """Analyze correlations between binding site and overall dynamics.

        Args:
            site_residues: List of residues in binding site

        Returns:
            Dictionary containing correlation analysis including:
            - Site-specific correlations
            - Correlation matrix
            - Correlation statistics
        """
        try:
            # Get correlation matrix
            correlation_matrix = self._get_correlation_matrix()
            if correlation_matrix is None:
                # Need to calculate normal modes first
                self.analyze_normal_modes()
                correlation_matrix = self._get_correlation_matrix()
                if correlation_matrix is None:
                    return {}

            # Extract site-specific correlations
            site_correlations = {}
            for res in site_residues:
                res_id = res.get_id()[1]
                correlations = correlation_matrix[res_id]
                site_correlations[res_id] = {
                    "mean": float(np.mean(correlations)),
                    "std": float(np.std(correlations)),
                    "max": float(np.max(correlations)),
                    "min": float(np.min(correlations)),
                }

            return {
                "site_correlations": site_correlations,
                "correlation_matrix": correlation_matrix.tolist(),
            }

        except Exception as e:
            logger.error(f"Error analyzing site correlations: {str(e)}")
            return {}

    def _get_correlation_matrix(self) -> Optional[np.ndarray]:
        """Get correlation matrix from cached normal modes analysis.

        Returns:
            Correlation matrix if available, None otherwise
        """
        return self._cache.get("correlation_matrix")

    def _get_domain_assignments(self) -> Dict[str, List[int]]:
        """Get domain assignments for protein structure.

        Returns:
            Dictionary mapping domain IDs to lists of residue numbers
        """
        # TODO: Implement domain assignment algorithm
        # For now, return empty domains
        return {}

    def _get_residue_domain(self, residue: Residue, domains: Dict[str, List[int]]) -> Optional[str]:
        """Get domain assignment for a residue.

        Args:
            residue: Residue to get domain for
            domains: Domain assignments

        Returns:
            Domain ID if residue is in a domain, None otherwise
        """
        res_id = residue.get_id()[1]
        for domain_id, domain_residues in domains.items():
            if res_id in domain_residues:
                return domain_id
        return None

    def _analyze_domain_motions(self, domains: Dict[str, List[int]]) -> Dict[str, Any]:
        """Analyze domain motions from normal modes.

        Args:
            domains: Domain assignments

        Returns:
            Dictionary containing domain motion analysis
        """
        # TODO: Implement domain motion analysis
        # For now, return empty results
        return {}
