"""Domain motion analysis functionality."""

import logging
from typing import Dict, List, Optional, Any
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue

logger = logging.getLogger(__name__)


class DomainAnalyzer:
    """Analyzes protein domain motions and interactions."""

    def __init__(self):
        """Initialize domain analyzer."""
        self.logger = logging.getLogger(__name__)

    def analyze_domain_motions(
        self,
        structure: Structure,
        contact_data: Optional[Dict[str, Any]] = None,
        flexibility_data: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Analyze domain motions and interactions.

        Args:
            structure: BioPython Structure object
            contact_data: Optional contact network data
            flexibility_data: Optional flexibility analysis data

        Returns:
            Dictionary containing domain analysis results
        """
        try:
            # Identify domains using contact patterns
            domains = self._identify_domains(structure, contact_data)
            if not domains:
                return {}

            # Analyze domain interfaces
            interfaces = self._analyze_interfaces(domains, structure)

            # Analyze domain motions if flexibility data available
            motions = {}
            if flexibility_data:
                motions = self._analyze_domain_motions(domains, flexibility_data)

            return {
                "domains": domains,
                "interfaces": interfaces,
                "motions": motions,
                "metrics": self._calculate_domain_metrics(domains, interfaces, motions),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing domain motions: {str(e)}")
            return {}

    def analyze_site_domain_context(
        self,
        site_residues: List[Residue],
        structure: Structure,
    ) -> Dict[str, Any]:
        """Analyze domain context of binding site residues.

        Args:
            site_residues: List of binding site residues
            structure: Full structure for context

        Returns:
            Dictionary containing site domain analysis
        """
        try:
            # Get domain assignments
            domains = self._identify_domains(structure)
            if not domains:
                return {}

            # Map residues to domains
            site_domains = {}
            for residue in site_residues:
                res_id = residue.get_id()[1]
                for domain_id, domain_residues in domains.items():
                    if res_id in domain_residues:
                        site_domains[res_id] = domain_id

            # Analyze domain distribution
            domain_stats = self._analyze_domain_distribution(site_domains, domains)

            # Analyze interface proximity
            interfaces = self._analyze_interfaces(domains, structure)
            interface_proximity = self._analyze_interface_proximity(
                site_residues,
                interfaces,
                structure,
            )

            return {
                "site_domains": site_domains,
                "domain_stats": domain_stats,
                "interface_proximity": interface_proximity,
                "metrics": {
                    "num_domains": len(set(site_domains.values())),
                    "interface_residues": len(interface_proximity["interface_residues"]),
                    "max_interface_distance": interface_proximity["max_distance"],
                },
            }

        except Exception as e:
            self.logger.error(f"Error analyzing site domain context: {str(e)}")
            return {}

    def _identify_domains(
        self,
        structure: Structure,
        contact_data: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, List[int]]:
        """Identify protein domains using contact patterns.

        Args:
            structure: BioPython Structure object
            contact_data: Optional contact network data

        Returns:
            Dictionary mapping domain IDs to lists of residue numbers
        """
        try:
            # Use contact data if available, otherwise calculate contacts
            if not contact_data:
                from .contacts import ContactAnalyzer

                contact_analyzer = ContactAnalyzer()
                contact_data = contact_analyzer.analyze_contacts(structure)

            network = contact_data.get("contact_network", {})
            if not network:
                return {}

            # Convert network to adjacency matrix
            residues = sorted(network.keys())
            n_res = len(residues)
            adj_matrix = np.zeros((n_res, n_res))
            res_to_idx = {res: i for i, res in enumerate(residues)}

            for res1, neighbors in network.items():
                for res2 in neighbors:
                    i, j = res_to_idx[res1], res_to_idx[res2]
                    adj_matrix[i, j] = adj_matrix[j, i] = 1

            # Use spectral clustering to identify domains
            from sklearn.cluster import SpectralClustering

            n_clusters = max(2, int(np.sqrt(n_res) / 2))
            clustering = SpectralClustering(
                n_clusters=n_clusters,
                affinity="precomputed",
                random_state=42,
            )
            labels = clustering.fit_predict(adj_matrix)

            # Convert clusters to domain assignments
            domains = {}
            for i, label in enumerate(labels):
                domain_id = f"domain_{label + 1}"
                if domain_id not in domains:
                    domains[domain_id] = []
                domains[domain_id].append(residues[i])

            return domains

        except Exception as e:
            self.logger.error(f"Error identifying domains: {str(e)}")
            return {}

    def _analyze_interfaces(
        self,
        domains: Dict[str, List[int]],
        structure: Structure,
    ) -> Dict[str, Any]:
        """Analyze domain interfaces.

        Args:
            domains: Dictionary mapping domain IDs to residue lists
            structure: BioPython Structure object

        Returns:
            Dictionary containing interface analysis
        """
        try:
            interfaces = {}
            structure_residues = {res.get_id()[1]: res for res in structure.get_residues()}

            for domain1 in domains:
                for domain2 in domains:
                    if domain1 >= domain2:
                        continue

                    interface_key = f"{domain1}_{domain2}"
                    interface_residues = set()

                    # Find residues in contact between domains
                    for res1_id in domains[domain1]:
                        if res1_id not in structure_residues:
                            continue
                        res1 = structure_residues[res1_id]

                        for res2_id in domains[domain2]:
                            if res2_id not in structure_residues:
                                continue
                            res2 = structure_residues[res2_id]

                            # Check if residues are in contact
                            min_distance = float("inf")
                            for atom1 in res1:
                                for atom2 in res2:
                                    distance = atom1 - atom2
                                    min_distance = min(min_distance, distance)

                            if min_distance < 8.0:  # Å
                                interface_residues.add(res1_id)
                                interface_residues.add(res2_id)

                    if interface_residues:
                        interfaces[interface_key] = {
                            "residues": sorted(interface_residues),
                            "size": len(interface_residues),
                        }

            return interfaces

        except Exception as e:
            self.logger.error(f"Error analyzing interfaces: {str(e)}")
            return {}

    def _analyze_domain_motions(
        self,
        domains: Dict[str, List[int]],
        flexibility_data: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Analyze domain motions using flexibility data.

        Args:
            domains: Dictionary mapping domain IDs to residue lists
            flexibility_data: Flexibility analysis data

        Returns:
            Dictionary containing domain motion analysis
        """
        try:
            motions = {}
            residue_flexibility = flexibility_data.get("residue_flexibility", {})

            for domain_id, residues in domains.items():
                # Calculate domain flexibility statistics
                domain_flex = [residue_flexibility.get(res, 0.0) for res in residues]

                if domain_flex:
                    motions[domain_id] = {
                        "mean_flexibility": float(np.mean(domain_flex)),
                        "flexibility_std": float(np.std(domain_flex)),
                        "max_flexibility": float(np.max(domain_flex)),
                        "min_flexibility": float(np.min(domain_flex)),
                    }

            return motions

        except Exception as e:
            self.logger.error(f"Error analyzing domain motions: {str(e)}")
            return {}

    def _analyze_domain_distribution(
        self,
        site_domains: Dict[int, str],
        domains: Dict[str, List[int]],
    ) -> Dict[str, Any]:
        """Analyze distribution of binding site residues across domains.

        Args:
            site_domains: Dictionary mapping residue IDs to domain IDs
            domains: Dictionary mapping domain IDs to residue lists

        Returns:
            Dictionary containing domain distribution analysis
        """
        try:
            # Count residues per domain
            domain_counts = {}
            for domain_id in site_domains.values():
                if domain_id not in domain_counts:
                    domain_counts[domain_id] = 0
                domain_counts[domain_id] += 1

            # Calculate domain fractions
            total_residues = len(site_domains)
            domain_fractions = {domain_id: count / total_residues for domain_id, count in domain_counts.items()}

            # Identify primary domain
            primary_domain = max(
                domain_counts.items(),
                key=lambda x: x[1],
            )[0]

            return {
                "counts": domain_counts,
                "fractions": domain_fractions,
                "primary_domain": primary_domain,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing domain distribution: {str(e)}")
            return {}

    def _analyze_interface_proximity(
        self,
        site_residues: List[Residue],
        interfaces: Dict[str, Any],
        structure: Structure,
    ) -> Dict[str, Any]:
        """Analyze proximity of binding site to domain interfaces.

        Args:
            site_residues: List of binding site residues
            interfaces: Dictionary of domain interface analysis
            structure: BioPython Structure object

        Returns:
            Dictionary containing interface proximity analysis
        """
        try:
            interface_residues = set()
            min_distances = {}
            max_distance = 0.0

            # Get all interface residues
            for interface_data in interfaces.values():
                interface_residues.update(interface_data["residues"])

            # Calculate distances to interfaces
            structure_residues = {res.get_id()[1]: res for res in structure.get_residues()}
            for site_res in site_residues:
                site_id = site_res.get_id()[1]
                min_dist = float("inf")

                for interface_id in interface_residues:
                    if interface_id not in structure_residues:
                        continue
                    interface_res = structure_residues[interface_id]

                    # Calculate minimum atomic distance
                    for atom1 in site_res:
                        for atom2 in interface_res:
                            distance = atom1 - atom2
                            min_dist = min(min_dist, distance)

                if min_dist < float("inf"):
                    min_distances[site_id] = min_dist
                    max_distance = max(max_distance, min_dist)

            return {
                "interface_residues": sorted(interface_residues),
                "min_distances": min_distances,
                "max_distance": max_distance,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing interface proximity: {str(e)}")
            return {}

    def _calculate_domain_metrics(
        self,
        domains: Dict[str, List[int]],
        interfaces: Dict[str, Any],
        motions: Dict[str, Any],
    ) -> Dict[str, float]:
        """Calculate overall domain analysis metrics.

        Args:
            domains: Dictionary mapping domain IDs to residue lists
            interfaces: Dictionary of interface analysis
            motions: Dictionary of domain motion analysis

        Returns:
            Dictionary of domain metrics
        """
        try:
            metrics = {
                "num_domains": len(domains),
                "num_interfaces": len(interfaces),
                "total_interface_residues": sum(data["size"] for data in interfaces.values()),
            }

            # Add motion metrics if available
            if motions:
                domain_flex = [data["mean_flexibility"] for data in motions.values()]
                if domain_flex:
                    metrics.update(
                        {
                            "mean_domain_flexibility": float(np.mean(domain_flex)),
                            "max_domain_flexibility": float(np.max(domain_flex)),
                            "domain_flexibility_variation": float(np.std(domain_flex)),
                        }
                    )

            return metrics

        except Exception as e:
            self.logger.error(f"Error calculating domain metrics: {str(e)}")
            return {}
