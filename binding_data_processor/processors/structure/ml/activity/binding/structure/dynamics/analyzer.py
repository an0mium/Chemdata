"""Main dynamics analysis coordinator."""

import logging
from typing import Dict, List, Optional, Any
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

from .modes import NormalModeAnalyzer
from .contacts import ContactAnalyzer
from .flexibility import FlexibilityAnalyzer
from .domains import DomainAnalyzer
from .bfactors import BFactorAnalyzer

logger = logging.getLogger(__name__)


class DynamicsAnalyzer:
    """Comprehensive protein dynamics analysis with enhanced capabilities."""

    def __init__(self):
        """Initialize dynamics analyzer with all components."""
        self.logger = logging.getLogger(__name__)
        self.mode_analyzer = NormalModeAnalyzer()
        self.contact_analyzer = ContactAnalyzer()
        self.flexibility_analyzer = FlexibilityAnalyzer()
        self.domain_analyzer = DomainAnalyzer()
        self.bfactor_analyzer = BFactorAnalyzer()

        # Enhanced settings
        self.settings = {
            "min_region_size": 3,  # Minimum residues for a region
            "contact_cutoff": 8.0,  # Å, distance cutoff for contacts
            "flexibility_threshold": 1.0,  # Standard deviations above mean
            "rigidity_threshold": -1.0,  # Standard deviations below mean
            "secondary_structure_weights": {
                "H": 0.3,  # Alpha helix - more rigid
                "G": 0.4,  # 3-10 helix
                "I": 0.3,  # Pi helix
                "E": 0.3,  # Beta strand - more rigid
                "B": 0.5,  # Beta bridge
                "T": 0.7,  # Turn - more flexible
                "S": 0.8,  # Bend - more flexible
                "-": 1.0,  # Coil - most flexible
            },
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
                bfactor_props = self.bfactor_analyzer.analyze_bfactors(structure)
                dynamics.update(bfactor_props)

            # Normal mode analysis
            if use_normal_modes:
                mode_props = self.mode_analyzer.analyze_normal_modes(structure)
                dynamics.update(mode_props)

            # Contact network analysis
            if use_contacts:
                contact_props = self.contact_analyzer.analyze_contacts(structure)
                dynamics.update(contact_props)

            # Flexibility analysis
            if use_flexibility:
                flex_props = self.flexibility_analyzer.analyze_flexibility(structure)
                dynamics.update(flex_props)

            # Domain motion analysis
            if use_domain_motions:
                motion_props = self.domain_analyzer.analyze_domain_motions(
                    structure,
                    contact_props if use_contacts else None,
                    flex_props if use_flexibility else None,
                )
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
        include_correlations: bool = True,
    ) -> Dict[str, Any]:
        """Analyze dynamics of binding site residues.

        Args:
            residues: List of residues in binding site
            structure: Full structure for context
            include_domain_context: Whether to analyze domain context
            include_correlations: Whether to analyze correlations

        Returns:
            Dictionary of dynamics metrics
        """
        try:
            analysis = {}

            # Get residue IDs
            residue_ids = [res.get_id()[1] for res in residues]

            # B-factor analysis
            bfactors = self.bfactor_analyzer.get_residue_bfactors(residues)
            analysis["bfactors"] = bfactors

            # Flexibility analysis
            flexibility = self.flexibility_analyzer.analyze_site_flexibility(
                residue_ids,
                structure,
                include_context=True,
            )
            analysis["flexibility"] = flexibility

            # Contact network analysis
            contacts = self.contact_analyzer.analyze_contacts(structure)
            if contacts:
                site_contacts = self.contact_analyzer.analyze_site_contacts(
                    residue_ids,
                    contacts.get("contact_network", {}),
                )
                analysis["contacts"] = site_contacts

            # Domain context analysis
            if include_domain_context:
                domain_context = self.domain_analyzer.analyze_site_domain_context(
                    residues,
                    structure,
                )
                analysis["domain_context"] = domain_context

            # Context analysis
            context = self._analyze_site_context(
                residue_ids,
                structure,
                contacts,
                flexibility,
            )
            analysis["context"] = context

            # Correlation analysis
            if include_correlations:
                correlations = self._analyze_site_correlations(
                    residue_ids,
                    structure,
                    contacts,
                    flexibility,
                )
                analysis["correlations"] = correlations

            # Calculate site metrics
            metrics = {
                "average_bfactor": float(bfactors["mean"]),
                "bfactor_std": float(bfactors["std"]),
                "relative_bfactor": self.bfactor_analyzer.calculate_relative_bfactor(
                    bfactors["values"],
                    structure,
                ),
                "flexibility_score": float(flexibility["mean"]),
                "rigidity_score": 1.0 - float(flexibility["mean"]),
                "variability": float(flexibility["std"]),
                "contact_density": site_contacts.get("contact_density", 0.0),
                "surface_exposure": site_contacts.get("surface_exposure", 0.0),
            }

            # Add domain context metrics
            if include_domain_context:
                metrics.update(domain_context.get("metrics", {}))

            analysis["metrics"] = metrics

            return analysis

        except Exception as e:
            self.logger.error(f"Error analyzing site dynamics: {str(e)}")
            return {}

    def _analyze_correlations(self, dynamics: Dict[str, Any]) -> Dict[str, Any]:
        """Analyze correlations between different dynamics properties.

        Args:
            dynamics: Dictionary of dynamics data

        Returns:
            Dictionary of correlation properties
        """
        try:
            correlations = {}

            # Extract relevant properties
            properties = {
                "bfactors": dynamics.get("bfactors", {}),
                "flexibility": dynamics.get("residue_flexibility", {}),
                "contacts": dynamics.get("degrees", {}),
                "modes": dynamics.get("mode_fluctuations", {}),
            }

            # Calculate correlations between properties
            for prop1 in properties:
                for prop2 in properties:
                    if prop1 < prop2:
                        corr = self._calculate_property_correlation(
                            properties[prop1],
                            properties[prop2],
                        )
                        if corr is not None:
                            correlations[f"{prop1}_{prop2}_correlation"] = float(corr)

            # Add secondary structure correlations if available
            if "secondary_structure" in dynamics:
                ss_corr = self._analyze_secondary_structure_correlations(
                    dynamics["secondary_structure"],
                    properties.get("flexibility", {}),
                    properties.get("contacts", {}),
                )
                correlations["secondary_structure"] = ss_corr

            return correlations

        except Exception as e:
            self.logger.error(f"Error analyzing correlations: {str(e)}")
            return {}

    def _calculate_property_correlation(
        self,
        prop1: Dict[int, float],
        prop2: Dict[int, float],
    ) -> Optional[float]:
        """Calculate correlation between two property dictionaries.

        Args:
            prop1: First property dictionary
            prop2: Second property dictionary

        Returns:
            Correlation coefficient or None if calculation fails
        """
        try:
            # Get common residues
            common_residues = set(prop1.keys()) & set(prop2.keys())
            if len(common_residues) < 3:
                return None

            # Extract values
            values1 = [prop1[res] for res in common_residues]
            values2 = [prop2[res] for res in common_residues]

            # Calculate correlation
            correlation = np.corrcoef(values1, values2)[0, 1]
            return float(correlation)

        except Exception as e:
            self.logger.error(f"Error calculating correlation: {str(e)}")
            return None

    def _analyze_site_context(
        self,
        site_residues: List[int],
        structure: Structure,
        contacts: Dict[str, Any],
        flexibility: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Analyze dynamics context around binding site.

        Args:
            site_residues: List of residue numbers in site
            structure: Full structure
            contacts: Contact network analysis results
            flexibility: Flexibility analysis results

        Returns:
            Dictionary of context analysis results
        """
        try:
            # Get surrounding residues
            network = contacts.get("contact_network", {})
            if not network:
                return {}

            context_residues = set()
            for res_id in site_residues:
                if res_id in network:
                    context_residues.update(network[res_id])
            context_residues = context_residues - set(site_residues)

            if not context_residues:
                return {}

            # Analyze context properties
            context = {
                "residues": sorted(list(context_residues)),
                "size": len(context_residues),
            }

            # Contact analysis
            context["contacts"] = self._analyze_context_contacts(
                site_residues,
                context_residues,
                network,
            )

            # Flexibility analysis
            flex_scores = flexibility.get("residue_flexibility", {})
            if flex_scores:
                context["flexibility"] = self._analyze_context_flexibility(
                    site_residues,
                    context_residues,
                    flex_scores,
                )

            # Secondary structure analysis
            try:
                dssp = dssp_dict_from_pdb_file(structure)
                context["secondary_structure"] = self._get_secondary_structure_composition(
                    list(context_residues),
                    dssp,
                )
            except Exception:
                pass

            return context

        except Exception as e:
            self.logger.error(f"Error analyzing site context: {str(e)}")
            return {}

    def _analyze_site_correlations(
        self,
        site_residues: List[int],
        structure: Structure,
        contacts: Dict[str, Any],
        flexibility: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Analyze correlations between site and surrounding residues.

        Args:
            site_residues: List of residue numbers in site
            structure: Full structure
            contacts: Contact network analysis results
            flexibility: Flexibility analysis results

        Returns:
            Dictionary of correlation analysis results
        """
        try:
            correlations = {}

            # Get required data
            network = contacts.get("contact_network", {})
            flex_scores = flexibility.get("residue_flexibility", {})
            if not network or not flex_scores:
                return {}

            # Get site and context residues
            context_residues = set()
            for res_id in site_residues:
                if res_id in network:
                    context_residues.update(network[res_id])
            context_residues = context_residues - set(site_residues)

            if not context_residues:
                return {}

            # Calculate correlations
            site_flex = [flex_scores[res] for res in site_residues if res in flex_scores]
            context_flex = [flex_scores[res] for res in context_residues if res in flex_scores]

            if site_flex and context_flex:
                correlations["site_context_flexibility"] = float(np.corrcoef(site_flex, context_flex)[0, 1])

            # Add secondary structure correlations
            try:
                dssp = dssp_dict_from_pdb_file(structure)
                ss_corr = self._analyze_secondary_structure_correlations(
                    dssp,
                    {res: flex_scores[res] for res in site_residues if res in flex_scores},
                    {res: len(network[res]) for res in site_residues if res in network},
                )
                correlations["secondary_structure"] = ss_corr
            except Exception:
                pass

            return correlations

        except Exception as e:
            self.logger.error(f"Error analyzing site correlations: {str(e)}")
            return {}

    def _analyze_context_contacts(
        self,
        site_residues: List[int],
        context_residues: set,
        network: Dict[int, set],
    ) -> Dict[str, Any]:
        """Analyze contacts between site and context.

        Args:
            site_residues: List of site residue numbers
            context_residues: Set of context residue numbers
            network: Contact network

        Returns:
            Dictionary of contact properties
        """
        try:
            site_context_contacts = 0
            context_contacts = 0

            for res_id in context_residues:
                if res_id in network:
                    neighbors = network[res_id]
                    site_context_contacts += len(neighbors & set(site_residues))
                    context_contacts += len(neighbors & context_residues)

            # Adjust for double counting
            context_contacts //= 2

            return {
                "site_context_contacts": site_context_contacts,
                "context_contacts": context_contacts,
                "contact_density": float(context_contacts / (len(context_residues) * (len(context_residues) - 1) / 2)) if len(context_residues) > 1 else 0.0,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing context contacts: {str(e)}")
            return {}

    def _analyze_context_flexibility(
        self,
        site_residues: List[int],
        context_residues: set,
        flexibility_scores: Dict[int, float],
    ) -> Dict[str, float]:
        """Analyze flexibility of context residues.

        Args:
            site_residues: List of site residue numbers
            context_residues: Set of context residue numbers
            flexibility_scores: Dictionary mapping residues to flexibility scores

        Returns:
            Dictionary of flexibility properties
        """
        try:
            site_scores = [flexibility_scores[res] for res in site_residues if res in flexibility_scores]
            context_scores = [flexibility_scores[res] for res in context_residues if res in flexibility_scores]

            if not site_scores or not context_scores:
                return {}

            return {
                "mean_context_flexibility": float(np.mean(context_scores)),
                "context_flexibility_std": float(np.std(context_scores)),
                "flexibility_difference": float(np.mean(site_scores) - np.mean(context_scores)),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing context flexibility: {str(e)}")
            return {}

    def _get_secondary_structure_composition(
        self,
        residues: List[int],
        dssp: Dict,
    ) -> Dict[str, float]:
        """Calculate secondary structure composition.

        Args:
            residues: List of residue numbers
            dssp: DSSP dictionary

        Returns:
            Dictionary mapping SS types to fractions
        """
        try:
            ss_counts = {"H": 0, "E": 0, "C": 0}  # Helix, Sheet, Other
            total = 0

            for res_id in residues:
                if res_id in dssp:
                    ss = dssp[res_id][2]
                    if ss in ["H", "G", "I"]:  # All helices
                        ss_counts["H"] += 1
                    elif ss in ["E", "B"]:  # All sheets
                        ss_counts["E"] += 1
                    else:  # Everything else
                        ss_counts["C"] += 1
                    total += 1

            if total == 0:
                return {}

            return {ss: count / total for ss, count in ss_counts.items()}

        except Exception as e:
            self.logger.error(f"Error calculating SS composition: {str(e)}")
            return {}

    def _analyze_secondary_structure_correlations(
        self,
        dssp: Dict,
        flexibility_scores: Dict[int, float],
        contact_counts: Dict[int, int],
    ) -> Dict[str, float]:
        """Analyze correlations with secondary structure.

        Args:
            dssp: DSSP dictionary
            flexibility_scores: Residue flexibility scores
            contact_counts: Residue contact counts

        Returns:
            Dictionary of correlation coefficients
        """
        try:
            # Group scores by secondary structure
            ss_flex = {"H": [], "E": [], "C": []}
            ss_contacts = {"H": [], "E": [], "C": []}

            for res_id in dssp:
                ss = dssp[res_id][2]
                if ss in ["H", "G", "I"]:
                    ss_type = "H"
                elif ss in ["E", "B"]:
                    ss_type = "E"
                else:
                    ss_type = "C"

                if res_id in flexibility_scores:
                    ss_flex[ss_type].append(flexibility_scores[res_id])
                if res_id in contact_counts:
                    ss_contacts[ss_type].append(contact_counts[res_id])

            # Calculate statistics
            correlations = {}
            for ss_type in ["H", "E", "C"]:
                if ss_flex[ss_type] and ss_contacts[ss_type]:
                    corr = float(np.corrcoef(ss_flex[ss_type], ss_contacts[ss_type])[0, 1])
                    correlations[f"{ss_type}_correlation"] = corr

            return correlations

        except Exception as e:
            self.logger.error(f"Error analyzing SS correlations: {str(e)}")
            return {}

    def _calculate_dynamics_metrics(self, dynamics: Dict[str, Any]) -> Dict[str, float]:
        """Calculate overall dynamics metrics.

        Args:
            dynamics: Complete dynamics analysis results

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

            # Domain based metrics
            if "domains" in dynamics:
                domains = dynamics["domains"]
                if isinstance(domains, dict):
                    interfaces = domains.get("interfaces", [])
                    metrics["num_domains"] = len(domains) - 1  # Subtract 1 for interfaces key
                    metrics["num_interfaces"] = len(interfaces)

            return metrics

        except Exception as e:
            self.logger.error(f"Error calculating dynamics metrics: {str(e)}")
            return {}
