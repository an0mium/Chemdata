"""Protein dynamics analysis functionality."""

from typing import Dict, List, Optional, Any
from Bio.PDB.Structure import Structure

from .contacts import ContactAnalyzer
from .flexibility import FlexibilityAnalyzer
from .modes import NormalModeAnalyzer

__all__ = [
    "DynamicsAnalyzer",
    "analyze_normal_modes",
    "analyze_contacts",
    "analyze_flexibility",
]


def analyze_normal_modes(structure: Structure) -> Dict[str, Any]:
    """Analyze protein normal modes.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary containing normal mode analysis results
    """
    analyzer = NormalModeAnalyzer()
    return analyzer.analyze_normal_modes(structure)


def analyze_contacts(structure: Structure) -> Dict[str, Any]:
    """Analyze protein contact network.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary containing contact analysis results
    """
    analyzer = ContactAnalyzer()
    return analyzer.analyze_contacts(structure)


def analyze_flexibility(structure: Structure) -> Dict[str, Any]:
    """Analyze protein flexibility.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary containing flexibility analysis results
    """
    analyzer = FlexibilityAnalyzer()
    return analyzer.analyze_flexibility(structure)


class DynamicsAnalyzer:
    """Comprehensive protein dynamics analysis."""

    def __init__(self):
        """Initialize dynamics analyzer with all subcomponents."""
        self.contact_analyzer = ContactAnalyzer()
        self.flexibility_analyzer = FlexibilityAnalyzer()

    def analyze_dynamics(
        self,
        structure: Structure,
        include_contacts: bool = True,
        include_flexibility: bool = True,
    ) -> Dict[str, Any]:
        """Analyze protein dynamics using multiple methods.

        Args:
            structure: BioPython Structure object
            include_contacts: Whether to analyze contact networks
            include_flexibility: Whether to analyze flexibility

        Returns:
            Dictionary containing comprehensive dynamics analysis
        """
        analysis = {}

        # Contact network analysis
        if include_contacts:
            contacts = self.contact_analyzer.analyze_contacts(structure)
            analysis["contacts"] = contacts

        # Flexibility analysis
        if include_flexibility:
            flexibility = self.flexibility_analyzer.analyze_flexibility(structure)
            analysis["flexibility"] = flexibility

        return analysis

    def analyze_site_dynamics(
        self,
        site_residues: List[int],
        structure: Structure,
        include_context: bool = True,
    ) -> Dict[str, Any]:
        """Analyze dynamics of specific binding site residues.

        Args:
            site_residues: List of residue numbers in site
            structure: Full structure for context
            include_context: Whether to analyze surrounding context

        Returns:
            Dictionary containing site dynamics analysis
        """
        analysis = {}

        # Contact network analysis
        contacts = self.contact_analyzer.analyze_contacts(structure)
        if contacts:
            site_contacts = self.contact_analyzer.analyze_site_contacts(
                site_residues,
                contacts.get("contact_network", {}),
            )
            analysis["contacts"] = site_contacts

        # Flexibility analysis
        site_flexibility = self.flexibility_analyzer.analyze_site_flexibility(
            site_residues,
            structure,
            include_context=include_context,
        )
        analysis["flexibility"] = site_flexibility

        return analysis

    def get_residue_dynamics(
        self,
        residue_id: int,
        structure: Structure,
    ) -> Dict[str, Any]:
        """Get dynamics properties for single residue.

        Args:
            residue_id: Residue number
            structure: Full structure for context

        Returns:
            Dictionary of residue dynamics properties
        """
        # Get contact properties
        contacts = self.contact_analyzer.analyze_contacts(structure)
        contact_props = {}
        if contacts:
            network = contacts.get("contact_network", {})
            if residue_id in network:
                contact_props = {
                    "num_contacts": len(network[residue_id]),
                    "contacting_residues": sorted(list(network[residue_id])),
                }

        # Get flexibility properties
        flex_analysis = self.flexibility_analyzer.analyze_flexibility(structure)
        flex_props = {}
        if flex_analysis:
            scores = flex_analysis.get("residue_flexibility", {})
            if residue_id in scores:
                flex_props = {
                    "flexibility_score": scores[residue_id],
                }

                # Add region type if available
                for region_type, regions in flex_analysis.get("regions", {}).items():
                    for region in regions:
                        if residue_id in region:
                            flex_props["region_type"] = region_type
                            break

        return {
            "contacts": contact_props,
            "flexibility": flex_props,
        }
