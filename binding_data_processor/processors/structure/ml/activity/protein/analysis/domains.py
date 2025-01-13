"""Domain analysis functionality for protein structure analysis."""

from typing import Dict, Any
from Bio.PDB.Structure import Structure

from ...binding.structure.dynamics.domains import DomainAnalyzer

# Re-export DomainAnalyzer for protein analysis module
__all__ = ["DomainAnalyzer", "analyze_domains"]


def analyze_domains(
    structure: Structure,
    include_motions: bool = True,
    include_interfaces: bool = True,
) -> Dict[str, Any]:
    """Analyze protein domains using the DomainAnalyzer.

    This is a convenience wrapper around DomainAnalyzer that provides
    a simpler interface for basic domain analysis in the protein
    analysis module.

    Args:
        structure: BioPython Structure object
        include_motions: Whether to analyze domain motions
        include_interfaces: Whether to analyze domain interfaces

    Returns:
        Dictionary containing domain analysis results
    """
    analyzer = DomainAnalyzer()
    return analyzer.analyze_domain_motions(
        structure,
        contact_data=None if not include_interfaces else {},
        flexibility_data=None if not include_motions else {},
    )
