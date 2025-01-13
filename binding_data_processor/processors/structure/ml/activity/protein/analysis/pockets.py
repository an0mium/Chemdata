"""Protein binding pocket analysis functionality."""

from typing import Dict, List, Any, Optional
from Bio.PDB.Structure import Structure

from ...binding.structure.pockets import PocketDetector

# Initialize pocket detector
_detector = PocketDetector()


def find_binding_pockets(
    structure: Structure,
    min_volume: float = 100.0,
    probe_radius: float = 1.4,
    properties: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Find potential binding pockets in protein structure.

    Args:
        structure: BioPython Structure object
        min_volume: Minimum pocket volume in Å³
        probe_radius: Probe radius for surface calculation
        properties: Pre-calculated structure properties

    Returns:
        List of dictionaries containing pocket properties
    """
    return _detector.find_pockets(
        structure,
        min_volume=min_volume,
        probe_radius=probe_radius,
        properties=properties,
    )


def analyze_pocket_properties(
    pockets: List[Dict[str, Any]],
    structure: Structure,
    properties: Dict[str, Any],
    include_scores: bool = True,
) -> List[Dict[str, Any]]:
    """Analyze and score binding pocket properties.

    Args:
        pockets: List of pocket dictionaries
        structure: BioPython Structure object
        properties: Structure properties
        include_scores: Whether to include detailed scoring

    Returns:
        List of scored pocket dictionaries
    """
    return _detector.score_pockets(
        pockets,
        structure,
        properties,
        include_scores=include_scores,
    )
