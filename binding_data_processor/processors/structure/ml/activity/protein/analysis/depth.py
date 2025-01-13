"""Residue depth calculation utilities."""

import logging
import numpy as np
from Bio.PDB import Residue

logger = logging.getLogger(__name__)


def calculate_residue_depth(residue: Residue, surface_points: np.ndarray) -> float:
    """Calculate depth of residue from protein surface.

    Args:
        residue: BioPython Residue object
        surface_points: Array of surface point coordinates (N x 3)

    Returns:
        Depth from surface in Angstroms
    """
    try:
        # Calculate residue center
        residue_coords = []
        for atom in residue:
            residue_coords.append(atom.get_coord())
        if not residue_coords:
            return 0.0

        residue_center = np.mean(residue_coords, axis=0)

        # Calculate minimum distance to any surface point
        min_distance = float("inf")
        for point in surface_points:
            dist = np.linalg.norm(point - residue_center)
            min_distance = min(min_distance, dist)

        return min_distance

    except Exception as e:
        logger.error(f"Error calculating residue depth: {str(e)}")
        return 0.0
