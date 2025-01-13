"""Domain analysis module for protein dynamics."""

import logging
from typing import Dict, List, Any, Optional
import numpy as np
from Bio.PDB import Structure, Residue

from binding_data_processor.core.config import ProteinAnalysisConfig
from ..utils import get_ca_atoms

logger = logging.getLogger(__name__)


def analyze_domain_motions(
    structure: Structure,
    domains: Dict[str, List[int]],
    config: Optional[ProteinAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze domain motions using normal modes.

    Args:
        structure: BioPython Structure object
        domains: Dictionary mapping domain IDs to lists of residue numbers
        config: Optional configuration object

    Returns:
        Dictionary containing domain motion analysis including:
        - Domain-specific motions
        - Inter-domain correlations
        - Hinge regions
    """
    try:
        if config is None:
            config = ProteinAnalysisConfig()

        # Get CA atoms and coordinates
        ca_atoms = []
        ca_coords = []
        residue_indices = {}  # Map residue numbers to indices

        for i, residue in enumerate(structure.get_residues()):
            if "CA" in residue:
                ca_atoms.append(residue["CA"])
                ca_coords.append(residue["CA"].get_coord())
                residue_indices[residue.get_id()[1]] = i

        ca_coords = np.array(ca_coords)

        # Map domains to coordinate indices
        domain_indices = {}
        for domain_id, residues in domains.items():
            indices = [residue_indices[res_id] for res_id in residues if res_id in residue_indices]
            if indices:
                domain_indices[domain_id] = indices

        if not domain_indices:
            return {}

        # Calculate domain centers
        domain_centers = {}
        for domain_id, indices in domain_indices.items():
            center = np.mean(ca_coords[indices], axis=0)
            domain_centers[domain_id] = center

        # Calculate domain motions from normal modes
        domain_motions = {}
        for domain_id, indices in domain_indices.items():
            # Calculate domain displacement vectors
            displacements = ca_coords[indices] - domain_centers[domain_id]

            # Calculate principal axes of motion
            try:
                U, S, Vt = np.linalg.svd(displacements)
                principal_axes = Vt
                variances = S**2 / len(indices)

                domain_motions[domain_id] = {
                    "principal_axes": principal_axes.tolist(),
                    "variances": variances.tolist(),
                    "center": domain_centers[domain_id].tolist(),
                }
            except np.linalg.LinAlgError:
                logger.warning(f"Could not calculate SVD for domain {domain_id}")
                continue

        # Calculate inter-domain correlations
        correlations = {}
        for d1 in domain_indices:
            for d2 in domain_indices:
                if d1 < d2:
                    corr = _calculate_domain_correlation(
                        ca_coords[domain_indices[d1]],
                        ca_coords[domain_indices[d2]],
                    )
                    correlations[f"{d1}-{d2}"] = float(corr)

        # Identify hinge regions
        hinge_regions = _identify_hinge_regions(
            ca_coords,
            domain_indices,
            domain_centers,
        )

        return {
            "domain_motions": domain_motions,
            "correlations": correlations,
            "hinge_regions": hinge_regions,
        }

    except Exception as e:
        logger.error(f"Error analyzing domain motions: {str(e)}")
        return {}


def _calculate_domain_correlation(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Calculate correlation between domain motions.

    Args:
        coords1: Coordinates of first domain
        coords2: Coordinates of second domain

    Returns:
        Correlation coefficient
    """
    try:
        # Calculate displacement vectors
        disp1 = coords1 - np.mean(coords1, axis=0)
        disp2 = coords2 - np.mean(coords2, axis=0)

        # Calculate correlation matrix
        corr = np.corrcoef(disp1.flatten(), disp2.flatten())[0, 1]
        return float(corr)

    except Exception as e:
        logger.error(f"Error calculating domain correlation: {str(e)}")
        return 0.0


def _identify_hinge_regions(
    coords: np.ndarray,
    domain_indices: Dict[str, List[int]],
    domain_centers: Dict[str, np.ndarray],
) -> List[Dict[str, Any]]:
    """Identify hinge regions between domains.

    Args:
        coords: CA coordinates
        domain_indices: Mapping of domain IDs to coordinate indices
        domain_centers: Domain center coordinates

    Returns:
        List of dictionaries containing hinge region information
    """
    try:
        hinge_regions = []

        # For each pair of domains
        domain_ids = list(domain_indices.keys())
        for i, d1 in enumerate(domain_ids):
            for d2 in domain_ids[i + 1 :]:
                # Get residues at domain interface
                interface_residues = _get_interface_residues(
                    coords,
                    domain_indices[d1],
                    domain_indices[d2],
                )

                if interface_residues:
                    # Calculate hinge score based on distance to domain centers
                    scores = []
                    for res_idx in interface_residues:
                        dist1 = np.linalg.norm(coords[res_idx] - domain_centers[d1])
                        dist2 = np.linalg.norm(coords[res_idx] - domain_centers[d2])
                        scores.append(min(dist1, dist2))

                    hinge_regions.append(
                        {
                            "domains": [d1, d2],
                            "residues": interface_residues,
                            "scores": [float(s) for s in scores],
                        }
                    )

        return hinge_regions

    except Exception as e:
        logger.error(f"Error identifying hinge regions: {str(e)}")
        return []


def _get_interface_residues(
    coords: np.ndarray,
    indices1: List[int],
    indices2: List[int],
    cutoff: float = 8.0,
) -> List[int]:
    """Get residues at domain interface.

    Args:
        coords: CA coordinates
        indices1: Indices of first domain
        indices2: Indices of second domain
        cutoff: Distance cutoff for interface residues

    Returns:
        List of residue indices at interface
    """
    try:
        interface = set()

        # Calculate distances between all pairs
        for i in indices1:
            for j in indices2:
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < cutoff:
                    interface.add(i)
                    interface.add(j)

        return sorted(list(interface))

    except Exception as e:
        logger.error(f"Error getting interface residues: {str(e)}")
        return []
