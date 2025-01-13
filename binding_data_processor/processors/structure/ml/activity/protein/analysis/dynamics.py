"""Protein dynamics analysis module.

This module provides comprehensive protein dynamics analysis capabilities including:
- Normal modes analysis using elastic network model
- Flexibility analysis combining B-factors and normal modes
- Domain motion analysis with AlphaFold confidence integration
- Site-specific dynamics analysis
- Contact network and correlation analysis
"""

import logging
from typing import Dict, List, Any, Optional, Tuple
import numpy as np
from Bio.PDB import Structure, Residue

from binding_data_processor.core.config import ProteinAnalysisConfig, DynamicsAnalysisConfig
from .base import BaseProteinAnalyzer
from .utils import get_ca_atoms
from .dynamics.analyzer import DynamicsAnalyzer
from .dynamics.domains import analyze_domain_motions

logger = logging.getLogger(__name__)


def analyze_dynamics(
    structure: Structure,
    config: Optional[DynamicsAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze protein dynamics using multiple methods.

    Performs comprehensive dynamics analysis including:
    - Normal modes calculation
    - Flexibility analysis
    - Domain motion analysis
    - Contact network analysis
    - AlphaFold confidence integration

    Args:
        structure: BioPython Structure object
        config: Optional dynamics analysis configuration

    Returns:
        Dictionary containing:
        - flexibility: Per-residue flexibility scores
        - correlations: Residue-residue correlation matrix
        - modes: Top normal modes
        - domain_motions: Domain motion analysis results
        - contact_network: Residue contact analysis
        - site_dynamics: Site-specific dynamics (if sites defined)
    """
    try:
        if config is None:
            config = DynamicsAnalysisConfig()

        # Initialize analyzer
        analyzer = DynamicsAnalyzer(structure, config)

        # Get structure coordinates
        coords = analyzer.get_coordinates()

        # Calculate normal modes if enabled
        modes = None
        if config.use_normal_modes:
            modes_results = analyzer.analyze_normal_modes(n_modes=config.protein_config.n_modes)
            modes = np.array(modes_results["modes"])

        # Calculate flexibility profile
        flexibility_results = analyzer.analyze_flexibility()
        flexibility = np.array([flexibility_results["residue_b_factors"].get(i, 0.0) for i in range(len(list(structure.get_residues())))])

        # Calculate correlations if enabled
        correlations = None
        if config.use_correlations and modes is not None:
            correlations = analyzer._get_correlation_matrix()

        # Analyze domain motions if enabled
        domain_motions = None
        if config.use_domain_motions:
            domain_motions = _analyze_domain_motions(
                structure,
                modes=modes,
                config=config,
            )

        # Weight results by AlphaFold confidence if available
        if config.protein_config.use_alphafold_confidence:
            flexibility = _apply_confidence_weights(
                flexibility,
                structure,
                config.advanced_params["min_plddt"],
            )

        results = {
            "flexibility": flexibility,
            "correlations": correlations,
            "modes": modes,
            "domain_motions": domain_motions,
        }

        # Add contact network analysis
        if config.use_contacts:
            results["contact_network"] = analyze_contacts(
                structure,
                cutoff=config.protein_config.contact_cutoff,
                config=config,
            )

        return results

    except Exception as e:
        logger.error(f"Error analyzing dynamics: {str(e)}")
        return {}


def analyze_site_dynamics(
    residues: List[int],
    structure: Structure,
    config: Optional[DynamicsAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze dynamics of specific binding site residues.

    Args:
        residues: List of residue numbers in binding site
        structure: Full structure for context
        config: Optional dynamics analysis configuration

    Returns:
        Dictionary containing site dynamics analysis including:
        - Site-specific flexibility
        - Local correlations
        - Contact network
        - Domain context
        - AlphaFold confidence metrics
    """
    try:
        if config is None:
            config = DynamicsAnalysisConfig()

        # Initialize analyzer
        analyzer = DynamicsAnalyzer(structure, config)

        # Get residue objects
        structure_residues = {res.get_id()[1]: res for res in structure.get_residues()}
        site_residues = [structure_residues[res_id] for res_id in residues if res_id in structure_residues]

        if not site_residues:
            return {}

        # Get site-specific dynamics
        results = analyzer.analyze_site_dynamics(
            site_residues,
            structure,
            include_domain_context=config.use_domain_motions,
            include_correlations=config.use_correlations,
        )

        # Add AlphaFold confidence analysis if available
        if config.protein_config.use_alphafold_confidence:
            results["confidence"] = _analyze_site_confidence(
                site_residues,
                structure,
                min_plddt=config.advanced_params["min_plddt"],
            )

        return results

    except Exception as e:
        logger.error(f"Error analyzing site dynamics: {str(e)}")
        return {}


def analyze_flexibility(
    structure: Structure,
    config: Optional[DynamicsAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze protein flexibility using multiple methods.

    Combines information from:
    - B-factors
    - Normal modes
    - Contact topology
    - AlphaFold confidence scores

    Args:
        structure: BioPython Structure object
        config: Optional dynamics analysis configuration

    Returns:
        Dictionary containing:
        - per_residue: Per-residue flexibility scores
        - regions: Identified flexible and rigid regions
        - statistics: Overall flexibility metrics
        - confidence: AlphaFold-based confidence metrics
    """
    try:
        if config is None:
            config = DynamicsAnalysisConfig()

        # Get full dynamics analysis
        results = analyze_dynamics(structure, config)

        # Extract and process flexibility information
        flexibility_results = {
            "per_residue": results["flexibility"],
            "regions": _identify_flexibility_regions(
                results["flexibility"],
                window=config.advanced_params["flexibility_window"],
            ),
            "statistics": _calculate_flexibility_statistics(results["flexibility"]),
        }

        # Add confidence metrics if available
        if config.protein_config.use_alphafold_confidence:
            flexibility_results["confidence"] = _get_confidence_metrics(structure)

        return flexibility_results

    except Exception as e:
        logger.error(f"Error analyzing flexibility: {str(e)}")
        return {}


def analyze_contacts(
    structure: Structure,
    cutoff: float = 8.0,
    config: Optional[DynamicsAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze protein contact network.

    Args:
        structure: BioPython Structure object
        cutoff: Distance cutoff for contacts in Angstroms
        config: Optional dynamics analysis configuration

    Returns:
        Dictionary containing:
        - contact_map: Residue contact matrix
        - statistics: Contact network metrics
        - key_contacts: Important contact residues
        - communities: Contact-based communities
    """
    try:
        if config is None:
            config = DynamicsAnalysisConfig()

        # Get C-alpha atoms
        ca_atoms = get_ca_atoms(structure)
        if not ca_atoms:
            return {}

        # Calculate contact matrix
        n_res = len(ca_atoms)
        contacts = np.zeros((n_res, n_res))

        for i in range(n_res):
            for j in range(i + 1, n_res):
                diff = ca_atoms[i].get_coord() - ca_atoms[j].get_coord()
                dist = np.sqrt(np.sum(diff * diff))
                if dist < cutoff:
                    contacts[i, j] = contacts[j, i] = 1

        # Weight contacts by AlphaFold PAE if available
        if config.protein_config.use_alphafold_pae:
            contacts = _weight_contacts_by_pae(
                contacts,
                structure,
                config.advanced_params["pae_power"],
            )

        return {
            "contact_map": contacts,
            "statistics": _calculate_network_statistics(contacts),
            "key_contacts": _identify_key_contacts(contacts),
            "communities": _identify_contact_communities(contacts),
        }

    except Exception as e:
        logger.error(f"Error analyzing contacts: {str(e)}")
        return {}


def _calculate_normal_modes(
    coords: np.ndarray,
    n_modes: int = 10,
    cutoff: float = 0.1,
) -> np.ndarray:
    """Calculate protein normal modes using elastic network model."""
    # Calculate distance matrix
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff * diff, axis=-1))

    # Build Hessian matrix
    hessian = _build_hessian(distances, cutoff)

    # Calculate eigenvalues and eigenvectors
    eigenvals, eigenvecs = np.linalg.eigh(hessian)

    # Sort and filter modes
    sorted_idx = np.argsort(eigenvals)[6 : n_modes + 6]  # Skip first 6 trivial modes
    modes = eigenvecs[:, sorted_idx]

    return modes


def _calculate_flexibility(
    structure: Structure,
    modes: Optional[np.ndarray] = None,
    window: int = 5,
    use_bfactors: bool = True,
) -> np.ndarray:
    """Calculate per-residue flexibility scores."""
    n_residues = len(list(structure.get_residues()))
    flexibility = np.zeros(n_residues)

    # Use B-factors if available and enabled
    if use_bfactors:
        for i, residue in enumerate(structure.get_residues()):
            b_factors = [atom.get_bfactor() for atom in residue]
            flexibility[i] = np.mean(b_factors)

    # Add contribution from normal modes
    if modes is not None:
        mode_flex = np.sum(modes * modes, axis=1)
        mode_flex = mode_flex.reshape(n_residues, -1).mean(axis=1)

        if use_bfactors:
            # Combine with B-factors
            flexibility = 0.5 * (flexibility + mode_flex)
        else:
            flexibility = mode_flex

    # Smooth using sliding window
    if window > 1:
        flexibility = _sliding_window_smooth(flexibility, window)

    return flexibility


def _calculate_correlations(
    modes: np.ndarray,
    cutoff: float = 0.5,
) -> np.ndarray:
    """Calculate residue-residue correlations from normal modes."""
    correlations = np.zeros((modes.shape[0], modes.shape[0]))

    for i in range(modes.shape[1]):
        mode = modes[:, i].reshape(-1, 3)
        corr = np.outer(mode[:, 0], mode[:, 0])
        corr += np.outer(mode[:, 1], mode[:, 1])
        corr += np.outer(mode[:, 2], mode[:, 2])
        correlations += corr

    # Normalize
    correlations /= modes.shape[1]

    # Apply cutoff
    correlations[np.abs(correlations) < cutoff] = 0

    return correlations


def _analyze_domain_motions(
    structure: Structure,
    modes: Optional[np.ndarray],
    config: DynamicsAnalysisConfig,
) -> Dict:
    """Analyze domain motions using normal modes."""
    # Get domain definitions from structure
    domains = _identify_domains(
        structure,
        contact_cutoff=config.protein_config.domain_analysis_params["contact_cutoff"],
        min_size=config.protein_config.domain_analysis_params["min_domain_size"],
    )

    if modes is None or not domains:
        return None

    # Analyze domain-specific mode components
    domain_motions = {}
    for domain_id, domain_residues in domains.items():
        # Get mode components for domain
        domain_modes = modes[domain_residues]

        # Calculate domain motion magnitude and direction
        magnitude = np.sqrt(np.sum(domain_modes * domain_modes))
        direction = domain_modes.mean(axis=0)
        direction /= np.linalg.norm(direction)

        domain_motions[domain_id] = {
            "magnitude": magnitude,
            "direction": direction,
        }

    return domain_motions


def _build_hessian(distances: np.ndarray, cutoff: float) -> np.ndarray:
    """Build Hessian matrix for elastic network model."""
    n_atoms = distances.shape[0]
    hessian = np.zeros((3 * n_atoms, 3 * n_atoms))

    # Build spring network
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if distances[i, j] < cutoff:
                diff = distances[i, j]
                force_const = 1.0  # Uniform spring constant

                # Fill Hessian blocks
                for di in range(3):
                    for dj in range(3):
                        if di == dj:
                            val = force_const * diff * diff
                        else:
                            val = 0

                        hessian[3 * i + di, 3 * j + dj] = -val
                        hessian[3 * j + dj, 3 * i + di] = -val
                        hessian[3 * i + di, 3 * i + di] += val
                        hessian[3 * j + dj, 3 * j + dj] += val

    return hessian


def _identify_domains(
    structure: Structure,
    contact_cutoff: float = 8.0,
    min_size: int = 30,
) -> Dict[str, List[int]]:
    """Identify protein domains based on contact topology."""
    # Get C-alpha atoms
    ca_atoms = get_ca_atoms(structure)
    if not ca_atoms:
        return {}

    # Calculate contact matrix
    n_res = len(ca_atoms)
    contacts = np.zeros((n_res, n_res))

    for i in range(n_res):
        for j in range(i + 1, n_res):
            diff = ca_atoms[i].get_coord() - ca_atoms[j].get_coord()
            dist = np.sqrt(np.sum(diff * diff))
            if dist < contact_cutoff:
                contacts[i, j] = contacts[j, i] = 1

    # Identify domains using contact clustering
    domains = {}
    visited = set()
    domain_id = 1

    for i in range(n_res):
        if i in visited:
            continue

        # Find connected residues
        domain = set([i])
        stack = [i]

        while stack:
            current = stack.pop()
            neighbors = np.where(contacts[current] > 0)[0]

            for neighbor in neighbors:
                if neighbor not in visited:
                    domain.add(neighbor)
                    visited.add(neighbor)
                    stack.append(neighbor)

        # Add domain if large enough
        if len(domain) >= min_size:
            domains[f"D{domain_id}"] = sorted(list(domain))
            domain_id += 1

    return domains


def _sliding_window_smooth(
    array: np.ndarray,
    window: int,
) -> np.ndarray:
    """Smooth array using sliding window average."""
    if window == 1:
        return array

    window = min(window, len(array))
    if window % 2 == 0:
        window += 1

    half = window // 2
    smoothed = np.zeros_like(array)

    for i in range(len(array)):
        start = max(0, i - half)
        end = min(len(array), i + half + 1)
        smoothed[i] = array[start:end].mean()

    return smoothed


def _apply_confidence_weights(
    values: np.ndarray,
    structure: Structure,
    min_plddt: float,
) -> np.ndarray:
    """Weight values by AlphaFold pLDDT confidence scores."""
    weighted = values.copy()

    for i, residue in enumerate(structure.get_residues()):
        # Get pLDDT score from B-factor field
        plddt = residue["CA"].get_bfactor() if "CA" in residue else 0

        # Apply confidence weighting
        if plddt < min_plddt:
            weighted[i] = 0
        else:
            weighted[i] *= plddt / 100.0

    return weighted


def _analyze_site_confidence(
    site_residues: List[Residue],
    structure: Structure,
    min_plddt: float,
) -> Dict[str, float]:
    """Analyze AlphaFold confidence scores for binding site residues."""
    site_plddt = []
    for residue in site_residues:
        if "CA" in residue:
            site_plddt.append(residue["CA"].get_bfactor())

    if not site_plddt:
        return {}

    return {
        "mean_plddt": np.mean(site_plddt),
        "min_plddt": np.min(site_plddt),
        "max_plddt": np.max(site_plddt),
        "fraction_confident": np.mean(np.array(site_plddt) >= min_plddt),
    }


def _identify_flexibility_regions(
    flexibility: np.ndarray,
    window: int = 5,
) -> Dict[str, List[Tuple[int, int]]]:
    """Identify contiguous regions of high/low flexibility."""
    # Smooth flexibility profile
    smoothed = _sliding_window_smooth(flexibility, window)

    # Calculate mean and standard deviation
    mean_flex = np.mean(smoothed)
    std_flex = np.std(smoothed)

    # Identify regions
    flexible_regions = []
    rigid_regions = []

    current_start = 0
    current_type = smoothed[0] > mean_flex

    for i in range(1, len(smoothed)):
        is_flexible = smoothed[i] > mean_flex

        if is_flexible != current_type:
            if current_type:
                flexible_regions.append((current_start, i - 1))
            else:
                rigid_regions.append((current_start, i - 1))
            current_start = i
            current_type = is_flexible

    # Add final region
    if current_type:
        flexible_regions.append((current_start, len(smoothed) - 1))
    else:
        rigid_regions.append((current_start, len(smoothed) - 1))

    return {
        "flexible": flexible_regions,
        "rigid": rigid_regions,
    }


def _calculate_flexibility_statistics(
    flexibility: np.ndarray,
) -> Dict[str, float]:
    """Calculate overall flexibility statistics."""
    return {
        "mean": np.mean(flexibility),
        "std": np.std(flexibility),
        "min": np.min(flexibility),
        "max": np.max(flexibility),
        "median": np.median(flexibility),
    }


def _get_confidence_metrics(
    structure: Structure,
) -> Dict[str, float]:
    """Get AlphaFold confidence metrics for the structure."""
    plddt_scores = []
    for residue in structure.get_residues():
        if "CA" in residue:
            plddt_scores.append(residue["CA"].get_bfactor())

    if not plddt_scores:
        return {}

    plddt_scores = np.array(plddt_scores)

    return {
        "mean_plddt": np.mean(plddt_scores),
        "median_plddt": np.median(plddt_scores),
        "std_plddt": np.std(plddt_scores),
        "min_plddt": np.min(plddt_scores),
        "max_plddt": np.max(plddt_scores),
        "very_high_confidence": np.mean(plddt_scores >= 90),
        "high_confidence": np.mean(plddt_scores >= 70),
        "low_confidence": np.mean(plddt_scores < 50),
    }


def _calculate_network_statistics(
    contacts: np.ndarray,
) -> Dict[str, float]:
    """Calculate statistics for contact network."""
    return {
        "mean_contacts": np.mean(np.sum(contacts, axis=1)),
        "max_contacts": np.max(np.sum(contacts, axis=1)),
        "min_contacts": np.min(np.sum(contacts, axis=1)),
        "density": np.sum(contacts) / (contacts.shape[0] * contacts.shape[1]),
    }


def _identify_key_contacts(
    contacts: np.ndarray,
    percentile: float = 90,
) -> List[int]:
    """Identify residues with high number of contacts."""
    contact_counts = np.sum(contacts, axis=1)
    threshold = np.percentile(contact_counts, percentile)
    return list(np.where(contact_counts >= threshold)[0])


def _identify_contact_communities(
    contacts: np.ndarray,
    min_size: int = 5,
) -> List[List[int]]:
    """Identify communities in contact network using simple clustering."""
    communities = []
    visited = set()

    for i in range(contacts.shape[0]):
        if i in visited:
            continue

        # Find connected residues
        community = set([i])
        stack = [i]

        while stack:
            current = stack.pop()
            neighbors = np.where(contacts[current] > 0)[0]

            for neighbor in neighbors:
                if neighbor not in visited:
                    community.add(neighbor)
                    visited.add(neighbor)
                    stack.append(neighbor)

        # Add community if large enough
        if len(community) >= min_size:
            communities.append(sorted(list(community)))

    return communities


def _weight_contacts_by_pae(
    contacts: np.ndarray,
    structure: Structure,
    pae_power: float = 1.0,
) -> np.ndarray:
    """Weight contact matrix by AlphaFold PAE scores."""
    weighted = contacts.copy()

    # Get PAE scores from structure (implementation depends on how PAE is stored)
    # This is a placeholder - actual implementation would depend on PAE data format
    pae = np.ones_like(contacts)  # Replace with actual PAE matrix

    # Apply PAE weighting
    weighted *= (1 - pae) ** pae_power

    return weighted
