"""Protein dynamics analysis module."""

import logging
from typing import Dict, List, Optional, Any
import numpy as np
from Bio.PDB import Structure

from binding_data_processor.core.config import ProteinAnalysisConfig
from ..utils import get_ca_atoms
from .analyzer import DynamicsAnalyzer

logger = logging.getLogger(__name__)


def analyze_dynamics(
    structure: Structure,
    config: Optional[ProteinAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze protein dynamics using multiple methods.

    Args:
        structure: BioPython Structure object
        config: Optional configuration object

    Returns:
        Dictionary containing dynamics analysis results including:
        - Basic dynamics metrics
        - Flexibility analysis
        - Normal mode analysis
        - Contact network analysis
        - Overall dynamics metrics
    """
    try:
        if config is None:
            config = ProteinAnalysisConfig()

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
            domain_motions = analyzer.analyze_domain_motions(
                modes=modes,
                config=config,
            )

        # Weight results by AlphaFold confidence if available
        if config.protein_config.use_alphafold_confidence:
            flexibility = analyzer.apply_confidence_weights(
                flexibility,
                min_plddt=config.advanced_params["min_plddt"],
            )

        results = {
            "flexibility": flexibility,
            "correlations": correlations,
            "modes": modes,
            "domain_motions": domain_motions,
        }

        # Add contact network analysis
        if config.use_contacts:
            results["contact_network"] = analyzer.analyze_contacts(
                cutoff=config.protein_config.contact_cutoff,
                config=config,
            )

        return results

    except Exception as e:
        logger.error(f"Error analyzing dynamics: {str(e)}")
        return {}


def analyze_flexibility(
    structure: Structure,
    config: Optional[ProteinAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze protein flexibility using multiple methods.

    Args:
        structure: BioPython Structure object
        config: Optional configuration object

    Returns:
        Dictionary containing flexibility analysis results including:
        - B-factor analysis
        - Flexibility scores
        - Flexible and rigid regions
        - Overall flexibility metrics
    """
    try:
        if config is None:
            config = ProteinAnalysisConfig()

        # Initialize analyzer
        analyzer = DynamicsAnalyzer(structure, config)

        # Get basic flexibility analysis
        basic_results = analyzer.analyze_flexibility()

        # Get enhanced flexibility analysis
        enhanced_results = analyzer.analyze_enhanced_flexibility()

        # Combine results
        results = {
            "basic": basic_results,
            "detailed": enhanced_results,
        }

        return results

    except Exception as e:
        logger.error(f"Error analyzing flexibility: {str(e)}")
        return {}


def analyze_site_dynamics(
    residues: List[int],
    structure: Structure,
    config: Optional[ProteinAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze dynamics of specific binding site residues.

    Args:
        residues: List of residue numbers in binding site
        structure: Full structure for context
        config: Optional configuration object

    Returns:
        Dictionary containing site dynamics analysis including:
        - Site-specific dynamics
        - Local flexibility
        - Contact network
        - Correlation with overall dynamics
    """
    try:
        if config is None:
            config = ProteinAnalysisConfig()

        # Initialize analyzer
        analyzer = DynamicsAnalyzer(structure, config)

        # Get residue objects
        structure_residues = {res.get_id()[1]: res for res in structure.get_residues()}
        site_residues = [structure_residues[res_id] for res_id in residues if res_id in structure_residues]

        if not site_residues:
            return {}

        # Analyze site dynamics
        site_dynamics = analyzer.analyze_site_dynamics(
            site_residues,
            structure,
            include_domain_context=True,
            include_correlations=True,
        )

        # Analyze site flexibility
        site_flexibility = analyzer.analyze_site_flexibility(
            residues,
            structure,
            include_context=True,
        )

        # Analyze site contacts
        contact_network = analyzer.analyze_contacts(structure)["contact_network"]
        site_contacts = analyzer.analyze_site_contacts(residues, contact_network)

        # Calculate site metrics
        metrics = {
            "average_flexibility": float(np.mean([site_flexibility["scores"].get(res, 0.0) for res in residues])),
            "flexibility_variation": float(np.std([site_flexibility["scores"].get(res, 0.0) for res in residues])),
            "contact_density": site_contacts.get("contact_density", 0.0),
            "surface_exposure": site_contacts.get("surface_exposure", 0.0),
            "average_bfactor": float(np.mean([res.get_bfactor() for res in site_residues])),
        }

        results = {
            "dynamics": site_dynamics,
            "flexibility": site_flexibility,
            "contacts": site_contacts,
            "metrics": metrics,
        }

        return results

    except Exception as e:
        logger.error(f"Error analyzing site dynamics: {str(e)}")
        return {}


def analyze_normal_modes(
    structure: Structure,
    n_modes: int = 10,
    config: Optional[ProteinAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze protein normal modes using multiple methods.

    Args:
        structure: BioPython Structure object
        n_modes: Number of normal modes to calculate
        config: Optional configuration object

    Returns:
        Dictionary containing normal modes analysis results including:
        - Basic normal modes
        - Enhanced mode analysis
        - Mode collectivity
        - Mode frequencies
    """
    try:
        if config is None:
            config = ProteinAnalysisConfig()

        # Initialize analyzer
        analyzer = DynamicsAnalyzer(structure, config)

        # Get basic normal modes analysis
        basic_results = analyzer.analyze_normal_modes(n_modes=n_modes)

        # Get enhanced normal modes analysis
        enhanced_results = analyzer.analyze_enhanced_modes(
            n_modes=n_modes,
            include_collectivity=True,
            include_correlations=True,
        )

        # Combine results
        results = {
            "basic": basic_results,
            "detailed": enhanced_results,
        }

        return results

    except Exception as e:
        logger.error(f"Error analyzing normal modes: {str(e)}")
        return {}


def analyze_contacts(
    structure: Structure,
    cutoff: float = 8.0,
    config: Optional[ProteinAnalysisConfig] = None,
) -> Dict[str, Any]:
    """Analyze protein contact network using multiple methods.

    Args:
        structure: BioPython Structure object
        cutoff: Distance cutoff for contacts in Angstroms
        config: Optional configuration object

    Returns:
        Dictionary containing contact analysis results including:
        - Basic contacts
        - Contact network
        - Interface analysis
        - Network metrics
    """
    try:
        if config is None:
            config = ProteinAnalysisConfig()

        # Initialize analyzer
        analyzer = DynamicsAnalyzer(structure, config)

        # Get basic contact analysis
        basic_results = analyzer.analyze_contacts(cutoff=cutoff)

        # Get enhanced contact analysis
        enhanced_results = analyzer.analyze_enhanced_contacts(
            cutoff=cutoff,
            include_interfaces=True,
            include_network=True,
        )

        # Combine results
        results = {
            "basic": basic_results,
            "detailed": enhanced_results,
        }

        return results

    except Exception as e:
        logger.error(f"Error analyzing contacts: {str(e)}")
        return {}
