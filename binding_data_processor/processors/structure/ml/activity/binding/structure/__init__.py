"""Structure analysis module for protein binding sites with enhanced capabilities.

This module provides comprehensive protein structure analysis for binding site detection,
including:
1. Geometric analysis and shape complementarity
2. Conservation analysis and evolutionary features  
3. Dynamics analysis and flexibility prediction
4. Surface analysis and pocket detection
5. Integration with AlphaFold for structure prediction
6. Integration with pharmacophore detection
7. Machine learning-based scoring and prediction
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union, Any

import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue

from .analyzer import StructureAnalyzer
from .geometry import GeometryCalculator
from .conservation import ConservationAnalyzer
from .dynamics import DynamicsAnalyzer
from .surface import SurfaceAnalyzer
from .pockets import PocketDetector
from ..pharmacophore import PharmacophoreGenerator
from ..alphafold import AlphaFoldIntegrator

logger = logging.getLogger(__name__)

# Version info
__version__ = "0.3.0"
__author__ = "Cline"
__description__ = "Advanced protein structure analysis for binding site detection"

# Feature flags
FEATURES = {
    "use_alphafold": True,  # Use AlphaFold integration when available
    "use_ml": True,  # Use machine learning models when available
    "use_conservation": True,  # Include conservation analysis
    "use_dynamics": True,  # Include dynamics analysis
    "use_shape_analysis": True,  # Include molecular shape analysis
    "use_pharmacophore": True,  # Include pharmacophore detection
}

# Default parameters
DEFAULTS = {
    "probe_radius": 1.4,  # Probe radius for surface calculations (Å)
    "min_pocket_volume": 100.0,  # Minimum pocket volume (Å³)
    "max_pockets": 5,  # Maximum number of pockets to return
    "distance_cutoff": 4.0,  # Distance cutoff for residue contacts (Å)
    "conservation_cutoff": 0.7,  # Minimum conservation score
    "druggability_threshold": 0.5,  # Minimum druggability score
    "shape_weight": 0.3,  # Weight for shape complementarity
    "energy_weight": 0.4,  # Weight for interaction energy
    "conservation_weight": 0.3,  # Weight for conservation score
}

# Analysis components
COMPONENTS = {
    "geometry": {
        "enabled": True,
        "description": "Basic geometric analysis",
        "class": GeometryCalculator,
        "dependencies": [],
    },
    "conservation": {
        "enabled": True,
        "description": "Sequence conservation analysis",
        "class": ConservationAnalyzer,
        "dependencies": ["geometry"],
    },
    "dynamics": {
        "enabled": True,
        "description": "Protein dynamics analysis",
        "class": DynamicsAnalyzer,
        "dependencies": ["geometry"],
    },
    "surface": {
        "enabled": True,
        "description": "Surface property analysis",
        "class": SurfaceAnalyzer,
        "dependencies": ["geometry"],
    },
    "pockets": {
        "enabled": True,
        "description": "Binding pocket detection",
        "class": PocketDetector,
        "dependencies": ["geometry", "surface"],
    },
}

# Enhanced capabilities
CAPABILITIES = {
    "structure_prediction": {
        "alphafold": True,
        "rosetta": False,
    },
    "pocket_detection": {
        "geometric": True,
        "energy_based": True,
        "ml_based": True,
        "consensus": True,
    },
    "analysis": {
        "conservation": True,
        "dynamics": True,
        "electrostatics": True,
        "hydrophobicity": True,
        "shape": True,
        "pharmacophore": True,
    },
}


class StructureManager:
    """Manages protein structure analysis with enhanced integration."""

    def __init__(
        self,
        use_alphafold: bool = FEATURES["use_alphafold"],
        use_ml: bool = FEATURES["use_ml"],
        use_conservation: bool = FEATURES["use_conservation"],
        use_dynamics: bool = FEATURES["use_dynamics"],
        use_shape_analysis: bool = FEATURES["use_shape_analysis"],
        use_pharmacophore: bool = FEATURES["use_pharmacophore"],
        **kwargs,
    ):
        """Initialize structure manager with specified capabilities.

        Args:
            use_alphafold: Whether to use AlphaFold integration
            use_ml: Whether to use machine learning models
            use_conservation: Whether to include conservation analysis
            use_dynamics: Whether to include dynamics analysis
            use_shape_analysis: Whether to include shape analysis
            use_pharmacophore: Whether to include pharmacophore detection
            **kwargs: Additional arguments passed to analyzers
        """
        # Update feature flags
        self.features = FEATURES.copy()
        self.features.update(
            {
                "use_alphafold": use_alphafold,
                "use_ml": use_ml,
                "use_conservation": use_conservation,
                "use_dynamics": use_dynamics,
                "use_shape_analysis": use_shape_analysis,
                "use_pharmacophore": use_pharmacophore,
            }
        )

        # Initialize core analyzers
        self.structure_analyzer = StructureAnalyzer(**kwargs)
        self.geometry_calculator = GeometryCalculator(**kwargs)
        self.surface_analyzer = SurfaceAnalyzer(**kwargs)
        self.pocket_detector = PocketDetector(**kwargs)

        # Optional analyzers based on features
        self.conservation_analyzer = ConservationAnalyzer(**kwargs) if use_conservation else None
        self.dynamics_analyzer = DynamicsAnalyzer(**kwargs) if use_dynamics else None

        # Optional integrations
        self.pharmacophore_generator = PharmacophoreGenerator() if use_pharmacophore else None
        self.alphafold_integrator = AlphaFoldIntegrator() if use_alphafold else None

    def analyze_structure(
        self,
        structure: Structure,
        include_dynamics: bool = True,
        include_conservation: bool = True,
        include_pharmacophore: bool = True,
    ) -> Dict[str, Any]:
        """Analyze protein structure with all enabled components.

        Args:
            structure: BioPython Structure object
            include_dynamics: Whether to include dynamics analysis
            include_conservation: Whether to include conservation analysis
            include_pharmacophore: Whether to include pharmacophore analysis

        Returns:
            Dictionary of structure properties
        """
        try:
            # Basic structure properties
            properties = self.structure_analyzer.get_structure_properties(structure)

            # Core analyses
            properties.update(
                {
                    "geometry": self.geometry_calculator.analyze_geometry(structure),
                    "surface": self.surface_analyzer.analyze_surface(structure),
                    "pockets": self.pocket_detector.detect_pockets(structure),
                }
            )

            # Optional analyses based on features
            if include_conservation and self.conservation_analyzer:
                properties["conservation"] = self.conservation_analyzer.analyze_conservation(structure)

            if include_dynamics and self.dynamics_analyzer:
                properties["dynamics"] = self.dynamics_analyzer.analyze_dynamics(structure)

            if include_pharmacophore and self.pharmacophore_generator:
                properties["pharmacophore"] = self.pharmacophore_generator.get_structure_features(structure)

            return properties

        except Exception as e:
            logger.error(f"Error analyzing structure: {str(e)}")
            return {}

    def analyze_binding_site(
        self,
        structure: Structure,
        site_residues: List[int],
        include_pharmacophore: bool = True,
    ) -> Dict[str, Any]:
        """Analyze binding site properties with all enabled components.

        Args:
            structure: BioPython Structure object
            site_residues: List of binding site residue numbers
            include_pharmacophore: Whether to include pharmacophore analysis

        Returns:
            Dictionary of binding site properties
        """
        try:
            # Get site residues
            residues = [res for res in structure.get_residues() if res.get_id()[1] in site_residues]
            if not residues:
                return {}

            # Basic properties
            properties = {
                "residues": site_residues,
                "center": self.geometry_calculator.get_center(residues),
                "volume": self.geometry_calculator.get_volume(residues),
            }

            # Core analyses
            properties.update(
                {
                    "geometry": self.geometry_calculator.analyze_site_geometry(residues),
                    "surface": self.surface_analyzer.analyze_site_surface(structure, residues),
                }
            )

            # Optional analyses based on features
            if self.conservation_analyzer:
                properties["conservation"] = self.conservation_analyzer.get_site_conservation(site_residues)

            if self.dynamics_analyzer:
                properties["dynamics"] = self.dynamics_analyzer.analyze_site_dynamics(structure, residues)

            if include_pharmacophore and self.pharmacophore_generator:
                properties["pharmacophore"] = self.pharmacophore_generator.get_site_features(structure, site_residues)

            return properties

        except Exception as e:
            logger.error(f"Error analyzing binding site: {str(e)}")
            return {}

    def predict_structure(
        self,
        sequence: str,
        **kwargs,
    ) -> Optional[Structure]:
        """Predict protein structure using AlphaFold integration.

        Args:
            sequence: Protein sequence
            **kwargs: Additional arguments passed to AlphaFold

        Returns:
            Predicted structure or None if prediction fails
        """
        try:
            if not self.alphafold_integrator:
                logger.warning("AlphaFold integration not enabled")
                return None

            return self.alphafold_integrator.predict_structure(sequence, **kwargs)

        except Exception as e:
            logger.error(f"Error predicting structure: {str(e)}")
            return None

    def get_enabled_components(self) -> Dict[str, str]:
        """Get currently enabled analysis components.

        Returns:
            Dictionary mapping component names to descriptions
        """
        return {name: info["description"] for name, info in COMPONENTS.items() if info["enabled"]}

    def get_version_info(self) -> Dict[str, Any]:
        """Get version and configuration information.

        Returns:
            Dictionary containing version info and enabled features
        """
        versions = {
            "version": __version__,
            "author": __author__,
            "description": __description__,
            "features": self.features,
            "components": self.get_enabled_components(),
            "capabilities": CAPABILITIES,
        }

        # Add component versions
        component_versions = {
            "structure_analyzer": self.structure_analyzer.version,
            "geometry_calculator": self.geometry_calculator.version,
            "surface_analyzer": self.surface_analyzer.version,
            "pocket_detector": self.pocket_detector.version,
        }

        if self.conservation_analyzer:
            component_versions["conservation_analyzer"] = self.conservation_analyzer.version
        if self.dynamics_analyzer:
            component_versions["dynamics_analyzer"] = self.dynamics_analyzer.version
        if self.pharmacophore_generator:
            component_versions["pharmacophore"] = self.pharmacophore_generator.version
        if self.alphafold_integrator:
            component_versions["alphafold"] = self.alphafold_integrator.version

        versions["component_versions"] = component_versions
        return versions


def create_manager(**kwargs) -> StructureManager:
    """Create a configured StructureManager instance.

    Args:
        **kwargs: Arguments passed to StructureManager constructor

    Returns:
        Configured StructureManager instance
    """
    return StructureManager(**kwargs)


__all__ = [
    # Main classes
    "StructureManager",
    "StructureAnalyzer",
    "GeometryCalculator",
    "ConservationAnalyzer",
    "DynamicsAnalyzer",
    "SurfaceAnalyzer",
    "PocketDetector",
    # Factory function
    "create_manager",
    # Constants and info
    "FEATURES",
    "DEFAULTS",
    "COMPONENTS",
    "CAPABILITIES",
    "__version__",
    "__author__",
    "__description__",
]
