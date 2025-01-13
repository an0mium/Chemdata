"""Coordinator module for integrating structure analysis components.

This module provides:
1. High-level coordination between analysis components
2. Unified interfaces for structure analysis
3. Efficient data sharing and caching
4. Error handling and validation
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union, Any

import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from .analyzer import StructureAnalyzer
from .integration import StructureIntegrator
from ..pharmacophore import PharmacophoreGenerator
from ..alphafold import AlphaFoldIntegrator
from ..constants import BINDING_SITE_PARAMS

logger = logging.getLogger(__name__)


class StructureCoordinator:
    """High-level coordinator for structure analysis components."""

    def __init__(
        self,
        use_alphafold: bool = True,
        use_ml: bool = True,
        use_conservation: bool = True,
        use_dynamics: bool = True,
        use_shape_analysis: bool = True,
        use_pharmacophore: bool = True,
        **kwargs,
    ):
        """Initialize structure coordinator.

        Args:
            use_alphafold: Whether to use AlphaFold integration
            use_ml: Whether to use machine learning models
            use_conservation: Whether to include conservation analysis
            use_dynamics: Whether to include dynamics analysis
            use_shape_analysis: Whether to include shape analysis
            use_pharmacophore: Whether to include pharmacophore detection
            **kwargs: Additional configuration parameters
        """
        # Initialize integrator
        self.integrator = StructureIntegrator(
            use_alphafold=use_alphafold,
            use_ml=use_ml,
            use_conservation=use_conservation,
            use_dynamics=use_dynamics,
            use_shape_analysis=use_shape_analysis,
            use_pharmacophore=use_pharmacophore,
            **kwargs,
        )

        # Store settings
        self.settings = {
            "use_alphafold": use_alphafold,
            "use_ml": use_ml,
            "use_conservation": use_conservation,
            "use_dynamics": use_dynamics,
            "use_shape_analysis": use_shape_analysis,
            "use_pharmacophore": use_pharmacophore,
        }

        # Initialize cache
        self._cache = {}

    def analyze_structure(
        self,
        structure: Structure,
        include_conservation: bool = True,
        include_dynamics: bool = True,
        include_pharmacophore: bool = True,
        cache_results: bool = True,
    ) -> Dict[str, Any]:
        """Analyze protein structure with all enabled components.

        Args:
            structure: BioPython Structure object
            include_conservation: Whether to include conservation analysis
            include_dynamics: Whether to include dynamics analysis
            include_pharmacophore: Whether to include pharmacophore features
            cache_results: Whether to cache results for reuse

        Returns:
            Dictionary of structure properties
        """
        try:
            # Check cache first
            cache_key = self._get_cache_key(structure)
            if cache_key in self._cache and cache_results:
                return self._cache[cache_key]

            # Perform analysis
            properties = self.integrator.analyze_structure(
                structure,
                include_conservation=include_conservation,
                include_dynamics=include_dynamics,
                include_pharmacophore=include_pharmacophore,
            )

            # Cache results if requested
            if cache_results:
                self._cache[cache_key] = properties

            return properties

        except Exception as e:
            logger.error(f"Error in structure analysis: {str(e)}")
            return {}

    def find_binding_sites(
        self,
        structure: Structure,
        properties: Optional[Dict[str, Any]] = None,
        min_volume: float = 100.0,
        max_sites: int = 5,
        include_scores: bool = True,
        include_pharmacophore: bool = True,
        reuse_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """Find and analyze binding sites.

        Args:
            structure: BioPython Structure object
            properties: Pre-calculated structure properties
            min_volume: Minimum pocket volume in Å³
            max_sites: Maximum number of sites to return
            include_scores: Whether to include detailed scoring
            include_pharmacophore: Whether to include pharmacophore analysis
            reuse_cache: Whether to reuse cached results

        Returns:
            List of dictionaries containing binding site properties
        """
        try:
            # Get or calculate properties
            if properties is None:
                cache_key = self._get_cache_key(structure)
                if reuse_cache and cache_key in self._cache:
                    properties = self._cache[cache_key]
                else:
                    properties = self.analyze_structure(structure)

            # Find binding sites
            sites = self.integrator.find_binding_sites(
                structure,
                properties=properties,
                min_volume=min_volume,
                max_sites=max_sites,
                include_scores=include_scores,
                include_pharmacophore=include_pharmacophore,
            )

            return sites

        except Exception as e:
            logger.error(f"Error finding binding sites: {str(e)}")
            return []

    def predict_structure(
        self,
        sequence: str,
        analyze: bool = True,
        **kwargs,
    ) -> Tuple[Optional[Structure], Optional[Dict[str, Any]]]:
        """Predict and analyze protein structure.

        Args:
            sequence: Protein sequence
            analyze: Whether to analyze predicted structure
            **kwargs: Additional arguments passed to AlphaFold

        Returns:
            Tuple of (predicted structure, analysis results)
        """
        try:
            # Predict structure
            structure = self.integrator.predict_structure(sequence, **kwargs)
            if structure is None:
                return None, None

            # Analyze if requested
            properties = None
            if analyze:
                properties = self.analyze_structure(structure)

            return structure, properties

        except Exception as e:
            logger.error(f"Error predicting structure: {str(e)}")
            return None, None

    def get_enabled_components(self) -> Dict[str, str]:
        """Get currently enabled analysis components.

        Returns:
            Dictionary mapping component names to descriptions
        """
        return self.integrator.get_enabled_components()

    def get_version_info(self) -> Dict[str, Any]:
        """Get version and configuration information.

        Returns:
            Dictionary containing version info and enabled features
        """
        return {
            "coordinator_version": "0.4.0",
            "settings": self.settings,
            "integrator": self.integrator.get_version_info(),
        }

    def clear_cache(self) -> None:
        """Clear cached results."""
        self._cache.clear()

    def _get_cache_key(self, structure: Structure) -> str:
        """Generate cache key for structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Cache key string
        """
        try:
            # Use structure ID and number of atoms as key
            structure_id = structure.get_id()
            num_atoms = len(list(structure.get_atoms()))
            return f"{structure_id}_{num_atoms}"

        except Exception as e:
            logger.error(f"Error generating cache key: {str(e)}")
            return str(id(structure))  # Fallback to object ID
