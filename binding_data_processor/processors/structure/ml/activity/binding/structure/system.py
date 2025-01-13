"""System module for initializing and managing the structure analysis system.

This module provides:
1. High-level interface for users
2. Configuration management
3. System lifecycle management
4. Error handling and logging
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

from Bio.PDB.Structure import Structure

from .analyzer import StructureAnalyzer
from .integration import StructureIntegrator
from .coordinator import StructureCoordinator
from .config import SystemConfig, create_default_config, load_config
from .factory import ComponentFactory

logger = logging.getLogger(__name__)


class StructureAnalysisSystem:
    """Main interface for structure analysis system."""

    def __init__(
        self,
        config: Optional[Union[SystemConfig, Dict[str, Any], str, Path]] = None,
        debug: bool = False,
    ):
        """Initialize structure analysis system.

        Args:
            config: Configuration as SystemConfig object, dict, or path to JSON file
            debug: Whether to enable debug mode
        """
        self.config = self._load_config(config)
        if debug:
            self.config.debug_mode = True
            logging.getLogger().setLevel(logging.DEBUG)

        self.factory = ComponentFactory(self.config)
        self.coordinator = None
        self._initialized = False

    def initialize(self) -> bool:
        """Initialize the system.

        Returns:
            True if initialization successful
        """
        try:
            if self._initialized:
                logger.warning("System already initialized")
                return True

            # Create system through factory
            self.coordinator = self.factory.create_system()
            self._initialized = True
            logger.info("Structure analysis system initialized successfully")
            return True

        except Exception as e:
            logger.error(f"Error initializing system: {str(e)}")
            self.cleanup()
            return False

    def cleanup(self):
        """Clean up system resources."""
        try:
            if self.factory:
                self.factory.cleanup()
            self.coordinator = None
            self._initialized = False
            logger.info("System cleanup completed")

        except Exception as e:
            logger.error(f"Error during cleanup: {str(e)}")

    def analyze_structure(
        self,
        structure: Structure,
        include_conservation: bool = True,
        include_dynamics: bool = True,
        include_pharmacophore: bool = True,
    ) -> Dict[str, Any]:
        """Analyze protein structure.

        Args:
            structure: BioPython Structure object
            include_conservation: Whether to include conservation analysis
            include_dynamics: Whether to include dynamics analysis
            include_pharmacophore: Whether to include pharmacophore features

        Returns:
            Dictionary of structure properties
        """
        self._ensure_initialized()
        return self.coordinator.analyze_structure(
            structure,
            include_conservation=include_conservation,
            include_dynamics=include_dynamics,
            include_pharmacophore=include_pharmacophore,
        )

    def find_binding_sites(
        self,
        structure: Structure,
        properties: Optional[Dict[str, Any]] = None,
        min_volume: float = 100.0,
        max_sites: int = 5,
        include_scores: bool = True,
        include_pharmacophore: bool = True,
    ) -> List[Dict[str, Any]]:
        """Find and analyze binding sites.

        Args:
            structure: BioPython Structure object
            properties: Pre-calculated structure properties
            min_volume: Minimum pocket volume in Å³
            max_sites: Maximum number of sites to return
            include_scores: Whether to include detailed scoring
            include_pharmacophore: Whether to include pharmacophore analysis

        Returns:
            List of dictionaries containing binding site properties
        """
        self._ensure_initialized()
        return self.coordinator.find_binding_sites(
            structure,
            properties=properties,
            min_volume=min_volume,
            max_sites=max_sites,
            include_scores=include_scores,
            include_pharmacophore=include_pharmacophore,
        )

    def predict_structure(
        self,
        sequence: str,
        analyze: bool = True,
        **kwargs,
    ) -> Structure:
        """Predict protein structure using AlphaFold.

        Args:
            sequence: Protein sequence
            analyze: Whether to analyze predicted structure
            **kwargs: Additional arguments passed to AlphaFold

        Returns:
            Predicted BioPython Structure object
        """
        self._ensure_initialized()
        return self.coordinator.predict_structure(
            sequence,
            analyze=analyze,
            **kwargs,
        )

    def get_enabled_components(self) -> Dict[str, str]:
        """Get currently enabled analysis components.

        Returns:
            Dictionary mapping component names to descriptions
        """
        self._ensure_initialized()
        return self.coordinator.get_enabled_components()

    def get_version_info(self) -> Dict[str, Any]:
        """Get version and configuration information.

        Returns:
            Dictionary containing version info and enabled features
        """
        info = {
            "system_version": "0.4.0",
            "initialized": self._initialized,
            "debug_mode": self.config.debug_mode,
        }
        if self._initialized:
            info.update(self.coordinator.get_version_info())
        return info

    def _load_config(
        self,
        config: Optional[Union[SystemConfig, Dict[str, Any], str, Path]],
    ) -> SystemConfig:
        """Load configuration from various sources.

        Args:
            config: Configuration source

        Returns:
            SystemConfig instance
        """
        try:
            if config is None:
                return create_default_config()

            if isinstance(config, SystemConfig):
                return config

            if isinstance(config, (str, Path)):
                with open(config) as f:
                    config_dict = json.load(f)
                return load_config(config_dict)

            if isinstance(config, dict):
                return load_config(config)

            raise ValueError(f"Invalid config type: {type(config)}")

        except Exception as e:
            logger.error(f"Error loading configuration: {str(e)}")
            return create_default_config()

    def _ensure_initialized(self):
        """Ensure system is initialized.

        Raises:
            RuntimeError: If system not initialized
        """
        if not self._initialized:
            raise RuntimeError("System not initialized. Call initialize() first.")


def create_analysis_system(
    config: Optional[Union[SystemConfig, Dict[str, Any], str, Path]] = None,
    debug: bool = False,
) -> StructureAnalysisSystem:
    """Create and initialize structure analysis system.

    Args:
        config: Optional configuration
        debug: Whether to enable debug mode

    Returns:
        Initialized StructureAnalysisSystem instance
    """
    system = StructureAnalysisSystem(config, debug=debug)
    if not system.initialize():
        raise RuntimeError("Failed to initialize structure analysis system")
    return system
