"""Configuration module for structure analysis system.

This module provides:
1. Centralized configuration management
2. Component initialization and dependency handling
3. System-wide settings and defaults
4. Configuration validation
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Any

from .constants import BINDING_SITE_PARAMS

logger = logging.getLogger(__name__)


@dataclass
class ComponentConfig:
    """Configuration for individual analysis components."""

    enabled: bool = True
    params: Dict[str, Any] = field(default_factory=dict)
    dependencies: Set[str] = field(default_factory=set)
    version: str = "0.1.0"


@dataclass
class SystemConfig:
    """System-wide configuration for structure analysis."""

    # Core components
    analyzer: ComponentConfig = field(default_factory=ComponentConfig)
    geometry: ComponentConfig = field(default_factory=ComponentConfig)
    surface: ComponentConfig = field(default_factory=ComponentConfig)
    pockets: ComponentConfig = field(default_factory=ComponentConfig)

    # Optional components
    conservation: Optional[ComponentConfig] = None
    dynamics: Optional[ComponentConfig] = None
    pharmacophore: Optional[ComponentConfig] = None
    alphafold: Optional[ComponentConfig] = None
    ml: Optional[ComponentConfig] = None

    # System settings
    cache_enabled: bool = True
    debug_mode: bool = False
    log_level: str = "INFO"

    def __post_init__(self):
        """Initialize with default configurations."""
        # Set default analyzer config
        if not self.analyzer.params:
            self.analyzer.params = {
                "version": "0.4.0",
                **BINDING_SITE_PARAMS,
            }

        # Set default geometry config
        if not self.geometry.params:
            self.geometry.params = {
                "distance_cutoff": 6.0,
                "neighbor_cutoff": 4,
            }

        # Set default surface config
        if not self.surface.params:
            self.surface.params = {
                "probe_radius": 1.4,
                "min_neighbors": 4,
            }

        # Set default pockets config
        if not self.pockets.params:
            self.pockets.params = {
                "min_volume": 100.0,
                "max_sites": 5,
                "score_weights": {
                    "geometry": 0.4,
                    "conservation": 0.2,
                    "dynamics": 0.2,
                    "pharmacophore": 0.2,
                },
            }

        # Initialize optional components if enabled
        if self.conservation and not self.conservation.params:
            self.conservation.params = {
                "window_size": 3,
                "gap_penalty": -1,
            }

        if self.dynamics and not self.dynamics.params:
            self.dynamics.params = {
                "modes": 10,
                "cutoff": 15.0,
            }

        if self.pharmacophore and not self.pharmacophore.params:
            self.pharmacophore.params = {
                "feature_types": [
                    "hydrophobic",
                    "aromatic",
                    "hbond_donor",
                    "hbond_acceptor",
                    "positive",
                    "negative",
                ],
            }

        if self.alphafold and not self.alphafold.params:
            self.alphafold.params = {
                "max_length": 1000,
                "num_recycles": 3,
            }

        if self.ml and not self.ml.params:
            self.ml.params = {
                "model_type": "ensemble",
                "batch_size": 32,
            }

        # Set component dependencies
        self._set_dependencies()

    def _set_dependencies(self):
        """Set up component dependencies."""
        # Analyzer depends on geometry and surface
        self.analyzer.dependencies.update({"geometry", "surface"})

        # Pockets depend on geometry and surface
        self.pockets.dependencies.update({"geometry", "surface"})

        # ML depends on all enabled components
        if self.ml and self.ml.enabled:
            self.ml.dependencies.update(
                {
                    "analyzer",
                    "geometry",
                    "surface",
                    "pockets",
                }
            )
            if self.conservation and self.conservation.enabled:
                self.ml.dependencies.add("conservation")
            if self.dynamics and self.dynamics.enabled:
                self.ml.dependencies.add("dynamics")
            if self.pharmacophore and self.pharmacophore.enabled:
                self.ml.dependencies.add("pharmacophore")

    def validate(self) -> List[str]:
        """Validate configuration.

        Returns:
            List of validation error messages
        """
        errors = []

        # Check required components
        if not self.analyzer.enabled:
            errors.append("Analyzer component must be enabled")
        if not self.geometry.enabled:
            errors.append("Geometry component must be enabled")
        if not self.surface.enabled:
            errors.append("Surface component must be enabled")
        if not self.pockets.enabled:
            errors.append("Pockets component must be enabled")

        # Check dependencies
        enabled_components = self._get_enabled_components()
        for component in enabled_components:
            missing_deps = component.dependencies - {c.name for c in enabled_components}
            if missing_deps:
                errors.append(f"Component {component.name} missing dependencies: {missing_deps}")

        # Validate component-specific parameters
        self._validate_analyzer_params(errors)
        self._validate_geometry_params(errors)
        self._validate_surface_params(errors)
        self._validate_pockets_params(errors)
        self._validate_optional_components(errors)

        return errors

    def _get_enabled_components(self) -> List[ComponentConfig]:
        """Get list of enabled components.

        Returns:
            List of enabled component configurations
        """
        components = [
            self.analyzer,
            self.geometry,
            self.surface,
            self.pockets,
        ]
        if self.conservation and self.conservation.enabled:
            components.append(self.conservation)
        if self.dynamics and self.dynamics.enabled:
            components.append(self.dynamics)
        if self.pharmacophore and self.pharmacophore.enabled:
            components.append(self.pharmacophore)
        if self.alphafold and self.alphafold.enabled:
            components.append(self.alphafold)
        if self.ml and self.ml.enabled:
            components.append(self.ml)
        return components

    def _validate_analyzer_params(self, errors: List[str]):
        """Validate analyzer parameters."""
        params = self.analyzer.params
        if "version" not in params:
            errors.append("Analyzer version must be specified")
        if not isinstance(params.get("version", ""), str):
            errors.append("Analyzer version must be a string")

    def _validate_geometry_params(self, errors: List[str]):
        """Validate geometry parameters."""
        params = self.geometry.params
        if "distance_cutoff" not in params:
            errors.append("Geometry distance_cutoff must be specified")
        if not isinstance(params.get("distance_cutoff", 0), (int, float)):
            errors.append("Geometry distance_cutoff must be a number")
        if params.get("distance_cutoff", 0) <= 0:
            errors.append("Geometry distance_cutoff must be positive")

    def _validate_surface_params(self, errors: List[str]):
        """Validate surface parameters."""
        params = self.surface.params
        if "probe_radius" not in params:
            errors.append("Surface probe_radius must be specified")
        if not isinstance(params.get("probe_radius", 0), (int, float)):
            errors.append("Surface probe_radius must be a number")
        if params.get("probe_radius", 0) <= 0:
            errors.append("Surface probe_radius must be positive")

    def _validate_pockets_params(self, errors: List[str]):
        """Validate pockets parameters."""
        params = self.pockets.params
        if "min_volume" not in params:
            errors.append("Pockets min_volume must be specified")
        if not isinstance(params.get("min_volume", 0), (int, float)):
            errors.append("Pockets min_volume must be a number")
        if params.get("min_volume", 0) <= 0:
            errors.append("Pockets min_volume must be positive")

        weights = params.get("score_weights", {})
        if not isinstance(weights, dict):
            errors.append("Pockets score_weights must be a dictionary")
        if sum(weights.values()) != 1.0:
            errors.append("Pockets score weights must sum to 1.0")

    def _validate_optional_components(self, errors: List[str]):
        """Validate optional component parameters."""
        # Validate conservation
        if self.conservation and self.conservation.enabled:
            params = self.conservation.params
            if "window_size" not in params:
                errors.append("Conservation window_size must be specified")
            if not isinstance(params.get("window_size", 0), int):
                errors.append("Conservation window_size must be an integer")
            if params.get("window_size", 0) <= 0:
                errors.append("Conservation window_size must be positive")

        # Validate dynamics
        if self.dynamics and self.dynamics.enabled:
            params = self.dynamics.params
            if "modes" not in params:
                errors.append("Dynamics modes must be specified")
            if not isinstance(params.get("modes", 0), int):
                errors.append("Dynamics modes must be an integer")
            if params.get("modes", 0) <= 0:
                errors.append("Dynamics modes must be positive")

        # Validate pharmacophore
        if self.pharmacophore and self.pharmacophore.enabled:
            params = self.pharmacophore.params
            if "feature_types" not in params:
                errors.append("Pharmacophore feature_types must be specified")
            if not isinstance(params.get("feature_types", []), list):
                errors.append("Pharmacophore feature_types must be a list")
            if not params.get("feature_types", []):
                errors.append("Pharmacophore feature_types cannot be empty")

        # Validate AlphaFold
        if self.alphafold and self.alphafold.enabled:
            params = self.alphafold.params
            if "max_length" not in params:
                errors.append("AlphaFold max_length must be specified")
            if not isinstance(params.get("max_length", 0), int):
                errors.append("AlphaFold max_length must be an integer")
            if params.get("max_length", 0) <= 0:
                errors.append("AlphaFold max_length must be positive")

        # Validate ML
        if self.ml and self.ml.enabled:
            params = self.ml.params
            if "model_type" not in params:
                errors.append("ML model_type must be specified")
            if not isinstance(params.get("model_type", ""), str):
                errors.append("ML model_type must be a string")
            if "batch_size" not in params:
                errors.append("ML batch_size must be specified")
            if not isinstance(params.get("batch_size", 0), int):
                errors.append("ML batch_size must be an integer")
            if params.get("batch_size", 0) <= 0:
                errors.append("ML batch_size must be positive")


def create_default_config() -> SystemConfig:
    """Create default system configuration.

    Returns:
        Default SystemConfig instance
    """
    return SystemConfig(
        analyzer=ComponentConfig(
            enabled=True,
            params={
                "version": "0.4.0",
                **BINDING_SITE_PARAMS,
            },
        ),
        geometry=ComponentConfig(
            enabled=True,
            params={
                "distance_cutoff": 6.0,
                "neighbor_cutoff": 4,
            },
        ),
        surface=ComponentConfig(
            enabled=True,
            params={
                "probe_radius": 1.4,
                "min_neighbors": 4,
            },
        ),
        pockets=ComponentConfig(
            enabled=True,
            params={
                "min_volume": 100.0,
                "max_sites": 5,
                "score_weights": {
                    "geometry": 0.4,
                    "conservation": 0.2,
                    "dynamics": 0.2,
                    "pharmacophore": 0.2,
                },
            },
        ),
        conservation=ComponentConfig(
            enabled=True,
            params={
                "window_size": 3,
                "gap_penalty": -1,
            },
        ),
        dynamics=ComponentConfig(
            enabled=True,
            params={
                "modes": 10,
                "cutoff": 15.0,
            },
        ),
        pharmacophore=ComponentConfig(
            enabled=True,
            params={
                "feature_types": [
                    "hydrophobic",
                    "aromatic",
                    "hbond_donor",
                    "hbond_acceptor",
                    "positive",
                    "negative",
                ],
            },
        ),
        alphafold=ComponentConfig(
            enabled=True,
            params={
                "max_length": 1000,
                "num_recycles": 3,
            },
        ),
        ml=ComponentConfig(
            enabled=True,
            params={
                "model_type": "ensemble",
                "batch_size": 32,
            },
        ),
    )


def load_config(config_dict: Dict[str, Any]) -> SystemConfig:
    """Load configuration from dictionary.

    Args:
        config_dict: Configuration dictionary

    Returns:
        SystemConfig instance
    """
    try:
        config = SystemConfig()

        # Update component configs
        for component in [
            "analyzer",
            "geometry",
            "surface",
            "pockets",
            "conservation",
            "dynamics",
            "pharmacophore",
            "alphafold",
            "ml",
        ]:
            if component in config_dict:
                comp_config = config_dict[component]
                if hasattr(config, component):
                    setattr(
                        config,
                        component,
                        ComponentConfig(
                            enabled=comp_config.get("enabled", True),
                            params=comp_config.get("params", {}),
                            dependencies=set(comp_config.get("dependencies", [])),
                            version=comp_config.get("version", "0.1.0"),
                        ),
                    )

        # Update system settings
        if "settings" in config_dict:
            settings = config_dict["settings"]
            config.cache_enabled = settings.get("cache_enabled", True)
            config.debug_mode = settings.get("debug_mode", False)
            config.log_level = settings.get("log_level", "INFO")

        return config

    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return create_default_config()
