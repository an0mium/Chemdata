"""Factory module for creating and managing structure analysis components.

This module provides:
1. Component creation and initialization
2. Dependency management
3. Component lifecycle management
4. System setup and teardown
"""

import logging
from typing import Dict, List, Optional, Set, Any, Type

from .analyzer import StructureAnalyzer
from .integration import StructureIntegrator
from .coordinator import StructureCoordinator
from .config import SystemConfig, ComponentConfig, create_default_config
from .geometry import GeometryCalculator
from .conservation import ConservationAnalyzer
from .dynamics import DynamicsAnalyzer
from .surface import SurfaceAnalyzer
from .pockets import PocketDetector
from ..pharmacophore import PharmacophoreGenerator
from ..alphafold import AlphaFoldIntegrator

logger = logging.getLogger(__name__)


class ComponentFactory:
    """Factory for creating and managing analysis components."""

    def __init__(self, config: Optional[SystemConfig] = None):
        """Initialize component factory.

        Args:
            config: System configuration
        """
        self.config = config or create_default_config()
        self._components: Dict[str, Any] = {}
        self._initialized = False

    def create_system(self) -> StructureCoordinator:
        """Create complete analysis system.

        Returns:
            Configured StructureCoordinator instance
        """
        try:
            # Validate configuration
            errors = self.config.validate()
            if errors:
                error_msg = "\n".join(errors)
                raise ValueError(f"Invalid configuration:\n{error_msg}")

            # Create components in dependency order
            self._create_core_components()
            self._create_optional_components()
            self._initialize_components()

            # Create coordinator
            coordinator = self._create_coordinator()

            self._initialized = True
            return coordinator

        except Exception as e:
            logger.error(f"Error creating system: {str(e)}")
            self.cleanup()
            raise

    def cleanup(self):
        """Clean up all components."""
        try:
            for component in self._components.values():
                if hasattr(component, "cleanup"):
                    component.cleanup()
            self._components.clear()
            self._initialized = False

        except Exception as e:
            logger.error(f"Error during cleanup: {str(e)}")

    def _create_core_components(self):
        """Create required core components."""
        try:
            # Create in dependency order
            self._create_component(
                "geometry",
                GeometryCalculator,
                self.config.geometry,
            )
            self._create_component(
                "surface",
                SurfaceAnalyzer,
                self.config.surface,
            )
            self._create_component(
                "analyzer",
                StructureAnalyzer,
                self.config.analyzer,
            )
            self._create_component(
                "pockets",
                PocketDetector,
                self.config.pockets,
            )

        except Exception as e:
            logger.error(f"Error creating core components: {str(e)}")
            raise

    def _create_optional_components(self):
        """Create enabled optional components."""
        try:
            # Conservation analysis
            if self.config.conservation and self.config.conservation.enabled:
                self._create_component(
                    "conservation",
                    ConservationAnalyzer,
                    self.config.conservation,
                )

            # Dynamics analysis
            if self.config.dynamics and self.config.dynamics.enabled:
                self._create_component(
                    "dynamics",
                    DynamicsAnalyzer,
                    self.config.dynamics,
                )

            # Pharmacophore analysis
            if self.config.pharmacophore and self.config.pharmacophore.enabled:
                self._create_component(
                    "pharmacophore",
                    PharmacophoreGenerator,
                    self.config.pharmacophore,
                )

            # AlphaFold integration
            if self.config.alphafold and self.config.alphafold.enabled:
                self._create_component(
                    "alphafold",
                    AlphaFoldIntegrator,
                    self.config.alphafold,
                )

        except Exception as e:
            logger.error(f"Error creating optional components: {str(e)}")
            raise

    def _create_component(
        self,
        name: str,
        component_class: Type,
        config: ComponentConfig,
    ):
        """Create individual component.

        Args:
            name: Component name
            component_class: Component class
            config: Component configuration
        """
        try:
            # Check dependencies
            self._check_dependencies(name, config.dependencies)

            # Create component
            component = component_class(**config.params)
            self._components[name] = component

            logger.debug(f"Created component: {name}")

        except Exception as e:
            logger.error(f"Error creating component {name}: {str(e)}")
            raise

    def _check_dependencies(self, name: str, dependencies: Set[str]):
        """Check if dependencies are satisfied.

        Args:
            name: Component name
            dependencies: Required dependencies

        Raises:
            ValueError: If dependencies are not met
        """
        missing = dependencies - set(self._components.keys())
        if missing:
            raise ValueError(f"Component {name} missing dependencies: {missing}")

    def _initialize_components(self):
        """Initialize all created components."""
        try:
            for name, component in self._components.items():
                if hasattr(component, "initialize"):
                    component.initialize()
                logger.debug(f"Initialized component: {name}")

        except Exception as e:
            logger.error(f"Error initializing components: {str(e)}")
            raise

    def _create_coordinator(self) -> StructureCoordinator:
        """Create system coordinator.

        Returns:
            Configured StructureCoordinator instance
        """
        try:
            # Create integrator first
            integrator = StructureIntegrator(
                analyzer=self._components.get("analyzer"),
                geometry=self._components.get("geometry"),
                conservation=self._components.get("conservation"),
                dynamics=self._components.get("dynamics"),
                surface=self._components.get("surface"),
                pockets=self._components.get("pockets"),
                pharmacophore=self._components.get("pharmacophore"),
                alphafold=self._components.get("alphafold"),
            )

            # Create coordinator
            coordinator = StructureCoordinator(
                integrator=integrator,
                cache_enabled=self.config.cache_enabled,
                debug_mode=self.config.debug_mode,
            )

            return coordinator

        except Exception as e:
            logger.error(f"Error creating coordinator: {str(e)}")
            raise

    @property
    def initialized(self) -> bool:
        """Check if system is initialized.

        Returns:
            True if system is initialized
        """
        return self._initialized

    def get_component(self, name: str) -> Optional[Any]:
        """Get component by name.

        Args:
            name: Component name

        Returns:
            Component instance or None if not found
        """
        return self._components.get(name)

    def get_enabled_components(self) -> Dict[str, Any]:
        """Get all enabled components.

        Returns:
            Dictionary mapping component names to instances
        """
        return self._components.copy()


def create_analysis_system(
    config: Optional[SystemConfig] = None,
) -> StructureCoordinator:
    """Create complete analysis system.

    Args:
        config: Optional system configuration

    Returns:
        Configured StructureCoordinator instance
    """
    factory = ComponentFactory(config)
    return factory.create_system()
