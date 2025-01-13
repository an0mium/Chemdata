"""Enhanced base implementation for responsive system with scientific data support."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any


@dataclass
class ViewportConfig:
    """Viewport configuration settings."""

    width: str = "device-width"
    initial_scale: float = 1.0
    maximum_scale: float = 5.0
    user_scalable: bool = True
    safe_area_inset: bool = True
    color_scheme: str = "light dark"
    orientation: str = "portrait landscape"


@dataclass
class ComponentConfig:
    """Enhanced configuration for responsive components."""

    id: str = ""
    classes: List[str] = field(default_factory=list)
    styles: Dict[str, str] = field(default_factory=dict)
    attributes: Dict[str, str] = field(default_factory=dict)
    data_attributes: Dict[str, str] = field(default_factory=dict)
    aria_attributes: Dict[str, str] = field(default_factory=dict)
    scientific_attributes: Dict[str, Any] = field(default_factory=dict)
    performance_attributes: Dict[str, str] = field(default_factory=dict)


@dataclass
class Breakpoint:
    """Responsive breakpoint configuration."""

    name: str
    min_width: int
    max_width: Optional[int] = None
    columns: int = 12
    gap: int = 16
    margin: int = 16


class ViewportManager:
    """Enhanced viewport configuration manager."""

    def __init__(self):
        self.config = ViewportConfig()
        self._media_features: Dict[str, str] = {}

    async def configure(self, **kwargs) -> None:
        """Configure viewport settings.

        Args:
            **kwargs: Viewport configuration parameters
        """
        for key, value in kwargs.items():
            if hasattr(self.config, key):
                setattr(self.config, key, value)

    def get_meta_tag(self) -> str:
        """Get viewport meta tag with enhanced features.

        Returns:
            Complete viewport meta tag string
        """
        settings = [
            f"width={self.config.width}",
            f"initial-scale={self.config.initial_scale}",
            f"maximum-scale={self.config.maximum_scale}",
            f"user-scalable={'yes' if self.config.user_scalable else 'no'}",
            f"viewport-fit={'cover' if self.config.safe_area_inset else 'contain'}",
        ]
        return f'<meta name="viewport" content="{", ".join(settings)}">'

    def add_media_feature(self, name: str, value: str) -> None:
        """Add custom media feature.

        Args:
            name: Media feature name
            value: Media feature value
        """
        self._media_features[name] = value


class FlexibleGrid:
    """Enhanced responsive grid system."""

    def __init__(self):
        self.columns = 12
        self.gap = 16
        self.margin = 16
        self.container_padding = 16

    async def adjust(self, breakpoint: Breakpoint) -> None:
        """Adjust grid based on breakpoint.

        Args:
            breakpoint: Current breakpoint configuration
        """
        self.columns = breakpoint.columns
        self.gap = breakpoint.gap
        self.margin = breakpoint.margin

        # Adjust container padding based on breakpoint
        if breakpoint.name in ["xs", "sm"]:
            self.container_padding = 16
        else:
            self.container_padding = 24


class EnhancedResponsiveBase:
    """Enhanced base class for responsive components with scientific data support."""

    def __init__(self):
        self.viewport = ViewportManager()
        self.grid = FlexibleGrid()
        self.config = ComponentConfig()
        self._breakpoints: Dict[str, Breakpoint] = {}
        self._current_breakpoint: Optional[str] = None

    async def setup(self) -> None:
        """Set up responsive support with enhanced features."""
        # Configure viewport
        await self.viewport.configure(
            width="device-width",
            initial_scale=1.0,
            maximum_scale=5.0,
            user_scalable=True,
            safe_area_inset=True,
        )

        # Set up default breakpoints
        self._breakpoints = {
            "xs": Breakpoint("xs", 0, 575, columns=4, gap=8, margin=8),
            "sm": Breakpoint("sm", 576, 767, columns=8, gap=16, margin=16),
            "md": Breakpoint("md", 768, 991, columns=12, gap=24, margin=24),
            "lg": Breakpoint("lg", 992, 1199, columns=12, gap=32, margin=32),
            "xl": Breakpoint("xl", 1200, None, columns=12, gap=32, margin=32),
        }

        # Initialize with default breakpoint
        self._current_breakpoint = "md"
        await self.grid.adjust(self._breakpoints[self._current_breakpoint])

    async def handle_resize(self, width: int) -> None:
        """Handle viewport resize with performance optimization.

        Args:
            width: New viewport width
        """
        for breakpoint in self._breakpoints.values():
            if breakpoint.min_width <= width and (breakpoint.max_width is None or width <= breakpoint.max_width):
                if breakpoint.name != self._current_breakpoint:
                    self._current_breakpoint = breakpoint.name
                    await self.grid.adjust(breakpoint)
                break

    def get_component_classes(self) -> List[str]:
        """Get enhanced component classes.

        Returns:
            List of CSS classes
        """
        classes = self.config.classes.copy()

        # Add responsive classes
        if self._current_breakpoint and self._current_breakpoint != "xs":
            classes.append(f"{self._current_breakpoint}:block")

        # Add scientific data classes
        if self.config.scientific_attributes:
            classes.extend(
                [
                    "data-scientific",
                    "data-monospace",
                    "data-tabular-nums",
                ]
            )

        # Add performance classes
        if self.config.performance_attributes.get("gpu_accelerated"):
            classes.append("gpu-accelerated")

        return classes

    def get_component_attributes(self) -> Dict[str, str]:
        """Get enhanced component attributes.

        Returns:
            Dictionary of HTML attributes
        """
        attributes = self.config.attributes.copy()

        # Add ID if specified
        if self.config.id:
            attributes["id"] = self.config.id

        # Add data attributes
        for key, value in self.config.data_attributes.items():
            attributes[f"data-{key}"] = str(value)

        # Add ARIA attributes
        for key, value in self.config.aria_attributes.items():
            attributes[f"aria-{key}"] = str(value)

        # Add scientific data attributes
        if self.config.scientific_attributes:
            for key, value in self.config.scientific_attributes.items():
                attributes[f"data-scientific-{key}"] = str(value)

        # Add performance attributes
        if self.config.performance_attributes:
            attributes.update({"loading": "lazy", "decoding": "async", "fetchpriority": "high" if self.config.performance_attributes.get("high_priority") else "auto"})

        return attributes

    def get_component_styles(self) -> Dict[str, str]:
        """Get enhanced component styles.

        Returns:
            Dictionary of CSS properties
        """
        styles = self.config.styles.copy()

        # Add grid styles if component uses grid
        if "grid" in self.config.classes:
            styles.update(
                {
                    "grid-template-columns": f"repeat({self.grid.columns}, 1fr)",
                    "gap": f"{self.grid.gap}px",
                    "margin": f"{self.grid.margin}px",
                }
            )

        # Add scientific visualization styles
        if self.config.scientific_attributes:
            styles.update(
                {
                    "font-family": "var(--data-font)",
                    "font-variant-numeric": "tabular-nums",
                }
            )

        # Add performance styles
        if self.config.performance_attributes.get("gpu_accelerated"):
            styles.update(
                {
                    "transform": "translateZ(0)",
                    "backface-visibility": "hidden",
                    "perspective": "1000px",
                }
            )

        return styles

    def get_media_query(self, breakpoint_name: str) -> str:
        """Get media query for breakpoint.

        Args:
            breakpoint_name: Name of breakpoint

        Returns:
            Media query string

        Raises:
            ValueError: If breakpoint name is invalid
        """
        if breakpoint_name not in self._breakpoints:
            raise ValueError(f"Invalid breakpoint: {breakpoint_name}")

        breakpoint = self._breakpoints[breakpoint_name]
        if breakpoint.max_width is None:
            return f"@media (min-width: {breakpoint.min_width}px)"
        return f"@media (min-width: {breakpoint.min_width}px) and (max-width: {breakpoint.max_width}px)"

    def get_container_width(self, breakpoint_name: str) -> str:
        """Get container width for breakpoint.

        Args:
            breakpoint_name: Name of breakpoint

        Returns:
            Container width value

        Raises:
            ValueError: If breakpoint name is invalid
        """
        widths = {
            "xs": "100%",
            "sm": "540px",
            "md": "720px",
            "lg": "960px",
            "xl": "1140px",
        }
        if breakpoint_name not in widths:
            raise ValueError(f"Invalid breakpoint: {breakpoint_name}")
        return widths[breakpoint_name]
