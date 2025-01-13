"""Enhanced grid system optimized for scientific data visualization and compound displays."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set, Union

from .viewport import EnhancedViewportManager


class GridType(Enum):
    """Grid layout types for different visualization needs."""

    STANDARD = "standard"  # Default grid
    COMPOUND = "compound"  # Optimized for compound displays
    SCIENTIFIC = "scientific"  # For scientific data visualization
    ANALYSIS = "analysis"  # For data analysis views
    DASHBOARD = "dashboard"  # For metric dashboards


@dataclass
class GridConfig:
    """Enhanced grid configuration with scientific display options."""

    columns: int = 12
    gap: int = 16
    margin: int = 16
    container_padding: int = 16
    row_height: int = 40
    min_column_width: int = 60
    max_column_width: int = 120
    scientific_data_columns: int = 16  # Higher precision for data displays
    compound_display_columns: int = 24  # Ultra-wide for compound visualization
    data_precision_columns: int = 20  # For numerical data display
    grid_type: GridType = GridType.STANDARD


@dataclass
class GridArea:
    """Grid area specification with enhanced positioning."""

    start_col: int
    end_col: int
    start_row: int
    end_row: int
    min_height: Optional[int] = None
    max_height: Optional[int] = None
    priority: int = 0  # For responsive reordering
    sticky: bool = False  # For fixed positioning


@dataclass
class CompoundGridLayout:
    """Specialized grid layout for compound visualization."""

    structure_area: GridArea
    properties_area: GridArea
    metrics_area: GridArea
    analysis_area: GridArea
    visualization_area: GridArea
    annotations: Dict[str, GridArea] = field(default_factory=dict)


class EnhancedGridSystem:
    """Advanced grid system with scientific and compound visualization optimization."""

    def __init__(self, viewport_manager: EnhancedViewportManager):
        """Initialize enhanced grid system.

        Args:
            viewport_manager: Viewport manager instance
        """
        self.viewport = viewport_manager
        self.config = GridConfig()
        self._breakpoint_configs: Dict[str, Dict[GridType, GridConfig]] = {
            "xs": {
                GridType.STANDARD: GridConfig(columns=4, gap=8, margin=8),
                GridType.COMPOUND: GridConfig(columns=8, gap=12, margin=12, compound_display_columns=12),
                GridType.SCIENTIFIC: GridConfig(columns=8, gap=12, margin=12, scientific_data_columns=12),
                GridType.ANALYSIS: GridConfig(columns=6, gap=10, margin=10, data_precision_columns=10),
                GridType.DASHBOARD: GridConfig(columns=4, gap=8, margin=8),
            },
            "sm": {
                GridType.STANDARD: GridConfig(columns=8, gap=16, margin=16),
                GridType.COMPOUND: GridConfig(columns=12, gap=20, margin=16, compound_display_columns=16),
                GridType.SCIENTIFIC: GridConfig(columns=12, gap=20, margin=16, scientific_data_columns=16),
                GridType.ANALYSIS: GridConfig(columns=10, gap=16, margin=16, data_precision_columns=14),
                GridType.DASHBOARD: GridConfig(columns=8, gap=16, margin=16),
            },
            "md": {
                GridType.STANDARD: GridConfig(columns=12, gap=24, margin=24),
                GridType.COMPOUND: GridConfig(columns=16, gap=24, margin=24, compound_display_columns=20),
                GridType.SCIENTIFIC: GridConfig(columns=16, gap=24, margin=24, scientific_data_columns=20),
                GridType.ANALYSIS: GridConfig(columns=14, gap=20, margin=20, data_precision_columns=18),
                GridType.DASHBOARD: GridConfig(columns=12, gap=24, margin=24),
            },
            "lg": {
                GridType.STANDARD: GridConfig(columns=12, gap=32, margin=32),
                GridType.COMPOUND: GridConfig(columns=20, gap=32, margin=32, compound_display_columns=24),
                GridType.SCIENTIFIC: GridConfig(columns=20, gap=32, margin=32, scientific_data_columns=24),
                GridType.ANALYSIS: GridConfig(columns=16, gap=24, margin=24, data_precision_columns=20),
                GridType.DASHBOARD: GridConfig(columns=12, gap=32, margin=32),
            },
            "xl": {
                GridType.STANDARD: GridConfig(columns=12, gap=32, margin=32),
                GridType.COMPOUND: GridConfig(columns=24, gap=32, margin=32, compound_display_columns=28),
                GridType.SCIENTIFIC: GridConfig(columns=24, gap=32, margin=32, scientific_data_columns=28),
                GridType.ANALYSIS: GridConfig(columns=20, gap=28, margin=28, data_precision_columns=24),
                GridType.DASHBOARD: GridConfig(columns=12, gap=32, margin=32),
            },
        }
        self._compound_layouts: Dict[str, CompoundGridLayout] = {}
        self._active_optimizations: Set[str] = set()

    def get_grid_classes(self, grid_type: GridType = GridType.STANDARD) -> str:
        """Get grid classes with optimization flags.

        Args:
            grid_type: Type of grid layout

        Returns:
            Space-separated class string
        """
        config = self._get_config(grid_type)
        classes = [
            "grid",
            f"grid-cols-{config.columns}",
            f"gap-{config.gap}",
            "gpu-accelerated" if "gpu_optimization" in self._active_optimizations else "",
            "high-precision" if grid_type in [GridType.SCIENTIFIC, GridType.ANALYSIS] else "",
            "compound-optimized" if grid_type == GridType.COMPOUND else "",
        ]
        return " ".join(filter(None, classes))

    def get_container_classes(self, grid_type: GridType = GridType.STANDARD) -> str:
        """Get container classes with layout optimizations.

        Args:
            grid_type: Type of grid layout

        Returns:
            Space-separated class string
        """
        config = self._get_config(grid_type)
        classes = [
            "container",
            "mx-auto",
            f"px-{config.container_padding}",
            "scientific-container" if grid_type == GridType.SCIENTIFIC else "",
            "compound-container" if grid_type == GridType.COMPOUND else "",
            "analysis-container" if grid_type == GridType.ANALYSIS else "",
        ]
        return " ".join(filter(None, classes))

    def get_column_styles(
        self,
        span: int = 1,
        offset: int = 0,
        align: str = "stretch",
        grid_type: GridType = GridType.STANDARD,
    ) -> Dict[str, str]:
        """Get column styles with scientific optimizations.

        Args:
            span: Number of columns to span
            offset: Number of columns to offset
            align: Vertical alignment
            grid_type: Type of grid layout

        Returns:
            Dictionary of CSS properties
        """
        config = self._get_config(grid_type)
        width = self._calculate_column_width(span, config.columns)
        margin_left = self._calculate_column_width(offset, config.columns)

        styles = {
            "flex": f"0 0 {width}%",
            "max-width": f"{width}%",
            "margin-left": f"{margin_left}%",
            "padding-left": f"{config.gap / 2}px",
            "padding-right": f"{config.gap / 2}px",
            "align-self": align,
        }

        # Add scientific data optimizations
        if grid_type in [GridType.SCIENTIFIC, GridType.ANALYSIS]:
            styles.update(self._get_scientific_data_styles())

        # Add compound visualization optimizations
        if grid_type == GridType.COMPOUND:
            styles.update(self._get_compound_visualization_styles())

        return styles

    def get_compound_grid_styles(self, area_name: str) -> Dict[str, str]:
        """Get grid styles optimized for compound visualization.

        Args:
            area_name: Name of grid area

        Returns:
            Dictionary of CSS properties
        """
        layout = self._get_compound_layout()
        area = getattr(layout, f"{area_name}_area")

        styles = {
            "grid-column-start": str(area.start_col),
            "grid-column-end": str(area.end_col),
            "grid-row-start": str(area.start_row),
            "grid-row-end": str(area.end_row),
            "order": str(area.priority),
        }

        if area.sticky:
            styles.update(
                {
                    "position": "sticky",
                    "top": "0",
                    "z-index": "10",
                }
            )

        # Add area-specific optimizations
        if area_name == "structure":
            styles.update(self._get_structure_visualization_styles())
        elif area_name == "properties":
            styles.update(self._get_properties_display_styles())
        elif area_name == "metrics":
            styles.update(self._get_metrics_display_styles())
        elif area_name == "analysis":
            styles.update(self._get_analysis_display_styles())

        return styles

    def get_responsive_grid_template(
        self,
        grid_type: GridType = GridType.STANDARD,
        row_heights: Optional[List[str]] = None,
    ) -> str:
        """Get responsive grid template CSS.

        Args:
            grid_type: Type of grid layout
            row_heights: Optional list of row height values

        Returns:
            CSS grid template string
        """
        config = self._get_config(grid_type)
        template = [
            "display: grid;",
            f"grid-template-columns: repeat({config.columns}, 1fr);",
            f"gap: {config.gap}px;",
            f"margin: {config.margin}px;",
        ]

        if row_heights:
            template.append(f"grid-template-rows: {' '.join(row_heights)};")

        # Add layout-specific optimizations
        if grid_type == GridType.COMPOUND:
            template.extend(self._get_compound_grid_optimizations())
        elif grid_type in [GridType.SCIENTIFIC, GridType.ANALYSIS]:
            template.extend(self._get_scientific_grid_optimizations())

        return " ".join(template)

    def _get_config(self, grid_type: GridType) -> GridConfig:
        """Get grid configuration for current breakpoint and type."""
        return self._breakpoint_configs[self.viewport._current_breakpoint][grid_type]

    def _get_compound_layout(self) -> CompoundGridLayout:
        """Get compound layout for current breakpoint."""
        breakpoint = self.viewport._current_breakpoint
        if breakpoint not in self._compound_layouts:
            self._compound_layouts[breakpoint] = self._create_compound_layout(breakpoint)
        return self._compound_layouts[breakpoint]

    def _create_compound_layout(self, breakpoint: str) -> CompoundGridLayout:
        """Create compound layout for breakpoint.

        Args:
            breakpoint: Breakpoint name

        Returns:
            Compound grid layout configuration
        """
        config = self._get_config(GridType.COMPOUND)
        is_mobile = breakpoint in ["xs", "sm"]

        if is_mobile:
            # Stack vertically on mobile
            return CompoundGridLayout(
                structure_area=GridArea(1, config.compound_display_columns + 1, 1, 7, priority=1, sticky=True),
                properties_area=GridArea(1, config.compound_display_columns + 1, 7, 13, priority=2),
                metrics_area=GridArea(1, config.compound_display_columns + 1, 13, 17, priority=3),
                analysis_area=GridArea(1, config.compound_display_columns + 1, 17, 23, priority=4),
                visualization_area=GridArea(1, config.compound_display_columns + 1, 23, 29, priority=5),
            )
        else:
            # Side-by-side layout on larger screens
            return CompoundGridLayout(
                structure_area=GridArea(1, config.compound_display_columns // 2 + 1, 1, 13, priority=1, sticky=True),
                properties_area=GridArea(config.compound_display_columns // 2 + 1, config.compound_display_columns + 1, 1, 7, priority=2),
                metrics_area=GridArea(config.compound_display_columns // 2 + 1, config.compound_display_columns + 1, 7, 13, priority=3),
                analysis_area=GridArea(1, config.compound_display_columns // 2 + 1, 13, 23, priority=4),
                visualization_area=GridArea(config.compound_display_columns // 2 + 1, config.compound_display_columns + 1, 13, 23, priority=5),
            )

    def _calculate_column_width(self, columns: int, total_columns: int) -> float:
        """Calculate percentage width for number of columns."""
        return (columns * 100) / total_columns

    def _get_scientific_data_styles(self) -> Dict[str, str]:
        """Get styles optimized for scientific data display."""
        return {
            "font-family": "var(--data-font)",
            "font-variant-numeric": "tabular-nums",
            "line-height": "1.5",
            "letter-spacing": "0.01em",
            "white-space": "pre",
            "overflow-x": "auto",
        }

    def _get_compound_visualization_styles(self) -> Dict[str, str]:
        """Get styles optimized for compound visualization."""
        return {
            "transform": "translateZ(0)",
            "backface-visibility": "hidden",
            "perspective": "1000px",
            "transform-style": "preserve-3d",
            "will-change": "transform" if self.viewport.metrics.capabilities.supports_3d else "auto",
        }

    def _get_structure_visualization_styles(self) -> Dict[str, str]:
        """Get styles for compound structure visualization."""
        return {
            "background": "var(--structure-bg)",
            "border-radius": "8px",
            "box-shadow": "var(--structure-shadow)",
            "transition": "transform 0.2s ease-in-out",
        }

    def _get_properties_display_styles(self) -> Dict[str, str]:
        """Get styles for compound properties display."""
        return {
            "background": "var(--properties-bg)",
            "border-radius": "8px",
            "padding": "16px",
            "font-family": "var(--data-font)",
        }

    def _get_metrics_display_styles(self) -> Dict[str, str]:
        """Get styles for metrics display."""
        return {
            "background": "var(--metrics-bg)",
            "border-radius": "8px",
            "padding": "16px",
            "font-variant-numeric": "tabular-nums",
        }

    def _get_analysis_display_styles(self) -> Dict[str, str]:
        """Get styles for analysis display."""
        return {
            "background": "var(--analysis-bg)",
            "border-radius": "8px",
            "padding": "16px",
            "overflow": "auto",
        }

    def _get_compound_grid_optimizations(self) -> List[str]:
        """Get compound grid optimization styles."""
        return [
            "contain: layout style paint;",
            "content-visibility: auto;",
            "contain-intrinsic-size: 1000px;",
        ]

    def _get_scientific_grid_optimizations(self) -> List[str]:
        """Get scientific grid optimization styles."""
        return [
            "contain: layout style;",
            "text-rendering: optimizeLegibility;",
            "-webkit-font-smoothing: antialiased;",
        ]

    async def update_layout(self, breakpoint: str) -> None:
        """Update grid configuration for breakpoint.

        Args:
            breakpoint: New breakpoint name
        """
        self._compound_layouts = {}  # Clear cached layouts
        self._update_optimizations(breakpoint)

    def _update_optimizations(self, breakpoint: str) -> None:
        """Update active optimizations for breakpoint."""
        self._active_optimizations.clear()

        # Add GPU optimization for larger screens
        if breakpoint in ["lg", "xl"]:
            self._active_optimizations.add("gpu_optimization")

        # Add performance optimizations for mobile
        if breakpoint in ["xs", "sm"]:
            self._active_optimizations.add("mobile_optimization")
