"""Enhanced responsive layout system optimized for scientific data visualization and compound displays."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set, Tuple, Union, Any

from .base import EnhancedResponsiveBase
from .viewport import ViewportManager, EnhancedViewportManager
from .grid import FlexibleGrid, EnhancedGridSystem, GridType


class LayoutType(Enum):
    """Layout types for different content displays."""

    STANDARD = "standard"  # Default layout
    COMPOUND = "compound"  # For compound visualization
    SCIENTIFIC = "scientific"  # For scientific data
    ANALYSIS = "analysis"  # For data analysis
    DASHBOARD = "dashboard"  # For metrics display
    QUANTUM = "quantum"  # For quantum visualization
    MOLECULAR = "molecular"  # For molecular structure
    COMPARISON = "comparison"  # For compound comparison


@dataclass
class LayoutConfig:
    """Enhanced layout configuration with scientific display options."""

    # Grid settings
    columns: int = 12
    gap: int = 16
    margin: int = 16
    container_padding: int = 16

    # Breakpoint-specific settings
    xs_columns: int = 4
    sm_columns: int = 8
    md_columns: int = 12
    lg_columns: int = 12
    xl_columns: int = 12

    # Component sizes
    structure_width: int = 200
    structure_height: int = 200
    chart_width: int = 400
    chart_height: int = 300

    # Touch targets
    min_touch_target: int = 44
    touch_spacing: int = 16

    # Typography
    base_font_size: int = 16
    data_font: str = "var(--data-font)"
    heading_scale: float = 1.2

    # Layout type and dimensions
    layout_type: LayoutType = LayoutType.STANDARD
    padding: Tuple[int, int, int, int] = (16, 16, 16, 16)  # top, right, bottom, left
    sidebar_width: int = 280
    content_max_width: int = 1200
    data_panel_height: int = 600
    visualization_height: int = 480
    controls_height: int = 80
    min_height: Optional[int] = None

    # Enhanced features
    enable_animations: bool = True
    enable_transitions: bool = True
    high_performance: bool = False
    scientific_precision: bool = False
    quantum_enabled: bool = False
    molecular_detail: str = "high"  # low, medium, high
    comparison_mode: str = "side-by-side"  # side-by-side, overlay, split


@dataclass
class LayoutArea:
    """Layout area specification with enhanced positioning."""

    name: str
    grid_area: str
    min_height: Optional[int] = None
    max_height: Optional[int] = None
    priority: int = 0
    sticky: bool = False
    collapsible: bool = False
    scrollable: bool = True
    resizable: bool = False
    draggable: bool = False
    quantum_ready: bool = False
    molecular_view: bool = False
    scientific_data: bool = False


class EnhancedResponsiveLayout(EnhancedResponsiveBase):
    """Enhanced responsive layout manager with scientific optimization."""

    def __init__(self):
        """Initialize enhanced layout manager."""
        super().__init__()
        self.config = LayoutConfig()
        self.viewport = EnhancedViewportManager()
        self.grid = EnhancedGridSystem(self.viewport)
        self._active_areas: Dict[str, LayoutArea] = {}
        self._active_optimizations: Set[str] = set()
        self._breakpoint_layouts: Dict[str, LayoutConfig] = self._initialize_breakpoint_layouts()

    def _initialize_breakpoint_layouts(self) -> Dict[str, LayoutConfig]:
        """Initialize breakpoint-specific layout configurations."""
        return {
            "xs": LayoutConfig(
                columns=4,
                gap=8,
                margin=8,
                container_padding=16,
                content_max_width=None,
                min_height=300,
                visualization_height=320,
                controls_height=60,
                padding=(8, 8, 8, 8),
            ),
            "sm": LayoutConfig(
                columns=8,
                gap=16,
                margin=16,
                container_padding=16,
                content_max_width=540,
                min_height=400,
                visualization_height=400,
                controls_height=70,
                padding=(12, 12, 12, 12),
            ),
            "md": LayoutConfig(
                columns=12,
                gap=24,
                margin=24,
                container_padding=24,
                content_max_width=720,
                min_height=500,
                visualization_height=480,
                controls_height=80,
                padding=(16, 16, 16, 16),
            ),
            "lg": LayoutConfig(
                columns=12,
                gap=32,
                margin=32,
                container_padding=24,
                content_max_width=960,
                min_height=600,
                visualization_height=540,
                controls_height=80,
                padding=(24, 24, 24, 24),
            ),
            "xl": LayoutConfig(
                columns=12,
                gap=32,
                margin=32,
                container_padding=24,
                content_max_width=1140,
                min_height=700,
                visualization_height=600,
                controls_height=80,
                padding=(32, 32, 32, 32),
            ),
        }

    async def setup(self) -> None:
        """Set up responsive layout system."""
        await super().setup()

        # Configure viewport
        await self.viewport.configure(
            width="device-width",
            initial_scale=1.0,
            maximum_scale=5.0,
            user_scalable=True,
            safe_area_inset=True,
        )

        await self.adjust_to_breakpoint(self.viewport._current_breakpoint)
        await self._initialize_optimizations()

    async def handle_resize(self, width: int) -> None:
        """Handle viewport resize.

        Args:
            width: New viewport width
        """
        await super().handle_resize(width)

        # Update grid configuration
        if width < 576:  # xs
            self.config.columns = self.config.xs_columns
            self.config.structure_width = width - 32
            self.config.chart_width = width - 32
        elif width < 768:  # sm
            self.config.columns = self.config.sm_columns
            self.config.structure_width = (width - 48) // 2
            self.config.chart_width = width - 32
        elif width < 992:  # md
            self.config.columns = self.config.md_columns
            self.config.structure_width = (width - 96) // 3
            self.config.chart_width = (width - 72) // 2
        elif width < 1200:  # lg
            self.config.columns = self.config.lg_columns
            self.config.structure_width = (width - 128) // 4
            self.config.chart_width = (width - 96) // 2
        else:  # xl
            self.config.columns = self.config.xl_columns
            self.config.structure_width = (width - 160) // 5
            self.config.chart_width = (width - 96) // 2

        # Update grid
        await self.grid.adjust(self._current_breakpoint)

    async def adjust_to_breakpoint(self, breakpoint: str) -> None:
        """Adjust layout system to breakpoint."""
        if breakpoint in self._breakpoint_layouts:
            self.config = self._breakpoint_layouts[breakpoint]
            await self._update_layout_optimizations()

    def get_layout_styles(self, layout_type: Optional[LayoutType] = None) -> Dict[str, str]:
        """Get layout styles with optimizations."""
        layout = layout_type or self.config.layout_type
        grid_type = self._get_grid_type(layout)

        styles = {
            "display": "grid",
            "width": "100%",
            "max-width": f"{self.config.content_max_width}px" if self.config.content_max_width else "100%",
            "margin": " ".join(str(m) + "px" for m in self.config.margin),
            "padding": " ".join(str(p) + "px" for p in self.config.padding),
            "gap": f"{self.config.gap}px",
            "min-height": f"{self.config.min_height}px" if self.config.min_height else "auto",
            "grid-template-columns": f"repeat({self.config.columns}, 1fr)",
        }

        # Add layout-specific styles
        if layout == LayoutType.COMPOUND:
            styles.update(self._get_compound_layout_styles())
        elif layout == LayoutType.SCIENTIFIC:
            styles.update(self._get_scientific_layout_styles())
        elif layout == LayoutType.ANALYSIS:
            styles.update(self._get_analysis_layout_styles())
        elif layout == LayoutType.DASHBOARD:
            styles.update(self._get_dashboard_layout_styles())
        elif layout == LayoutType.QUANTUM:
            styles.update(self._get_quantum_layout_styles())
        elif layout == LayoutType.MOLECULAR:
            styles.update(self._get_molecular_layout_styles())
        elif layout == LayoutType.COMPARISON:
            styles.update(self._get_comparison_layout_styles())

        return styles

    def get_area_styles(self, area_name: str) -> Dict[str, str]:
        """Get styles for layout area."""
        area = self._active_areas.get(area_name)
        if not area:
            return {}

        styles = {
            "grid-area": area.grid_area,
            "overflow": "auto" if area.scrollable else "hidden",
        }

        if area.min_height:
            styles["min-height"] = f"{area.min_height}px"
        if area.max_height:
            styles["max-height"] = f"{area.max_height}px"

        # Add area-specific styles
        if area.quantum_ready:
            styles.update(self._get_quantum_area_styles())
        elif area.molecular_view:
            styles.update(self._get_molecular_area_styles())
        elif area.scientific_data:
            styles.update(self._get_scientific_data_styles())
        elif "structure" in area_name:
            styles.update(self._get_structure_area_styles())
        elif "data" in area_name:
            styles.update(self._get_data_area_styles())
        elif "visualization" in area_name:
            styles.update(self._get_visualization_area_styles())
        elif "controls" in area_name:
            styles.update(self._get_controls_area_styles())

        # Add interaction styles
        if area.resizable:
            styles.update(self._get_resizable_styles())
        if area.draggable:
            styles.update(self._get_draggable_styles())

        return styles

    def get_component_config(self, component_type: str) -> Dict[str, Any]:
        """Get component configuration.

        Args:
            component_type: Type of component

        Returns:
            Component configuration dictionary
        """
        base_config = {
            "breakpoint": self._current_breakpoint,
            "container_class": self.get_container_class(),
            "grid_style": self.get_grid_style(),
            "touch_style": self.get_touch_style(),
            "data_style": self.get_data_style(),
        }

        if component_type == "list":
            base_config.update(
                {
                    "structure_style": self.get_structure_style(),
                    "chart_style": self.get_chart_style(),
                }
            )

        elif component_type == "detail":
            base_config.update(
                {
                    "structure_style": {
                        "width": f"{self.config.structure_width * 2}px",
                        "height": f"{self.config.structure_height * 2}px",
                    },
                    "chart_style": {
                        "width": f"{self.config.chart_width * 1.5}px",
                        "height": f"{self.config.chart_height * 1.5}px",
                    },
                }
            )

        elif component_type == "search":
            base_config.update(
                {
                    "input_style": self.get_touch_style(),
                }
            )

        return base_config

    def get_container_class(self) -> str:
        """Get container class with current breakpoint."""
        return f"container {self._current_breakpoint}"

    def get_touch_style(self) -> Dict[str, str]:
        """Get touch-friendly styles."""
        return {
            "min-height": f"{self.config.min_touch_target}px",
            "min-width": f"{self.config.min_touch_target}px",
            "padding": f"{self.config.touch_spacing}px",
        }

    def get_structure_style(self) -> Dict[str, str]:
        """Get structure viewer styles."""
        return {
            "width": f"{self.config.structure_width}px",
            "height": f"{self.config.structure_height}px",
        }

    def get_chart_style(self) -> Dict[str, str]:
        """Get chart styles."""
        return {
            "width": f"{self.config.chart_width}px",
            "height": f"{self.config.chart_height}px",
        }

    def get_data_style(self) -> Dict[str, str]:
        """Get data display styles."""
        return {
            "font-family": self.config.data_font,
            "font-variant-numeric": "tabular-nums",
        }

    def _get_quantum_area_styles(self) -> Dict[str, str]:
        """Get styles for quantum visualization areas."""
        return {
            "background": "var(--quantum-bg)",
            "border-radius": "12px",
            "box-shadow": "var(--quantum-glow)",
            "backdrop-filter": "blur(8px)",
            "transition": "all 0.3s cubic-bezier(0.4, 0, 0.2, 1)",
        }

    def _get_molecular_area_styles(self) -> Dict[str, str]:
        """Get styles for molecular visualization areas."""
        detail_map = {
            "low": {"shadow": "sm", "blur": "4px"},
            "medium": {"shadow": "md", "blur": "8px"},
            "high": {"shadow": "lg", "blur": "12px"},
        }
        detail = detail_map[self.config.molecular_detail]

        return {
            "background": "var(--molecular-bg)",
            "border-radius": "12px",
            f"box-shadow": f"var(--molecular-shadow-{detail['shadow']})",
            "backdrop-filter": f"blur({detail['blur']})",
            "transition": "transform 0.2s ease-in-out",
        }

    def _get_comparison_layout_styles(self) -> Dict[str, str]:
        """Get styles for comparison layouts."""
        if self.config.comparison_mode == "overlay":
            return {
                "position": "relative",
                "isolation": "isolate",
                "mix-blend-mode": "difference",
            }
        elif self.config.comparison_mode == "split":
            return {
                "position": "relative",
                "clip-path": "polygon(0 0, 50% 0, 50% 100%, 0 100%)",
                "transition": "clip-path 0.3s ease-in-out",
            }
        return {}  # side-by-side is default

    def _get_resizable_styles(self) -> Dict[str, str]:
        """Get styles for resizable areas."""
        return {
            "resize": "both",
            "overflow": "auto",
            "min-width": "100px",
            "min-height": "100px",
            "position": "relative",
        }

    def _get_draggable_styles(self) -> Dict[str, str]:
        """Get styles for draggable areas."""
        return {
            "cursor": "move",
            "user-select": "none",
            "position": "relative",
            "z-index": "10",
            "touch-action": "none",
        }

    async def _initialize_optimizations(self) -> None:
        """Initialize layout optimizations."""
        # Performance optimizations
        if self.config.high_performance:
            self._active_optimizations.add("gpu_acceleration")
            self._active_optimizations.add("content_visibility")
            self._active_optimizations.add("will_change_transform")

        # Scientific optimizations
        if self.config.scientific_precision:
            self._active_optimizations.add("high_precision")
            self._active_optimizations.add("tabular_nums")
            self._active_optimizations.add("monospace_font")

        # Animation optimizations
        if self.config.enable_animations:
            self._active_optimizations.add("smooth_animations")
        if self.config.enable_transitions:
            self._active_optimizations.add("smooth_transitions")

        # Quantum optimizations
        if self.config.quantum_enabled:
            self._active_optimizations.add("quantum_effects")
            self._active_optimizations.add("quantum_transitions")

        # Update grid optimizations
        grid_type = self._get_grid_type(self.config.layout_type)
        self.grid.config.grid_type = grid_type

    def _get_grid_type(self, layout_type: LayoutType) -> GridType:
        """Get corresponding grid type for layout."""
        return {
            LayoutType.STANDARD: GridType.STANDARD,
            LayoutType.COMPOUND: GridType.COMPOUND,
            LayoutType.SCIENTIFIC: GridType.SCIENTIFIC,
            LayoutType.ANALYSIS: GridType.ANALYSIS,
            LayoutType.DASHBOARD: GridType.DASHBOARD,
            LayoutType.QUANTUM: GridType.SCIENTIFIC,
            LayoutType.MOLECULAR: GridType.COMPOUND,
            LayoutType.COMPARISON: GridType.ANALYSIS,
        }[layout_type]

    def register_area(self, area: LayoutArea) -> None:
        """Register a layout area."""
        self._active_areas[area.name] = area
        self._update_layout_optimizations()

    def remove_area(self, area_name: str) -> None:
        """Remove a layout area."""
        if area_name in self._active_areas:
            del self._active_areas[area_name]
            self._update_layout_optimizations()

    async def _update_layout_optimizations(self) -> None:
        """Update active layout optimizations."""
        await self._initialize_optimizations()

        # Add area-specific optimizations
        for area in self._active_areas.values():
            if area.quantum_ready:
                self._active_optimizations.add("quantum_ready")
            if area.molecular_view:
                self._active_optimizations.add("molecular_view")
            if area.scientific_data:
                self._active_optimizations.add("scientific_data")
