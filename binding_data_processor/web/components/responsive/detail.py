"""Enhanced responsive detail component optimized for scientific data visualization."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Any, Set

from .base import EnhancedResponsiveBase
from .layout import EnhancedResponsiveLayout, LayoutType, LayoutArea
from .list import ScientificListItem, ScientificMetadata


class DetailDisplayMode(Enum):
    """Display modes for compound details."""

    OVERVIEW = "overview"  # Basic compound information
    STRUCTURE = "structure"  # 3D structure visualization
    QUANTUM = "quantum"  # Quantum properties visualization
    SPECTRAL = "spectral"  # Spectroscopic data display
    CRYSTAL = "crystal"  # Crystallographic data display
    ANALYSIS = "analysis"  # Analysis results view
    BINDING = "binding"  # Binding data visualization
    SAFETY = "safety"  # Safety data display
    LITERATURE = "literature"  # Literature references
    COMPARISON = "comparison"  # Side-by-side comparison


@dataclass
class DetailConfig:
    """Configuration for detail component."""

    # Display options
    display_mode: DetailDisplayMode = DetailDisplayMode.OVERVIEW
    show_3d_structure: bool = True
    show_quantum_viz: bool = True
    show_spectral_data: bool = True
    show_crystal_data: bool = True
    show_analysis: bool = True
    show_binding: bool = True
    show_safety: bool = True
    show_literature: bool = True
    enable_comparison: bool = True

    # Layout options
    layout_type: LayoutType = LayoutType.MOLECULAR
    grid_columns: int = 12
    structure_columns: int = 6
    data_columns: int = 6
    gap: int = 24
    margin: int = 24

    # Component sizes
    structure_width: int = 600
    structure_height: int = 600
    chart_width: int = 800
    chart_height: int = 400
    sidebar_width: int = 320

    # Scientific display
    scientific_precision: bool = True
    molecular_detail: str = "high"  # low, medium, high
    quantum_enabled: bool = False
    visualization_quality: str = "high"  # low, medium, high
    data_precision_level: str = "high"  # low, medium, high

    # Interaction options
    enable_structure_interaction: bool = True
    enable_data_export: bool = True
    enable_sharing: bool = True
    enable_annotations: bool = True
    enable_history: bool = True
    enable_accessibility: bool = True
    enable_keyboard_nav: bool = True
    enable_touch_gestures: bool = True

    # Performance
    high_performance: bool = True
    enable_animations: bool = True
    enable_transitions: bool = True

    # Export options
    export_formats: List[str] = field(default_factory=lambda: ["mol", "sdf", "pdb", "json", "csv"])


class EnhancedResponsiveDetail(EnhancedResponsiveBase):
    """Enhanced responsive detail component with scientific visualization."""

    def __init__(self, config: Optional[DetailConfig] = None):
        """Initialize detail component.

        Args:
            config: Optional detail configuration
        """
        super().__init__()
        self.config = config or DetailConfig()
        self.layout = EnhancedResponsiveLayout()
        self._active_optimizations: Set[str] = set()
        self._current_item: Optional[ScientificListItem] = None
        self._comparison_item: Optional[ScientificListItem] = None

    async def setup(self) -> None:
        """Set up responsive detail component."""
        await super().setup()
        await self.layout.setup()
        await self._initialize_optimizations()

        # Register layout areas
        self.layout.register_area(
            LayoutArea(
                name="structure_viewer",
                grid_area=f"1 / 1 / 2 / {self.config.structure_columns + 1}",
                molecular_view=True,
                quantum_ready=self.config.quantum_enabled,
            )
        )

        self.layout.register_area(
            LayoutArea(
                name="data_panel",
                grid_area=f"1 / {self.config.structure_columns + 1} / 2 / {self.config.grid_columns + 1}",
                scientific_data=True,
                scrollable=True,
            )
        )

        if self.config.show_spectral_data:
            self.layout.register_area(
                LayoutArea(
                    name="spectral_panel",
                    grid_area=f"2 / 1 / 3 / {self.config.grid_columns + 1}",
                    scientific_data=True,
                )
            )

        # Configure layout
        self.layout.config.layout_type = self.config.layout_type
        self.layout.config.scientific_precision = self.config.scientific_precision
        self.layout.config.molecular_detail = self.config.molecular_detail
        self.layout.config.quantum_enabled = self.config.quantum_enabled
        self.layout.config.high_performance = self.config.high_performance

    async def set_item(self, item: ScientificListItem) -> None:
        """Set current item for display.

        Args:
            item: Item to display
        """
        self._current_item = item
        await self._process_scientific_data()
        await self._update_display()

    async def set_comparison_item(self, item: Optional[ScientificListItem]) -> None:
        """Set comparison item.

        Args:
            item: Item to compare with current item
        """
        self._comparison_item = item
        if item:
            await self._process_scientific_data(item)
        await self._update_display()

    def get_layout_styles(self) -> Dict[str, str]:
        """Get layout styles.

        Returns:
            Dictionary of CSS properties
        """
        return self.layout.get_layout_styles()

    def get_area_styles(self, area_name: str) -> Dict[str, str]:
        """Get styles for layout area.

        Args:
            area_name: Name of area

        Returns:
            Dictionary of CSS properties
        """
        return self.layout.get_area_styles(area_name)

    async def handle_resize(self, width: int) -> None:
        """Handle viewport resize.

        Args:
            width: New viewport width
        """
        await super().handle_resize(width)
        await self.layout.handle_resize(width)

        # Update component sizes
        if width < 576:  # xs
            self.config.structure_columns = self.config.grid_columns
            self.config.structure_width = width - 32
            self.config.chart_width = width - 32
        elif width < 768:  # sm
            self.config.structure_columns = self.config.grid_columns
            self.config.structure_width = width - 48
            self.config.chart_width = width - 48
        elif width < 992:  # md
            self.config.structure_columns = 6
            self.config.structure_width = (width - 72) // 2
            self.config.chart_width = (width - 72) // 2
        elif width < 1200:  # lg
            self.config.structure_columns = 6
            self.config.structure_width = (width - 96) // 2
            self.config.chart_width = (width - 96) // 2
        else:  # xl
            self.config.structure_columns = 6
            self.config.structure_width = (width - 120) // 2
            self.config.chart_width = (width - 120) // 2

        await self._update_display()

    async def _initialize_optimizations(self) -> None:
        """Initialize detail optimizations."""
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

        # Visualization optimizations
        if self.config.show_3d_structure:
            self._active_optimizations.add("3d_ready")
            if self.config.visualization_quality == "high":
                self._active_optimizations.add("high_quality_3d")

        # Data precision optimizations
        if self.config.data_precision_level == "high":
            self._active_optimizations.add("high_precision_data")
            self._active_optimizations.add("uncertainty_visualization")

        # Interaction optimizations
        if self.config.enable_structure_interaction:
            self._active_optimizations.add("structure_interaction")
        if self.config.enable_touch_gestures:
            self._active_optimizations.add("touch_optimized")

        # Accessibility optimizations
        if self.config.enable_accessibility:
            self._active_optimizations.add("screen_reader_support")
            self._active_optimizations.add("keyboard_navigation")
            self._active_optimizations.add("high_contrast")

    async def _process_scientific_data(self, item: Optional[ScientificListItem] = None) -> None:
        """Process and optimize scientific data for display.

        Args:
            item: Optional item to process, defaults to current item
        """
        target = item or self._current_item
        if not target:
            return

        # Process structure data
        if target.structure_data and self.config.show_3d_structure:
            await self._optimize_structure_data(target)

        # Process quantum data
        if target.quantum_data and self.config.quantum_enabled:
            await self._optimize_quantum_data(target)

        # Process spectral data
        if target.spectral_data and self.config.show_spectral_data:
            await self._optimize_spectral_data(target)

        # Process crystal data
        if target.crystal_data and self.config.show_crystal_data:
            await self._optimize_crystal_data(target)

    async def _update_display(self) -> None:
        """Update display based on current state."""
        if not self._current_item:
            return

        # Update layout areas based on display mode
        if self.config.display_mode == DetailDisplayMode.STRUCTURE:
            self.layout.config.layout_type = LayoutType.MOLECULAR
        elif self.config.display_mode == DetailDisplayMode.QUANTUM:
            self.layout.config.layout_type = LayoutType.QUANTUM
        elif self.config.display_mode == DetailDisplayMode.SPECTRAL:
            self.layout.config.layout_type = LayoutType.SCIENTIFIC
        elif self.config.display_mode == DetailDisplayMode.CRYSTAL:
            self.layout.config.layout_type = LayoutType.SCIENTIFIC
        elif self.config.display_mode == DetailDisplayMode.ANALYSIS:
            self.layout.config.layout_type = LayoutType.ANALYSIS
        elif self.config.display_mode == DetailDisplayMode.COMPARISON:
            self.layout.config.layout_type = LayoutType.COMPARISON

        # Update optimizations based on data
        if self._current_item.quantum_data:
            self._active_optimizations.add("quantum_ready")
        if self._current_item.spectral_data:
            self._active_optimizations.add("spectral_ready")
        if self._current_item.crystal_data:
            self._active_optimizations.add("crystal_ready")

    async def _optimize_structure_data(self, item: ScientificListItem) -> None:
        """Optimize structure data for visualization.

        Args:
            item: Item to optimize structure data for
        """
        if not item.structure_data:
            return

        # Add structure optimization flags
        self._active_optimizations.add("structure_ready")
        if self.config.molecular_detail == "high":
            self._active_optimizations.add("high_quality_structure")

    async def _optimize_quantum_data(self, item: ScientificListItem) -> None:
        """Optimize quantum data for visualization.

        Args:
            item: Item to optimize quantum data for
        """
        if not item.quantum_data:
            return

        # Add quantum optimization flags
        self._active_optimizations.add("quantum_ready")
        if self.config.visualization_quality == "high":
            self._active_optimizations.add("high_quality_quantum")

    async def _optimize_spectral_data(self, item: ScientificListItem) -> None:
        """Optimize spectral data for visualization.

        Args:
            item: Item to optimize spectral data for
        """
        if not item.spectral_data:
            return

        # Add spectral optimization flags
        self._active_optimizations.add("spectral_ready")
        if self.config.visualization_quality == "high":
            self._active_optimizations.add("high_quality_spectral")

    async def _optimize_crystal_data(self, item: ScientificListItem) -> None:
        """Optimize crystal data for visualization.

        Args:
            item: Item to optimize crystal data for
        """
        if not item.crystal_data:
            return

        # Add crystal optimization flags
        self._active_optimizations.add("crystal_ready")
        if self.config.visualization_quality == "high":
            self._active_optimizations.add("high_quality_crystal")
