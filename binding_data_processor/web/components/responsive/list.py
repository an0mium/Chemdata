"""Enhanced responsive list component optimized for scientific data visualization and compound displays."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Set, Union, Callable

from .base import EnhancedResponsiveBase
from .layout import EnhancedResponsiveLayout, LayoutType, LayoutArea


class ListDisplayMode(Enum):
    """Display modes for compound lists."""

    GRID = "grid"  # Grid layout with cards
    TABLE = "table"  # Tabular data display
    COMPACT = "compact"  # Condensed list view
    DETAILED = "detailed"  # Expanded list with details
    MOLECULAR = "molecular"  # Molecular structure focused
    QUANTUM = "quantum"  # Quantum properties display
    COMPARISON = "comparison"  # Side-by-side comparison
    ANALYSIS = "analysis"  # Analysis results view
    SPECTRAL = "spectral"  # Spectroscopic data view
    CRYSTALLOGRAPHIC = "crystallographic"  # Crystal structure view
    NMR = "nmr"  # NMR spectroscopy view
    MASS_SPEC = "mass_spec"  # Mass spectrometry view
    CHROMATOGRAPHY = "chromatography"  # Chromatography data view
    PHARMACOLOGY = "pharmacology"  # Pharmacological data view
    TOXICOLOGY = "toxicology"  # Toxicological data view


@dataclass
class ScientificMetadata:
    """Scientific metadata for enhanced data handling."""

    data_quality: float = 1.0  # Quality score (0-1)
    confidence_level: float = 1.0  # Confidence score (0-1)
    validation_status: str = "verified"  # Data validation status
    last_verified: Optional[str] = None  # ISO timestamp
    measurement_method: Optional[str] = None
    instrumentation: Optional[str] = None
    resolution: Optional[float] = None
    precision: Optional[float] = None
    uncertainty: Optional[float] = None
    units: Optional[str] = None
    conditions: Optional[Dict[str, Any]] = None


@dataclass
class ScientificListItem:
    """Enhanced list item with scientific data capabilities."""

    id: str
    title: str
    subtitle: Optional[str] = None
    description: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None
    scientific_metadata: Optional[ScientificMetadata] = None
    image_url: Optional[str] = None
    structure_data: Optional[Dict[str, Any]] = None
    quantum_data: Optional[Dict[str, Any]] = None
    spectral_data: Optional[Dict[str, Any]] = None
    crystal_data: Optional[Dict[str, Any]] = None
    nmr_data: Optional[Dict[str, Any]] = None
    mass_spec_data: Optional[Dict[str, Any]] = None
    chromatography_data: Optional[Dict[str, Any]] = None
    pharmacology_data: Optional[Dict[str, Any]] = None
    toxicology_data: Optional[Dict[str, Any]] = None
    analysis_results: Optional[Dict[str, Any]] = None
    molecular_properties: Optional[Dict[str, Any]] = None
    binding_data: Optional[Dict[str, Any]] = None
    safety_data: Optional[Dict[str, Any]] = None
    citations: Optional[List[Dict[str, str]]] = None
    tags: Optional[List[str]] = None
    priority: int = 0
    visibility: float = 1.0


@dataclass
class ListConfig:
    """Enhanced list configuration with scientific display options."""

    # Display options
    display_mode: ListDisplayMode = ListDisplayMode.GRID
    page_size: int = 10
    show_structures: bool = True
    show_predictions: bool = True
    show_analysis: bool = True
    show_visualizations: bool = True
    enable_filtering: bool = True
    enable_sorting: bool = True
    enable_selection: bool = True
    enable_batch_ops: bool = True
    enable_export: bool = True
    enable_accessibility: bool = True
    enable_keyboard_nav: bool = True
    enable_touch_gestures: bool = True
    enable_data_export: bool = True
    enable_citations: bool = True
    enable_3d_visualization: bool = True
    enable_data_validation: bool = True
    enable_quality_indicators: bool = True
    enable_uncertainty_visualization: bool = True

    # Layout options
    layout_type: LayoutType = LayoutType.COMPOUND
    grid_columns: int = 12
    structure_columns: int = 4
    data_columns: int = 8
    gap: int = 16
    margin: int = 16
    min_card_width: int = 280
    max_card_width: int = 400
    card_aspect_ratio: float = 1.4

    # Component sizes
    structure_width: int = 200
    structure_height: int = 200
    chart_width: int = 400
    chart_height: int = 300

    # Scientific display
    scientific_precision: bool = True
    molecular_detail: str = "high"  # low, medium, high
    quantum_enabled: bool = False
    visualization_quality: str = "high"  # low, medium, high
    data_precision_level: str = "high"  # low, medium, high

    # Performance
    high_performance: bool = True
    enable_virtualization: bool = True
    batch_size: int = 50
    prefetch_count: int = 10
    enable_animations: bool = True
    enable_transitions: bool = True

    # Export options
    export_formats: List[str] = field(default_factory=lambda: ["tsv", "csv", "json", "sdf"])
    default_export_columns: List[str] = field(
        default_factory=lambda: [
            "name",
            "smiles",
            "cas",
            "binding_affinities",
            "toxicity_score",
            "abuse_potential",
            "bbb_permeability",
        ]
    )

    # Interaction
    interaction_mode: str = "standard"  # standard, expert, touch
    pagination_style: str = "numbered"  # numbered, infinite-scroll, load-more
    comparison_mode: str = "side-by-side"  # side-by-side, overlay, split


class EnhancedResponsiveList(EnhancedResponsiveBase):
    """Enhanced responsive list component with scientific optimization."""

    def __init__(
        self,
        config: Optional[ListConfig] = None,
        on_select: Optional[Callable[[str], None]] = None,
        on_export: Optional[Callable[[List[str], str], None]] = None,
    ):
        """Initialize list component.

        Args:
            config: Optional list configuration
            on_select: Optional callback when item is selected
            on_export: Optional callback when items are exported
        """
        super().__init__()
        self.config = config or ListConfig()
        self.layout = EnhancedResponsiveLayout()
        self.on_select = on_select
        self.on_export = on_export
        self._active_optimizations: Set[str] = set()
        self._selected_items: List[str] = []
        self._current_page = 1
        self._total_pages = 1
        self._items: List[ScientificListItem] = []
        self._filtered_items: List[ScientificListItem] = []
        self._sort_by: Optional[str] = None
        self._sort_ascending: bool = True
        self._filters: Dict[str, Any] = {}

    async def setup(self) -> None:
        """Set up responsive list component."""
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

        if self.config.show_visualizations:
            self.layout.register_area(
                LayoutArea(
                    name="visualization_panel",
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

    async def update_items(self, items: List[ScientificListItem]) -> None:
        """Update list items.

        Args:
            items: List of items to display
        """
        self._items = items
        await self._process_scientific_data()
        self._apply_filters()
        self._apply_sort()
        self._update_pagination()

    def set_filters(self, filters: Dict[str, Any]) -> None:
        """Set active filters.

        Args:
            filters: Dictionary of filter settings
        """
        self._filters = filters
        self._apply_filters()
        self._update_pagination()
        self._current_page = 1

    def set_sort(self, field: str, ascending: bool = True) -> None:
        """Set sort settings.

        Args:
            field: Field to sort by
            ascending: Sort direction
        """
        self._sort_by = field
        self._sort_ascending = ascending
        self._apply_sort()

    def get_page_items(self, page: int) -> List[ScientificListItem]:
        """Get items for specified page.

        Args:
            page: Page number (1-based)

        Returns:
            List of items for page
        """
        if not self._filtered_items:
            return []

        start = (page - 1) * self.config.page_size
        end = start + self.config.page_size
        return self._filtered_items[start:end]

    def get_selected_items(self) -> List[ScientificListItem]:
        """Get selected items.

        Returns:
            List of selected items
        """
        return [item for item in self._filtered_items if item.id in self._selected_items]

    def select_item(self, item_id: str) -> None:
        """Select an item.

        Args:
            item_id: ID of item to select
        """
        if item_id not in self._selected_items:
            self._selected_items.append(item_id)
            if self.on_select:
                self.on_select(item_id)

    def deselect_item(self, item_id: str) -> None:
        """Deselect an item.

        Args:
            item_id: ID of item to deselect
        """
        if item_id in self._selected_items:
            self._selected_items.remove(item_id)
            if self.on_select:
                self.on_select(item_id)

    def clear_selection(self) -> None:
        """Clear all selected items."""
        self._selected_items = []

    def export_items(
        self,
        format: str = "tsv",
        columns: Optional[List[str]] = None,
    ) -> None:
        """Export selected items.

        Args:
            format: Export format
            columns: Optional list of columns to export
        """
        if format not in self.config.export_formats:
            raise ValueError(f"Invalid export format: {format}")

        items = self.get_selected_items() or self._filtered_items
        if self.on_export:
            self.on_export([item.id for item in items], format)

    async def _process_scientific_data(self) -> None:
        """Process and optimize scientific data for display."""
        for item in self._items:
            # Process each data type
            if item.structure_data:
                await self._optimize_structure_data(item)
            if item.quantum_data:
                await self._optimize_quantum_data(item)
            if item.spectral_data:
                await self._optimize_spectral_data(item)
            if item.crystal_data:
                await self._optimize_crystal_data(item)
            if item.nmr_data:
                await self._optimize_nmr_data(item)
            if item.mass_spec_data:
                await self._optimize_mass_spec_data(item)
            if item.chromatography_data:
                await self._optimize_chromatography_data(item)

            # Validate data quality
            if self.config.enable_data_validation and item.scientific_metadata:
                await self._validate_scientific_data(item)

    async def _validate_scientific_data(self, item: ScientificListItem) -> None:
        """Validate scientific data quality and completeness."""
        if not item.scientific_metadata:
            return

        # Add validation flags
        self._active_optimizations.add("data_validated")

        # Update validation status based on checks
        if self._check_data_completeness(item) and self._check_data_consistency(item):
            item.scientific_metadata.validation_status = "verified"
        else:
            item.scientific_metadata.validation_status = "needs_review"

    def _check_data_completeness(self, item: ScientificListItem) -> bool:
        """Check if all required data fields are present."""
        required_fields = self._get_required_fields(item)
        return all(self._has_field(item, field) for field in required_fields)

    def _check_data_consistency(self, item: ScientificListItem) -> bool:
        """Check internal data consistency."""
        # Add specific consistency checks based on data type
        if item.structure_data and item.molecular_properties:
            return self._check_structure_property_consistency(item)
        return True

    def _get_required_fields(self, item: ScientificListItem) -> Set[str]:
        """Get required fields based on display mode."""
        required = {"id", "title"}

        if self.config.display_mode == ListDisplayMode.MOLECULAR:
            required.add("structure_data")
        elif self.config.display_mode == ListDisplayMode.SPECTRAL:
            required.add("spectral_data")

        return required

    def _has_field(self, item: ScientificListItem, field: str) -> bool:
        """Check if item has required field with valid data."""
        return hasattr(item, field) and getattr(item, field) is not None

    async def _initialize_optimizations(self) -> None:
        """Initialize list optimizations."""
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
        if self.config.enable_3d_visualization:
            self._active_optimizations.add("3d_ready")
            if self.config.visualization_quality == "high":
                self._active_optimizations.add("high_quality_3d")

        # Data precision optimizations
        if self.config.data_precision_level == "high":
            self._active_optimizations.add("high_precision_data")
            self._active_optimizations.add("uncertainty_visualization")

        # Interaction optimizations
        if self.config.interaction_mode == "expert":
            self._active_optimizations.add("expert_controls")
        elif self.config.interaction_mode == "touch":
            self._active_optimizations.add("touch_optimized")

        # Accessibility optimizations
        if self.config.enable_accessibility:
            self._active_optimizations.add("screen_reader_support")
            self._active_optimizations.add("keyboard_navigation")
            self._active_optimizations.add("high_contrast")

    def _apply_filters(self) -> None:
        """Apply active filters to items."""
        filtered = self._items

        for field, value in self._filters.items():
            if not value:
                continue

            if isinstance(value, dict):
                # Range filter
                min_val = value.get("min")
                max_val = value.get("max")
                filtered = [
                    item
                    for item in filtered
                    if (hasattr(item, field) and (min_val is None or getattr(item, field) >= min_val) and (max_val is None or getattr(item, field) <= max_val))
                ]
            else:
                # Exact match filter
                filtered = [item for item in filtered if hasattr(item, field) and getattr(item, field) == value]

        self._filtered_items = filtered

    def _apply_sort(self) -> None:
        """Apply sort to filtered items."""
        if not self._sort_by:
            return

        def get_sort_value(item: ScientificListItem) -> Any:
            if hasattr(item, self._sort_by):
                return getattr(item, self._sort_by) or ""
            return ""

        self._filtered_items.sort(
            key=get_sort_value,
            reverse=not self._sort_ascending,
        )

    def _update_pagination(self) -> None:
        """Update pagination state."""
        if not self._filtered_items:
            self._total_pages = 1
            return

        self._total_pages = (len(self._filtered_items) + self.config.page_size - 1) // self.config.page_size

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
