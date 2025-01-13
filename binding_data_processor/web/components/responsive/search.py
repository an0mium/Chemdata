"""Enhanced responsive search component optimized for scientific data."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Any, Set, Callable

from .base import EnhancedResponsiveBase
from .layout import EnhancedResponsiveLayout, LayoutType, LayoutArea
from .list import ScientificListItem


class SearchMode(Enum):
    """Search modes for compound data."""

    TEXT = "text"  # Basic text search
    STRUCTURE = "structure"  # Structure/substructure search
    SIMILARITY = "similarity"  # Similarity search
    PROPERTY = "property"  # Property-based search
    ADVANCED = "advanced"  # Advanced query builder
    QUANTUM = "quantum"  # Quantum property search
    BINDING = "binding"  # Binding data search
    SAFETY = "safety"  # Safety data search
    LITERATURE = "literature"  # Literature search


@dataclass
class SearchConfig:
    """Configuration for search component."""

    # Search options
    search_mode: SearchMode = SearchMode.TEXT
    enable_structure_search: bool = True
    enable_similarity_search: bool = True
    enable_property_search: bool = True
    enable_advanced_search: bool = True
    enable_quantum_search: bool = True
    enable_binding_search: bool = True
    enable_safety_search: bool = True
    enable_literature_search: bool = True

    # Layout options
    layout_type: LayoutType = LayoutType.STANDARD
    grid_columns: int = 12
    search_columns: int = 4
    results_columns: int = 8
    gap: int = 16
    margin: int = 16

    # Component sizes
    search_width: int = 400
    structure_editor_width: int = 600
    structure_editor_height: int = 400
    results_width: int = 800

    # Scientific display
    scientific_precision: bool = True
    molecular_detail: str = "high"  # low, medium, high
    quantum_enabled: bool = False
    visualization_quality: str = "high"  # low, medium, high
    data_precision_level: str = "high"  # low, medium, high

    # Interaction options
    enable_structure_drawing: bool = True
    enable_property_ranges: bool = True
    enable_query_builder: bool = True
    enable_history: bool = True
    enable_suggestions: bool = True
    enable_autocomplete: bool = True
    enable_accessibility: bool = True
    enable_keyboard_nav: bool = True
    enable_touch_gestures: bool = True

    # Performance
    high_performance: bool = True
    enable_search_debounce: bool = True
    debounce_delay: int = 300
    enable_results_virtualization: bool = True
    batch_size: int = 50
    prefetch_count: int = 10

    # Search fields
    default_fields: List[str] = field(
        default_factory=lambda: [
            "name",
            "smiles",
            "inchi",
            "cas",
            "description",
            "tags",
        ]
    )
    property_fields: List[str] = field(
        default_factory=lambda: [
            "molecular_weight",
            "logP",
            "PSA",
            "HBA",
            "HBD",
        ]
    )
    binding_fields: List[str] = field(
        default_factory=lambda: [
            "5HT2A",
            "D2",
            "binding_affinities",
        ]
    )
    safety_fields: List[str] = field(
        default_factory=lambda: [
            "ld50",
            "toxicity_score",
            "abuse_potential",
        ]
    )


class EnhancedResponsiveSearch(EnhancedResponsiveBase):
    """Enhanced responsive search component with scientific capabilities."""

    def __init__(
        self,
        config: Optional[SearchConfig] = None,
        on_search: Optional[Callable[[Dict[str, Any]], None]] = None,
        on_select: Optional[Callable[[str], None]] = None,
    ):
        """Initialize search component.

        Args:
            config: Optional search configuration
            on_search: Optional callback when search is performed
            on_select: Optional callback when result is selected
        """
        super().__init__()
        self.config = config or SearchConfig()
        self.layout = EnhancedResponsiveLayout()
        self.on_search = on_search
        self.on_select = on_select
        self._active_optimizations: Set[str] = set()
        self._current_query: Dict[str, Any] = {}
        self._search_history: List[Dict[str, Any]] = []
        self._suggestions: List[str] = []
        self._results: List[ScientificListItem] = []

    async def setup(self) -> None:
        """Set up responsive search component."""
        await super().setup()
        await self.layout.setup()
        await self._initialize_optimizations()

        # Register layout areas
        self.layout.register_area(
            LayoutArea(
                name="search_panel",
                grid_area=f"1 / 1 / 2 / {self.config.search_columns + 1}",
                scrollable=True,
            )
        )

        if self.config.enable_structure_drawing:
            self.layout.register_area(
                LayoutArea(
                    name="structure_editor",
                    grid_area=f"2 / 1 / 3 / {self.config.search_columns + 1}",
                    molecular_view=True,
                )
            )

        self.layout.register_area(
            LayoutArea(
                name="results_panel",
                grid_area=f"1 / {self.config.search_columns + 1} / 3 / {self.config.grid_columns + 1}",
                scientific_data=True,
                scrollable=True,
            )
        )

        # Configure layout
        self.layout.config.layout_type = self.config.layout_type
        self.layout.config.scientific_precision = self.config.scientific_precision
        self.layout.config.molecular_detail = self.config.molecular_detail
        self.layout.config.quantum_enabled = self.config.quantum_enabled
        self.layout.config.high_performance = self.config.high_performance

    async def search(self, query: Dict[str, Any]) -> None:
        """Perform search with given query.

        Args:
            query: Search query parameters
        """
        self._current_query = query
        self._search_history.append(query)

        if self.on_search:
            self.on_search(query)

    async def select_result(self, result_id: str) -> None:
        """Select search result.

        Args:
            result_id: ID of result to select
        """
        if self.on_select:
            self.on_select(result_id)

    async def update_results(self, results: List[ScientificListItem]) -> None:
        """Update search results.

        Args:
            results: List of search results
        """
        self._results = results
        await self._process_scientific_data()

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
            self.config.search_columns = self.config.grid_columns
            self.config.search_width = width - 32
            self.config.results_width = width - 32
        elif width < 768:  # sm
            self.config.search_columns = self.config.grid_columns
            self.config.search_width = width - 48
            self.config.results_width = width - 48
        elif width < 992:  # md
            self.config.search_columns = 4
            self.config.search_width = (width - 72) // 3
            self.config.results_width = (width - 72) * 2 // 3
        elif width < 1200:  # lg
            self.config.search_columns = 4
            self.config.search_width = (width - 96) // 3
            self.config.results_width = (width - 96) * 2 // 3
        else:  # xl
            self.config.search_columns = 4
            self.config.search_width = (width - 120) // 3
            self.config.results_width = (width - 120) * 2 // 3

    async def _initialize_optimizations(self) -> None:
        """Initialize search optimizations."""
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

        # Search optimizations
        if self.config.enable_search_debounce:
            self._active_optimizations.add("search_debounce")
        if self.config.enable_results_virtualization:
            self._active_optimizations.add("results_virtualization")

        # Visualization optimizations
        if self.config.enable_structure_drawing:
            self._active_optimizations.add("structure_editor")
            if self.config.visualization_quality == "high":
                self._active_optimizations.add("high_quality_structure")

        # Interaction optimizations
        if self.config.enable_suggestions:
            self._active_optimizations.add("search_suggestions")
        if self.config.enable_autocomplete:
            self._active_optimizations.add("search_autocomplete")
        if self.config.enable_touch_gestures:
            self._active_optimizations.add("touch_optimized")

        # Accessibility optimizations
        if self.config.enable_accessibility:
            self._active_optimizations.add("screen_reader_support")
            self._active_optimizations.add("keyboard_navigation")
            self._active_optimizations.add("high_contrast")

    async def _process_scientific_data(self) -> None:
        """Process and optimize scientific data in results."""
        for result in self._results:
            # Process structure data
            if result.structure_data and self.config.enable_structure_drawing:
                await self._optimize_structure_data(result)

            # Process quantum data
            if result.quantum_data and self.config.quantum_enabled:
                await self._optimize_quantum_data(result)

            # Process spectral data
            if result.spectral_data:
                await self._optimize_spectral_data(result)

            # Process binding data
            if result.binding_data and self.config.enable_binding_search:
                await self._optimize_binding_data(result)

            # Process safety data
            if result.safety_data and self.config.enable_safety_search:
                await self._optimize_safety_data(result)

    async def _optimize_structure_data(self, result: ScientificListItem) -> None:
        """Optimize structure data for visualization.

        Args:
            result: Result to optimize structure data for
        """
        if not result.structure_data:
            return

        # Add structure optimization flags
        self._active_optimizations.add("structure_ready")
        if self.config.molecular_detail == "high":
            self._active_optimizations.add("high_quality_structure")

    async def _optimize_quantum_data(self, result: ScientificListItem) -> None:
        """Optimize quantum data for visualization.

        Args:
            result: Result to optimize quantum data for
        """
        if not result.quantum_data:
            return

        # Add quantum optimization flags
        self._active_optimizations.add("quantum_ready")
        if self.config.visualization_quality == "high":
            self._active_optimizations.add("high_quality_quantum")

    async def _optimize_spectral_data(self, result: ScientificListItem) -> None:
        """Optimize spectral data for visualization.

        Args:
            result: Result to optimize spectral data for
        """
        if not result.spectral_data:
            return

        # Add spectral optimization flags
        self._active_optimizations.add("spectral_ready")
        if self.config.visualization_quality == "high":
            self._active_optimizations.add("high_quality_spectral")

    async def _optimize_binding_data(self, result: ScientificListItem) -> None:
        """Optimize binding data for visualization.

        Args:
            result: Result to optimize binding data for
        """
        if not result.binding_data:
            return

        # Add binding optimization flags
        self._active_optimizations.add("binding_ready")
        if self.config.data_precision_level == "high":
            self._active_optimizations.add("high_precision_binding")

    async def _optimize_safety_data(self, result: ScientificListItem) -> None:
        """Optimize safety data for visualization.

        Args:
            result: Result to optimize safety data for
        """
        if not result.safety_data:
            return

        # Add safety optimization flags
        self._active_optimizations.add("safety_ready")
        if self.config.data_precision_level == "high":
            self._active_optimizations.add("high_precision_safety")
