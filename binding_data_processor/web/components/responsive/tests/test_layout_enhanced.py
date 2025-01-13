"""Tests for enhanced layout system."""

import pytest
from ..layout import EnhancedLayout, LayoutConfig
from ..breakpoints import BreakpointManager
from ..base import Breakpoint


@pytest.fixture
def breakpoint_manager():
    """Create a BreakpointManager instance."""
    return BreakpointManager()


@pytest.fixture
async def layout(breakpoint_manager):
    """Create an EnhancedLayout instance."""
    layout = EnhancedLayout(breakpoint_manager)
    await layout.setup()
    return layout


@pytest.mark.asyncio
async def test_layout_initialization(layout):
    """Test layout initialization."""
    assert isinstance(layout.current_grid, GridConfig)
    assert isinstance(layout.current_layout, LayoutConfig)
    assert layout.current_grid.columns == 12
    assert layout.current_grid.gap == 16
    assert layout.current_layout.columns == 12
    assert layout.current_layout.gap == 16


@pytest.mark.asyncio
async def test_breakpoint_adjustment(layout):
    """Test layout adjustment to breakpoints."""
    # Test xs breakpoint
    await layout.adjust_to_breakpoint(Breakpoint("xs", 0, 575))
    assert layout.current_grid.columns == 4
    assert layout.current_grid.gap == 8
    assert layout.current_layout.columns == 4
    assert layout.current_layout.gap == 8

    # Test xl breakpoint
    await layout.adjust_to_breakpoint(Breakpoint("xl", 1200, None))
    assert layout.current_grid.columns == 12
    assert layout.current_grid.gap == 32
    assert layout.current_layout.columns == 12
    assert layout.current_layout.gap == 32


def test_container_styles(layout):
    """Test container style generation."""
    styles = layout.get_container_style()
    assert "padding" in styles
    assert "margin" in styles
    assert "max-width" in styles
    assert "min-height" in styles


def test_grid_styles(layout):
    """Test grid style generation."""
    # Test default columns
    styles = layout.get_grid_style()
    assert "display: grid" in styles["display"]
    assert "repeat(1, 1fr)" in styles["grid-template-columns"]
    assert "16px" in styles["gap"]

    # Test custom columns
    styles = layout.get_grid_style(columns=3)
    assert "repeat(3, 1fr)" in styles["grid-template-columns"]


def test_scientific_data_grid(layout):
    """Test scientific data grid configuration."""
    # Test data grid
    data_grid = layout.get_data_grid_style()
    assert "font-family: var(--data-font)" in data_grid
    assert "gap: 1.5rem" in data_grid

    # Test chart grid
    chart_grid = layout.get_chart_grid_style()
    assert "aspect-ratio: var(--chart-aspect-ratio)" in chart_grid
    assert "padding: 1rem" in chart_grid


def test_responsive_layout_generation(layout):
    """Test responsive layout class generation."""
    breakpoints = {
        "xs": {"p": 4, "gap": 2},
        "md": {"p": 6, "gap": 4},
        "lg": {"p": "2rem", "gap": "1.5rem"},
    }

    classes = layout.get_responsive_layout(breakpoints)
    assert "xs:p-4" in classes
    assert "xs:gap-2" in classes
    assert "md:p-6" in classes
    assert "md:gap-4" in classes
    assert "lg:p-[2rem]" in classes
    assert "lg:gap-[1.5rem]" in classes


def test_data_visualization_support(layout):
    """Test data visualization style generation."""
    # Test chart container
    chart_styles = layout.get_chart_container_style()
    assert "position: relative" in chart_styles
    assert "width: 100%" in chart_styles
    assert "contain: layout size style" in chart_styles

    # Test data table
    table_styles = layout.get_data_table_style()
    assert "font-family: var(--data-font)" in table_styles
    assert "font-variant-numeric: tabular-nums" in table_styles


def test_performance_optimizations(layout):
    """Test performance optimization features."""
    # Test GPU acceleration
    gpu_styles = layout.get_gpu_optimized_style()
    assert "transform: translateZ(0)" in gpu_styles
    assert "backface-visibility: hidden" in gpu_styles

    # Test content visibility
    visibility_styles = layout.get_content_visibility_style()
    assert "content-visibility: auto" in visibility_styles
    assert "contain-intrinsic-size: 0 500px" in visibility_styles


def test_print_optimization(layout):
    """Test print style generation."""
    print_styles = layout.get_print_styles()
    assert "max-width: none" in print_styles
    assert "page-break-inside: avoid" in print_styles
    assert "break-inside: avoid" in print_styles


def test_accessibility_features(layout):
    """Test accessibility style generation."""
    # Test reduced motion
    motion_styles = layout.get_reduced_motion_styles()
    assert "animation: none !important" in motion_styles
    assert "transition: none !important" in motion_styles

    # Test high contrast
    contrast_styles = layout.get_high_contrast_styles()
    assert "--color-data-primary: #000000" in contrast_styles
    assert "--color-graph-grid: #000000" in contrast_styles


def test_touch_device_optimization(layout):
    """Test touch device style generation."""
    touch_styles = layout.get_touch_optimized_styles()
    assert "min-height: 44px" in touch_styles
    assert "min-width: 44px" in touch_styles
    assert "padding: 0.75rem" in touch_styles


def test_scientific_visualization_enhancements(layout):
    """Test scientific visualization style generation."""
    # Test plot line styles
    plot_styles = layout.get_plot_line_style()
    assert "stroke-width: var(--plot-line-width)" in plot_styles
    assert "vector-effect: non-scaling-stroke" in plot_styles

    # Test axis line styles
    axis_styles = layout.get_axis_line_style()
    assert "stroke-width: var(--axis-line-width)" in axis_styles
    assert "shape-rendering: crispEdges" in axis_styles

    # Test grid line styles
    grid_styles = layout.get_grid_line_style()
    assert "stroke-width: var(--grid-line-width)" in grid_styles
    assert "stroke: var(--color-graph-grid)" in grid_styles


def test_high_dpi_screen_support(layout):
    """Test high DPI screen optimizations."""
    dpi_styles = layout.get_high_dpi_styles()
    assert "transform: translateZ(0)" in dpi_styles
    assert "image-rendering: crisp-edges" in dpi_styles


def test_data_annotation_styles(layout):
    """Test data annotation style generation."""
    annotation_styles = layout.get_data_annotation_style()
    assert "font-family: var(--data-font)" in annotation_styles
    assert "fill: var(--color-annotation)" in annotation_styles
    assert "font-size: 0.875rem" in annotation_styles


def test_responsive_container_widths(layout):
    """Test responsive container width calculations."""
    # Test fluid container
    assert layout.get_container_width("xs") == "100%"

    # Test fixed containers
    assert layout.get_container_width("sm") == "540px"
    assert layout.get_container_width("md") == "720px"
    assert layout.get_container_width("lg") == "960px"
    assert layout.get_container_width("xl") == "1140px"


def test_grid_column_calculations(layout):
    """Test grid column width calculations."""
    # Test different column spans
    assert layout.get_column_width(1) > 0
    assert layout.get_column_width(2) > layout.get_column_width(1)
    assert layout.get_column_width(12) > layout.get_column_width(6)


def test_spacing_calculations(layout):
    """Test spacing system calculations."""
    # Test gap sizes
    assert layout.get_gap_size("xs") == 8
    assert layout.get_gap_size("md") == 16
    assert layout.get_gap_size("xl") == 32

    # Test margin sizes
    assert layout.get_margin_size("xs") == 8
    assert layout.get_margin_size("md") == 16
    assert layout.get_margin_size("xl") == 32


def test_layout_error_handling(layout):
    """Test layout error handling."""
    # Test invalid breakpoint
    with pytest.raises(ValueError):
        layout.get_container_width("invalid")

    # Test invalid column count
    with pytest.raises(ValueError):
        layout.get_column_width(13)

    # Test invalid gap size
    with pytest.raises(ValueError):
        layout.get_gap_size("invalid")
