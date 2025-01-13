"""Tests for enhanced layout system."""

import pytest
from ..layout import EnhancedLayout, GridConfig, LayoutConfig
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


def test_container_classes(layout):
    """Test container class generation."""
    # Test fluid container
    classes = layout.get_container_classes("md", fluid=True)
    assert "container" in classes
    assert "mx-auto" in classes
    assert "max-w-" not in classes

    # Test fixed container
    classes = layout.get_container_classes("md", fluid=False)
    assert "container" in classes
    assert "mx-auto" in classes
    assert "max-w-720" in classes

    # Test non-centered container
    classes = layout.get_container_classes("md", centered=False)
    assert "mx-auto" not in classes


def test_container_style(layout):
    """Test container style generation."""
    styles = layout.get_container_style()
    assert "padding" in styles
    assert "margin" in styles
    assert "max-width" in styles
    assert "min-height" in styles


def test_grid_style(layout):
    """Test grid style generation."""
    # Test default columns
    styles = layout.get_grid_style()
    assert "display: grid" in styles["display"]
    assert "repeat(1, 1fr)" in styles["grid-template-columns"]
    assert "16px" in styles["gap"]

    # Test custom columns
    styles = layout.get_grid_style(columns=3)
    assert "repeat(3, 1fr)" in styles["grid-template-columns"]


def test_column_width_calculation(layout):
    """Test column width calculations."""
    # Test single column
    width = layout.get_column_width(1)
    assert isinstance(width, int)
    assert width > 0

    # Test multiple columns
    width_3 = layout.get_column_width(3)
    width_1 = layout.get_column_width(1)
    assert width_3 > width_1


def test_stack_classes(layout):
    """Test stack class generation."""
    # Test default stack
    classes = layout.get_stack_classes("md")
    assert "flex" in classes
    assert "flex-col" in classes
    assert "gap-16" in classes

    # Test reversed stack
    classes = layout.get_stack_classes("md", reverse=True)
    assert "flex-col-reverse" in classes

    # Test custom spacing
    classes = layout.get_stack_classes("md", spacing=24)
    assert "gap-24" in classes


def test_row_classes(layout):
    """Test row class generation."""
    # Test default row
    classes = layout.get_row_classes("md")
    assert "flex" in classes
    assert "flex-wrap" in classes
    assert "items-start" in classes
    assert "justify-start" in classes
    assert "gap-16" in classes

    # Test custom alignment
    classes = layout.get_row_classes("md", align="center", justify="between")
    assert "items-center" in classes
    assert "justify-between" in classes

    # Test no wrap
    classes = layout.get_row_classes("md", wrap=False)
    assert "flex-wrap" not in classes


def test_responsive_layout(layout):
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


def test_sidebar_layout(layout):
    """Test sidebar layout class generation."""
    # Test numeric width
    classes = layout.get_sidebar_layout("md", sidebar_width=3)
    assert "grid" in classes
    assert "sidebar 9fr" in classes
    assert "gap-16" in classes

    # Test string width
    classes = layout.get_sidebar_layout("md", sidebar_width="300px", sidebar_position="right")
    assert "1fr 300px" in classes

    # Test custom spacing
    classes = layout.get_sidebar_layout("md", sidebar_width=3, spacing=24)
    assert "gap-24" in classes


def test_responsive_sidebar(layout):
    """Test responsive sidebar class generation."""
    breakpoints = {
        "xs": {"width": "100%", "position": "top"},
        "md": {"width": 3, "position": "left"},
        "lg": {"width": "300px", "position": "right", "gap": 24},
    }

    classes = layout.get_responsive_sidebar(breakpoints)
    assert "grid" in classes
    assert "md:grid-template-columns: sidebar 9fr" in classes
    assert "lg:grid-template-columns: 1fr 300px" in classes
    assert "lg:gap-24" in classes


def test_base_classes_integration(layout):
    """Test integration of base classes."""
    base = ["custom-class", "another-class"]

    # Test with responsive layout
    classes = layout.get_responsive_layout({"md": {"p": 4}}, base_classes=base)
    assert "custom-class" in classes
    assert "another-class" in classes
    assert "md:p-4" in classes

    # Test with responsive sidebar
    classes = layout.get_responsive_sidebar({"md": {"width": 3}}, base_classes=base)
    assert "custom-class" in classes
    assert "another-class" in classes
    assert "md:grid-template-columns" in classes
