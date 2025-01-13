"""Tests for grid system."""

import pytest
from ..grid import Grid
from ..breakpoints import BreakpointManager


@pytest.fixture
def breakpoint_manager():
    """Create a BreakpointManager instance."""
    return BreakpointManager()


@pytest.fixture
def grid(breakpoint_manager):
    """Create a Grid instance."""
    return Grid(breakpoint_manager)


def test_grid_initialization(grid):
    """Test grid initialization."""
    assert grid.default_gap == 16
    assert grid.default_margin == 16


def test_grid_classes(grid):
    """Test grid class generation."""
    # Test default columns
    classes = grid.get_grid_classes("md")
    assert "grid" in classes
    assert "grid-cols-12" in classes
    assert "gap-16" in classes

    # Test custom columns
    classes = grid.get_grid_classes("md", columns=6)
    assert "grid-cols-6" in classes


def test_container_classes(grid):
    """Test container class generation."""
    # Test xs container (fluid)
    classes = grid.get_container_classes("xs")
    assert "container" in classes
    assert "mx-auto" in classes
    assert "px-8" in classes
    assert "max-w-" not in classes

    # Test md container (fixed width)
    classes = grid.get_container_classes("md")
    assert "container" in classes
    assert "mx-auto" in classes
    assert "px-16" in classes
    assert "max-w-720" in classes


def test_column_classes(grid):
    """Test column class generation."""
    # Test basic span
    classes = grid.get_column_classes("md", span=6)
    assert "col-span-6" in classes
    assert "col-start-" not in classes

    # Test span with offset
    classes = grid.get_column_classes("md", span=6, offset=2)
    assert "col-span-6" in classes
    assert "col-start-3" in classes  # offset + 1


def test_responsive_classes(grid):
    """Test responsive class generation."""
    breakpoints = {
        "xs": 12,
        "md": {"span": 6, "offset": 3},
        "lg": {"span": 4, "offset": 4},
    }

    classes = grid.get_responsive_classes(breakpoints)
    assert "xs:col-span-12" in classes
    assert "md:col-span-6" in classes
    assert "md:col-start-4" in classes
    assert "lg:col-span-4" in classes
    assert "lg:col-start-5" in classes


def test_grid_template(grid):
    """Test grid template generation."""
    areas = [
        ["header", "header", "header"],
        ["nav", "main", "aside"],
        ["footer", "footer", "footer"],
    ]
    row_heights = ["auto", "1fr", "auto"]

    styles = grid.get_grid_template("md", areas, row_heights)
    assert "display: grid;" in styles
    assert 'grid-template-areas: "header header header"' in styles
    assert "grid-template-rows: auto 1fr auto;" in styles
    assert "grid-template-columns: repeat(3, 1fr);" in styles
    assert "gap: 16px;" in styles


def test_area_style(grid):
    """Test grid area style generation."""
    style = grid.get_area_style("header")
    assert style == "grid-area: header;"


def test_responsive_grid_template(grid):
    """Test responsive grid template generation."""
    breakpoints = {
        "xs": [["main"], ["nav"], ["aside"]],  # Stack on mobile
        "md": [["nav", "main", "aside"]],  # Side by side on desktop
    }

    styles = grid.get_responsive_grid_template(breakpoints)
    assert "@media (min-width: 0px)" in styles
    assert "@media (min-width: 768px)" in styles
    assert 'grid-template-areas: "main"' in styles
    assert 'grid-template-areas: "nav main aside"' in styles


def test_auto_grid_classes(grid):
    """Test auto-responsive grid class generation."""
    # Test with default max width
    classes = grid.get_auto_grid_classes("md", min_width="200px")
    assert "grid" in classes
    assert "grid-auto-flow-dense" in classes
    assert "gap-16" in classes
    assert "grid-template-columns: repeat(auto-fit, minmax(200px))" in classes

    # Test with custom max width
    classes = grid.get_auto_grid_classes("md", min_width="200px", max_width="400px")
    assert "minmax(200px, 400px)" in classes


def test_masonry_grid_classes(grid):
    """Test masonry grid class generation."""
    # Test default responsive columns
    classes = grid.get_masonry_grid_classes("md")
    assert "columns" in classes
    assert "gap-16" in classes
    assert "columns-1" in classes
    assert "sm:columns-2" in classes
    assert "md:columns-3" in classes
    assert "lg:columns-4" in classes
    assert "xl:columns-5" in classes

    # Test custom columns
    custom_columns = {
        "sm": 3,
        "md": 4,
        "lg": 6,
    }
    classes = grid.get_masonry_grid_classes("md", columns=custom_columns)
    assert "sm:columns-3" in classes
    assert "md:columns-4" in classes
    assert "lg:columns-6" in classes


def test_flex_grid_classes(grid):
    """Test flex grid class generation."""
    # Test defaults
    classes = grid.get_flex_grid_classes("md")
    assert "flex" in classes
    assert "gap-16" in classes
    assert "items-start" in classes
    assert "justify-start" in classes
    assert "flex-wrap" in classes

    # Test custom alignment
    classes = grid.get_flex_grid_classes("md", align="center", justify="between", wrap=False)
    assert "items-center" in classes
    assert "justify-between" in classes
    assert "flex-wrap" not in classes


def test_responsive_spacing(grid):
    """Test responsive spacing class generation."""
    spacing = {
        "xs": {"p": 4, "gap": 2},
        "md": {"p": 6, "gap": 4},
        "lg": {"p": "2rem", "gap": "1.5rem"},
    }

    classes = grid.get_responsive_spacing(spacing)
    assert "xs:p-4" in classes
    assert "xs:gap-2" in classes
    assert "md:p-6" in classes
    assert "md:gap-4" in classes
    assert "lg:p-[2rem]" in classes
    assert "lg:gap-[1.5rem]" in classes
