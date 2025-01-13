"""Tests for breakpoint manager."""

import pytest
from ..breakpoints import BreakpointManager, BreakpointConfig


@pytest.fixture
def breakpoint_manager():
    """Create a BreakpointManager instance."""
    return BreakpointManager()


def test_breakpoint_initialization(breakpoint_manager):
    """Test breakpoint manager initialization."""
    assert len(breakpoint_manager.breakpoints) == 5
    assert set(breakpoint_manager.breakpoints.keys()) == {"xs", "sm", "md", "lg", "xl"}
    assert breakpoint_manager.ordered_breakpoints == ["xs", "sm", "md", "lg", "xl"]


def test_breakpoint_configs(breakpoint_manager):
    """Test breakpoint configurations."""
    # Test xs config
    xs_config = breakpoint_manager.breakpoints["xs"]
    assert xs_config.name == "xs"
    assert xs_config.min_width == 0
    assert xs_config.max_width == 575
    assert xs_config.columns == 4
    assert xs_config.gap == 8
    assert xs_config.container_width is None

    # Test md config
    md_config = breakpoint_manager.breakpoints["md"]
    assert md_config.name == "md"
    assert md_config.min_width == 768
    assert md_config.max_width == 991
    assert md_config.columns == 12
    assert md_config.gap == 16
    assert md_config.container_width == 720

    # Test xl config
    xl_config = breakpoint_manager.breakpoints["xl"]
    assert xl_config.name == "xl"
    assert xl_config.min_width == 1200
    assert xl_config.max_width is None
    assert xl_config.columns == 12
    assert xl_config.gap == 32
    assert xl_config.container_width == 1140


def test_get_breakpoint(breakpoint_manager):
    """Test breakpoint detection."""
    # Test boundary conditions
    assert breakpoint_manager.get_breakpoint(0) == "xs"
    assert breakpoint_manager.get_breakpoint(575) == "xs"
    assert breakpoint_manager.get_breakpoint(576) == "sm"
    assert breakpoint_manager.get_breakpoint(767) == "sm"
    assert breakpoint_manager.get_breakpoint(768) == "md"
    assert breakpoint_manager.get_breakpoint(991) == "md"
    assert breakpoint_manager.get_breakpoint(992) == "lg"
    assert breakpoint_manager.get_breakpoint(1199) == "lg"
    assert breakpoint_manager.get_breakpoint(1200) == "xl"
    assert breakpoint_manager.get_breakpoint(9999) == "xl"


def test_get_config(breakpoint_manager):
    """Test getting breakpoint configuration."""
    config = breakpoint_manager.get_config("md")
    assert isinstance(config, BreakpointConfig)
    assert config.name == "md"
    assert config.min_width == 768
    assert config.max_width == 991


def test_container_widths(breakpoint_manager):
    """Test container width calculations."""
    assert breakpoint_manager.get_container_width("xs") == 100
    assert breakpoint_manager.get_container_width("sm") == 540
    assert breakpoint_manager.get_container_width("md") == 720
    assert breakpoint_manager.get_container_width("lg") == 960
    assert breakpoint_manager.get_container_width("xl") == 1140


def test_grid_columns(breakpoint_manager):
    """Test grid column calculations."""
    assert breakpoint_manager.get_grid_columns("xs") == 4
    assert breakpoint_manager.get_grid_columns("sm") == 8
    assert breakpoint_manager.get_grid_columns("md") == 12
    assert breakpoint_manager.get_grid_columns("lg") == 12
    assert breakpoint_manager.get_grid_columns("xl") == 12


def test_spacing(breakpoint_manager):
    """Test spacing calculations."""
    assert breakpoint_manager.get_spacing("xs") == 8
    assert breakpoint_manager.get_spacing("sm") == 12
    assert breakpoint_manager.get_spacing("md") == 16
    assert breakpoint_manager.get_spacing("lg") == 24
    assert breakpoint_manager.get_spacing("xl") == 32


def test_margins(breakpoint_manager):
    """Test margin calculations."""
    assert breakpoint_manager.get_margin("xs") == 8
    assert breakpoint_manager.get_margin("sm") == 12
    assert breakpoint_manager.get_margin("md") == 16
    assert breakpoint_manager.get_margin("lg") == 24
    assert breakpoint_manager.get_margin("xl") == 32


def test_device_type_checks(breakpoint_manager):
    """Test device type detection."""
    # Mobile
    assert breakpoint_manager.is_mobile("xs")
    assert breakpoint_manager.is_mobile("sm")
    assert not breakpoint_manager.is_mobile("md")
    assert not breakpoint_manager.is_mobile("lg")
    assert not breakpoint_manager.is_mobile("xl")

    # Tablet
    assert not breakpoint_manager.is_tablet("xs")
    assert not breakpoint_manager.is_tablet("sm")
    assert breakpoint_manager.is_tablet("md")
    assert not breakpoint_manager.is_tablet("lg")
    assert not breakpoint_manager.is_tablet("xl")

    # Desktop
    assert not breakpoint_manager.is_desktop("xs")
    assert not breakpoint_manager.is_desktop("sm")
    assert not breakpoint_manager.is_desktop("md")
    assert breakpoint_manager.is_desktop("lg")
    assert breakpoint_manager.is_desktop("xl")


def test_media_queries(breakpoint_manager):
    """Test media query generation."""
    # Min-width queries
    assert breakpoint_manager.get_media_query("sm") == "@media (min-width: 576px)"
    assert breakpoint_manager.get_media_query("md") == "@media (min-width: 768px)"

    # Max-width queries
    assert breakpoint_manager.get_media_query("sm", "max") == "@media (max-width: 767px)"
    assert breakpoint_manager.get_media_query("md", "max") == "@media (max-width: 991px)"

    # XL breakpoint (no max-width)
    assert breakpoint_manager.get_media_query("xl", "max") == "@media (min-width: 1200px)"


def test_breakpoint_ranges(breakpoint_manager):
    """Test breakpoint range queries."""
    # Test ranges
    assert breakpoint_manager.get_breakpoint_range("sm", "md") == ("@media (min-width: 576px) and (max-width: 991px)")
    assert breakpoint_manager.get_breakpoint_range("md", "lg") == ("@media (min-width: 768px) and (max-width: 1199px)")
    assert breakpoint_manager.get_breakpoint_range("lg", "xl") == ("@media (min-width: 992px) and (max-width: 9999px)")


def test_breakpoint_navigation(breakpoint_manager):
    """Test breakpoint navigation methods."""
    # Next breakpoint
    assert breakpoint_manager.get_next_breakpoint("xs") == "sm"
    assert breakpoint_manager.get_next_breakpoint("sm") == "md"
    assert breakpoint_manager.get_next_breakpoint("md") == "lg"
    assert breakpoint_manager.get_next_breakpoint("lg") == "xl"
    assert breakpoint_manager.get_next_breakpoint("xl") is None

    # Previous breakpoint
    assert breakpoint_manager.get_prev_breakpoint("xl") == "lg"
    assert breakpoint_manager.get_prev_breakpoint("lg") == "md"
    assert breakpoint_manager.get_prev_breakpoint("md") == "sm"
    assert breakpoint_manager.get_prev_breakpoint("sm") == "xs"
    assert breakpoint_manager.get_prev_breakpoint("xs") is None

    # Breakpoints up
    assert breakpoint_manager.get_breakpoints_up("md") == ["md", "lg", "xl"]
    assert breakpoint_manager.get_breakpoints_up("xl") == ["xl"]
    assert breakpoint_manager.get_breakpoints_up("invalid") == []

    # Breakpoints down
    assert breakpoint_manager.get_breakpoints_down("md") == ["xs", "sm", "md"]
    assert breakpoint_manager.get_breakpoints_down("xs") == ["xs"]
    assert breakpoint_manager.get_breakpoints_down("invalid") == []


def test_breakpoint_bounds(breakpoint_manager):
    """Test getting breakpoint bounds."""
    # Test min/max width pairs
    assert breakpoint_manager.get_breakpoint_bounds("xs") == (0, 575)
    assert breakpoint_manager.get_breakpoint_bounds("sm") == (576, 767)
    assert breakpoint_manager.get_breakpoint_bounds("md") == (768, 991)
    assert breakpoint_manager.get_breakpoint_bounds("lg") == (992, 1199)
    assert breakpoint_manager.get_breakpoint_bounds("xl") == (1200, None)
