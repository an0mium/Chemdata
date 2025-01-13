"""Tests for viewport manager."""

import pytest
from ..viewport import ViewportManager, ViewportState


@pytest.fixture
async def viewport_manager():
    """Create a ViewportManager instance."""
    manager = ViewportManager()
    await manager.setup()
    return manager


@pytest.mark.asyncio
async def test_viewport_initialization(viewport_manager):
    """Test viewport manager initialization."""
    assert viewport_manager.state.width == 1024
    assert viewport_manager.state.height == 768
    assert viewport_manager.state.breakpoint == "md"
    assert not viewport_manager.state.is_mobile
    assert not viewport_manager.state.is_touch


@pytest.mark.asyncio
async def test_breakpoint_detection(viewport_manager):
    """Test breakpoint detection for different widths."""
    # Test xs breakpoint
    assert viewport_manager.get_breakpoint(400) == "xs"

    # Test sm breakpoint
    assert viewport_manager.get_breakpoint(600) == "sm"

    # Test md breakpoint
    assert viewport_manager.get_breakpoint(800) == "md"

    # Test lg breakpoint
    assert viewport_manager.get_breakpoint(1100) == "lg"

    # Test xl breakpoint
    assert viewport_manager.get_breakpoint(1300) == "xl"


@pytest.mark.asyncio
async def test_mobile_detection(viewport_manager):
    """Test mobile width detection."""
    assert viewport_manager.is_mobile_width(400)  # xs
    assert viewport_manager.is_mobile_width(576)  # sm
    assert not viewport_manager.is_mobile_width(800)  # md
    assert not viewport_manager.is_mobile_width(1200)  # xl


@pytest.mark.asyncio
async def test_touch_detection(viewport_manager):
    """Test touch device detection."""
    # Update viewport to mobile size
    await viewport_manager.update_viewport(400, 800)
    assert viewport_manager.state.is_touch

    # Update viewport to desktop size
    await viewport_manager.update_viewport(1200, 800)
    assert not viewport_manager.state.is_touch


@pytest.mark.asyncio
async def test_viewport_update(viewport_manager):
    """Test viewport state updates."""
    # Update to mobile size
    await viewport_manager.update_viewport(400, 800)
    assert viewport_manager.state.width == 400
    assert viewport_manager.state.height == 800
    assert viewport_manager.state.breakpoint == "xs"
    assert viewport_manager.state.is_mobile
    assert viewport_manager.state.is_touch

    # Update to desktop size
    await viewport_manager.update_viewport(1200, 800)
    assert viewport_manager.state.width == 1200
    assert viewport_manager.state.height == 800
    assert viewport_manager.state.breakpoint == "xl"
    assert not viewport_manager.state.is_mobile
    assert not viewport_manager.state.is_touch


@pytest.mark.asyncio
async def test_viewport_listeners(viewport_manager):
    """Test viewport state change listeners."""
    states = []

    def listener(state: ViewportState):
        states.append(state)

    # Add listener
    viewport_manager.add_listener(listener)

    # Update viewport
    await viewport_manager.update_viewport(400, 800)
    assert len(states) == 1
    assert states[0].width == 400
    assert states[0].breakpoint == "xs"

    # Update again
    await viewport_manager.update_viewport(1200, 800)
    assert len(states) == 2
    assert states[1].width == 1200
    assert states[1].breakpoint == "xl"

    # Remove listener
    viewport_manager.remove_listener(listener)

    # Update should not trigger listener
    await viewport_manager.update_viewport(800, 600)
    assert len(states) == 2


@pytest.mark.asyncio
async def test_container_widths(viewport_manager):
    """Test container width calculations."""
    assert viewport_manager.get_container_width("xs") == viewport_manager.state.width
    assert viewport_manager.get_container_width("sm") == 540
    assert viewport_manager.get_container_width("md") == 720
    assert viewport_manager.get_container_width("lg") == 960
    assert viewport_manager.get_container_width("xl") == 1140


@pytest.mark.asyncio
async def test_grid_columns(viewport_manager):
    """Test grid column calculations."""
    assert viewport_manager.get_grid_columns("xs") == 4
    assert viewport_manager.get_grid_columns("sm") == 8
    assert viewport_manager.get_grid_columns("md") == 12
    assert viewport_manager.get_grid_columns("lg") == 12
    assert viewport_manager.get_grid_columns("xl") == 12


@pytest.mark.asyncio
async def test_spacing_units(viewport_manager):
    """Test spacing unit calculations."""
    assert viewport_manager.get_spacing("xs") == 8
    assert viewport_manager.get_spacing("sm") == 8
    assert viewport_manager.get_spacing("md") == 16
    assert viewport_manager.get_spacing("lg") == 24
    assert viewport_manager.get_spacing("xl") == 24


def test_viewport_meta_tag(viewport_manager):
    """Test viewport meta tag generation."""
    meta_tag = viewport_manager.get_viewport_meta_tag()
    assert 'name="viewport"' in meta_tag
    assert "width=device-width" in meta_tag
    assert "initial-scale=1.0" in meta_tag
    assert "maximum-scale=5.0" in meta_tag
    assert "user-scalable=yes" in meta_tag


def test_media_queries(viewport_manager):
    """Test media query generation."""
    # Test min-width queries
    assert viewport_manager.get_media_query("sm") == "@media (min-width: 576px)"
    assert viewport_manager.get_media_query("md") == "@media (min-width: 768px)"

    # Test max-width queries
    assert viewport_manager.get_media_query("sm", "max") == "@media (max-width: 767px)"
    assert viewport_manager.get_media_query("md", "max") == "@media (max-width: 991px)"

    # Test xl breakpoint (no max-width)
    assert viewport_manager.get_media_query("xl", "max") == "@media (min-width: 1200px)"


def test_breakpoint_ranges(viewport_manager):
    """Test breakpoint range queries."""
    # Test sm to md range
    assert viewport_manager.get_breakpoint_range("sm", "md") == ("@media (min-width: 576px) and (max-width: 991px)")

    # Test md to lg range
    assert viewport_manager.get_breakpoint_range("md", "lg") == ("@media (min-width: 768px) and (max-width: 1199px)")

    # Test range including xl (no max-width)
    assert viewport_manager.get_breakpoint_range("lg", "xl") == ("@media (min-width: 992px) and (max-width: 9999px)")
