"""Tests for enhanced responsive list component with scientific data capabilities."""

import pytest
from unittest.mock import Mock, patch
from typing import List, Dict, Any

from ..list import (
    ListDisplayMode,
    ScientificMetadata,
    ScientificListItem,
    ListConfig,
    EnhancedResponsiveList,
)
from ..viewport import EnhancedViewportManager


@pytest.fixture
def viewport_manager():
    """Create a viewport manager instance."""
    return EnhancedViewportManager()


@pytest.fixture
def scientific_metadata() -> ScientificMetadata:
    """Create test scientific metadata."""
    return ScientificMetadata(
        data_quality=0.95,
        confidence_level=0.9,
        validation_status="verified",
        measurement_method="HPLC-MS",
        instrumentation="Agilent 7890B",
        resolution=0.001,
        precision=0.0005,
        uncertainty=0.002,
        units="g/mol",
        conditions={"temperature": "25C", "pressure": "1 atm"},
    )


@pytest.fixture
def test_items(scientific_metadata: ScientificMetadata) -> List[ScientificListItem]:
    """Create test list items with comprehensive scientific data."""
    return [
        ScientificListItem(
            id="CMPD-001",
            title="Test Compound 1",
            subtitle="Test Series A",
            description="A test compound for validation",
            scientific_metadata=scientific_metadata,
            structure_data={
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-4-2-3-5-9(8)7(11)12/h2-5H,1H3,(H,11,12)",
                "molecular_weight": 180.159,
            },
            quantum_data={
                "energy": -512.345,
                "dipole": 2.34,
                "orbital_energies": [-10.2, -9.8, -8.5],
            },
            spectral_data={
                "nmr_peaks": [(2.1, "s", 3), (7.2, "m", 5)],
                "ir_bands": [1720, 1680, 1450],
            },
            crystal_data={
                "space_group": "P21/c",
                "cell_parameters": {"a": 10.2, "b": 15.3, "c": 8.7},
            },
            analysis_results={
                "purity": 99.5,
                "solubility": 125.3,
                "stability": "high",
            },
            molecular_properties={
                "logP": 2.34,
                "PSA": 45.6,
                "HBA": 3,
                "HBD": 1,
            },
            binding_data={"5HT2A": 8.2, "D2": 7.1},
            safety_data={"ld50": 1000},
            citations=[
                {"text": "Smith et al., J. Med. Chem. 2023", "doi": "10.1021/example"},
            ],
            tags=["validated", "high-quality", "reference"],
        ),
        ScientificListItem(
            id="CMPD-002",
            title="Test Compound 2",
            subtitle="Test Series B",
            description="Another test compound",
            scientific_metadata=scientific_metadata,
            structure_data={
                "smiles": "CC1=CC=C(C=C1)C(=O)O",
                "molecular_weight": 136.15,
            },
            molecular_properties={
                "logP": 1.89,
                "PSA": 37.3,
            },
            binding_data={"5HT2A": 7.8, "D2": 6.9},
            safety_data={"ld50": 800},
            tags=["preliminary"],
        ),
    ]


@pytest.fixture
def list_config() -> ListConfig:
    """Create test list configuration."""
    return ListConfig(
        display_mode=ListDisplayMode.GRID,
        page_size=10,
        show_structures=True,
        show_predictions=True,
        show_analysis=True,
        show_visualizations=True,
        enable_filtering=True,
        enable_sorting=True,
        scientific_precision=True,
        molecular_detail="high",
        quantum_enabled=True,
        enable_accessibility=True,
        enable_3d_visualization=True,
        high_performance=True,
    )


@pytest.fixture
async def list_component(list_config: ListConfig) -> EnhancedResponsiveList:
    """Create test list component."""
    component = EnhancedResponsiveList(config=list_config)
    await component.setup()
    return component


@pytest.mark.asyncio
async def test_initialization(list_component: EnhancedResponsiveList):
    """Test component initialization."""
    assert list_component.config is not None
    assert list_component.layout is not None
    assert len(list_component._active_optimizations) > 0


@pytest.mark.asyncio
async def test_update_items(list_component: EnhancedResponsiveList, test_items: List[ScientificListItem]):
    """Test updating list items."""
    await list_component.update_items(test_items)
    assert len(list_component._items) == 2
    assert len(list_component._filtered_items) == 2
    assert list_component._total_pages == 1


@pytest.mark.asyncio
async def test_filtering(list_component: EnhancedResponsiveList, test_items: List[ScientificListItem]):
    """Test item filtering."""
    await list_component.update_items(test_items)

    # Test exact match filter
    list_component.set_filters({"id": "CMPD-001"})
    assert len(list_component._filtered_items) == 1
    assert list_component._filtered_items[0].id == "CMPD-001"

    # Test range filter
    list_component.set_filters({"molecular_properties.logP": {"min": 2.0, "max": 3.0}})
    assert len(list_component._filtered_items) == 1
    assert list_component._filtered_items[0].id == "CMPD-001"

    # Test binding data filter
    list_component.set_filters({"binding_data.5HT2A": {"min": 8.0}})
    assert len(list_component._filtered_items) == 1
    assert list_component._filtered_items[0].id == "CMPD-001"


@pytest.mark.asyncio
async def test_sorting(list_component: EnhancedResponsiveList, test_items: List[ScientificListItem]):
    """Test item sorting."""
    await list_component.update_items(test_items)

    # Test ascending sort
    list_component.set_sort("title", ascending=True)
    assert list_component._filtered_items[0].title < list_component._filtered_items[1].title

    # Test descending sort
    list_component.set_sort("title", ascending=False)
    assert list_component._filtered_items[0].title > list_component._filtered_items[1].title

    # Test scientific data sorting
    list_component.set_sort("molecular_properties.logP", ascending=True)
    assert list_component._filtered_items[0].molecular_properties["logP"] < list_component._filtered_items[1].molecular_properties["logP"]


@pytest.mark.asyncio
async def test_pagination(list_component: EnhancedResponsiveList, test_items: List[ScientificListItem]):
    """Test pagination."""
    # Set small page size
    list_component.config.page_size = 1
    await list_component.update_items(test_items)

    # Check total pages
    assert list_component._total_pages == 2

    # Check page items
    page1 = list_component.get_page_items(1)
    assert len(page1) == 1
    assert page1[0].id == "CMPD-001"

    page2 = list_component.get_page_items(2)
    assert len(page2) == 1
    assert page2[0].id == "CMPD-002"


@pytest.mark.asyncio
async def test_selection(list_component: EnhancedResponsiveList, test_items: List[ScientificListItem]):
    """Test item selection."""
    await list_component.update_items(test_items)

    # Test select
    list_component.select_item("CMPD-001")
    assert "CMPD-001" in list_component._selected_items
    selected = list_component.get_selected_items()
    assert len(selected) == 1
    assert selected[0].id == "CMPD-001"

    # Test deselect
    list_component.deselect_item("CMPD-001")
    assert "CMPD-001" not in list_component._selected_items
    assert len(list_component.get_selected_items()) == 0

    # Test clear
    list_component.select_item("CMPD-001")
    list_component.select_item("CMPD-002")
    list_component.clear_selection()
    assert len(list_component._selected_items) == 0


@pytest.mark.asyncio
async def test_scientific_data_validation(list_component: EnhancedResponsiveList, test_items: List[ScientificListItem]):
    """Test scientific data validation."""
    await list_component.update_items(test_items)

    # Check validation flags
    assert "data_validated" in list_component._active_optimizations

    # Check validation status
    for item in list_component._items:
        assert item.scientific_metadata is not None
        assert item.scientific_metadata.validation_status == "verified"

    # Test data completeness validation
    assert list_component._check_data_completeness(test_items[0]) is True

    # Test with missing required data
    incomplete_item = test_items[0]
    incomplete_item.structure_data = None
    assert list_component._check_data_completeness(incomplete_item) is False


@pytest.mark.asyncio
async def test_layout_configuration(list_component: EnhancedResponsiveList):
    """Test layout configuration."""
    # Check layout areas
    assert "structure_viewer" in list_component.layout._active_areas
    assert "data_panel" in list_component.layout._active_areas
    if list_component.config.show_visualizations:
        assert "visualization_panel" in list_component.layout._active_areas

    # Check layout styles
    styles = list_component.get_layout_styles()
    assert "display" in styles
    assert "grid-template-columns" in styles

    # Check area styles
    structure_styles = list_component.get_area_styles("structure_viewer")
    assert "grid-area" in structure_styles
    assert "overflow" in structure_styles


@pytest.mark.asyncio
async def test_responsive_behavior(list_component: EnhancedResponsiveList):
    """Test responsive behavior."""
    viewport_sizes = [(320, 480), (768, 1024), (1920, 1080)]

    for width, height in viewport_sizes:
        await list_component.handle_resize(width)
        styles = list_component.get_layout_styles()

        if width <= 576:  # Mobile
            assert list_component.layout.config.columns <= 4
        elif width <= 768:  # Tablet
            assert list_component.layout.config.columns <= 8
        else:  # Desktop
            assert list_component.layout.config.columns <= 12


@pytest.mark.asyncio
async def test_optimizations(list_component: EnhancedResponsiveList):
    """Test optimization flags."""
    # Check performance optimizations
    if list_component.config.high_performance:
        assert "gpu_acceleration" in list_component._active_optimizations
        assert "content_visibility" in list_component._active_optimizations

    # Check scientific optimizations
    if list_component.config.scientific_precision:
        assert "high_precision" in list_component._active_optimizations
        assert "tabular_nums" in list_component._active_optimizations

    # Check visualization optimizations
    if list_component.config.enable_3d_visualization:
        assert "3d_ready" in list_component._active_optimizations

    # Check accessibility optimizations
    if list_component.config.enable_accessibility:
        assert "screen_reader_support" in list_component._active_optimizations
        assert "keyboard_navigation" in list_component._active_optimizations


@pytest.mark.asyncio
async def test_scientific_data_processing(list_component: EnhancedResponsiveList, test_items: List[ScientificListItem]):
    """Test scientific data processing pipeline."""
    await list_component.update_items(test_items)

    # Verify quantum data processing
    if test_items[0].quantum_data:
        assert "quantum_ready" in list_component._active_optimizations

    # Verify spectral data processing
    if test_items[0].spectral_data:
        assert "spectral_ready" in list_component._active_optimizations

    # Verify crystal data processing
    if test_items[0].crystal_data:
        assert "crystal_ready" in list_component._active_optimizations


@pytest.mark.asyncio
async def test_display_modes(list_component: EnhancedResponsiveList, test_items: List[ScientificListItem]):
    """Test different display modes."""
    await list_component.update_items(test_items)

    display_modes = [
        (ListDisplayMode.MOLECULAR, "molecular-view"),
        (ListDisplayMode.QUANTUM, "quantum-view"),
        (ListDisplayMode.SPECTRAL, "spectral-view"),
        (ListDisplayMode.CRYSTALLOGRAPHIC, "crystal-view"),
        (ListDisplayMode.ANALYSIS, "analysis-section"),
    ]

    for mode, expected_class in display_modes:
        list_component.config.display_mode = mode
        styles = list_component.get_layout_styles()
        assert styles is not None
