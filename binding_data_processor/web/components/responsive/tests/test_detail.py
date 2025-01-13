"""Tests for enhanced responsive detail component."""

import pytest
from unittest.mock import Mock, patch

from ..detail import (
    DetailDisplayMode,
    DetailConfig,
    EnhancedResponsiveDetail,
    ScientificListItem,
    ScientificMetadata,
)


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
def test_item(scientific_metadata: ScientificMetadata) -> ScientificListItem:
    """Create test item with comprehensive scientific data."""
    return ScientificListItem(
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
    )


@pytest.fixture
def comparison_item(scientific_metadata: ScientificMetadata) -> ScientificListItem:
    """Create comparison test item."""
    return ScientificListItem(
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
    )


@pytest.fixture
def detail_config() -> DetailConfig:
    """Create test detail configuration."""
    return DetailConfig(
        display_mode=DetailDisplayMode.OVERVIEW,
        show_3d_structure=True,
        show_quantum_viz=True,
        show_spectral_data=True,
        show_crystal_data=True,
        show_analysis=True,
        show_binding=True,
        show_safety=True,
        show_literature=True,
        enable_comparison=True,
        scientific_precision=True,
        molecular_detail="high",
        quantum_enabled=True,
        visualization_quality="high",
        data_precision_level="high",
        enable_structure_interaction=True,
        enable_accessibility=True,
        high_performance=True,
    )


@pytest.fixture
async def detail_component(detail_config: DetailConfig) -> EnhancedResponsiveDetail:
    """Create test detail component."""
    component = EnhancedResponsiveDetail(config=detail_config)
    await component.setup()
    return component


@pytest.mark.asyncio
async def test_initialization(detail_component: EnhancedResponsiveDetail):
    """Test component initialization."""
    assert detail_component.config is not None
    assert detail_component.layout is not None
    assert len(detail_component._active_optimizations) > 0
    assert detail_component._current_item is None
    assert detail_component._comparison_item is None


@pytest.mark.asyncio
async def test_set_item(detail_component: EnhancedResponsiveDetail, test_item: ScientificListItem):
    """Test setting current item."""
    await detail_component.set_item(test_item)
    assert detail_component._current_item == test_item

    # Check optimization flags
    assert "structure_ready" in detail_component._active_optimizations
    assert "quantum_ready" in detail_component._active_optimizations
    assert "spectral_ready" in detail_component._active_optimizations
    assert "crystal_ready" in detail_component._active_optimizations


@pytest.mark.asyncio
async def test_set_comparison_item(
    detail_component: EnhancedResponsiveDetail,
    test_item: ScientificListItem,
    comparison_item: ScientificListItem,
):
    """Test comparison functionality."""
    await detail_component.set_item(test_item)
    await detail_component.set_comparison_item(comparison_item)

    assert detail_component._current_item == test_item
    assert detail_component._comparison_item == comparison_item
    assert detail_component.layout.config.layout_type == detail_component.config.layout_type


@pytest.mark.asyncio
async def test_display_modes(detail_component: EnhancedResponsiveDetail, test_item: ScientificListItem):
    """Test different display modes."""
    await detail_component.set_item(test_item)

    # Test structure mode
    detail_component.config.display_mode = DetailDisplayMode.STRUCTURE
    await detail_component._update_display()
    assert "structure_ready" in detail_component._active_optimizations
    assert "high_quality_structure" in detail_component._active_optimizations

    # Test quantum mode
    detail_component.config.display_mode = DetailDisplayMode.QUANTUM
    await detail_component._update_display()
    assert "quantum_ready" in detail_component._active_optimizations
    assert "high_quality_quantum" in detail_component._active_optimizations

    # Test spectral mode
    detail_component.config.display_mode = DetailDisplayMode.SPECTRAL
    await detail_component._update_display()
    assert "spectral_ready" in detail_component._active_optimizations
    assert "high_quality_spectral" in detail_component._active_optimizations


@pytest.mark.asyncio
async def test_responsive_behavior(detail_component: EnhancedResponsiveDetail):
    """Test responsive layout behavior."""
    viewport_sizes = [(320, 480), (768, 1024), (1920, 1080)]

    for width, height in viewport_sizes:
        await detail_component.handle_resize(width)

        if width < 576:  # Mobile
            assert detail_component.config.structure_columns == detail_component.config.grid_columns
            assert detail_component.config.structure_width == width - 32
        elif width < 768:  # Tablet
            assert detail_component.config.structure_columns == detail_component.config.grid_columns
            assert detail_component.config.structure_width == width - 48
        else:  # Desktop
            assert detail_component.config.structure_columns == 6
            assert detail_component.config.structure_width == (width - 72) // 2


@pytest.mark.asyncio
async def test_optimizations(detail_component: EnhancedResponsiveDetail):
    """Test optimization flags."""
    # Check performance optimizations
    if detail_component.config.high_performance:
        assert "gpu_acceleration" in detail_component._active_optimizations
        assert "content_visibility" in detail_component._active_optimizations
        assert "will_change_transform" in detail_component._active_optimizations

    # Check scientific optimizations
    if detail_component.config.scientific_precision:
        assert "high_precision" in detail_component._active_optimizations
        assert "tabular_nums" in detail_component._active_optimizations
        assert "monospace_font" in detail_component._active_optimizations

    # Check visualization optimizations
    if detail_component.config.show_3d_structure:
        assert "3d_ready" in detail_component._active_optimizations
        if detail_component.config.visualization_quality == "high":
            assert "high_quality_3d" in detail_component._active_optimizations

    # Check accessibility optimizations
    if detail_component.config.enable_accessibility:
        assert "screen_reader_support" in detail_component._active_optimizations
        assert "keyboard_navigation" in detail_component._active_optimizations
        assert "high_contrast" in detail_component._active_optimizations


@pytest.mark.asyncio
async def test_layout_areas(detail_component: EnhancedResponsiveDetail):
    """Test layout area configuration."""
    # Check structure viewer area
    structure_area = detail_component.layout._active_areas.get("structure_viewer")
    assert structure_area is not None
    assert structure_area.molecular_view
    assert structure_area.quantum_ready == detail_component.config.quantum_enabled

    # Check data panel area
    data_area = detail_component.layout._active_areas.get("data_panel")
    assert data_area is not None
    assert data_area.scientific_data
    assert data_area.scrollable

    # Check spectral panel area if enabled
    if detail_component.config.show_spectral_data:
        spectral_area = detail_component.layout._active_areas.get("spectral_panel")
        assert spectral_area is not None
        assert spectral_area.scientific_data


@pytest.mark.asyncio
async def test_scientific_data_processing(detail_component: EnhancedResponsiveDetail, test_item: ScientificListItem):
    """Test scientific data processing."""
    await detail_component.set_item(test_item)

    # Check structure data processing
    if test_item.structure_data:
        assert "structure_ready" in detail_component._active_optimizations
        if detail_component.config.molecular_detail == "high":
            assert "high_quality_structure" in detail_component._active_optimizations

    # Check quantum data processing
    if test_item.quantum_data:
        assert "quantum_ready" in detail_component._active_optimizations
        if detail_component.config.visualization_quality == "high":
            assert "high_quality_quantum" in detail_component._active_optimizations

    # Check spectral data processing
    if test_item.spectral_data:
        assert "spectral_ready" in detail_component._active_optimizations
        if detail_component.config.visualization_quality == "high":
            assert "high_quality_spectral" in detail_component._active_optimizations

    # Check crystal data processing
    if test_item.crystal_data:
        assert "crystal_ready" in detail_component._active_optimizations
        if detail_component.config.visualization_quality == "high":
            assert "high_quality_crystal" in detail_component._active_optimizations
