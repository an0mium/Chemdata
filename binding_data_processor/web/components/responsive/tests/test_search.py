"""Tests for enhanced responsive search component."""

import pytest
from unittest.mock import Mock, patch

from ..search import (
    SearchMode,
    SearchConfig,
    EnhancedResponsiveSearch,
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
def test_results(scientific_metadata: ScientificMetadata) -> list[ScientificListItem]:
    """Create test search results."""
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
            binding_data={"5HT2A": 8.2, "D2": 7.1},
            safety_data={"ld50": 1000},
            tags=["validated", "high-quality"],
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
            binding_data={"5HT2A": 7.8, "D2": 6.9},
            safety_data={"ld50": 800},
            tags=["preliminary"],
        ),
    ]


@pytest.fixture
def search_config() -> SearchConfig:
    """Create test search configuration."""
    return SearchConfig(
        search_mode=SearchMode.TEXT,
        enable_structure_search=True,
        enable_similarity_search=True,
        enable_property_search=True,
        enable_advanced_search=True,
        enable_quantum_search=True,
        enable_binding_search=True,
        enable_safety_search=True,
        enable_literature_search=True,
        scientific_precision=True,
        molecular_detail="high",
        quantum_enabled=True,
        visualization_quality="high",
        data_precision_level="high",
        enable_structure_drawing=True,
        enable_accessibility=True,
        high_performance=True,
    )


@pytest.fixture
def mock_search_callback() -> Mock:
    """Create mock search callback."""
    return Mock()


@pytest.fixture
def mock_select_callback() -> Mock:
    """Create mock select callback."""
    return Mock()


@pytest.fixture
async def search_component(
    search_config: SearchConfig,
    mock_search_callback: Mock,
    mock_select_callback: Mock,
) -> EnhancedResponsiveSearch:
    """Create test search component."""
    component = EnhancedResponsiveSearch(
        config=search_config,
        on_search=mock_search_callback,
        on_select=mock_select_callback,
    )
    await component.setup()
    return component


@pytest.mark.asyncio
async def test_initialization(search_component: EnhancedResponsiveSearch):
    """Test component initialization."""
    assert search_component.config is not None
    assert search_component.layout is not None
    assert len(search_component._active_optimizations) > 0
    assert search_component._current_query == {}
    assert search_component._search_history == []
    assert search_component._suggestions == []
    assert search_component._results == []


@pytest.mark.asyncio
async def test_search(
    search_component: EnhancedResponsiveSearch,
    mock_search_callback: Mock,
):
    """Test search functionality."""
    query = {
        "text": "test compound",
        "mode": SearchMode.TEXT,
        "filters": {"tags": ["validated"]},
    }

    await search_component.search(query)

    assert search_component._current_query == query
    assert len(search_component._search_history) == 1
    assert search_component._search_history[0] == query
    mock_search_callback.assert_called_once_with(query)


@pytest.mark.asyncio
async def test_select_result(
    search_component: EnhancedResponsiveSearch,
    mock_select_callback: Mock,
):
    """Test result selection."""
    result_id = "CMPD-001"
    await search_component.select_result(result_id)
    mock_select_callback.assert_called_once_with(result_id)


@pytest.mark.asyncio
async def test_update_results(
    search_component: EnhancedResponsiveSearch,
    test_results: list[ScientificListItem],
):
    """Test updating search results."""
    await search_component.update_results(test_results)

    assert search_component._results == test_results

    # Check optimization flags for result data
    assert "structure_ready" in search_component._active_optimizations
    assert "quantum_ready" in search_component._active_optimizations
    assert "spectral_ready" in search_component._active_optimizations
    assert "binding_ready" in search_component._active_optimizations
    assert "safety_ready" in search_component._active_optimizations


@pytest.mark.asyncio
async def test_responsive_behavior(search_component: EnhancedResponsiveSearch):
    """Test responsive layout behavior."""
    viewport_sizes = [(320, 480), (768, 1024), (1920, 1080)]

    for width, height in viewport_sizes:
        await search_component.handle_resize(width)

        if width < 576:  # xs
            assert search_component.config.search_columns == search_component.config.grid_columns
            assert search_component.config.search_width == width - 32
            assert search_component.config.results_width == width - 32
        elif width < 768:  # sm
            assert search_component.config.search_columns == search_component.config.grid_columns
            assert search_component.config.search_width == width - 48
            assert search_component.config.results_width == width - 48
        else:  # md+
            assert search_component.config.search_columns == 4
            assert search_component.config.search_width == (width - 72) // 3
            assert search_component.config.results_width == (width - 72) * 2 // 3


@pytest.mark.asyncio
async def test_optimizations(search_component: EnhancedResponsiveSearch):
    """Test optimization flags."""
    # Check performance optimizations
    if search_component.config.high_performance:
        assert "gpu_acceleration" in search_component._active_optimizations
        assert "content_visibility" in search_component._active_optimizations
        assert "will_change_transform" in search_component._active_optimizations

    # Check scientific optimizations
    if search_component.config.scientific_precision:
        assert "high_precision" in search_component._active_optimizations
        assert "tabular_nums" in search_component._active_optimizations
        assert "monospace_font" in search_component._active_optimizations

    # Check search optimizations
    if search_component.config.enable_search_debounce:
        assert "search_debounce" in search_component._active_optimizations
    if search_component.config.enable_results_virtualization:
        assert "results_virtualization" in search_component._active_optimizations

    # Check visualization optimizations
    if search_component.config.enable_structure_drawing:
        assert "structure_editor" in search_component._active_optimizations
        if search_component.config.visualization_quality == "high":
            assert "high_quality_structure" in search_component._active_optimizations

    # Check interaction optimizations
    if search_component.config.enable_suggestions:
        assert "search_suggestions" in search_component._active_optimizations
    if search_component.config.enable_autocomplete:
        assert "search_autocomplete" in search_component._active_optimizations
    if search_component.config.enable_touch_gestures:
        assert "touch_optimized" in search_component._active_optimizations

    # Check accessibility optimizations
    if search_component.config.enable_accessibility:
        assert "screen_reader_support" in search_component._active_optimizations
        assert "keyboard_navigation" in search_component._active_optimizations
        assert "high_contrast" in search_component._active_optimizations


@pytest.mark.asyncio
async def test_layout_areas(search_component: EnhancedResponsiveSearch):
    """Test layout area configuration."""
    # Check search panel area
    search_area = search_component.layout._active_areas.get("search_panel")
    assert search_area is not None
    assert search_area.scrollable

    # Check structure editor area if enabled
    if search_component.config.enable_structure_drawing:
        structure_area = search_component.layout._active_areas.get("structure_editor")
        assert structure_area is not None
        assert structure_area.molecular_view

    # Check results panel area
    results_area = search_component.layout._active_areas.get("results_panel")
    assert results_area is not None
    assert results_area.scientific_data
    assert results_area.scrollable


@pytest.mark.asyncio
async def test_scientific_data_processing(
    search_component: EnhancedResponsiveSearch,
    test_results: list[ScientificListItem],
):
    """Test scientific data processing."""
    await search_component.update_results(test_results)

    for result in test_results:
        # Check structure data processing
        if result.structure_data and search_component.config.enable_structure_drawing:
            assert "structure_ready" in search_component._active_optimizations
            if search_component.config.molecular_detail == "high":
                assert "high_quality_structure" in search_component._active_optimizations

        # Check quantum data processing
        if result.quantum_data and search_component.config.quantum_enabled:
            assert "quantum_ready" in search_component._active_optimizations
            if search_component.config.visualization_quality == "high":
                assert "high_quality_quantum" in search_component._active_optimizations

        # Check spectral data processing
        if result.spectral_data:
            assert "spectral_ready" in search_component._active_optimizations
            if search_component.config.visualization_quality == "high":
                assert "high_quality_spectral" in search_component._active_optimizations

        # Check binding data processing
        if result.binding_data and search_component.config.enable_binding_search:
            assert "binding_ready" in search_component._active_optimizations
            if search_component.config.data_precision_level == "high":
                assert "high_precision_binding" in search_component._active_optimizations

        # Check safety data processing
        if result.safety_data and search_component.config.enable_safety_search:
            assert "safety_ready" in search_component._active_optimizations
            if search_component.config.data_precision_level == "high":
                assert "high_precision_safety" in search_component._active_optimizations
