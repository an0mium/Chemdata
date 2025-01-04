"""Tests for web visualization functionality."""

import pytest
from selenium.webdriver import Chrome
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.action_chains import ActionChains
from pathlib import Path

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..web import (
    WebApp,
    WebConfig,
)
from ..web_visualization import (
    WebVisualizer,
    PlotlyVisualizer,
    D3Visualizer,
    ChartGenerator,
)


@pytest.fixture
def test_compounds():
    """Create test compounds fixture."""
    compounds = []
    
    # Caffeine
    caffeine = PsychoactiveCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    caffeine.psychoactive_class = PsychoactiveClass.STIMULANT
    caffeine.add_receptor_binding(
        "A2A",
        affinity=0.8,
        confidence=0.95,
        activity="antagonist"
    )
    caffeine.add_receptor_binding(
        "A1",
        affinity=0.7,
        confidence=0.9,
        activity="antagonist"
    )
    caffeine.effect_profile = {
        "stimulation": (0.8, 0.9),
        "focus": (0.7, 0.8),
        "wakefulness": (0.9, 0.95),
    }
    caffeine.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
        "tachycardia": RiskLevel.LOW,
    }
    compounds.append(caffeine)
    
    # Amphetamine
    amphetamine = PsychoactiveCompound(
        name="Amphetamine",
        smiles="CC(N)CC1=CC=CC=C1",
        cas_number="300-62-9",
    )
    amphetamine.psychoactive_class = PsychoactiveClass.STIMULANT
    amphetamine.add_receptor_binding(
        "DAT",
        affinity=0.05,
        confidence=0.95,
        activity="inhibitor"
    )
    amphetamine.add_receptor_binding(
        "NET",
        affinity=0.1,
        confidence=0.9,
        activity="inhibitor"
    )
    amphetamine.effect_profile = {
        "stimulation": (0.9, 0.95),
        "euphoria": (0.8, 0.9),
        "focus": (0.85, 0.9),
    }
    amphetamine.safety_alerts = {
        "addiction": RiskLevel.HIGH,
        "cardiovascular": RiskLevel.HIGH,
        "neurotoxicity": RiskLevel.MODERATE,
    }
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def chrome_options():
    """Create Chrome options fixture."""
    options = Options()
    options.add_argument("--headless")  # Run in headless mode
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--window-size=1920,1080")  # Set window size for consistent testing
    return options


@pytest.fixture
def browser(chrome_options):
    """Create browser fixture."""
    driver = Chrome(options=chrome_options)
    yield driver
    driver.quit()


@pytest.fixture
def web_app():
    """Create test web application fixture."""
    config = WebConfig(
        host="localhost",
        port=8000,
        debug=True,
        static_dir=Path("static"),
        template_dir=Path("templates"),
        cache_dir=Path("cache")
    )
    app = WebApp(config=config)
    app.start()
    yield app
    app.stop()


@pytest.fixture
def web_visualizer():
    """Create test web visualizer fixture."""
    return WebVisualizer()


class TestBackendVisualization:
    """Tests for backend visualization functionality."""

    def test_initialization(self, web_visualizer):
        """Test initialization of WebVisualizer."""
        assert isinstance(web_visualizer.plotly_visualizer, PlotlyVisualizer)
        assert isinstance(web_visualizer.d3_visualizer, D3Visualizer)
        assert isinstance(web_visualizer.chart_generator, ChartGenerator)
        assert web_visualizer.stats == {}

    def test_binding_heatmap_plotly(self, test_compounds, web_visualizer):
        """Test Plotly binding heatmap generation."""
        # Generate plot
        plot_data = web_visualizer.generate_binding_heatmap(
            compounds=test_compounds,
            engine="plotly"
        )
        
        # Check plot data
        assert isinstance(plot_data, dict)
        assert "data" in plot_data
        assert "layout" in plot_data
        assert any("heatmap" in trace["type"].lower() for trace in plot_data["data"])
        assert "Receptor Binding Profile" in plot_data["layout"]["title"]["text"]

    def test_activity_plot_d3(self, test_compounds, web_visualizer):
        """Test D3.js activity plot generation."""
        # Generate plot
        plot_data = web_visualizer.generate_activity_plot(
            compounds=test_compounds,
            engine="d3"
        )
        
        # Check plot data
        assert isinstance(plot_data, dict)
        assert "nodes" in plot_data
        assert "links" in plot_data
        assert any(node["id"] == "A2A" for node in plot_data["nodes"])
        assert len(plot_data["links"]) > 0

    def test_safety_radar_chart(self, test_compounds, web_visualizer):
        """Test radar chart generation."""
        # Generate chart
        chart_data = web_visualizer.generate_safety_radar(
            compound=test_compounds[0]
        )
        
        # Check chart data
        assert isinstance(chart_data, dict)
        assert "labels" in chart_data
        assert "datasets" in chart_data
        assert "anxiety" in chart_data["labels"]
        assert "insomnia" in chart_data["labels"]
        assert len(chart_data["datasets"]) == 1

    def test_interactive_network(self, test_compounds, web_visualizer):
        """Test interactive network visualization."""
        # Generate network
        network_data = web_visualizer.generate_interaction_network(
            compounds=test_compounds,
            include_receptors=True,
            include_effects=True
        )
        
        # Check network data
        assert isinstance(network_data, dict)
        assert "nodes" in network_data
        assert "edges" in network_data
        assert any(n["type"] == "compound" for n in network_data["nodes"])
        assert any(n["type"] == "receptor" for n in network_data["nodes"])
        assert len(network_data["edges"]) > 0

    def test_plot_export(self, test_compounds, web_visualizer, tmp_path):
        """Test plot export functionality."""
        # Generate and export plot
        output_file = tmp_path / "web_plot.html"
        html = web_visualizer.export_interactive_plot(
            compounds=test_compounds,
            plot_type="binding",
            output_file=output_file
        )
        
        # Check output
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        assert isinstance(html, str)
        assert "plotly" in html.lower()
        assert test_compounds[0].name.lower() in html.lower()

    def test_visualization_customization(self, test_compounds, web_visualizer):
        """Test visualization customization options."""
        # Generate customized plot
        plot_data = web_visualizer.generate_binding_heatmap(
            compounds=test_compounds,
            engine="plotly",
            colorscale="Viridis",
            title="Custom Binding Profile",
            width=800,
            height=600
        )
        
        # Check customizations
        assert plot_data["layout"]["title"]["text"] == "Custom Binding Profile"
        assert plot_data["layout"]["width"] == 800
        assert plot_data["layout"]["height"] == 600
        assert plot_data["data"][0]["colorscale"] == "Viridis"

    def test_batch_visualization(self, test_compounds, web_visualizer, tmp_path):
        """Test batch visualization functionality."""
        # Generate multiple visualizations
        outputs = web_visualizer.generate_visualization_batch(
            compounds=test_compounds,
            output_dir=tmp_path,
            plot_types=["binding", "activity", "safety"]
        )
        
        # Check outputs
        assert len(outputs) == 3
        assert all(isinstance(data, dict) for data in outputs)
        assert len(list(tmp_path.glob("*.html"))) == 3

    def test_responsive_design(self, test_compounds, web_visualizer):
        """Test responsive design features."""
        # Generate responsive plot
        plot_data = web_visualizer.generate_binding_heatmap(
            compounds=test_compounds,
            engine="plotly",
            responsive=True
        )
        
        # Check responsive layout
        assert plot_data["layout"]["autosize"] is True
        assert "responsive" in plot_data["layout"]
        assert plot_data["layout"]["responsive"] is True

    def test_error_handling(self, test_compounds, web_visualizer):
        """Test error handling during visualization."""
        # Test invalid engine
        with pytest.raises(ValueError):
            web_visualizer.generate_binding_heatmap(
                compounds=test_compounds,
                engine="invalid_engine"
            )
        
        # Test invalid plot type
        with pytest.raises(ValueError):
            web_visualizer.generate_visualization_batch(
                compounds=test_compounds,
                plot_types=["invalid_type"]
            )

    def test_performance_monitoring(self, test_compounds, web_visualizer):
        """Test performance monitoring."""
        # Generate plot with performance monitoring
        web_visualizer.generate_binding_heatmap(
            compounds=test_compounds,
            engine="plotly"
        )
        
        # Check performance stats
        assert "visualization_time" in web_visualizer.stats
        assert isinstance(web_visualizer.stats["visualization_time"], float)
        assert web_visualizer.stats["visualization_time"] >= 0


class TestFrontendVisualization:
    """Tests for frontend visualization functionality."""

    def test_structure_viewer(self, browser, web_app):
        """Test structure viewer component."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Check structure viewer
        viewer = browser.find_element(By.ID, "structure-viewer")
        assert viewer.is_displayed()
        
        # Check SVG content
        svg = viewer.find_element(By.TAG_NAME, "svg")
        assert svg.is_displayed()
        assert len(svg.find_elements(By.TAG_NAME, "path")) > 0  # Has bonds
        assert len(svg.find_elements(By.TAG_NAME, "text")) > 0  # Has atoms

    def test_binding_plot(self, browser, web_app):
        """Test binding plot component."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Check plot
        plot = browser.find_element(By.ID, "binding-plot")
        assert plot.is_displayed()
        
        # Check heatmap cells
        cells = plot.find_elements(By.CLASS_NAME, "heatmap-cell")
        assert len(cells) > 0
        
        # Check legend
        legend = plot.find_element(By.CLASS_NAME, "plot-legend")
        assert legend.is_displayed()
        assert "Binding Affinity" in legend.text

    def test_safety_plot(self, browser, web_app):
        """Test safety plot component."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Check plot
        plot = browser.find_element(By.ID, "safety-plot")
        assert plot.is_displayed()
        
        # Check risk indicators
        indicators = plot.find_elements(By.CLASS_NAME, "risk-indicator")
        assert len(indicators) > 0
        
        # Check legend
        legend = plot.find_element(By.CLASS_NAME, "plot-legend")
        assert legend.is_displayed()
        assert "Risk Level" in legend.text

    def test_plot_interactions(self, browser, web_app):
        """Test plot interactions."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Get plot
        plot = browser.find_element(By.ID, "binding-plot")
        
        # Test hover interaction
        cell = plot.find_element(By.CLASS_NAME, "heatmap-cell")
        ActionChains(browser).move_to_element(cell).perform()
        
        # Check tooltip
        tooltip = WebDriverWait(browser, 10).until(
            EC.visibility_of_element_located((By.CLASS_NAME, "plot-tooltip"))
        )
        assert tooltip.is_displayed()
        assert "A2A" in tooltip.text
        assert "0.8" in tooltip.text

    def test_plot_customization(self, browser, web_app):
        """Test plot customization features."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Open customization panel
        customize_button = browser.find_element(By.ID, "customize-plots-button")
        customize_button.click()
        
        # Select color scheme
        scheme_select = WebDriverWait(browser, 10).until(
            EC.element_to_be_clickable((By.ID, "color-scheme-select"))
        )
        scheme_select.click()
        option = browser.find_element(
            By.XPATH,
            "//option[text()='Viridis']"
        )
        option.click()
        
        # Check plot colors
        plot = browser.find_element(By.ID, "binding-plot")
        cells = plot.find_elements(By.CLASS_NAME, "heatmap-cell")
        cell_colors = [cell.value_of_css_property("fill") for cell in cells]
        assert len(set(cell_colors)) > 1  # Multiple colors used

    def test_plot_export(self, browser, web_app):
        """Test plot export functionality."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Open export panel
        export_button = browser.find_element(By.ID, "export-plots-button")
        export_button.click()
        
        # Check export options
        panel = WebDriverWait(browser, 10).until(
            EC.presence_of_element_located((By.ID, "plot-export-panel"))
        )
        
        # Check format options
        format_select = panel.find_element(By.ID, "export-format-select")
        options = format_select.find_elements(By.TAG_NAME, "option")
        formats = [opt.text for opt in options]
        assert "PNG" in formats
        assert "SVG" in formats
        assert "PDF" in formats
        
        # Check resolution options
        resolution_select = panel.find_element(By.ID, "export-resolution-select")
        options = resolution_select.find_elements(By.TAG_NAME, "option")
        resolutions = [opt.text for opt in options]
        assert "72 DPI" in resolutions
        assert "300 DPI" in resolutions


if __name__ == "__main__":
    pytest.main([__file__])
