"""Tests for data visualization functionality."""

import pytest
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from bs4 import BeautifulSoup

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..visualization import (
    DataVisualizer,
    StructureVisualizer,
    PlotGenerator,
    BindingVisualizer,
    ActivityVisualizer,
    SafetyVisualizer,
    VisualizationConfig,
    PlotType,
    VisualizationResult,
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
    caffeine.effect_profile = {
        "stimulation": (0.8, 0.9),
        "focus": (0.7, 0.8),
    }
    caffeine.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
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
        affinity=0.07,
        confidence=0.90,
        activity="inhibitor"
    )
    amphetamine.effect_profile = {
        "stimulation": (0.9, 0.95),
        "euphoria": (0.8, 0.9),
    }
    amphetamine.safety_alerts = {
        "addiction": RiskLevel.HIGH,
        "cardiovascular": RiskLevel.HIGH,
    }
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def data_visualizer():
    """Create test data visualizer fixture."""
    return DataVisualizer()


@pytest.fixture
def structure_visualizer():
    """Create test structure visualizer fixture."""
    return StructureVisualizer()


@pytest.fixture
def plot_generator():
    """Create test plot generator fixture."""
    return PlotGenerator()


class TestDataVisualizer:
    """Tests for DataVisualizer class."""

    def test_initialization(self, data_visualizer):
        """Test initialization of DataVisualizer."""
        assert isinstance(data_visualizer.binding_visualizer, BindingVisualizer)
        assert isinstance(data_visualizer.activity_visualizer, ActivityVisualizer)
        assert isinstance(data_visualizer.safety_visualizer, SafetyVisualizer)
        assert isinstance(data_visualizer.structure_visualizer, StructureVisualizer)
        assert isinstance(data_visualizer.plot_generator, PlotGenerator)
        assert data_visualizer.stats == {}

    def test_generate_compound_view(self, data_visualizer, test_compounds):
        """Test generating compound detail view."""
        # Generate view
        result = data_visualizer.generate_compound_view(test_compounds[0])
        
        # Check result
        assert isinstance(result, VisualizationResult)
        assert result.success
        assert "structure_svg" in result.data
        assert "plots" in result.data
        assert len(result.data["plots"]) > 0
        
        # Check HTML content
        html = result.data["html"]
        soup = BeautifulSoup(html, "html.parser")
        assert soup.find(id="structure-viewer") is not None
        assert soup.find(id="binding-plot") is not None
        assert soup.find(id="safety-plot") is not None

    def test_binding_heatmap(self, data_visualizer, test_compounds, tmp_path):
        """Test binding profile heatmap generation."""
        # Generate static heatmap
        output_file = tmp_path / "binding_heatmap.png"
        fig = data_visualizer.plot_binding_heatmap(
            compounds=test_compounds,
            output_file=output_file
        )
        
        # Check static plot
        assert isinstance(fig, plt.Figure)
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        
        # Check plot data
        ax = fig.axes[0]
        assert len(ax.get_xticklabels()) == len(test_compounds)
        assert len(ax.get_yticklabels()) >= 3  # A2A, DAT, NET
        plt.close(fig)
        
        # Generate interactive heatmap
        result = data_visualizer.generate_binding_heatmap(test_compounds)
        
        # Check interactive plot
        assert isinstance(result, VisualizationResult)
        assert result.success
        assert isinstance(result.data["figure"], go.Figure)
        assert "A2A" in result.data["figure"].data[0].text

    def test_activity_visualization(self, data_visualizer, test_compounds, tmp_path):
        """Test activity visualization generation."""
        # Generate static plot
        output_file = tmp_path / "activity_plot.png"
        fig = data_visualizer.plot_activity_patterns(
            compounds=test_compounds,
            output_file=output_file
        )
        
        # Check static plot
        assert isinstance(fig, plt.Figure)
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        plt.close(fig)
        
        # Generate interactive plot
        result = data_visualizer.generate_activity_plot(test_compounds[0])
        
        # Check interactive plot
        assert isinstance(result, VisualizationResult)
        assert result.success
        assert isinstance(result.data["figure"], go.Figure)
        assert "stimulation" in result.data["figure"].data[0].text

    def test_safety_visualization(self, data_visualizer, test_compounds, tmp_path):
        """Test safety visualization generation."""
        # Generate static plot
        output_file = tmp_path / "safety_radar.png"
        fig = data_visualizer.plot_safety_radar(
            compounds=test_compounds,
            output_file=output_file
        )
        
        # Check static plot
        assert isinstance(fig, plt.Figure)
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        plt.close(fig)
        
        # Generate interactive plot
        result = data_visualizer.generate_safety_plot(test_compounds[0])
        
        # Check interactive plot
        assert isinstance(result, VisualizationResult)
        assert result.success
        assert isinstance(result.data["figure"], go.Figure)
        assert "anxiety" in result.data["figure"].data[0].text

    def test_structure_similarity_network(self, data_visualizer, test_compounds, tmp_path):
        """Test structure similarity network plot generation."""
        # Generate plot
        output_file = tmp_path / "similarity_network.png"
        fig = data_visualizer.plot_similarity_network(
            compounds=test_compounds,
            output_file=output_file,
            similarity_threshold=0.5
        )
        
        # Check plot
        assert isinstance(fig, plt.Figure)
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        
        # Check network elements
        ax = fig.axes[0]
        assert len(ax.collections) > 0  # Network nodes
        assert len(ax.lines) > 0  # Network edges
        plt.close(fig)

    def test_property_distribution(self, data_visualizer, test_compounds, tmp_path):
        """Test property distribution plot generation."""
        # Generate plot
        output_file = tmp_path / "property_dist.png"
        fig = data_visualizer.plot_property_distribution(
            compounds=test_compounds,
            property_name="affinity",
            output_file=output_file
        )
        
        # Check plot
        assert isinstance(fig, plt.Figure)
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        
        # Check distribution elements
        ax = fig.axes[0]
        assert len(ax.patches) > 0  # Histogram bars
        plt.close(fig)


class TestStructureVisualizer:
    """Tests for StructureVisualizer class."""

    def test_initialization(self, structure_visualizer):
        """Test initialization of StructureVisualizer."""
        assert structure_visualizer.stats == {}

    def test_generate_2d_structure(self, structure_visualizer, test_compounds):
        """Test generating 2D structure visualization."""
        # Generate structure
        result = structure_visualizer.generate_2d_structure(test_compounds[0])
        
        # Check result
        assert isinstance(result, VisualizationResult)
        assert result.success
        assert "svg" in result.data
        assert "<svg" in result.data["svg"]
        assert result.stats["structure_type"] == "2D"

    def test_generate_3d_structure(self, structure_visualizer, test_compounds):
        """Test generating 3D structure visualization."""
        # Generate structure
        result = structure_visualizer.generate_3d_structure(test_compounds[0])
        
        # Check result
        assert isinstance(result, VisualizationResult)
        assert result.success
        assert "html" in result.data
        assert "3Dmol.js" in result.data["html"]
        assert result.stats["structure_type"] == "3D"

    def test_highlight_substructure(self, structure_visualizer, test_compounds):
        """Test substructure highlighting."""
        # Generate highlighted structure
        result = structure_visualizer.highlight_substructure(
            test_compounds[0],
            substructure="C1=NC2=C1C(=O)N"
        )
        
        # Check result
        assert isinstance(result, VisualizationResult)
        assert result.success
        assert "svg" in result.data
        assert "highlight" in result.data["svg"]
        assert result.stats["highlighted_atoms"] > 0


class TestPlotGenerator:
    """Tests for PlotGenerator class."""

    def test_initialization(self, plot_generator):
        """Test initialization of PlotGenerator."""
        assert plot_generator.stats == {}

    def test_plot_customization(self, plot_generator, test_compounds):
        """Test plot customization options."""
        # Configure plot
        config = VisualizationConfig(
            title="Custom Plot",
            width=800,
            height=600,
            colormap="viridis",
            show_legend=True,
            interactive=True,
        )
        
        # Generate customized plot
        result = plot_generator.generate_binding_plot(
            test_compounds[0],
            config=config
        )
        
        # Check customization
        assert isinstance(result, VisualizationResult)
        assert result.success
        figure = result.data["figure"]
        assert figure.layout.title.text == "Custom Plot"
        assert figure.layout.width == 800
        assert figure.layout.height == 600
        assert figure.layout.showlegend is True

    def test_batch_visualization(self, plot_generator, test_compounds, tmp_path):
        """Test batch visualization functionality."""
        # Generate multiple plots
        results = plot_generator.generate_visualization_batch(
            compounds=test_compounds,
            plot_types=[PlotType.BINDING, PlotType.EFFECT, PlotType.SAFETY]
        )
        
        # Check results
        assert len(results) == 3
        assert all(isinstance(r, VisualizationResult) for r in results)
        assert all(r.success for r in results)
        assert all(isinstance(r.data["figure"], go.Figure) for r in results)

    def test_error_handling(self, plot_generator):
        """Test error handling during plot generation."""
        # Create invalid compound
        invalid_compound = PsychoactiveCompound(
            name="Invalid",
            smiles="INVALID",
            cas_number="invalid",
        )
        
        # Attempt plot generation
        result = plot_generator.generate_binding_plot(invalid_compound)
        
        # Check error handling
        assert not result.success
        assert "Error" in str(result.error)
        assert result.stats["failed_plots"] == 1

    def test_visualization_performance(self, plot_generator, test_compounds):
        """Test visualization performance monitoring."""
        # Generate plot with performance monitoring
        result = plot_generator.generate_binding_plot(test_compounds[0])
        
        # Check performance stats
        assert "plot_time" in result.stats
        assert isinstance(result.stats["plot_time"], float)
        assert result.stats["plot_time"] >= 0
        assert "memory_usage" in result.stats


if __name__ == "__main__":
    pytest.main([__file__])
