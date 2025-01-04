"""Tests for web interface functionality."""

import pytest
import pandas as pd
import json
from pathlib import Path
from bs4 import BeautifulSoup
from unittest.mock import patch, Mock

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..web import (
    WebInterface,
    WebServer,
    CompoundListView,
    CompoundDetailView,
    SearchInterface,
    FilterInterface,
    ExportInterface,
    APIEndpoints,
    SearchEngine,
    DataExporter,
    ComponentManager,
    VisualizationManager,
    WebConfig,
    ViewResult,
    WebResult,
    Dashboard,
    DataProcessor,
    DetailView,
    ListView,
    PlotManager,
    MLIntegration,
    AnalysisComponent,
    ExportComponent,
    InputComponent,
    VisualizationComponent,
    PatentSearchComponent,
    WebScrapingComponent,
    DataEnrichmentComponent
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
def web_config():
    """Create test web configuration fixture."""
    return WebConfig(
        host="localhost",
        port=8000,
        debug=True,
        static_dir=Path("static"),
        template_dir=Path("templates"),
        api_prefix="/api/v1",
        cors_origins=["http://localhost:3000"],
        cache_dir=Path("cache"),
        export_dir=Path("exports")
    )


@pytest.fixture
def web_interface():
    """Create test web interface fixture."""
    return WebInterface()


@pytest.fixture
def web_server(web_config, test_compounds):
    """Create test web server fixture."""
    server = WebServer(config=web_config)
    server.load_compounds(test_compounds)
    return server


@pytest.fixture
def dashboard(test_compounds):
    """Create test dashboard fixture."""
    return Dashboard(compounds=test_compounds)


@pytest.fixture
def list_view():
    """Create test compound list view fixture."""
    return CompoundListView()


@pytest.fixture
def detail_view():
    """Create test compound detail view fixture."""
    return CompoundDetailView()


@pytest.fixture
def search_interface():
    """Create test search interface fixture."""
    return SearchInterface()


@pytest.fixture
def filter_interface():
    """Create test filter interface fixture."""
    return FilterInterface()


@pytest.fixture
def export_interface():
    """Create test export interface fixture."""
    return ExportInterface()


@pytest.fixture
def visualization_manager():
    """Create test visualization manager fixture."""
    return VisualizationManager()


@pytest.fixture
def patent_search():
    """Create test patent search component fixture."""
    return PatentSearchComponent()


@pytest.fixture
def web_scraping():
    """Create test web scraping component fixture."""
    return WebScrapingComponent()


@pytest.fixture
def data_enrichment():
    """Create test data enrichment component fixture."""
    return DataEnrichmentComponent()


class TestWebInterface:
    """Tests for WebInterface class."""

    def test_initialization(self, web_interface):
        """Test initialization of WebInterface."""
        assert isinstance(web_interface.list_view, CompoundListView)
        assert isinstance(web_interface.detail_view, CompoundDetailView)
        assert isinstance(web_interface.search_interface, SearchInterface)
        assert isinstance(web_interface.filter_interface, FilterInterface)
        assert isinstance(web_interface.export_interface, ExportInterface)
        assert isinstance(web_interface.visualization_manager, VisualizationManager)
        assert web_interface.stats == {}

    def test_render_main_page(self, web_interface, test_compounds):
        """Test rendering main page."""
        # Render page
        result = web_interface.render_main_page(test_compounds)
        
        # Check result
        assert isinstance(result, ViewResult)
        assert result.success
        
        # Check HTML content
        html = result.data["html"]
        soup = BeautifulSoup(html, "html.parser")
        assert soup.find(id="compound-list") is not None
        assert soup.find(id="search-bar") is not None
        assert soup.find(id="filter-panel") is not None
        assert soup.find(id="export-panel") is not None
        
        # Check compound list
        compound_list = soup.find(id="compound-list")
        assert len(compound_list.find_all("tr")) == len(test_compounds) + 1  # +1 for header

    def test_handle_search(self, web_interface, test_compounds):
        """Test search functionality."""
        # Perform search
        result = web_interface.handle_search(
            compounds=test_compounds,
            query="caffeine",
            search_type="name"
        )
        
        # Check result
        assert isinstance(result, ViewResult)
        assert result.success
        assert len(result.data["compounds"]) == 1
        assert result.data["compounds"][0].name == "Caffeine"

    def test_handle_filtering(self, web_interface, test_compounds):
        """Test filtering functionality."""
        # Apply filters
        result = web_interface.handle_filtering(
            compounds=test_compounds,
            filters={
                "psychoactive_class": [PsychoactiveClass.STIMULANT],
                "min_binding_affinity": 0.7,
                "min_risk_level": RiskLevel.HIGH,
            }
        )
        
        # Check result
        assert isinstance(result, ViewResult)
        assert result.success
        assert len(result.data["compounds"]) == 1
        assert result.data["compounds"][0].name == "Amphetamine"

    def test_handle_export(self, web_interface, test_compounds, tmp_path):
        """Test export functionality."""
        # Configure export
        config = WebConfig(
            export_format="tsv",
            selected_columns=["name", "cas_number", "smiles"],
            include_metadata=True,
        )
        
        # Perform export
        output_file = tmp_path / "export.tsv"
        result = web_interface.handle_export(
            compounds=test_compounds,
            output_file=output_file,
            config=config
        )
        
        # Check result
        assert isinstance(result, ViewResult)
        assert result.success
        assert output_file.exists()
        
        # Verify export data
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == len(test_compounds)
        assert all(col in df.columns for col in config.selected_columns)

    def test_handle_visualization(self, web_interface, test_compounds):
        """Test visualization functionality."""
        # Configure visualization
        config = WebConfig(
            plot_type="2d_structure",
            width=400,
            height=400,
            interactive=True
        )
        
        # Generate visualization
        result = web_interface.handle_visualization(
            compound=test_compounds[0],
            config=config
        )
        
        # Check result
        assert isinstance(result, ViewResult)
        assert result.success
        assert "svg" in result.data
        assert result.data["width"] == config.width
        assert result.data["height"] == config.height


class TestWebServer:
    """Tests for WebServer class."""

    def test_initialization(self, web_server, web_config):
        """Test initialization of WebServer."""
        assert isinstance(web_server.api, APIEndpoints)
        assert isinstance(web_server.search, SearchEngine)
        assert isinstance(web_server.exporter, DataExporter)
        assert isinstance(web_server.components, ComponentManager)
        assert isinstance(web_server.visualizations, VisualizationManager)
        assert web_server.config == web_config
        assert web_server.stats == {}

    def test_load_compounds(self, web_server, test_compounds):
        """Test compound loading."""
        # Load compounds
        result = web_server.load_compounds(test_compounds)
        
        # Check result
        assert isinstance(result, WebResult)
        assert result.success
        assert len(result.compounds) == 2
        assert all(isinstance(c, PsychoactiveCompound) for c in result.compounds)
        assert "data_loading" in web_server.stats

    def test_api_endpoints(self, web_server, test_compounds):
        """Test API endpoints."""
        # Test compound list endpoint
        response = web_server.api.get_compounds()
        assert response.status_code == 200
        data = json.loads(response.data)
        assert len(data["compounds"]) == 2
        assert data["compounds"][0]["name"] == "Caffeine"
        
        # Test compound detail endpoint
        response = web_server.api.get_compound("58-08-2")  # Caffeine CAS
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data["compound"]["name"] == "Caffeine"
        assert data["compound"]["psychoactive_class"] == "STIMULANT"
        
        # Test search endpoint
        response = web_server.api.search_compounds(query="stimulant")
        assert response.status_code == 200
        data = json.loads(response.data)
        assert len(data["results"]) == 2  # Both are stimulants
        
        # Test export endpoint
        response = web_server.api.export_compounds(format="tsv")
        assert response.status_code == 200
        assert response.headers["Content-Type"] == "text/tab-separated-values"

    def test_search_functionality(self, web_server, test_compounds):
        """Test search functionality."""
        # Test text search
        results = web_server.search.text_search("caffeine")
        assert len(results) == 1
        assert results[0].name == "Caffeine"
        
        # Test structure search
        results = web_server.search.structure_search(
            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            threshold=0.7
        )
        assert len(results) == 1
        assert results[0].name == "Caffeine"
        
        # Test property filter
        results = web_server.search.filter_compounds(
            filters={
                "psychoactive_class": "STIMULANT",
                "risk_level": "HIGH"
            }
        )
        assert len(results) == 1
        assert results[0].name == "Amphetamine"

    def test_data_export(self, web_server, test_compounds, tmp_path):
        """Test data export functionality."""
        # Configure export
        export_file = tmp_path / "export.tsv"
        columns = ["name", "smiles", "cas_number", "psychoactive_class"]
        
        # Export data
        result = web_server.exporter.export_compounds(
            compounds=test_compounds,
            output_file=export_file,
            columns=columns,
            format="tsv"
        )
        
        # Check result
        assert result.success
        assert export_file.exists()
        content = export_file.read_text()
        assert all(col in content for col in columns)
        assert "Caffeine" in content
        assert "Amphetamine" in content

    def test_component_rendering(self, web_server, test_compounds):
        """Test component rendering."""
        # Test compound list component
        html = web_server.components.render_compound_list(test_compounds)
        assert isinstance(html, str)
        assert "Caffeine" in html
        assert "Amphetamine" in html
        
        # Test compound detail component
        html = web_server.components.render_compound_detail(test_compounds[0])
        assert isinstance(html, str)
        assert "Caffeine" in html
        assert "A2A" in html
        assert "antagonist" in html
        
        # Test structure viewer component
        html = web_server.components.render_structure_viewer(test_compounds[0])
        assert isinstance(html, str)
        assert test_compounds[0].smiles in html

    def test_visualization_generation(self, web_server, test_compounds):
        """Test visualization generation."""
        # Test receptor plot
        plot = web_server.visualizations.create_receptor_plot(test_compounds[0])
        assert isinstance(plot, dict)
        assert "data" in plot
        assert "layout" in plot
        
        # Test effect profile plot
        plot = web_server.visualizations.create_effect_plot(test_compounds[0])
        assert isinstance(plot, dict)
        assert "data" in plot
        assert "layout" in plot
        
        # Test safety plot
        plot = web_server.visualizations.create_safety_plot(test_compounds[0])
        assert isinstance(plot, dict)
        assert "data" in plot
        assert "layout" in plot

    def test_error_handling(self, web_server):
        """Test error handling."""
        # Test invalid compound request
        response = web_server.api.get_compound("invalid-cas")
        assert response.status_code == 404
        data = json.loads(response.data)
        assert "error" in data
        
        # Test invalid search
        response = web_server.api.search_compounds(query="")
        assert response.status_code == 400
        data = json.loads(response.data)
        assert "error" in data
        
        # Test invalid export
        response = web_server.api.export_compounds(format="invalid")
        assert response.status_code == 400
        data = json.loads(response.data)
        assert "error" in data

    def test_caching(self, web_server, test_compounds):
        """Test response caching."""
        # Make initial request
        response1 = web_server.api.get_compounds()
        assert response1.status_code == 200
        assert "X-Cache" not in response1.headers
        
        # Make same request again
        response2 = web_server.api.get_compounds()
        assert response2.status_code == 200
        assert response2.headers["X-Cache"] == "HIT"
        
        # Check cache stats
        assert "cache_hits" in web_server.stats
        assert "cache_misses" in web_server.stats

    def test_cors_handling(self, web_server):
        """Test CORS handling."""
        # Test allowed origin
        headers = {"Origin": "http://localhost:3000"}
        response = web_server.api.get_compounds(headers=headers)
        assert response.status_code == 200
        assert response.headers["Access-Control-Allow-Origin"] == headers["Origin"]
        
        # Test disallowed origin
        headers = {"Origin": "http://evil.com"}
        response = web_server.api.get_compounds(headers=headers)
        assert response.status_code == 200
        assert "Access-Control-Allow-Origin" not in response.headers


class TestDashboard:
    """Tests for Dashboard class."""

    def test_initialization(self, dashboard, test_compounds):
        """Test initialization of Dashboard."""
        assert isinstance(dashboard.data_processor, DataProcessor)
        assert isinstance(dashboard.detail_view, DetailView)
        assert isinstance(dashboard.list_view, ListView)
        assert isinstance(dashboard.plot_manager, PlotManager)
        assert isinstance(dashboard.ml_integration, MLIntegration)
        assert len(dashboard.compounds) == len(test_compounds)

    def test_update_compounds(self, dashboard, test_compounds):
        """Test updating compounds data."""
        # Add new compound
        new_compound = PsychoactiveCompound(
            name="New Compound",
            smiles="CC1=CC=CC=C1",
            cas_number="123-45-6",
        )
        test_compounds.append(new_compound)
        
        # Update dashboard
        dashboard.update_compounds(test_compounds)
        
        # Check update
        assert len(dashboard.compounds) == len(test_compounds)
        assert dashboard.compounds[-1].name == "New Compound"
        assert "update_time" in dashboard.stats

    def test_filter_compounds(self, dashboard):
        """Test compound filtering."""
        # Filter by class
        filtered = dashboard.filter_compounds(
            psychoactive_class=PsychoactiveClass.STIMULANT
        )
        assert len(filtered) == 2
        assert all(c.psychoactive_class == PsychoactiveClass.STIMULANT for c in filtered)
        
        # Filter by receptor
        filtered = dashboard.filter_compounds(receptor="A2A")
        assert len(filtered) == 1
        assert filtered[0].name == "Caffeine"
        
        # Filter by risk level
        filtered = dashboard.filter_compounds(
            min_risk_level=RiskLevel.HIGH
        )
        assert len(filtered) == 1
        assert filtered[0].name == "Amphetamine"

    def test_search_compounds(self, dashboard):
        """Test compound search."""
        # Search by name
        results = dashboard.search_compounds(query="caffeine")
        assert len(results) == 1
        assert results[0].name == "Caffeine"
        
        # Search by SMILES fragment
        results = dashboard.search_compounds(query="CC1=CC=CC=C1")
        assert len(results) == 1
        assert results[0].name == "Amphetamine"
        
        # Search by effect
        results = dashboard.search_compounds(query="stimulation")
        assert len(results) == 2


class TestDetailView:
    """Tests for DetailView class."""

    def test_compound_details(self, dashboard):
        """Test compound detail view."""
        # Get compound details
        details = dashboard.detail_view.get_compound_details(
            compound=dashboard.compounds[0]
        )
        
        # Check details
        assert details["name"] == "Caffeine"
        assert details["smiles"] == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        assert details["cas_number"] == "58-08-2"
        assert "structure_svg" in details
        assert "receptor_plot" in details
        assert "effect_plot" in details
        assert "safety_plot" in details

    def test_receptor_profile(self, dashboard):
        """Test receptor profile visualization."""
        # Get receptor profile
        profile = dashboard.detail_view.get_receptor_profile(
            compound=dashboard.compounds[0]
        )
        
        # Check profile
        assert "A2A" in profile["receptors"]
        assert profile["affinities"]["A2A"] == 0.8
        assert profile["activities"]["A2A"] == "antagonist"
        assert "plot" in profile

    def test_effect_profile(self, dashboard):
        """Test effect profile visualization."""
        # Get effect profile
        profile = dashboard.detail_view.get_effect_profile(
            compound=dashboard.compounds[0]
        )
        
        # Check profile
        assert "stimulation" in profile["effects"]
        assert profile["scores"]["stimulation"] == (0.8, 0.9)
        assert "plot" in profile

    def test_safety_profile(self, dashboard):
        """Test safety profile visualization."""
        # Get safety profile
        profile = dashboard.detail_view.get_safety_profile(
            compound=dashboard.compounds[0]
        )
        
        # Check profile
        assert "anxiety" in profile["alerts"]
        assert profile["risk_levels"]["anxiety"] == RiskLevel.MODERATE
        assert "plot" in profile


class TestListView:
    """Tests for ListView class."""

    def test_compound_list(self, dashboard):
        """Test compound list view."""
        # Get compound list
        compounds = dashboard.list_view.get_compound_list()
        
        # Check list
        assert len(compounds) == 2
        assert all("name" in c for c in compounds)
        assert all("smiles" in c for c in compounds)
        assert all("cas_number" in c for c in compounds)
        assert all("structure_svg" in c for c in compounds)

    def test_sorting(self, dashboard):
        """Test compound sorting."""
        # Sort by name
        compounds = dashboard.list_view.sort_compounds(sort_by="name")
        assert compounds[0]["name"] == "Amphetamine"
        assert compounds[1]["name"] == "Caffeine"
        
        # Sort by binding affinity
        compounds = dashboard.list_view.sort_compounds(
            sort_by="binding",
            receptor="DAT"
        )
        assert compounds[0]["name"] == "Amphetamine"  # Stronger DAT binding

    def test_pagination(self, dashboard):
        """Test compound pagination."""
        # Get first page
        page = dashboard.list_view.get_page(page=1, per_page=1)
        assert len(page["compounds"]) == 1
        assert page["total_pages"] == 2
        assert page["current_page"] == 1
        
        # Get second page
        page = dashboard.list_view.get_page(page=2, per_page=1)
        assert len(page["compounds"]) == 1
        assert page["current_page"] == 2


class TestPlotManager:
    """Tests for PlotManager class."""

    def test_structure_plot(self, dashboard):
        """Test structure visualization."""
        # Get structure plot
        plot = dashboard.plot_manager.plot_structure(
            compound=dashboard.compounds[0]
        )
        
        # Check plot
        assert isinstance(plot, str)  # SVG string
        assert plot.startswith("<?xml")
        assert "svg" in plot

    def test_receptor_plot(self, dashboard):
        """Test receptor profile plot."""
        # Get receptor plot
        plot = dashboard.plot_manager.plot_receptor_profile(
            compound=dashboard.compounds[0]
        )
        
        # Check plot
        assert isinstance(plot, dict)
        assert "data" in plot
        assert "layout" in plot
        assert plot["data"][0]["type"] == "bar"

    def test_effect_plot(self, dashboard):
        """Test effect profile plot."""
        # Get effect plot
        plot = dashboard.plot_manager.plot_effect_profile(
            compound=dashboard.compounds[0]
        )
        
        # Check plot
        assert isinstance(plot, dict)
        assert "data" in plot
        assert "layout" in plot
        assert plot["data"][0]["type"] == "radar"

    def test_safety_plot(self, dashboard):
        """Test safety profile plot."""
        # Get safety plot
        plot = dashboard.plot_manager.plot_safety_profile(
            compound=dashboard.compounds[0]
        )
        
        # Check plot
        assert isinstance(plot, dict)
        assert "data" in plot
        assert "layout" in plot
        assert plot["data"][0]["type"] == "heatmap"


class TestMLIntegration:
    """Tests for MLIntegration class."""

    def test_binding_predictions(self, dashboard):
        """Test binding affinity predictions."""
        # Get predictions
        predictions = dashboard.ml_integration.predict_binding(
            compound=dashboard.compounds[0]
        )
        
        # Check predictions
        assert isinstance(predictions, dict)
        assert all(isinstance(v, tuple) for v in predictions.values())
        assert all(len(v) == 2 for v in predictions.values())  # (affinity, confidence)

    def test_activity_predictions(self, dashboard):
        """Test activity predictions."""
        # Get predictions
        predictions = dashboard.ml_integration.predict_activity(
            compound=dashboard.compounds[0]
        )
        
        # Check predictions
        assert isinstance(predictions, dict)
        assert "class" in predictions
        assert "confidence" in predictions
        assert "effects" in predictions

    def test_safety_predictions(self, dashboard):
        """Test safety predictions."""
        # Get predictions
        predictions = dashboard.ml_integration.predict_safety(
            compound=dashboard.compounds[0]
        )
        
        # Check predictions
        assert isinstance(predictions, dict)
        assert "risks" in predictions
        assert "confidence" in predictions
        assert "alerts" in predictions


class TestComponents:
    """Tests for dashboard components."""

    def test_analysis_component(self, dashboard):
        """Test analysis component."""
        component = dashboard.get_component("analysis")
        assert isinstance(component, AnalysisComponent)
        
        # Test analysis
        result = component.analyze_compound(dashboard.compounds[0])
        assert "receptor_analysis" in result
        assert "effect_analysis" in result
        assert "safety_analysis" in result

    def test_export_component(self, dashboard):
        """Test export component."""
        component = dashboard.get_component("export")
        assert isinstance(component, ExportComponent)
        
        # Test export
        data = component.export_compounds(
            compounds=dashboard.compounds,
            format="tsv",
            columns=["name", "smiles", "cas_number"]
        )
        assert isinstance(data, str)
        assert len(data.split("\n")) == 3  # Header + 2 compounds

    def test_input_component(self, dashboard):
        """Test input component."""
        component = dashboard.get_component("input")
        assert isinstance(component, InputComponent)
        
        # Test input validation
        valid = component.validate_input({
            "name": "Test Compound",
            "smiles": "CC1=CC=CC=C1",
            "cas_number": "123-45-6"
        })
        assert valid
        
        # Test structure parsing
        structure = component.parse_structure("CC1=CC=CC=C1")
        assert structure is not None

    def test_visualization_component(self, dashboard):
        """Test visualization component."""
        component = dashboard.get_component("visualization")
        assert isinstance(component, VisualizationComponent)
        
        # Test visualization
        viz = component.visualize_compound(dashboard.compounds[0])
        assert "structure" in viz
        assert "receptor_profile" in viz
        assert "effect_profile" in viz
        assert "safety_profile" in viz


class TestPatentSearch:
    """Tests for patent search functionality."""
    
    def test_search_by_structure(self, dashboard, patent_search):
        """Test structure-based patent search."""
        results = patent_search.search_by_structure(
            compound=dashboard.compounds[0],
            similarity_threshold=0.7
        )
        
        assert isinstance(results, list)
        assert len(results) > 0
        assert all("patent_number" in r for r in results)
        assert all("similarity_score" in r for r in results)
        assert all("title" in r for r in results)

    def test_search_by_activity(self, dashboard, patent_search):
        """Test activity-based patent search."""
        results = patent_search.search_by_activity(
            compound=dashboard.compounds[0],
            activity_type="antagonist",
            target="A2A"
        )
        
        assert isinstance(results, list)
        assert len(results) > 0
        assert all("relevance_score" in r for r in results)

    def test_extract_compounds(self, dashboard, patent_search):
        """Test compound extraction from patents."""
        compounds = patent_search.extract_compounds(
            patent_number="US1234567"
        )
        
        assert isinstance(compounds, list)
        assert all(isinstance(c, PsychoactiveCompound) for c in compounds)


class TestWebScraping:
    """Tests for web scraping functionality."""
    
    def test_scrape_community_data(self, dashboard, web_scraping):
        """Test community data scraping."""
        data = web_scraping.scrape_community_data(
            compound=dashboard.compounds[0]
        )
        
        assert isinstance(data, dict)
        assert "experiences" in data
        assert "effects" in data
        assert "dosage" in data

    def test_scrape_literature(self, dashboard, web_scraping):
        """Test literature scraping."""
        data = web_scraping.scrape_literature(
            compound=dashboard.compounds[0]
        )
        
        assert isinstance(data, dict)
        assert "papers" in data
        assert "citations" in data

    def test_monitor_social_media(self, dashboard, web_scraping):
        """Test social media monitoring."""
        data = web_scraping.monitor_social_media(
            compound=dashboard.compounds[0],
            platforms=["reddit", "twitter"]
        )
        
        assert isinstance(data, dict)
        assert all(p in data for p in ["reddit", "twitter"])


class TestDataEnrichment:
    """Tests for data enrichment functionality."""
    
    def test_enrich_compound_data(self, dashboard, data_enrichment):
        """Test compound data enrichment."""
        enriched = data_enrichment.enrich_compound_data(
            compound=dashboard.compounds[0]
        )
        
        assert isinstance(enriched, PsychoactiveCompound)
        assert hasattr(enriched, "community_data")
        assert hasattr(enriched, "literature_data")
        assert hasattr(enriched, "patent_data")

    def test_validate_enriched_data(self, dashboard, data_enrichment):
        """Test enriched data validation."""
        validation = data_enrichment.validate_enriched_data(
            compound=dashboard.compounds[0]
        )
        
        assert isinstance(validation, dict)
        assert "is_valid" in validation
        assert "errors" in validation
        assert "warnings" in validation

    def test_merge_data_sources(self, dashboard, data_enrichment):
        """Test data source merging."""
        merged = data_enrichment.merge_data_sources(
            compound=dashboard.compounds[0],
            sources=["community", "literature", "patents"]
        )
        
        assert isinstance(merged, dict)
        assert all(s in merged for s in ["community", "literature", "patents"])


class TestPerformance:
    """Tests for performance and caching."""
    
    def test_batch_processing(self, dashboard):
        """Test batch processing performance."""
        with patch.object(dashboard.ml_integration, 'process_batch') as mock:
            dashboard.process_compounds_batch(
                compounds=dashboard.compounds,
                batch_size=2
            )
            
            assert mock.call_count == 1
            assert "batch_processing_time" in dashboard.stats

    def test_cache_invalidation(self, web_server):
        """Test cache invalidation."""
        response1 = web_server.api.get_compounds()
        assert response1.status_code == 200
        
        web_server.invalidate_cache()
        
        response2 = web_server.api.get_compounds()
        assert response2.status_code == 200
        assert "X-Cache" not in response2.headers

    def test_concurrent_requests(self, web_server):
        """Test concurrent request handling."""
        import concurrent.futures
        
        def make_request():
            return web_server.api.get_compounds()
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            futures = [executor.submit(make_request) for _ in range(4)]
            responses = [f.result() for f in futures]
        
        assert all(r.status_code == 200 for r in responses)
        assert "concurrent_requests" in web_server.stats


if __name__ == "__main__":
    pytest.main([__file__])
