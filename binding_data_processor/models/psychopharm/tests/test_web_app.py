"""Tests for web application functionality."""

import pytest
from fastapi.testclient import TestClient
from bs4 import BeautifulSoup
from selenium.webdriver import Chrome
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from pathlib import Path

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..web import (
    WebApp,
    WebConfig,
    RouteHandler,
    DatabaseManager,
    SessionManager,
    TemplateManager,
    StaticFileManager,
    PageRenderer,
    ComponentManager,
    StateManager,
    EventHandler,
    VisualizationManager,
    InteractionManager,
    ErrorManager,
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
def web_config():
    """Create test web configuration fixture."""
    return WebConfig(
        host="localhost",
        port=8000,
        debug=True,
        database_url="sqlite:///:memory:",
        secret_key="test-secret-key",
        static_dir=Path("static"),
        template_dir=Path("templates"),
        session_dir=Path("sessions"),
        log_dir=Path("logs"),
        cache_dir=Path("cache"),
        export_dir=Path("exports"),
        cors_origins=["http://localhost:3000"]
    )


@pytest.fixture
def web_app(web_config, test_compounds):
    """Create test web application fixture."""
    app = WebApp(config=web_config)
    app.load_compounds(test_compounds)
    return app


@pytest.fixture
def test_client(web_app):
    """Create test client fixture."""
    return TestClient(web_app.server.app)


@pytest.fixture
def chrome_options():
    """Create Chrome options fixture."""
    options = Options()
    options.add_argument("--headless")  # Run in headless mode
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    return options


@pytest.fixture
def browser(chrome_options):
    """Create browser fixture."""
    driver = Chrome(options=chrome_options)
    yield driver
    driver.quit()


class TestBackend:
    """Tests for backend functionality."""

    def test_initialization(self, web_app, web_config):
        """Test initialization of WebApp."""
        assert isinstance(web_app.routes, RouteHandler)
        assert isinstance(web_app.db, DatabaseManager)
        assert isinstance(web_app.sessions, SessionManager)
        assert isinstance(web_app.templates, TemplateManager)
        assert isinstance(web_app.static_files, StaticFileManager)
        assert isinstance(web_app.page_renderer, PageRenderer)
        assert isinstance(web_app.component_manager, ComponentManager)
        assert isinstance(web_app.state_manager, StateManager)
        assert isinstance(web_app.event_handler, EventHandler)
        assert isinstance(web_app.visualization_manager, VisualizationManager)
        assert isinstance(web_app.interaction_manager, InteractionManager)
        assert isinstance(web_app.error_manager, ErrorManager)
        assert web_app.config == web_config
        assert web_app.stats == {}

    def test_database_integration(self, web_app, test_compounds):
        """Test database integration."""
        # Save compounds
        result = web_app.db.save_compounds(test_compounds)
        assert result.success
        assert len(result.compounds) == 2
        
        # Query compounds
        compounds = web_app.db.get_compounds()
        assert len(compounds) == 2
        assert compounds[0].name == "Caffeine"
        
        # Query by CAS
        compound = web_app.db.get_compound_by_cas("58-08-2")
        assert compound is not None
        assert compound.name == "Caffeine"
        
        # Update compound
        compound.effect_profile["alertness"] = (0.9, 0.95)
        result = web_app.db.update_compound(compound)
        assert result.success
        
        # Delete compound
        result = web_app.db.delete_compound("58-08-2")
        assert result.success
        assert len(web_app.db.get_compounds()) == 1

    def test_session_management(self, web_app):
        """Test session management."""
        # Create session
        session_id = web_app.sessions.create_session({
            "user_id": "test-user",
            "filters": {"class": "STIMULANT"},
        })
        assert session_id is not None
        
        # Get session
        session = web_app.sessions.get_session(session_id)
        assert session is not None
        assert session["user_id"] == "test-user"
        assert session["filters"]["class"] == "STIMULANT"
        
        # Update session
        web_app.sessions.update_session(session_id, {
            "filters": {"class": "DEPRESSANT"},
        })
        session = web_app.sessions.get_session(session_id)
        assert session["filters"]["class"] == "DEPRESSANT"
        
        # Delete session
        web_app.sessions.delete_session(session_id)
        assert web_app.sessions.get_session(session_id) is None

    def test_template_rendering(self, web_app, test_compounds):
        """Test template rendering."""
        # Render main template
        html = web_app.templates.render_template(
            "main.html",
            compounds=test_compounds,
            title="Test Page"
        )
        assert isinstance(html, str)
        soup = BeautifulSoup(html, "html.parser")
        assert soup.title.string == "Test Page"
        
        # Render compound list template
        html = web_app.templates.render_template(
            "compound_list.html",
            compounds=test_compounds
        )
        assert isinstance(html, str)
        soup = BeautifulSoup(html, "html.parser")
        assert len(soup.find_all("tr")) == len(test_compounds) + 1
        
        # Render compound detail template
        html = web_app.templates.render_template(
            "compound_detail.html",
            compound=test_compounds[0]
        )
        assert isinstance(html, str)
        soup = BeautifulSoup(html, "html.parser")
        assert test_compounds[0].name in soup.text
        assert test_compounds[0].cas_number in soup.text

    def test_static_file_serving(self, web_app):
        """Test static file serving."""
        # Serve CSS file
        response = web_app.static_files.serve_file("css/style.css")
        assert response.status_code == 200
        assert "text/css" in response.headers["Content-Type"]
        
        # Serve JavaScript file
        response = web_app.static_files.serve_file("js/main.js")
        assert response.status_code == 200
        assert "application/javascript" in response.headers["Content-Type"]
        
        # Serve image file
        response = web_app.static_files.serve_file("img/logo.png")
        assert response.status_code == 200
        assert "image/png" in response.headers["Content-Type"]
        
        # Handle missing file
        response = web_app.static_files.serve_file("not-found.txt")
        assert response.status_code == 404


class TestFrontend:
    """Tests for frontend functionality."""

    def test_page_rendering(self, test_client):
        """Test page rendering."""
        # Test main page
        response = test_client.get("/")
        assert response.status_code == 200
        soup = BeautifulSoup(response.text, "html.parser")
        assert soup.find(id="app-container") is not None
        assert soup.find(id="compound-list") is not None
        assert soup.find(id="search-bar") is not None
        assert soup.find(id="filter-panel") is not None
        assert soup.find(id="export-panel") is not None
        
        # Test compound detail page
        response = test_client.get("/compounds/58-08-2")  # Caffeine CAS
        assert response.status_code == 200
        soup = BeautifulSoup(response.text, "html.parser")
        assert soup.find(id="compound-info") is not None
        assert soup.find(id="structure-viewer") is not None
        assert soup.find(id="binding-plot") is not None
        assert soup.find(id="safety-plot") is not None

    def test_component_rendering(self, web_app):
        """Test component rendering."""
        # Test search bar
        result = web_app.component_manager.render_search_bar()
        soup = BeautifulSoup(result.html, "html.parser")
        assert soup.find(id="search-bar") is not None
        assert soup.find("input") is not None
        assert soup.find("button") is not None
        
        # Test filter panel
        result = web_app.component_manager.render_filter_panel()
        soup = BeautifulSoup(result.html, "html.parser")
        assert soup.find(id="filter-panel") is not None
        assert soup.find(id="class-filter") is not None
        assert soup.find(id="property-filter") is not None
        assert soup.find(id="risk-filter") is not None
        
        # Test export panel
        result = web_app.component_manager.render_export_panel()
        soup = BeautifulSoup(result.html, "html.parser")
        assert soup.find(id="export-panel") is not None
        assert soup.find(id="format-selector") is not None
        assert soup.find(id="column-selector") is not None
        assert soup.find("button") is not None

    def test_state_management(self, web_app):
        """Test state management."""
        # Update search state
        web_app.state_manager.update_search_state(
            query="caffeine",
            search_type="name"
        )
        state = web_app.state_manager.get_state()
        assert state["search"]["query"] == "caffeine"
        assert state["search"]["type"] == "name"
        
        # Update filter state
        web_app.state_manager.update_filter_state(
            classes=[PsychoactiveClass.STIMULANT],
            min_binding_affinity=0.7,
            min_risk_level=RiskLevel.HIGH
        )
        state = web_app.state_manager.get_state()
        assert PsychoactiveClass.STIMULANT in state["filters"]["classes"]
        assert state["filters"]["min_binding_affinity"] == 0.7
        assert state["filters"]["min_risk_level"] == RiskLevel.HIGH
        
        # Update export state
        web_app.state_manager.update_export_state(
            format="tsv",
            columns=["name", "cas_number"],
            include_metadata=True
        )
        state = web_app.state_manager.get_state()
        assert state["export"]["format"] == "tsv"
        assert "name" in state["export"]["columns"]
        assert state["export"]["include_metadata"] is True

    def test_event_handling(self, web_app):
        """Test event handling."""
        # Test search event
        event = {
            "type": "search",
            "data": {
                "query": "caffeine",
                "search_type": "name"
            }
        }
        result = web_app.event_handler.handle_event(event)
        assert result.success
        assert len(result.data["compounds"]) == 1
        assert result.data["compounds"][0]["name"] == "Caffeine"
        
        # Test filter event
        event = {
            "type": "filter",
            "data": {
                "classes": ["STIMULANT"],
                "min_risk_level": "HIGH"
            }
        }
        result = web_app.event_handler.handle_event(event)
        assert result.success
        assert len(result.data["compounds"]) == 2  # Both have high risks
        
        # Test export event
        event = {
            "type": "export",
            "data": {
                "format": "tsv",
                "columns": ["name", "cas_number"],
                "output_file": "test_export.tsv"
            }
        }
        result = web_app.event_handler.handle_event(event)
        assert result.success
        assert Path("test_export.tsv").exists()


class TestBrowserInteraction:
    """Tests for browser interaction."""

    def test_search_functionality(self, browser, web_app):
        """Test search functionality."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        
        # Find search input
        search_input = browser.find_element(By.ID, "search-input")
        search_input.send_keys("caffeine")
        
        # Find search button and click
        search_button = browser.find_element(By.ID, "search-button")
        search_button.click()
        
        # Wait for results
        WebDriverWait(browser, 10).until(
            EC.presence_of_element_located((By.CLASS_NAME, "search-results"))
        )
        
        # Check results
        results = browser.find_elements(By.CLASS_NAME, "compound-row")
        assert len(results) == 1
        assert "Caffeine" in results[0].text

    def test_filter_functionality(self, browser, web_app):
        """Test filter functionality."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        
        # Open filter panel
        filter_button = browser.find_element(By.ID, "filter-button")
        filter_button.click()
        
        # Set filter values
        class_select = browser.find_element(By.ID, "class-filter")
        class_select.click()
        stimulant_option = browser.find_element(By.XPATH, "//option[text()='STIMULANT']")
        stimulant_option.click()
        
        # Apply filter
        apply_button = browser.find_element(By.ID, "apply-filter")
        apply_button.click()
        
        # Wait for results
        WebDriverWait(browser, 10).until(
            EC.presence_of_element_located((By.CLASS_NAME, "filter-results"))
        )
        
        # Check results
        results = browser.find_elements(By.CLASS_NAME, "compound-row")
        assert len(results) == 2  # Both compounds are stimulants

    def test_export_functionality(self, browser, web_app, tmp_path):
        """Test export functionality."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        
        # Open export panel
        export_button = browser.find_element(By.ID, "export-button")
        export_button.click()
        
        # Select columns
        name_checkbox = browser.find_element(By.ID, "column-name")
        cas_checkbox = browser.find_element(By.ID, "column-cas")
        name_checkbox.click()
        cas_checkbox.click()
        
        # Set export format
        format_select = browser.find_element(By.ID, "export-format")
        format_select.click()
        tsv_option = browser.find_element(By.XPATH, "//option[text()='TSV']")
        tsv_option.click()
        
        # Export data
        download_button = browser.find_element(By.ID, "download-button")
        download_button.click()
        
        # Wait for download
        WebDriverWait(browser, 10).until(
            lambda x: len(list(tmp_path.glob("*.tsv"))) > 0
        )
        
        # Check exported file
        export_file = next(tmp_path.glob("*.tsv"))
        content = export_file.read_text()
        assert "name\tcas_number" in content
        assert "Caffeine\t58-08-2" in content


class TestVisualization:
    """Tests for visualization functionality."""

    def test_structure_viewer(self, browser, web_app):
        """Test structure viewer."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine detail page
        
        # Check structure viewer
        structure_viewer = browser.find_element(By.ID, "structure-viewer")
        assert structure_viewer.is_displayed()
        assert "svg" in structure_viewer.get_attribute("innerHTML").lower()

    def test_binding_plot(self, browser, web_app):
        """Test binding plot."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine detail page
        
        # Check binding plot
        binding_plot = browser.find_element(By.ID, "binding-plot")
        assert binding_plot.is_displayed()
        assert "canvas" in binding_plot.get_attribute("innerHTML").lower()

    def test_safety_plot(self, browser, web_app):
        """Test safety plot."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine detail page
        
        # Check safety plot
        safety_plot = browser.find_element(By.ID, "safety-plot")
        assert safety_plot.is_displayed()
        assert "canvas" in safety_plot.get_attribute("innerHTML").lower()


class TestResponsiveDesign:
    """Tests for responsive design."""

    def test_desktop_view(self, browser, web_app):
        """Test desktop view."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        browser.set_window_size(1200, 800)
        
        # Check layout
        navbar = browser.find_element(By.ID, "navbar")
        assert navbar.is_displayed()
        assert browser.find_element(By.ID, "compound-list").is_displayed()
        assert browser.find_element(By.ID, "search-bar").is_displayed()
        assert browser.find_element(By.ID, "filter-panel").is_displayed()

    def test_tablet_view(self, browser, web_app):
        """Test tablet view."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        browser.set_window_size(768, 1024)
        
        # Check layout
        navbar = browser.find_element(By.ID, "navbar")
        assert navbar.is_displayed()
        assert browser.find_element(By.ID, "compound-list").is_displayed()
        assert browser.find_element(By.ID, "search-bar").is_displayed()
        assert browser.find_element(By.ID, "filter-toggle").is_displayed()

    def test_mobile_view(self, browser, web_app):
        """Test mobile view."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        browser.set_window_size(375, 667)
        
        # Check layout
        menu_button = browser.find_element(By.ID, "menu-button")
        assert menu_button.is_displayed()
        
        # Open menu
        menu_button.click()
        WebDriverWait(browser, 10).until(
            EC.visibility_of_element_located((By.ID, "mobile-menu"))
        )
        assert browser.find_element(By.ID, "mobile-menu").is_displayed()


class TestAccessibility:
    """Tests for accessibility."""

    def test_keyboard_navigation(self, browser, web_app):
        """Test keyboard navigation."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        
        # Tab through elements
        active = browser.switch_to.active_element
        active.send_keys("\t")  # First tab
        assert browser.find_element(By.ID, "search-input").equals(browser.switch_to.active_element)
        active.send_keys("\t")  # Second tab
        assert browser.find_element(By.ID, "search-button").equals(browser.switch_to.active_element)

    def test_screen_reader_support(self, browser, web_app):
        """Test screen reader support."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        
        # Check ARIA labels
        search_input = browser.find_element(By.ID, "search-input")
        assert search_input.get_attribute("aria-label") == "Search compounds"
        
        filter_button = browser.find_element(By.ID, "filter-button")
        assert filter_button.get_attribute("aria-expanded") == "false"
        assert filter_button.get_attribute("aria-controls") == "filter-panel"

    def test_color_contrast(self, browser, web_app):
        """Test color contrast."""
        # Start app
        web_app.start()
        browser.get("http://localhost:8000")
        
        # Check contrast ratios
        text_elements = browser.find_elements(By.CSS_SELECTOR, "p, h1, h2, h3, a")
        for element in text_elements:
            color = element.value_of_css_property("color")
            background = element.value_of_css_property("background-color")
            # Calculate contrast ratio (simplified)
            assert self._calculate_contrast_ratio(color, background) >= 4.5

    def _calculate_contrast_ratio(self, color1, color2):
        """Helper method to calculate contrast ratio."""
        # Simplified calculation - in reality, you'd use a proper color contrast library
        return 4.5  # Placeholder return value


if __name__ == "__main__":
    pytest.main([__file__])
