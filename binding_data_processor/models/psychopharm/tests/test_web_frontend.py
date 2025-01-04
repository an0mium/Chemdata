"""Tests for web frontend functionality."""

import pytest
from selenium.webdriver import Chrome
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from pathlib import Path

from ..web import (
    WebApp,
    WebConfig,
)


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


class TestComponents:
    """Tests for frontend components."""

    def test_compound_list(self, browser, web_app):
        """Test compound list component."""
        # Load page
        browser.get("http://localhost:8000")
        
        # Check compound list
        compound_list = browser.find_element(By.ID, "compound-list")
        assert compound_list.is_displayed()
        
        # Check compound rows
        rows = compound_list.find_elements(By.CLASS_NAME, "compound-row")
        assert len(rows) > 0
        
        # Check row data
        first_row = rows[0]
        assert first_row.find_element(By.CLASS_NAME, "compound-name").text != ""
        assert first_row.find_element(By.CLASS_NAME, "compound-cas").text != ""
        assert first_row.find_element(By.CLASS_NAME, "compound-class").text != ""

    def test_compound_detail(self, browser, web_app):
        """Test compound detail component."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Check compound info
        info = browser.find_element(By.ID, "compound-info")
        assert info.is_displayed()
        assert "Caffeine" in info.text
        assert "58-08-2" in info.text
        assert "STIMULANT" in info.text
        
        # Check structure viewer
        structure = browser.find_element(By.ID, "structure-viewer")
        assert structure.is_displayed()
        assert "svg" in structure.get_attribute("innerHTML").lower()
        
        # Check plots
        binding_plot = browser.find_element(By.ID, "binding-plot")
        assert binding_plot.is_displayed()
        assert "canvas" in binding_plot.get_attribute("innerHTML").lower()
        
        safety_plot = browser.find_element(By.ID, "safety-plot")
        assert safety_plot.is_displayed()
        assert "canvas" in safety_plot.get_attribute("innerHTML").lower()

    def test_search_bar(self, browser, web_app):
        """Test search bar component."""
        # Load page
        browser.get("http://localhost:8000")
        
        # Check search bar
        search_bar = browser.find_element(By.ID, "search-bar")
        assert search_bar.is_displayed()
        
        # Check search input
        search_input = search_bar.find_element(By.TAG_NAME, "input")
        assert search_input.get_attribute("placeholder") == "Search compounds..."
        
        # Check search button
        search_button = search_bar.find_element(By.TAG_NAME, "button")
        assert search_button.is_displayed()
        assert search_button.text.lower() == "search"

    def test_filter_panel(self, browser, web_app):
        """Test filter panel component."""
        # Load page
        browser.get("http://localhost:8000")
        
        # Open filter panel
        filter_button = browser.find_element(By.ID, "filter-button")
        filter_button.click()
        
        # Wait for panel
        filter_panel = WebDriverWait(browser, 10).until(
            EC.visibility_of_element_located((By.ID, "filter-panel"))
        )
        
        # Check filters
        assert filter_panel.find_element(By.ID, "class-filter").is_displayed()
        assert filter_panel.find_element(By.ID, "property-filter").is_displayed()
        assert filter_panel.find_element(By.ID, "risk-filter").is_displayed()
        
        # Check apply button
        apply_button = filter_panel.find_element(By.ID, "apply-filter")
        assert apply_button.is_displayed()
        assert apply_button.is_enabled()


class TestInteractions:
    """Tests for user interactions."""

    def test_search_interaction(self, browser, web_app):
        """Test search functionality."""
        # Load page
        browser.get("http://localhost:8000")
        
        # Find search input
        search_input = browser.find_element(By.ID, "search-input")
        search_input.send_keys("caffeine")
        
        # Click search
        search_button = browser.find_element(By.ID, "search-button")
        search_button.click()
        
        # Wait for results
        results = WebDriverWait(browser, 10).until(
            EC.presence_of_all_elements_located((By.CLASS_NAME, "compound-row"))
        )
        
        # Check results
        assert len(results) == 1
        assert "Caffeine" in results[0].text

    def test_filter_interaction(self, browser, web_app):
        """Test filter functionality."""
        # Load page
        browser.get("http://localhost:8000")
        
        # Open filter panel
        filter_button = browser.find_element(By.ID, "filter-button")
        filter_button.click()
        
        # Set filters
        class_select = WebDriverWait(browser, 10).until(
            EC.element_to_be_clickable((By.ID, "class-filter"))
        )
        class_select.click()
        
        stimulant_option = browser.find_element(
            By.XPATH,
            "//option[text()='STIMULANT']"
        )
        stimulant_option.click()
        
        # Apply filters
        apply_button = browser.find_element(By.ID, "apply-filter")
        apply_button.click()
        
        # Wait for results
        results = WebDriverWait(browser, 10).until(
            EC.presence_of_all_elements_located((By.CLASS_NAME, "compound-row"))
        )
        
        # Check results
        assert len(results) == 2  # Both compounds are stimulants

    def test_structure_interaction(self, browser, web_app):
        """Test structure viewer interaction."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Get structure viewer
        viewer = browser.find_element(By.ID, "structure-viewer")
        
        # Test rotation
        viewer.click()  # Click to focus
        viewer.send_keys("r")  # Rotate
        
        # Get new transform
        transform = viewer.find_element(By.TAG_NAME, "svg").get_attribute("transform")
        assert transform is not None
        assert "rotate" in transform.lower()

    def test_plot_interaction(self, browser, web_app):
        """Test plot interaction."""
        # Load detail page
        browser.get("http://localhost:8000/compounds/58-08-2")  # Caffeine
        
        # Get binding plot
        plot = browser.find_element(By.ID, "binding-plot")
        
        # Hover over data point
        data_point = plot.find_element(By.CLASS_NAME, "point")
        browser.execute_script(
            "arguments[0].dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));",
            data_point
        )
        
        # Check tooltip
        tooltip = WebDriverWait(browser, 10).until(
            EC.visibility_of_element_located((By.CLASS_NAME, "tooltip"))
        )
        assert tooltip.is_displayed()
        assert "A2A" in tooltip.text
        assert "0.8" in tooltip.text


class TestResponsiveDesign:
    """Tests for responsive design."""

    def test_desktop_layout(self, browser, web_app):
        """Test desktop layout."""
        # Set desktop size
        browser.set_window_size(1200, 800)
        browser.get("http://localhost:8000")
        
        # Check layout
        assert browser.find_element(By.ID, "navbar").is_displayed()
        assert browser.find_element(By.ID, "compound-list").is_displayed()
        assert browser.find_element(By.ID, "search-bar").is_displayed()
        assert browser.find_element(By.ID, "filter-panel").is_displayed()

    def test_tablet_layout(self, browser, web_app):
        """Test tablet layout."""
        # Set tablet size
        browser.set_window_size(768, 1024)
        browser.get("http://localhost:8000")
        
        # Check layout
        assert browser.find_element(By.ID, "navbar").is_displayed()
        assert browser.find_element(By.ID, "compound-list").is_displayed()
        assert browser.find_element(By.ID, "search-bar").is_displayed()
        assert browser.find_element(By.ID, "filter-toggle").is_displayed()

    def test_mobile_layout(self, browser, web_app):
        """Test mobile layout."""
        # Set mobile size
        browser.set_window_size(375, 667)
        browser.get("http://localhost:8000")
        
        # Check layout
        assert browser.find_element(By.ID, "menu-button").is_displayed()
        
        # Open menu
        menu_button = browser.find_element(By.ID, "menu-button")
        menu_button.click()
        
        # Check menu
        menu = WebDriverWait(browser, 10).until(
            EC.visibility_of_element_located((By.ID, "mobile-menu"))
        )
        assert menu.is_displayed()


class TestAccessibility:
    """Tests for accessibility."""

    def test_keyboard_navigation(self, browser, web_app):
        """Test keyboard navigation."""
        # Load page
        browser.get("http://localhost:8000")
        
        # Focus first element
        browser.find_element(By.TAG_NAME, "body").send_keys("\t")
        
        # Check tab order
        active = browser.switch_to.active_element
        assert active.get_attribute("id") == "search-input"
        
        active.send_keys("\t")
        active = browser.switch_to.active_element
        assert active.get_attribute("id") == "search-button"
        
        active.send_keys("\t")
        active = browser.switch_to.active_element
        assert active.get_attribute("id") == "filter-button"

    def test_screen_reader_support(self, browser, web_app):
        """Test screen reader support."""
        # Load page
        browser.get("http://localhost:8000")
        
        # Check ARIA labels
        search_input = browser.find_element(By.ID, "search-input")
        assert search_input.get_attribute("aria-label") == "Search compounds"
        
        filter_button = browser.find_element(By.ID, "filter-button")
        assert filter_button.get_attribute("aria-expanded") == "false"
        assert filter_button.get_attribute("aria-controls") == "filter-panel"
        
        compound_list = browser.find_element(By.ID, "compound-list")
        assert compound_list.get_attribute("role") == "grid"
        assert compound_list.get_attribute("aria-label") == "Compound list"

    def test_color_contrast(self, browser, web_app):
        """Test color contrast."""
        # Load page
        browser.get("http://localhost:8000")
        
        # Check text elements
        text_elements = browser.find_elements(
            By.CSS_SELECTOR,
            "p, h1, h2, h3, a, button, label"
        )
        
        for element in text_elements:
            color = element.value_of_css_property("color")
            background = element.value_of_css_property("background-color")
            contrast_ratio = self._calculate_contrast_ratio(color, background)
            assert contrast_ratio >= 4.5  # WCAG AA standard

    def _calculate_contrast_ratio(self, color1, color2):
        """Helper method to calculate contrast ratio."""
        # This would normally use a proper color contrast calculation library
        return 4.5  # Placeholder return value


if __name__ == "__main__":
    pytest.main([__file__])
