"""Web interface functionality for compound data.

This module provides functionality to:
1. Create web interface components
2. Handle user interactions
3. Manage data updates
4. Support filtering and searching
5. Enable real-time visualization updates
"""

import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from pathlib import Path
import json
import datetime

from ....models.validation import ValidationResult
from .data_enrichment import EnrichedData
from .data_analysis import AnalysisResult
from .data_visualization import VisualizationResult
from .web_visualization import WebVisualizationResult


@dataclass
class WebInterfaceResult(ValidationResult):
    """Result of web interface generation."""
    
    components: Dict[str, Dict[str, Any]]
    templates: Dict[str, str]
    assets: Dict[str, Path]
    configs: Dict[str, Dict[str, Any]]
    stats: Dict[str, Any]
    issues: List[str]


class WebInterface:
    """Web interface for compound data."""

    # Default interface configurations
    CONFIGS = {
        "theme": {
            "primary": "#007bff",
            "secondary": "#6c757d",
            "success": "#28a745",
            "danger": "#dc3545",
            "warning": "#ffc107",
            "info": "#17a2b8",
        },
        "layout": {
            "sidebar_width": 250,
            "content_width": "calc(100% - 250px)",
            "header_height": 60,
            "footer_height": 40,
        },
        "features": {
            "dark_mode": True,
            "responsive": True,
            "animations": True,
            "tooltips": True,
            "keyboard_shortcuts": True,
        },
    }

    def __init__(
        self,
        output_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize web interface."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.output_dir = Path(output_dir) if output_dir else None

    def create_interface(
        self,
        compounds: List[EnrichedData],
        analysis: AnalysisResult,
        visualization: VisualizationResult,
        web_visualization: WebVisualizationResult,
    ) -> WebInterfaceResult:
        """Create web interface."""
        self.logger.debug("Creating web interface")
        
        try:
            components = {}
            templates = {}
            assets = {}
            configs = {}
            stats = {}
            
            # Create interface components
            components = self._create_components(
                compounds, analysis, visualization, web_visualization
            )
            
            # Generate templates
            templates = self._generate_templates(
                components
            )
            
            # Copy assets
            assets = self._copy_assets()
            
            # Generate configs
            configs = self._generate_configs(
                components
            )
            
            # Calculate interface stats
            stats = self._calculate_stats(
                components, templates, assets, configs
            )
            
            # Save interface files
            if self.output_dir:
                self._save_interface_files(
                    components, templates, assets, configs
                )
            
            return WebInterfaceResult(
                is_valid=True,
                components=components,
                templates=templates,
                assets=assets,
                configs=configs,
                stats=stats,
                issues=[],
            )
            
        except Exception as e:
            self.logger.error(
                f"Error creating web interface: {str(e)}",
                exc_info=True
            )
            return WebInterfaceResult(
                is_valid=False,
                components={},
                templates={},
                assets={},
                configs={},
                stats={},
                issues=[str(e)],
            )

    def _create_components(
        self,
        compounds: List[EnrichedData],
        analysis: AnalysisResult,
        visualization: VisualizationResult,
        web_visualization: WebVisualizationResult,
    ) -> Dict[str, Dict[str, Any]]:
        """Create interface components."""
        components = {
            "header": self._create_header_component(),
            "sidebar": self._create_sidebar_component(compounds),
            "content": self._create_content_component(
                compounds, analysis, visualization, web_visualization
            ),
            "footer": self._create_footer_component(),
        }
        return components

    def _create_header_component(self) -> Dict[str, Any]:
        """Create header component."""
        return {
            "type": "header",
            "data": {
                "title": "ChemData",
                "subtitle": "Psychopharmacological Compound Database",
                "menu": [
                    {"label": "Home", "route": "/"},
                    {"label": "Browse", "route": "/browse"},
                    {"label": "Analysis", "route": "/analysis"},
                    {"label": "Settings", "route": "/settings"},
                ],
            },
        }

    def _create_sidebar_component(
        self,
        compounds: List[EnrichedData],
    ) -> Dict[str, Any]:
        """Create sidebar component."""
        return {
            "type": "sidebar",
            "data": {
                "filters": self._create_filters(compounds),
                "stats": self._create_stats(compounds),
                "legend": self._create_legend(),
            },
        }

    def _create_content_component(
        self,
        compounds: List[EnrichedData],
        analysis: AnalysisResult,
        visualization: VisualizationResult,
        web_visualization: WebVisualizationResult,
    ) -> Dict[str, Any]:
        """Create content component."""
        return {
            "type": "content",
            "data": {
                "visualizations": web_visualization.components,
                "layouts": web_visualization.layouts,
                "configs": web_visualization.configs,
                "tables": self._create_tables(compounds),
                "details": self._create_details(compounds),
            },
        }

    def _create_footer_component(self) -> Dict[str, Any]:
        """Create footer component."""
        return {
            "type": "footer",
            "data": {
                "copyright": f"© {datetime.datetime.now().year} ChemData",
                "version": "1.0.0",
                "links": [
                    {"label": "About", "url": "/about"},
                    {"label": "API", "url": "/api"},
                    {"label": "Contact", "url": "/contact"},
                ],
            },
        }

    def _create_filters(
        self,
        compounds: List[EnrichedData],
    ) -> Dict[str, Any]:
        """Create filter configurations."""
        return {
            "text": {
                "type": "text",
                "label": "Search",
                "placeholder": "Search compounds...",
            },
            "target": {
                "type": "select",
                "label": "Target",
                "options": sorted(set(
                    c.targets[0] if c.targets else "Unknown"
                    for c in compounds
                )),
            },
            "activity": {
                "type": "select",
                "label": "Activity",
                "options": sorted(set(
                    c.activity_type if hasattr(c, "activity_type") else "Unknown"
                    for c in compounds
                )),
            },
            "property": {
                "type": "range",
                "label": "Property",
                "options": {
                    prop: {
                        "min": min(
                            float(c.properties.get(prop, 0))
                            for c in compounds
                            if c.properties.get(prop) is not None
                        ),
                        "max": max(
                            float(c.properties.get(prop, 0))
                            for c in compounds
                            if c.properties.get(prop) is not None
                        ),
                    }
                    for prop in {"molecular_weight", "logp", "psa"}
                },
            },
        }

    def _create_stats(
        self,
        compounds: List[EnrichedData],
    ) -> Dict[str, Any]:
        """Create statistics summary."""
        return {
            "total": len(compounds),
            "with_activity": sum(
                1 for c in compounds if c.activities
            ),
            "with_predictions": sum(
                1 for c in compounds
                if any(
                    getattr(c, f"{pred_type}_predictions", None)
                    for pred_type in {
                        "bbb", "activity", "toxicity", "abuse"
                    }
                )
            ),
            "with_web_data": sum(
                1 for c in compounds if c.web_data
            ),
        }

    def _create_legend(self) -> Dict[str, Any]:
        """Create visualization legend."""
        return {
            "colors": {
                "activity": {
                    "high": "#28a745",
                    "medium": "#ffc107",
                    "low": "#dc3545",
                },
                "predictions": {
                    "confident": "#007bff",
                    "uncertain": "#6c757d",
                },
                "safety": {
                    "safe": "#28a745",
                    "caution": "#ffc107",
                    "warning": "#dc3545",
                },
            },
            "symbols": {
                "experimental": "●",
                "predicted": "○",
                "literature": "■",
                "community": "▲",
            },
        }

    def _create_tables(
        self,
        compounds: List[EnrichedData],
    ) -> Dict[str, Any]:
        """Create table configurations."""
        return {
            "compounds": {
                "columns": [
                    {"key": "name", "label": "Name"},
                    {"key": "smiles", "label": "SMILES"},
                    {"key": "target", "label": "Target"},
                    {"key": "activity", "label": "Activity"},
                    {"key": "predictions", "label": "Predictions"},
                    {"key": "safety", "label": "Safety"},
                ],
                "data": [
                    self._format_compound_row(compound)
                    for compound in compounds
                ],
                "options": {
                    "pagination": True,
                    "sorting": True,
                    "filtering": True,
                    "export": True,
                },
            },
        }

    def _create_details(
        self,
        compounds: List[EnrichedData],
    ) -> Dict[str, Any]:
        """Create detail view configurations."""
        return {
            "tabs": [
                {"key": "overview", "label": "Overview"},
                {"key": "activity", "label": "Activity"},
                {"key": "predictions", "label": "Predictions"},
                {"key": "safety", "label": "Safety"},
                {"key": "community", "label": "Community"},
            ],
            "sections": {
                "overview": [
                    {"key": "structure", "label": "Structure"},
                    {"key": "properties", "label": "Properties"},
                    {"key": "targets", "label": "Targets"},
                ],
                "activity": [
                    {"key": "experimental", "label": "Experimental"},
                    {"key": "literature", "label": "Literature"},
                    {"key": "predictions", "label": "Predictions"},
                ],
                "predictions": [
                    {"key": "binding", "label": "Binding"},
                    {"key": "activity", "label": "Activity"},
                    {"key": "toxicity", "label": "Toxicity"},
                ],
                "safety": [
                    {"key": "warnings", "label": "Warnings"},
                    {"key": "risks", "label": "Risks"},
                    {"key": "contraindications", "label": "Contraindications"},
                ],
                "community": [
                    {"key": "reports", "label": "Reports"},
                    {"key": "discussions", "label": "Discussions"},
                    {"key": "references", "label": "References"},
                ],
            },
        }

    def _format_compound_row(
        self,
        compound: EnrichedData,
    ) -> Dict[str, Any]:
        """Format compound data for table row."""
        return {
            "name": compound.compound.name,
            "smiles": compound.compound.smiles,
            "target": compound.targets[0] if compound.targets else "Unknown",
            "activity": {
                "type": compound.activity_type if hasattr(
                    compound, "activity_type"
                ) else "Unknown",
                "value": compound.activity_value if hasattr(
                    compound, "activity_value"
                ) else None,
                "unit": compound.activity_unit if hasattr(
                    compound, "activity_unit"
                ) else None,
            },
            "predictions": {
                pred_type: getattr(compound, f"{pred_type}_predictions", {})
                for pred_type in {"bbb", "activity", "toxicity", "abuse"}
            },
            "safety": {
                "warnings": compound.warnings if hasattr(
                    compound, "warnings"
                ) else [],
                "risks": compound.risks if hasattr(
                    compound, "risks"
                ) else [],
            },
        }

    def _generate_templates(
        self,
        components: Dict[str, Dict[str, Any]],
    ) -> Dict[str, str]:
        """Generate HTML templates."""
        templates = {}
        
        # Generate base template
        templates["base"] = self._generate_base_template()
        
        # Generate component templates
        for name, component in components.items():
            templates[name] = self._generate_component_template(
                name, component
            )
        
        return templates

    def _generate_base_template(self) -> str:
        """Generate base HTML template."""
        return """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChemData</title>
    <link rel="stylesheet" href="/static/css/styles.css">
</head>
<body>
    <div id="app">
        {{ header }}
        <div class="main">
            {{ sidebar }}
            <div class="content">
                {{ content }}
            </div>
        </div>
        {{ footer }}
    </div>
    <script src="/static/js/app.js"></script>
</body>
</html>
"""

    def _generate_component_template(
        self,
        name: str,
        component: Dict[str, Any],
    ) -> str:
        """Generate component HTML template."""
        if component["type"] == "header":
            return self._generate_header_template(component)
        elif component["type"] == "sidebar":
            return self._generate_sidebar_template(component)
        elif component["type"] == "content":
            return self._generate_content_template(component)
        elif component["type"] == "footer":
            return self._generate_footer_template(component)
        return ""

    def _generate_header_template(
        self,
        component: Dict[str, Any],
    ) -> str:
        """Generate header HTML template."""
        return """
<header class="header">
    <div class="header-title">
        <h1>{{ title }}</h1>
        <p>{{ subtitle }}</p>
    </div>
    <nav class="header-menu">
        {% for item in menu %}
        <a href="{{ item.route }}" class="menu-item">{{ item.label }}</a>
        {% endfor %}
    </nav>
</header>
"""

    def _generate_sidebar_template(
        self,
        component: Dict[str, Any],
    ) -> str:
        """Generate sidebar HTML template."""
        return """
<aside class="sidebar">
    <div class="filters">
        {% for filter in filters %}
        <div class="filter">
            <label>{{ filter.label }}</label>
            {% if filter.type == "text" %}
            <input type="text" placeholder="{{ filter.placeholder }}">
            {% elif filter.type == "select" %}
            <select>
                {% for option in filter.options %}
                <option value="{{ option }}">{{ option }}</option>
                {% endfor %}
            </select>
            {% elif filter.type == "range" %}
            <div class="range">
                <input type="range" min="{{ filter.min }}" max="{{ filter.max }}">
                <span class="value"></span>
            </div>
            {% endif %}
        </div>
        {% endfor %}
    </div>
    <div class="stats">
        {% for stat in stats %}
        <div class="stat">
            <label>{{ stat.label }}</label>
            <span>{{ stat.value }}</span>
        </div>
        {% endfor %}
    </div>
    <div class="legend">
        {% for item in legend %}
        <div class="legend-item">
            <span class="symbol" style="color: {{ item.color }}">{{ item.symbol }}</span>
            <span class="label">{{ item.label }}</span>
        </div>
        {% endfor %}
    </div>
</aside>
"""

    def _generate_content_template(
        self,
        component: Dict[str, Any],
    ) -> str:
        """Generate content HTML template."""
        return """
<main class="content">
    <div class="visualizations">
        {% for viz in visualizations %}
        <div class="visualization" data-config="{{ viz.config }}">
            {{ viz.content }}
        </div>
        {% endfor %}
    </div>
    <div class="tables">
        {% for table in tables %}
        <div class="table" data-config="{{ table.config }}">
            {{ table.content }}
        </div>
        {% endfor %}
    </div>
    <div class="details">
        {% for detail in details %}
        <div class="detail" data-config="{{ detail.config }}">
            {{ detail.content }}
        </div>
        {% endfor %}
    </div>
</main>
"""

    def _generate_footer_template(
        self,
        component: Dict[str, Any],
    ) -> str:
        """Generate footer HTML template."""
        return """
<footer class="footer">
    <div class="copyright">{{ copyright }}</div>
    <div class="version">v{{ version }}</div>
    <nav class="footer-links">
        {% for link in links %}
        <a href="{{ link.url }}">{{ link.label }}</a>
        {% endfor %}
    </nav>
</footer>
"""

    def _copy_assets(self) -> Dict[str, Path]:
        """Copy static assets to output directory."""
        if not self.output_dir:
            return {}
        
        assets = {}
        static_dir = self.output_dir / "static"
        
        # Create static directories
        for subdir in ["css", "js", "img"]:
            (static_dir / subdir).mkdir(parents=True, exist_ok=True)
        
        # Copy CSS files
        css_file = static_dir / "css" / "styles.css"
        css_file.write_text(self._generate_css())
        assets["css"] = css_file
        
        # Copy JavaScript files
        js_file = static_dir / "js" / "app.js"
        js_file.write_text(self._generate_js())
        assets["js"] = js_file
        
        return assets

    def _generate_css(self) -> str:
        """Generate CSS styles."""
        return """
:root {
    --primary: #007bff;
    --secondary: #6c757d;
    --success: #28a745;
    --danger: #dc3545;
    --warning: #ffc107;
    --info: #17a2b8;
}

/* Layout */
.app {
    display: flex;
    flex-direction: column;
    min-height: 100vh;
}

.main {
    display: flex;
    flex: 1;
}

/* Components */
.header {
    height: 60px;
    padding: 0 20px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    background: var(--primary);
    color: white;
}

.sidebar {
    width: 250px;
    padding: 20px;
    background: #f8f9fa;
}

.content {
    flex: 1;
    padding: 20px;
}

.footer {
    height: 40px;
    padding: 0 20px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    background: #f8f9fa;
}

/* Visualizations */
.visualization {
    margin-bottom: 20px;
    padding: 20px;
    background: white;
    border-radius: 4px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}

/* Tables */
.table {
    margin-bottom: 20px;
    background: white;
    border-radius: 4px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}

/* Details */
.detail {
    margin-bottom: 20px;
    padding: 20px;
    background: white;
    border-radius: 4px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}

/* Responsive */
@media (max-width: 768px) {
    .main {
        flex-direction: column;
    }
    
    .sidebar {
        width: 100%;
    }
}
"""

    def _generate_js(self) -> str:
        """Generate JavaScript code."""
        return """
// Initialize application
document.addEventListener('DOMContentLoaded', function() {
    initializeFilters();
    initializeVisualizations();
    initializeTables();
    initializeDetails();
});

// Initialize filters
function initializeFilters() {
    const filters = document.querySelectorAll('.filter');
    filters.forEach(filter => {
        const input = filter.querySelector('input, select');
        if (input) {
            input.addEventListener('change', function() {
                updateResults();
            });
        }
    });
}

// Initialize visualizations
function initializeVisualizations() {
    const visualizations = document.querySelectorAll('.visualization');
    visualizations.forEach(viz => {
        const config = JSON.parse(viz.dataset.config);
        renderVisualization(viz, config);
    });
}

// Initialize tables
function initializeTables() {
    const tables = document.querySelectorAll('.table');
    tables.forEach(table => {
        const config = JSON.parse(table.dataset.config);
        renderTable(table, config);
    });
}

// Initialize details
function initializeDetails() {
    const details = document.querySelectorAll('.detail');
    details.forEach(detail => {
        const config = JSON.parse(detail.dataset.config);
        renderDetail(detail, config);
    });
}

// Update results based on filters
function updateResults() {
    const filters = getFilterValues();
    updateVisualizations(filters);
    updateTables(filters);
    updateDetails(filters);
}

// Get current filter values
function getFilterValues() {
    const filters = {};
    document.querySelectorAll('.filter').forEach(filter => {
        const input = filter.querySelector('input, select');
        if (input) {
            filters[input.name] = input.value;
        }
    });
    return filters;
}

// Render visualization
function renderVisualization(container, config) {
    // Implementation depends on visualization type
}

// Render table
function renderTable(container, config) {
    // Implementation depends on table type
}

// Render detail
function renderDetail(container, config) {
    // Implementation depends on detail type
}

// Update visualizations
function updateVisualizations(filters) {
    document.querySelectorAll('.visualization').forEach(viz => {
        const config = JSON.parse(viz.dataset.config);
        updateVisualization(viz, config, filters);
    });
}

// Update tables
function updateTables(filters) {
    document.querySelectorAll('.table').forEach(table => {
        const config = JSON.parse(table.dataset.config);
        updateTable(table, config, filters);
    });
}

// Update details
function updateDetails(filters) {
    document.querySelectorAll('.detail').forEach(detail => {
        const config = JSON.parse(detail.dataset.config);
        updateDetail(detail, config, filters);
    });
}

// Update visualization
function updateVisualization(container, config, filters) {
    // Implementation depends on visualization type
}

// Update table
function updateTable(container, config, filters) {
    // Implementation depends on table type
}

// Update detail
function updateDetail(container, config, filters) {
    // Implementation depends on detail type
}
"""

    def _generate_configs(
        self,
        components: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Dict[str, Any]]:
        """Generate interface configurations."""
        configs = self.CONFIGS.copy()
        
        # Add component-specific configs
        configs["components"] = {
            name: self._get_component_config(component)
            for name, component in components.items()
        }
        
        return configs

    def _get_component_config(
        self,
        component: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Get component-specific configuration."""
        config = {
            "type": component["type"],
            "visible": True,
            "enabled": True,
        }
        
        if component["type"] == "header":
            config.update({
                "fixed": True,
                "transparent": False,
            })
        elif component["type"] == "sidebar":
            config.update({
                "collapsible": True,
                "collapsed": False,
            })
        elif component["type"] == "content":
            config.update({
                "scrollable": True,
                "padding": 20,
            })
        elif component["type"] == "footer":
            config.update({
                "fixed": True,
                "transparent": True,
            })
        
        return config

    def _calculate_stats(
        self,
        components: Dict[str, Dict[str, Any]],
        templates: Dict[str, str],
        assets: Dict[str, Path],
        configs: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Calculate interface statistics."""
        stats = {
            "components": {
                "total": len(components),
                "types": {
                    ctype: sum(
                        1 for c in components.values()
                        if c["type"] == ctype
                    )
                    for ctype in {
                        "header", "sidebar", "content", "footer"
                    }
                },
            },
            "templates": {
                "total": len(templates),
                "size": sum(
                    len(template)
                    for template in templates.values()
                ),
            },
            "assets": {
                "total": len(assets),
                "size": sum(
                    asset.stat().st_size
                    for asset in assets.values()
                ),
            },
            "configs": {
                "total": len(configs),
                "features": len(configs.get("features", {})),
            },
        }
        return stats

    def _save_interface_files(
        self,
        components: Dict[str, Dict[str, Any]],
        templates: Dict[str, str],
        assets: Dict[str, Path],
        configs: Dict[str, Dict[str, Any]],
    ) -> None:
        """Save interface files."""
        if not self.output_dir:
            return
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save components
        components_file = self.output_dir / "components.json"
        with open(components_file, "w") as f:
            json.dump(components, f, indent=2)
        
        # Save templates
        templates_dir = self.output_dir / "templates"
        templates_dir.mkdir(exist_ok=True)
        for name, template in templates.items():
            template_file = templates_dir / f"{name}.html"
            template_file.write_text(template)
        
        # Save configs
        config_file = self.output_dir / "config.json"
        with open(config_file, "w") as f:
            json.dump(configs, f, indent=2)
