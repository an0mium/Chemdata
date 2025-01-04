"""Web visualization functionality for compound data.

This module provides functionality to:
1. Create interactive web visualizations
2. Generate visualization components
3. Handle visualization updates
4. Manage visualization layouts
5. Support real-time updates
"""

import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import json
from pathlib import Path

from ....models.validation import ValidationResult
from .data_enrichment import EnrichedData
from .data_analysis import AnalysisResult
from .data_visualization import VisualizationResult


@dataclass
class WebVisualizationResult(ValidationResult):
    """Result of web visualization."""
    
    components: Dict[str, Dict[str, Any]]
    layouts: Dict[str, Dict[str, Any]]
    configs: Dict[str, Dict[str, Any]]
    stats: Dict[str, Any]
    issues: List[str]


class WebVisualizer:
    """Web visualizer for compound data."""

    # Default visualization layouts
    LAYOUTS = {
        "overview": {
            "distributions": {"x": 0, "y": 0, "w": 12, "h": 6},
            "correlations": {"x": 0, "y": 6, "w": 6, "h": 6},
            "clusters": {"x": 6, "y": 6, "w": 6, "h": 6},
            "outliers": {"x": 0, "y": 12, "w": 12, "h": 6},
            "trends": {"x": 0, "y": 18, "w": 12, "h": 6},
        },
        "detail": {
            "structure": {"x": 0, "y": 0, "w": 4, "h": 4},
            "properties": {"x": 4, "y": 0, "w": 8, "h": 4},
            "predictions": {"x": 0, "y": 4, "w": 12, "h": 4},
            "activity": {"x": 0, "y": 8, "w": 6, "h": 6},
            "safety": {"x": 6, "y": 8, "w": 6, "h": 6},
        },
    }

    def __init__(
        self,
        output_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize web visualizer."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.output_dir = Path(output_dir) if output_dir else None

    def create_visualizations(
        self,
        compounds: List[EnrichedData],
        analysis: AnalysisResult,
        visualization: VisualizationResult,
        layout_type: str = "overview",
    ) -> WebVisualizationResult:
        """Create web visualizations."""
        self.logger.debug("Creating web visualizations")
        
        try:
            components = {}
            layouts = {}
            configs = {}
            stats = {}
            
            # Create visualization components
            components = self._create_components(
                compounds, analysis, visualization
            )
            
            # Get layout configuration
            layouts = self._get_layout(layout_type)
            
            # Generate component configs
            configs = self._generate_configs(
                components, layouts
            )
            
            # Calculate visualization stats
            stats = self._calculate_stats(
                components, layouts, configs
            )
            
            # Save visualization files
            if self.output_dir:
                self._save_visualization_files(
                    components, layouts, configs
                )
            
            return WebVisualizationResult(
                is_valid=True,
                components=components,
                layouts=layouts,
                configs=configs,
                stats=stats,
                issues=[],
            )
            
        except Exception as e:
            self.logger.error(
                f"Error creating web visualizations: {str(e)}",
                exc_info=True
            )
            return WebVisualizationResult(
                is_valid=False,
                components={},
                layouts={},
                configs={},
                stats={},
                issues=[str(e)],
            )

    def _create_components(
        self,
        compounds: List[EnrichedData],
        analysis: AnalysisResult,
        visualization: VisualizationResult,
    ) -> Dict[str, Dict[str, Any]]:
        """Create visualization components."""
        components = {}
        
        # Convert Plotly figures to JSON
        for name, fig in visualization.figures.items():
            if fig:
                components[name] = {
                    "type": "plotly",
                    "data": json.loads(fig.to_json()),
                }
        
        # Add analysis components
        if analysis.property_stats:
            components["property_table"] = {
                "type": "table",
                "data": self._format_property_stats(
                    analysis.property_stats
                ),
            }
        
        if analysis.clusters:
            components["cluster_summary"] = {
                "type": "summary",
                "data": self._format_cluster_summary(
                    analysis.clusters
                ),
            }
        
        if analysis.outliers:
            components["outlier_summary"] = {
                "type": "summary",
                "data": self._format_outlier_summary(
                    analysis.outliers
                ),
            }
        
        if analysis.trends:
            components["trend_summary"] = {
                "type": "summary",
                "data": self._format_trend_summary(
                    analysis.trends
                ),
            }
        
        return components

    def _get_layout(
        self,
        layout_type: str,
    ) -> Dict[str, Dict[str, Any]]:
        """Get layout configuration."""
        if layout_type in self.LAYOUTS:
            return self.LAYOUTS[layout_type]
        return self.LAYOUTS["overview"]

    def _generate_configs(
        self,
        components: Dict[str, Dict[str, Any]],
        layouts: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Dict[str, Any]]:
        """Generate component configurations."""
        configs = {}
        
        for name, component in components.items():
            if name in layouts:
                configs[name] = {
                    "layout": layouts[name],
                    "options": self._get_component_options(
                        name, component
                    ),
                }
        
        return configs

    def _get_component_options(
        self,
        name: str,
        component: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Get component-specific options."""
        options = {
            "title": name.replace("_", " ").title(),
            "showLegend": True,
            "responsive": True,
        }
        
        if component["type"] == "plotly":
            options.update({
                "displayModeBar": True,
                "scrollZoom": True,
                "editable": True,
            })
        elif component["type"] == "table":
            options.update({
                "pagination": True,
                "search": True,
                "sorting": True,
            })
        elif component["type"] == "summary":
            options.update({
                "collapsible": True,
                "expanded": True,
            })
        
        return options

    def _format_property_stats(
        self,
        property_stats: Dict[str, Dict[str, float]],
    ) -> List[Dict[str, Any]]:
        """Format property statistics for table display."""
        rows = []
        
        for prop, stats in property_stats.items():
            row = {"property": prop}
            row.update({
                k: f"{v:.2f}" if isinstance(v, float) else v
                for k, v in stats.items()
            })
            rows.append(row)
        
        return rows

    def _format_cluster_summary(
        self,
        clusters: Dict[str, List[str]],
    ) -> List[Dict[str, Any]]:
        """Format cluster summary."""
        summary = []
        
        for cluster_id, compounds in clusters.items():
            summary.append({
                "cluster": cluster_id,
                "size": len(compounds),
                "compounds": compounds,
            })
        
        return summary

    def _format_outlier_summary(
        self,
        outliers: Dict[str, List[str]],
    ) -> List[Dict[str, Any]]:
        """Format outlier summary."""
        summary = []
        
        for prop, compounds in outliers.items():
            summary.append({
                "property": prop,
                "count": len(compounds),
                "compounds": compounds,
            })
        
        return summary

    def _format_trend_summary(
        self,
        trends: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Format trend summary."""
        summary = []
        
        if "property_trends" in trends:
            for prop, trend in trends["property_trends"].items():
                summary.append({
                    "property": prop,
                    "distribution": trend["distribution"],
                    "modality": trend["modality"],
                    "trend": trend["trend"],
                })
        
        return summary

    def _calculate_stats(
        self,
        components: Dict[str, Dict[str, Any]],
        layouts: Dict[str, Dict[str, Any]],
        configs: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Calculate visualization statistics."""
        stats = {
            "total_components": len(components),
            "component_types": {
                ctype: sum(
                    1 for comp in components.values()
                    if comp["type"] == ctype
                )
                for ctype in {"plotly", "table", "summary"}
            },
            "layout_size": {
                "width": max(
                    layout["x"] + layout["w"]
                    for layout in layouts.values()
                ),
                "height": max(
                    layout["y"] + layout["h"]
                    for layout in layouts.values()
                ),
            },
        }
        return stats

    def _save_visualization_files(
        self,
        components: Dict[str, Dict[str, Any]],
        layouts: Dict[str, Dict[str, Any]],
        configs: Dict[str, Dict[str, Any]],
    ) -> None:
        """Save visualization files."""
        if not self.output_dir:
            return
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save components
        for name, component in components.items():
            filepath = self.output_dir / f"{name}.json"
            with open(filepath, "w") as f:
                json.dump(component, f, indent=2)
        
        # Save layouts
        layout_file = self.output_dir / "layouts.json"
        with open(layout_file, "w") as f:
            json.dump(layouts, f, indent=2)
        
        # Save configs
        config_file = self.output_dir / "configs.json"
        with open(config_file, "w") as f:
            json.dump(configs, f, indent=2)
