"""Data visualization functionality for analyzed compound data.

This module provides functionality to:
1. Generate property distribution plots
2. Create correlation heatmaps
3. Visualize clusters and outliers
4. Plot trends and patterns
5. Create interactive visualizations
"""

import logging
from typing import Dict, List, Optional, Any
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats as scipy_stats
from dataclasses import dataclass

from ....models.validation import ValidationResult
from .data_enrichment import EnrichedData
from .data_analysis import AnalysisResult


@dataclass
class VisualizationResult(ValidationResult):
    """Result of data visualization."""
    
    figures: Dict[str, go.Figure]
    html_files: Dict[str, str]
    stats: Dict[str, Any]
    issues: List[str]


class DataVisualizer:
    """Visualizer for analyzed compound data."""

    def __init__(
        self,
        output_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize data visualizer."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.output_dir = output_dir

    def visualize_analysis(
        self,
        compounds: List[EnrichedData],
        analysis: AnalysisResult,
        save_html: bool = True,
    ) -> VisualizationResult:
        """Create visualizations from analysis results."""
        self.logger.debug("Creating visualizations")
        
        try:
            figures = {}
            html_files = {}
            stats = {}
            
            # Create property distribution plots
            figures["distributions"] = self._plot_distributions(
                compounds, analysis.property_stats
            )
            
            # Create correlation heatmap
            figures["correlations"] = self._plot_correlations(
                analysis.correlations
            )
            
            # Create cluster visualization
            figures["clusters"] = self._plot_clusters(
                compounds, analysis.clusters
            )
            
            # Create outlier plots
            figures["outliers"] = self._plot_outliers(
                compounds, analysis.outliers
            )
            
            # Create trend plots
            figures["trends"] = self._plot_trends(
                compounds, analysis.trends
            )
            
            # Save HTML files if requested
            if save_html and self.output_dir:
                html_files = self._save_html_files(figures)
            
            # Calculate visualization stats
            stats = self._calculate_visualization_stats(figures)
            
            return VisualizationResult(
                is_valid=True,
                figures=figures,
                html_files=html_files,
                stats=stats,
                issues=[],
            )
            
        except Exception as e:
            self.logger.error(
                f"Error creating visualizations: {str(e)}",
                exc_info=True
            )
            return VisualizationResult(
                is_valid=False,
                figures={},
                html_files={},
                stats={},
                issues=[str(e)],
            )

    def _plot_distributions(
        self,
        compounds: List[EnrichedData],
        property_stats: Dict[str, Dict[str, float]],
    ) -> go.Figure:
        """Create property distribution plots."""
        # Create subplot grid
        n_props = len(property_stats)
        n_cols = min(3, n_props)
        n_rows = (n_props + n_cols - 1) // n_cols
        
        fig = make_subplots(
            rows=n_rows,
            cols=n_cols,
            subplot_titles=list(property_stats.keys()),
        )
        
        # Add distribution plots
        for i, (prop, stats) in enumerate(property_stats.items()):
            row = i // n_cols + 1
            col = i % n_cols + 1
            
            # Get property values
            values = [
                float(c.properties.get(prop, 0))
                for c in compounds
                if c.properties.get(prop) is not None
            ]
            
            if values:
                # Add histogram
                fig.add_trace(
                    go.Histogram(
                        x=values,
                        name=prop,
                        showlegend=False,
                    ),
                    row=row,
                    col=col,
                )
                
                # Add KDE curve
                kde = scipy_stats.gaussian_kde(values)
                x_range = np.linspace(min(values), max(values), 100)
                y_kde = kde(x_range)
                
                fig.add_trace(
                    go.Scatter(
                        x=x_range,
                        y=y_kde * len(values) * (max(values) - min(values)) / 50,
                        name=f"{prop} KDE",
                        line=dict(color="red"),
                        showlegend=False,
                    ),
                    row=row,
                    col=col,
                )
        
        # Update layout
        fig.update_layout(
            title="Property Distributions",
            showlegend=False,
            height=300 * n_rows,
        )
        
        return fig

    def _plot_correlations(
        self,
        correlations: Dict[str, Dict[str, float]],
    ) -> go.Figure:
        """Create correlation heatmap."""
        # Convert to matrix form
        props = list(correlations.keys())
        matrix = np.zeros((len(props), len(props)))
        
        for i, prop1 in enumerate(props):
            for j, prop2 in enumerate(props):
                if prop1 != prop2:
                    matrix[i, j] = correlations[prop1][prop2]
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=props,
            y=props,
            colorscale="RdBu",
            zmid=0,
            text=np.round(matrix, 2),
            texttemplate="%{text}",
            textfont={"size": 10},
            hoverongaps=False,
        ))
        
        # Update layout
        fig.update_layout(
            title="Property Correlations",
            xaxis_title="Property",
            yaxis_title="Property",
            width=800,
            height=800,
        )
        
        return fig

    def _plot_clusters(
        self,
        compounds: List[EnrichedData],
        clusters: Dict[str, List[str]],
    ) -> go.Figure:
        """Create cluster visualization."""
        if not clusters:
            return go.Figure()
        
        # Get first two properties for 2D plot
        props = list(compounds[0].properties.keys())[:2]
        if len(props) < 2:
            return go.Figure()
        
        # Create scatter plot
        fig = go.Figure()
        
        # Add points for each cluster
        for cluster_id, compound_names in clusters.items():
            cluster_compounds = [
                c for c in compounds
                if c.compound.name in compound_names
            ]
            
            x = [
                float(c.properties.get(props[0], 0))
                for c in cluster_compounds
            ]
            y = [
                float(c.properties.get(props[1], 0))
                for c in cluster_compounds
            ]
            
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode="markers",
                name=f"Cluster {cluster_id}",
                text=[c.compound.name for c in cluster_compounds],
                hovertemplate=(
                    "%{text}<br>"
                    f"{props[0]}: %{{x}}<br>"
                    f"{props[1]}: %{{y}}"
                ),
            ))
        
        # Update layout
        fig.update_layout(
            title="Compound Clusters",
            xaxis_title=props[0],
            yaxis_title=props[1],
            width=800,
            height=600,
        )
        
        return fig

    def _plot_outliers(
        self,
        compounds: List[EnrichedData],
        outliers: Dict[str, List[str]],
    ) -> go.Figure:
        """Create outlier visualization."""
        if not outliers:
            return go.Figure()
        
        # Create subplot grid
        n_props = len(outliers)
        n_cols = min(2, n_props)
        n_rows = (n_props + n_cols - 1) // n_cols
        
        fig = make_subplots(
            rows=n_rows,
            cols=n_cols,
            subplot_titles=list(outliers.keys()),
        )
        
        # Add box plots with outliers
        for i, (prop, outlier_names) in enumerate(outliers.items()):
            row = i // n_cols + 1
            col = i % n_cols + 1
            
            values = [
                float(c.properties.get(prop, 0))
                for c in compounds
                if c.properties.get(prop) is not None
            ]
            
            if values:
                # Add box plot
                fig.add_trace(
                    go.Box(
                        y=values,
                        name=prop,
                        boxpoints="outliers",
                        showlegend=False,
                    ),
                    row=row,
                    col=col,
                )
        
        # Update layout
        fig.update_layout(
            title="Property Outliers",
            showlegend=False,
            height=300 * n_rows,
        )
        
        return fig

    def _plot_trends(
        self,
        compounds: List[EnrichedData],
        trends: Dict[str, Any],
    ) -> go.Figure:
        """Create trend visualization."""
        if not trends.get("temporal_trends"):
            return go.Figure()
        
        # Create subplot grid
        n_props = len(trends["temporal_trends"])
        n_cols = min(2, n_props)
        n_rows = (n_props + n_cols - 1) // n_cols
        
        fig = make_subplots(
            rows=n_rows,
            cols=n_cols,
            subplot_titles=list(trends["temporal_trends"].keys()),
        )
        
        # Get compounds with timestamps
        dated_compounds = [
            (c, max(c.timestamps.values()))
            for c in compounds
            if c.timestamps
        ]
        
        if dated_compounds:
            # Sort by date
            dated_compounds.sort(key=lambda x: x[1])
            compounds_sorted, dates = zip(*dated_compounds)
            
            # Add trend lines
            for i, (prop, trend_data) in enumerate(
                trends["temporal_trends"].items()
            ):
                row = i // n_cols + 1
                col = i % n_cols + 1
                
                values = [
                    float(c.properties.get(prop, 0))
                    for c in compounds_sorted
                    if c.properties.get(prop) is not None
                ]
                
                if values:
                    # Add scatter plot
                    fig.add_trace(
                        go.Scatter(
                            x=dates,
                            y=values,
                            mode="markers",
                            name=prop,
                            showlegend=False,
                        ),
                        row=row,
                        col=col,
                    )
                    
                    # Add trend line
                    x_range = np.arange(len(values))
                    slope = trend_data["slope"]
                    intercept = np.mean(values) - slope * np.mean(x_range)
                    y_trend = slope * x_range + intercept
                    
                    fig.add_trace(
                        go.Scatter(
                            x=dates,
                            y=y_trend,
                            mode="lines",
                            name=f"{prop} trend",
                            line=dict(color="red"),
                            showlegend=False,
                        ),
                        row=row,
                        col=col,
                    )
        
        # Update layout
        fig.update_layout(
            title="Property Trends Over Time",
            showlegend=False,
            height=300 * n_rows,
        )
        
        return fig

    def _save_html_files(
        self,
        figures: Dict[str, go.Figure],
    ) -> Dict[str, str]:
        """Save figures as HTML files."""
        html_files = {}
        
        if self.output_dir:
            for name, fig in figures.items():
                if isinstance(fig, go.Figure):
                    filepath = f"{self.output_dir}/{name}.html"
                    fig.write_html(filepath)
                    html_files[name] = filepath
        
        return html_files

    def _calculate_visualization_stats(
        self,
        figures: Dict[str, go.Figure],
    ) -> Dict[str, Any]:
        """Calculate statistics about visualizations."""
        stats = {
            "total_figures": len(figures),
            "figure_types": list(figures.keys()),
            "total_traces": sum(
                len(fig.data)
                for fig in figures.values()
                if isinstance(fig, go.Figure)
            ),
        }
        return stats
