"""Advanced plot management and ML visualization.

This module provides:
1. Advanced plot creation and customization
2. ML visualization integration
3. Plot combination and layout management
4. Interactive visualization features
5. Color scheme management
"""

import logging
from typing import Dict, List, Optional, Tuple, Union
from functools import lru_cache

import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from dash import html, dcc
from wordcloud import WordCloud
import networkx as nx
from datetime import datetime, timedelta

from binding_data_processor.processors.structure.ml.visualization import MLVisualizer
from binding_data_processor.processors.structure.ml.advanced_visualization import (
    AdvancedVisualizer,
)


class PlotManager:
    """Advanced plot management with ML visualization integration."""

    # Color schemes for different plot types
    COLORSCALES = {
        "binding": "RdBu",
        "abuse": "Reds",
        "toxicity": "RdYlGn_r",
        "activity": "Viridis",
        "psychoactive": "Plasma",
        "nootropic": "Magma",
        "interactions": "Inferno",
        "similarity": "Blues",
        "uncertainty": "Greys",
        "timeline": "Viridis",
        "network": "Set3",
        "heatmap": "RdBu",
        "wordcloud": "hsv",
    }

    # Plot configurations
    PLOT_CONFIGS = {
        "binding": {
            "type": "heatmap",
            "title": "Binding Affinity Predictions",
            "xaxis_title": "Receptor",
            "yaxis_title": "Affinity (pKi)",
        },
        "abuse": {
            "type": "radar",
            "title": "Abuse Potential Assessment",
            "radial_title": "Score",
            "angular_title": "Category",
        },
        "toxicity": {
            "type": "bar",
            "title": "Toxicity Risk Assessment",
            "xaxis_title": "Endpoint",
            "yaxis_title": "Risk Level",
        },
        "activity": {
            "type": "scatter",
            "title": "Activity Profile",
            "xaxis_title": "Target",
            "yaxis_title": "Activity",
        },
        "psychoactive": {
            "type": "radar",
            "title": "Psychoactive Effects Profile",
            "radial_title": "Score",
            "angular_title": "Effect",
        },
        "nootropic": {
            "type": "radar",
            "title": "Nootropic Effects Profile",
            "radial_title": "Score",
            "angular_title": "Effect",
        },
        "interactions": {
            "type": "network",
            "title": "Drug Interaction Network",
            "node_title": "Compound",
            "edge_title": "Interaction",
        },
        "similarity": {
            "type": "heatmap",
            "title": "Structural Similarity Matrix",
            "xaxis_title": "Compound",
            "yaxis_title": "Compound",
        },
        "uncertainty": {
            "type": "violin",
            "title": "Prediction Uncertainty",
            "xaxis_title": "Model",
            "yaxis_title": "Uncertainty",
        },
    }

    def __init__(
        self,
        ml_visualizer: Optional[MLVisualizer] = None,
        advanced_visualizer: Optional[AdvancedVisualizer] = None,
    ):
        """Initialize plot manager.

        Args:
            ml_visualizer: ML visualization component
            advanced_visualizer: Advanced visualization component
        """
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        self.logger.addHandler(handler)

        self.ml_visualizer = ml_visualizer or MLVisualizer()
        self.advanced_visualizer = advanced_visualizer or AdvancedVisualizer()

    def create_combined_plot(
        self,
        predictions: Dict,
        selected_types: List[str],
        layout_type: str = "grid",
        height: int = 800,
    ) -> Union[go.Figure, html.Div]:
        """Create combined plot of selected prediction types.

        Args:
            predictions: Dictionary of predictions by type
            selected_types: List of prediction types to plot
            layout_type: Type of layout (grid, tabs, or combined)
            height: Plot height in pixels

        Returns:
            Combined plot as Plotly figure or Dash div
        """
        try:
            if not selected_types:
                return html.Div("No plot types selected")

            if layout_type == "tabs":
                return self._create_tabbed_layout(predictions, selected_types, height)
            elif layout_type == "combined":
                return self._create_combined_layout(predictions, selected_types, height)
            else:
                return self._create_grid_layout(predictions, selected_types, height)

        except Exception as e:
            self.logger.error(f"Error creating combined plot: {str(e)}")
            return html.Div(f"Error creating plot: {str(e)}")

    def create_web_visualization(
        self,
        web_data: Dict,
        selected_sources: List[str],
        viz_type: str = "timeline",
        height: int = 600,
    ) -> html.Div:
        """Create web data visualization.

        Args:
            web_data: Dictionary of web data by source
            selected_sources: List of data sources to visualize
            viz_type: Type of visualization
            height: Plot height in pixels

        Returns:
            Dash div with visualization
        """
        try:
            if not selected_sources:
                return html.Div("No data sources selected")

            # Filter data by selected sources
            filtered_data = {k: v for k, v in web_data.items() if k in selected_sources}

            # Create visualization based on type
            if viz_type == "timeline":
                return self._create_timeline_viz(filtered_data, height)
            elif viz_type == "network":
                return self._create_network_viz(filtered_data, height)
            elif viz_type == "heatmap":
                return self._create_heatmap_viz(filtered_data, height)
            elif viz_type == "wordcloud":
                return self._create_wordcloud_viz(filtered_data, height)
            else:
                return html.Div(f"Unknown visualization type: {viz_type}")

        except Exception as e:
            self.logger.error(f"Error creating web visualization: {str(e)}")
            return html.Div(f"Error creating visualization: {str(e)}")

    @lru_cache(maxsize=32)
    def _create_cached_plot(
        self,
        plot_type: str,
        data_key: str,
        data: Dict,
    ) -> go.Figure:
        """Create cached plot of specified type.

        Args:
            plot_type: Type of plot to create
            data_key: Key for cache invalidation
            data: Plot data

        Returns:
            Plotly figure
        """
        return self._create_plot(plot_type, data)

    def _create_grid_layout(
        self, predictions: Dict, selected_types: List[str], height: int
    ) -> go.Figure:
        """Create grid layout of plots.

        Args:
            predictions: Dictionary of predictions by type
            selected_types: List of prediction types to plot
            height: Plot height in pixels

        Returns:
            Combined Plotly figure
        """
        # Calculate grid dimensions
        n_plots = len(selected_types)
        n_rows = int(np.ceil(n_plots / 2))
        n_cols = min(2, n_plots)

        # Create figure with subplots
        fig = go.Figure()

        # Add each plot type
        for i, plot_type in enumerate(selected_types):
            if plot_type in predictions:
                row = (i // n_cols) + 1
                col = (i % n_cols) + 1

                subplot = self._create_cached_plot(
                    plot_type,
                    f"{plot_type}_{hash(str(predictions[plot_type]))}",
                    predictions[plot_type],
                )
                for trace in subplot.data:
                    trace.update(
                        xaxis=f"x{i + 1}" if i > 0 else "x",
                        yaxis=f"y{i + 1}" if i > 0 else "y",
                    )
                    fig.add_trace(trace)

        # Update layout
        fig.update_layout(
            title="ML Predictions Overview",
            showlegend=True,
            height=height,
            grid={"rows": n_rows, "columns": n_cols, "pattern": "independent"},
            template="plotly_white",
        )

        return fig

    def _create_tabbed_layout(
        self, predictions: Dict, selected_types: List[str], height: int
    ) -> html.Div:
        """Create tabbed layout of plots.

        Args:
            predictions: Dictionary of predictions by type
            selected_types: List of prediction types to plot
            height: Plot height in pixels

        Returns:
            Dash div with tabs
        """
        tabs = []
        for plot_type in selected_types:
            if plot_type in predictions:
                fig = self._create_cached_plot(
                    plot_type,
                    f"{plot_type}_{hash(str(predictions[plot_type]))}",
                    predictions[plot_type],
                )
                fig.update_layout(height=height)

                tabs.append(
                    dcc.Tab(
                        label=plot_type.title(),
                        children=[
                            dcc.Graph(
                                figure=fig,
                                config={"displayModeBar": True},
                            )
                        ],
                    )
                )

        return html.Div([dcc.Tabs(tabs)])

    def _create_combined_layout(
        self, predictions: Dict, selected_types: List[str], height: int
    ) -> html.Div:
        """Create combined layout of plots.

        Args:
            predictions: Dictionary of predictions by type
            selected_types: List of prediction types to plot
            height: Plot height in pixels

        Returns:
            Dash div with combined visualization
        """
        try:
            # Create main figure
            fig = go.Figure()

            # Add data from each prediction type
            for plot_type in selected_types:
                if plot_type in predictions:
                    subplot = self._create_cached_plot(
                        plot_type,
                        f"{plot_type}_{hash(str(predictions[plot_type]))}",
                        predictions[plot_type],
                    )
                    for trace in subplot.data:
                        fig.add_trace(trace)

            # Update layout
            fig.update_layout(
                title="Combined ML Predictions",
                showlegend=True,
                height=height,
                template="plotly_white",
                margin=dict(l=50, r=50, t=50, b=50),
            )

            return html.Div([dcc.Graph(figure=fig)])

        except Exception as e:
            self.logger.error(f"Error creating combined layout: {str(e)}")
            return html.Div(f"Error creating combined layout: {str(e)}")

    def _create_plot(self, plot_type: str, data: Dict) -> go.Figure:
        """Create plot of specified type.

        Args:
            plot_type: Type of plot to create
            data: Plot data

        Returns:
            Plotly figure
        """
        try:
            # Get plot creation function
            plot_func = getattr(self, f"_create_{plot_type}_plot", None)
            if plot_func is None:
                return self._create_default_plot(data)

            # Create plot
            fig = plot_func(data)

            # Update layout with common settings
            fig.update_layout(
                template="plotly_white",
                margin=dict(l=50, r=50, t=50, b=50),
                colorway=px.colors.qualitative.Set3,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating {plot_type} plot: {str(e)}")
            return go.Figure()

    def _create_binding_plot(self, data: Dict) -> go.Figure:
        """Create binding affinity prediction plot."""
        return self.ml_visualizer.create_binding_plot(
            data, colorscale=self.COLORSCALES["binding"]
        )

    def _create_abuse_plot(self, data: Dict) -> go.Figure:
        """Create abuse potential prediction plot."""
        return self.ml_visualizer.create_abuse_plot(
            data, colorscale=self.COLORSCALES["abuse"]
        )

    def _create_toxicity_plot(self, data: Dict) -> go.Figure:
        """Create toxicity prediction plot."""
        return self.ml_visualizer.create_toxicity_plot(
            data, colorscale=self.COLORSCALES["toxicity"]
        )

    def _create_activity_plot(self, data: Dict) -> go.Figure:
        """Create activity prediction plot."""
        return self.ml_visualizer.create_activity_plot(
            data, colorscale=self.COLORSCALES["activity"]
        )

    def _create_psychoactive_plot(self, data: Dict) -> go.Figure:
        """Create psychoactive effects prediction plot."""
        return self.ml_visualizer.create_psychoactive_plot(
            data, colorscale=self.COLORSCALES["psychoactive"]
        )

    def _create_nootropic_plot(self, data: Dict) -> go.Figure:
        """Create nootropic effects prediction plot."""
        return self.ml_visualizer.create_nootropic_plot(
            data, colorscale=self.COLORSCALES["nootropic"]
        )

    def _create_interaction_plot(self, data: Dict) -> go.Figure:
        """Create interaction network plot."""
        return self.advanced_visualizer.create_interaction_network(
            data, colorscale=self.COLORSCALES["interactions"]
        )

    def _create_similarity_plot(self, data: Dict) -> go.Figure:
        """Create structural similarity plot."""
        return self.advanced_visualizer.create_similarity_matrix(
            data, colorscale=self.COLORSCALES["similarity"]
        )

    def _create_uncertainty_plot(self, data: Dict) -> go.Figure:
        """Create prediction uncertainty plot."""
        return self.advanced_visualizer.create_uncertainty_plot(
            data, colorscale=self.COLORSCALES["uncertainty"]
        )

    def _create_timeline_viz(self, data: Dict, height: int) -> html.Div:
        """Create timeline visualization of web data.

        Args:
            data: Web data by source
            height: Plot height in pixels

        Returns:
            Dash div with timeline visualization
        """
        try:
            # Extract timeline data
            events = []
            for source, source_data in data.items():
                if isinstance(source_data, dict):
                    for item in source_data.get("timeline", []):
                        events.append(
                            {
                                "source": source,
                                "date": item.get("date"),
                                "type": item.get("type"),
                                "description": item.get("description"),
                            }
                        )

            # Create timeline plot
            fig = go.Figure()

            # Add events as scatter points
            for event in events:
                fig.add_trace(
                    go.Scatter(
                        x=[event["date"]],
                        y=[event["source"]],
                        mode="markers+text",
                        name=event["type"],
                        text=event["description"],
                        textposition="top center",
                    )
                )

            # Update layout
            fig.update_layout(
                title="Timeline of Web Mentions",
                xaxis_title="Date",
                yaxis_title="Source",
                height=height,
                showlegend=True,
                template="plotly_white",
            )

            return html.Div([dcc.Graph(figure=fig)])

        except Exception as e:
            self.logger.error(f"Error creating timeline visualization: {str(e)}")
            return html.Div(f"Error creating timeline: {str(e)}")

    def _create_network_viz(self, data: Dict, height: int) -> html.Div:
        """Create network visualization of web data.

        Args:
            data: Web data by source
            height: Plot height in pixels

        Returns:
            Dash div with network visualization
        """
        try:
            # Create network graph
            G = nx.Graph()

            # Add nodes and edges from data
            for source, source_data in data.items():
                if isinstance(source_data, dict):
                    for connection in source_data.get("connections", []):
                        G.add_edge(
                            connection["from"],
                            connection["to"],
                            weight=connection.get("weight", 1),
                            type=connection.get("type", "unknown"),
                        )

            # Create network plot
            pos = nx.spring_layout(G)
            edge_x = []
            edge_y = []
            for edge in G.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

            fig = go.Figure()

            # Add edges
            fig.add_trace(
                go.Scatter(
                    x=edge_x,
                    y=edge_y,
                    mode="lines",
                    line=dict(width=0.5, color="#888"),
                    hoverinfo="none",
                )
            )

            # Add nodes
            node_x = [pos[node][0] for node in G.nodes()]
            node_y = [pos[node][1] for node in G.nodes()]
            node_text = list(G.nodes())

            fig.add_trace(
                go.Scatter(
                    x=node_x,
                    y=node_y,
                    mode="markers+text",
                    text=node_text,
                    textposition="top center",
                    marker=dict(
                        size=10,
                        color=list(range(len(G.nodes()))),
                        colorscale=self.COLORSCALES["network"],
                        line_width=2,
                    ),
                )
            )

            # Update layout
            fig.update_layout(
                title="Network of Related Compounds",
                showlegend=False,
                height=height,
                template="plotly_white",
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            )

            return html.Div([dcc.Graph(figure=fig)])

        except Exception as e:
            self.logger.error(f"Error creating network visualization: {str(e)}")
            return html.Div(f"Error creating network: {str(e)}")

    def _create_heatmap_viz(self, data: Dict, height: int) -> html.Div:
        """Create heatmap visualization of web data.

        Args:
            data: Web data by source
            height: Plot height in pixels

        Returns:
            Dash div with heatmap visualization
        """
        try:
            # Extract heatmap data
            categories = set()
            sources = list(data.keys())
            values = {}

            for source, source_data in data.items():
                if isinstance(source_data, dict):
                    for category, value in source_data.get("metrics", {}).items():
                        categories.add(category)
                        values[(source, category)] = value

            categories = sorted(categories)

            # Create heatmap matrix
            matrix = []
            for source in sources:
                row = []
                for category in categories:
                    row.append(values.get((source, category), 0))
                matrix.append(row)

            # Create heatmap plot
            fig = go.Figure(
                data=go.Heatmap(
                    z=matrix,
                    x=categories,
                    y=sources,
                    colorscale=self.COLORSCALES["heatmap"],
                )
            )

            # Update layout
            fig.update_layout(
                title="Data Source Metrics",
                height=height,
                template="plotly_white",
            )

            return html.Div([dcc.Graph(figure=fig)])

        except Exception as e:
            self.logger.error(f"Error creating heatmap visualization: {str(e)}")
            return html.Div(f"Error creating heatmap: {str(e)}")

    def _create_wordcloud_viz(self, data: Dict, height: int) -> html.Div:
        """Create word cloud visualization of web data.

        Args:
            data: Web data by source
            height: Plot height in pixels

        Returns:
            Dash div with word cloud visualization
        """
        try:
            # Extract text data
            text = ""
            for source_data in data.values():
                if isinstance(source_data, dict):
                    # Add mentions
                    for mention in source_data.get("mentions", []):
                        text += f" {mention['text']}"
                    # Add descriptions
                    for desc in source_data.get("descriptions", []):
                        text += f" {desc}"
                    # Add effects
                    for effect in source_data.get("effects", []):
                        text += f" {effect}"

            # Generate word cloud
            if text:
                wordcloud = WordCloud(
                    width=800,
                    height=height,
                    background_color="white",
                    colormap=self.COLORSCALES["wordcloud"],
                ).generate(text)

                # Convert to base64 image
                import io
                import base64

                img = io.BytesIO()
                wordcloud.to_image().save(img, format="PNG")
                img_str = base64.b64encode(img.getvalue()).decode()

                return html.Div(
                    [
                        html.Img(
                            src=f"data:image/png;base64,{img_str}",
                            style={"height": f"{height}px", "width": "100%"},
                        )
                    ]
                )
            else:
                return html.Div("No text data available for word cloud")

        except Exception as e:
            self.logger.error(f"Error creating word cloud visualization: {str(e)}")
            return html.Div(f"Error creating word cloud: {str(e)}")

    def _create_default_plot(self, data: Dict) -> go.Figure:
        """Create default plot for unknown types."""
        fig = go.Figure()

        # Extract values and labels
        values = []
        labels = []
        for key, val in data.items():
            if isinstance(val, dict) and "value" in val:
                values.append(val["value"])
                labels.append(key)
            elif isinstance(val, (int, float)):
                values.append(val)
                labels.append(key)

        if values and labels:
            fig.add_trace(
                go.Scatter(
                    x=labels,
                    y=values,
                    mode="lines+markers",
                    name="Values",
                )
            )

            fig.update_layout(
                title="Data Values",
                xaxis_title="Category",
                yaxis_title="Value",
                showlegend=True,
            )

        return fig
