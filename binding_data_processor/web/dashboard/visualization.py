"""Visualization components for compound data.

This module provides interactive visualization components for:
1. 2D/3D molecular structure viewing
2. Activity and prediction visualization
3. Interactive compound data tables
4. ML prediction visualization
5. Structure-activity relationship plots
"""

import json
from typing import Dict, List, Optional, Union

import dash_bio as dashbio
import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, dash_table
import numpy as np

from binding_data_processor.models.compound import CompoundData
from binding_data_processor.processors.structure.ml.visualization import MLVisualizer
from binding_data_processor.processors.structure.ml.advanced_visualization import (
    AdvancedVisualizer,
)


class StructureViewer(html.Div):
    """2D/3D structure viewer component."""

    def __init__(
        self,
        id: str,
        smiles: str,
        view_type: str = "2d",
        height: int = 300,
        width: Optional[int] = None,
    ):
        """Initialize structure viewer.

        Args:
            id: Component ID
            smiles: SMILES string to display
            view_type: Type of view (2d or 3d)
            height: Viewer height in pixels
            width: Optional viewer width in pixels
        """
        styles = {"height": f"{height}px"}
        if width:
            styles["width"] = f"{width}px"

        if view_type == "3d":
            viewer = dashbio.Molecule3dViewer(
                id=f"{id}-3d",
                modelData={"smiles": smiles},
                styles=styles,
            )
        else:
            viewer = dashbio.Molecule2dViewer(
                id=f"{id}-2d",
                modelData={"smiles": smiles},
                styles=styles,
            )

        super().__init__(viewer)


class ActivityPlot(html.Div):
    """Activity prediction visualization component."""

    def __init__(
        self,
        id: str,
        data: Optional[Dict] = None,
        plot_type: str = "bar",
        height: int = 300,
    ):
        """Initialize activity plot.

        Args:
            id: Component ID
            data: Activity prediction data
            plot_type: Type of plot (bar, radar, or heatmap)
            height: Plot height in pixels
        """
        if not data:
            super().__init__("No activity predictions available")
            return

        # Create figure based on plot type
        if plot_type == "radar":
            fig = self._create_radar_plot(data)
        elif plot_type == "heatmap":
            fig = self._create_heatmap_plot(data)
        else:
            fig = self._create_bar_plot(data)

        # Update layout
        fig.update_layout(
            margin=dict(l=50, r=50, t=50, b=50),
            height=height,
            template="plotly_white",
        )

        super().__init__(
            dcc.Graph(
                id=id,
                figure=fig,
                config={"displayModeBar": False},
            )
        )

    def _create_bar_plot(self, data: Dict) -> go.Figure:
        """Create bar plot of activities."""
        activities = []
        probabilities = []
        for activity, prob in data.items():
            activities.append(activity)
            probabilities.append(prob)

        return go.Figure(
            data=[
                go.Bar(
                    x=activities,
                    y=probabilities,
                    marker_color=px.colors.qualitative.Set3,
                )
            ],
            layout=dict(
                title="Predicted Activities",
                xaxis_title="Activity Type",
                yaxis_title="Probability",
                yaxis_range=[0, 1],
                showlegend=False,
            ),
        )

    def _create_radar_plot(self, data: Dict) -> go.Figure:
        """Create radar plot of activities."""
        return go.Figure(
            data=[
                go.Scatterpolar(
                    r=list(data.values()),
                    theta=list(data.keys()),
                    fill="toself",
                )
            ],
            layout=dict(
                title="Activity Profile",
                polar=dict(radialaxis=dict(range=[0, 1])),
                showlegend=False,
            ),
        )

    def _create_heatmap_plot(self, data: Dict) -> go.Figure:
        """Create heatmap of activities."""
        return go.Figure(
            data=[
                go.Heatmap(
                    z=[list(data.values())],
                    x=list(data.keys()),
                    y=["Activity"],
                    colorscale="Viridis",
                )
            ],
            layout=dict(
                title="Activity Heatmap",
                xaxis_title="Activity Type",
                showlegend=False,
            ),
        )


class PredictionPlot(html.Div):
    """ML prediction visualization component."""

    def __init__(
        self,
        id: str,
        data: Optional[Dict] = None,
        plot_type: str = "gauge",
        height: int = 300,
        show_details: bool = True,
    ):
        """Initialize prediction plot.

        Args:
            id: Component ID
            data: Prediction data
            plot_type: Type of plot (gauge or bar)
            height: Plot height in pixels
            show_details: Whether to show detailed prediction data
        """
        if not data:
            super().__init__("No predictions available")
            return

        # Create figure based on plot type
        if plot_type == "bar":
            fig = self._create_bar_plot(data)
        else:
            fig = self._create_gauge_plot(data)

        # Update layout
        fig.update_layout(
            margin=dict(l=50, r=50, t=50, b=50),
            height=height,
            template="plotly_white",
        )

        # Create layout
        layout = [
            dcc.Graph(
                id=id,
                figure=fig,
                config={"displayModeBar": False},
            )
        ]

        # Add details if requested
        if show_details:
            layout.append(
                html.Div(
                    [
                        html.H5("Details"),
                        html.Pre(
                            json.dumps(
                                {k: v for k, v in data.items() if k not in ["score"]},
                                indent=2,
                            )
                        ),
                    ],
                    style={"marginTop": "1rem"},
                )
            )

        super().__init__(layout)

    def _create_gauge_plot(self, data: Dict) -> go.Figure:
        """Create gauge plot of prediction score."""
        return go.Figure(
            go.Indicator(
                mode="gauge+number",
                value=data.get("score", 0),
                domain={"x": [0, 1], "y": [0, 1]},
                gauge={
                    "axis": {"range": [0, 1]},
                    "bar": {"color": "rgb(55, 83, 109)"},
                    "steps": [
                        {"range": [0, 0.33], "color": "rgb(198, 246, 213)"},
                        {"range": [0.33, 0.67], "color": "rgb(254, 235, 200)"},
                        {"range": [0.67, 1], "color": "rgb(254, 215, 215)"},
                    ],
                },
            )
        )

    def _create_bar_plot(self, data: Dict) -> go.Figure:
        """Create bar plot of prediction scores."""
        return go.Figure(
            data=[
                go.Bar(
                    x=list(data.keys()),
                    y=list(data.values()),
                    marker_color=px.colors.qualitative.Set3,
                )
            ],
            layout=dict(
                title="Prediction Scores",
                xaxis_title="Category",
                yaxis_title="Score",
                yaxis_range=[0, 1],
                showlegend=False,
            ),
        )


class CompoundTable(html.Div):
    """Interactive compound data table component."""

    def __init__(
        self,
        id: str,
        data: Optional[List[Dict]] = None,
        columns: Optional[List[Dict]] = None,
        page_size: int = 20,
    ):
        """Initialize compound table.

        Args:
            id: Component ID
            data: List of compound data dictionaries
            columns: Optional custom column definitions
            page_size: Number of rows per page
        """
        if columns is None:
            columns = [
                {"name": "Name", "id": "name"},
                {"name": "SMILES", "id": "smiles"},
                {"name": "Source", "id": "source"},
                {"name": "Toxicity", "id": "toxicity_score"},
                {"name": "Abuse Potential", "id": "abuse_potential"},
                {"name": "Targets", "id": "binding_targets"},
                {"name": "Activities", "id": "activity_types"},
            ]

        super().__init__(
            dash_table.DataTable(
                id=id,
                data=data or [],
                columns=columns,
                style_table={"overflowX": "auto"},
                style_cell={
                    "textAlign": "left",
                    "padding": "15px",
                    "whiteSpace": "normal",
                    "height": "auto",
                },
                style_header={
                    "backgroundColor": "rgb(230, 230, 230)",
                    "fontWeight": "bold",
                },
                style_data_conditional=[
                    {
                        "if": {"row_index": "odd"},
                        "backgroundColor": "rgb(248, 248, 248)",
                    }
                ],
                sort_action="native",
                filter_action="native",
                page_size=page_size,
            )
        )
