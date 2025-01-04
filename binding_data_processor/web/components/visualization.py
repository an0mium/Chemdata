"""Visualization panel component for the web interface."""

import logging
from typing import Dict, List, Optional
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import pandas as pd
import numpy as np

from ...models.compound import CompoundData
from ...processors.structure.depiction import StructureDepiction
from ...processors.structure.ml.visualization import MLVisualizer


class VisualizationPanel:
    """Component for data visualization and plotting."""

    def __init__(
        self,
        parent_app: dash.Dash,
        structure_depiction: Optional[StructureDepiction] = None,
        ml_visualizer: Optional[MLVisualizer] = None,
    ):
        """Initialize visualization panel component.

        Args:
            parent_app: Parent Dash application
            structure_depiction: Optional structure depiction handler
            ml_visualizer: Optional ML visualization handler
        """
        self.app = parent_app
        self.structure_depiction = structure_depiction or StructureDepiction()
        self.ml_visualizer = ml_visualizer or MLVisualizer()
        self.logger = logging.getLogger(__name__)
        self._setup_callbacks()

    def get_layout(self) -> html.Div:
        """Get component layout."""
        return html.Div(
            [
                html.H3("Visualization"),
                dcc.Tabs(
                    [
                        dcc.Tab(
                            label="Activity Distribution",
                            children=[
                                dcc.Graph(id="activity-dist-plot"),
                                html.Label("Plot Type"),
                                dcc.Dropdown(
                                    id="activity-plot-type",
                                    options=[
                                        {"label": "Histogram", "value": "histogram"},
                                        {"label": "Box Plot", "value": "box"},
                                        {"label": "Violin Plot", "value": "violin"},
                                    ],
                                    value="histogram",
                                ),
                                html.Label("Group By"),
                                dcc.Dropdown(
                                    id="activity-group-by",
                                    options=[
                                        {
                                            "label": "Activity Type",
                                            "value": "activity_type",
                                        },
                                        {"label": "Target", "value": "target"},
                                        {"label": "Source", "value": "source"},
                                    ],
                                    value="activity_type",
                                ),
                            ],
                        ),
                        dcc.Tab(
                            label="Structure Analysis",
                            children=[
                                dcc.Graph(id="structure-plot"),
                                html.Label("Analysis Type"),
                                dcc.Dropdown(
                                    id="structure-analysis-type",
                                    options=[
                                        {
                                            "label": "Similarity Matrix",
                                            "value": "similarity",
                                        },
                                        {"label": "Clustering", "value": "clustering"},
                                        {"label": "MCS Analysis", "value": "mcs"},
                                    ],
                                    value="similarity",
                                ),
                                html.Label("Color By"),
                                dcc.Dropdown(
                                    id="structure-color-by",
                                    options=[
                                        {"label": "Activity", "value": "activity"},
                                        {"label": "Target", "value": "target"},
                                        {"label": "Cluster", "value": "cluster"},
                                    ],
                                    value="activity",
                                ),
                            ],
                        ),
                        dcc.Tab(
                            label="SAR Analysis",
                            children=[
                                dcc.Graph(id="sar-plot"),
                                html.Label("Plot Type"),
                                dcc.Dropdown(
                                    id="sar-plot-type",
                                    options=[
                                        {
                                            "label": "Activity vs Property",
                                            "value": "property",
                                        },
                                        {"label": "Activity Cliff", "value": "cliff"},
                                        {
                                            "label": "R-Group Analysis",
                                            "value": "rgroup",
                                        },
                                    ],
                                    value="property",
                                ),
                                html.Label("Property"),
                                dcc.Dropdown(
                                    id="sar-property",
                                    options=[
                                        {"label": "MW", "value": "MW"},
                                        {"label": "LogP", "value": "LogP"},
                                        {"label": "TPSA", "value": "TPSA"},
                                        {"label": "HBA", "value": "HBA"},
                                        {"label": "HBD", "value": "HBD"},
                                    ],
                                    value="LogP",
                                ),
                            ],
                        ),
                    ]
                ),
                html.Div(id="plot-error"),
            ],
            style={"padding": "20px"},
        )

    def _setup_callbacks(self):
        """Setup component callbacks."""

        @self.app.callback(
            Output("activity-dist-plot", "figure"),
            [
                Input("activity-plot-type", "value"),
                Input("activity-group-by", "value"),
            ],
            [State("stored-data", "data")],
        )
        def update_activity_plot(
            plot_type: str, group_by: str, data: List[Dict]
        ) -> go.Figure:
            """Update activity distribution plot."""
            if not data:
                return go.Figure()

            try:
                df = pd.DataFrame(data)

                if plot_type == "histogram":
                    fig = go.Figure()
                    for group in df[group_by].unique():
                        group_data = df[df[group_by] == group]
                        fig.add_trace(
                            go.Histogram(
                                x=group_data["activity_value"],
                                name=group,
                                opacity=0.7,
                            )
                        )
                    fig.update_layout(barmode="overlay")

                elif plot_type == "box":
                    fig = go.Figure(
                        go.Box(
                            x=df[group_by],
                            y=df["activity_value"],
                            boxpoints="all",
                        )
                    )

                elif plot_type == "violin":
                    fig = go.Figure(
                        go.Violin(
                            x=df[group_by],
                            y=df["activity_value"],
                            box_visible=True,
                            points="all",
                        )
                    )

                fig.update_layout(
                    title=f"Activity Distribution by {group_by}",
                    xaxis_title=group_by.replace("_", " ").title(),
                    yaxis_title="Activity (nM)",
                )
                return fig

            except Exception as e:
                self.logger.error(f"Error updating activity plot: {str(e)}")
                return go.Figure()

        @self.app.callback(
            Output("structure-plot", "figure"),
            [
                Input("structure-analysis-type", "value"),
                Input("structure-color-by", "value"),
            ],
            [State("stored-data", "data")],
        )
        def update_structure_plot(
            analysis_type: str, color_by: str, data: List[Dict]
        ) -> go.Figure:
            """Update structure analysis plot."""
            if not data:
                return go.Figure()

            try:
                df = pd.DataFrame(data)

                if analysis_type == "similarity":
                    # Create similarity matrix heatmap
                    similarity_matrix = (
                        self.structure_depiction.create_similarity_matrix(
                            [CompoundData(**d) for d in data]
                        )
                    )
                    fig = go.Figure(
                        go.Heatmap(
                            z=similarity_matrix,
                            colorscale="Viridis",
                        )
                    )
                    fig.update_layout(
                        title="Structural Similarity Matrix",
                        xaxis_title="Compound Index",
                        yaxis_title="Compound Index",
                    )

                elif analysis_type == "clustering":
                    # Create clustering plot
                    clusters = self.ml_visualizer.plot_clusters(
                        [CompoundData(**d) for d in data],
                        color_by=color_by,
                    )
                    fig = go.Figure(data=clusters)
                    fig.update_layout(
                        title=f"Structure Clustering (colored by {color_by})",
                        xaxis_title="UMAP 1",
                        yaxis_title="UMAP 2",
                    )

                elif analysis_type == "mcs":
                    # Create MCS analysis plot
                    mcs_data = self.structure_depiction.analyze_mcs(
                        [CompoundData(**d) for d in data]
                    )
                    fig = go.Figure(data=mcs_data)
                    fig.update_layout(
                        title="Maximum Common Substructure Analysis",
                        xaxis_title="Substructure Size",
                        yaxis_title="Frequency",
                    )

                return fig

            except Exception as e:
                self.logger.error(f"Error updating structure plot: {str(e)}")
                return go.Figure()

        @self.app.callback(
            Output("sar-plot", "figure"),
            [
                Input("sar-plot-type", "value"),
                Input("sar-property", "value"),
            ],
            [State("stored-data", "data")],
        )
        def update_sar_plot(
            plot_type: str, property_name: str, data: List[Dict]
        ) -> go.Figure:
            """Update SAR analysis plot."""
            if not data:
                return go.Figure()

            try:
                df = pd.DataFrame(data)

                if plot_type == "property":
                    # Create property vs activity scatter plot
                    fig = go.Figure(
                        go.Scatter(
                            x=df[f"descriptors.{property_name}"],
                            y=df["activity_value"],
                            mode="markers",
                            text=df["name"],
                            marker=dict(
                                size=10,
                                color=df["activity_value"],
                                colorscale="Viridis",
                                showscale=True,
                            ),
                        )
                    )
                    fig.update_layout(
                        title=f"Activity vs {property_name}",
                        xaxis_title=property_name,
                        yaxis_title="Activity (nM)",
                    )

                elif plot_type == "cliff":
                    # Create activity cliff plot
                    cliff_data = self.ml_visualizer.plot_activity_cliffs(
                        [CompoundData(**d) for d in data]
                    )
                    fig = go.Figure(data=cliff_data)
                    fig.update_layout(
                        title="Activity Cliff Analysis",
                        xaxis_title="Structural Similarity",
                        yaxis_title="Activity Ratio",
                    )

                elif plot_type == "rgroup":
                    # Create R-group analysis plot
                    rgroup_data = self.structure_depiction.analyze_rgroups(
                        [CompoundData(**d) for d in data]
                    )
                    fig = go.Figure(data=rgroup_data)
                    fig.update_layout(
                        title="R-Group Analysis",
                        xaxis_title="R-Group",
                        yaxis_title="Activity Change",
                    )

                return fig

            except Exception as e:
                self.logger.error(f"Error updating SAR plot: {str(e)}")
                return go.Figure()
