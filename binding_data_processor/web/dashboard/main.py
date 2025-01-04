"""Main dashboard component integrating all visualization components."""

import logging
from typing import Dict, List, Optional
from functools import lru_cache

import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from binding_data_processor.models.compound import CompoundData
from binding_data_processor.web.dashboard.list_view import CompoundListView
from binding_data_processor.web.dashboard.detail_view import CompoundDetailView
from binding_data_processor.web.dashboard.plot_manager import PlotManager
from binding_data_processor.web.dashboard.base import BaseDashboard


class MainDashboard(BaseDashboard):
    """Main dashboard integrating all visualization components."""

    def __init__(self):
        """Initialize main dashboard."""
        super().__init__()
        self.list_view = CompoundListView()
        self.detail_view = CompoundDetailView()
        self.plot_manager = PlotManager()
        self.compounds: List[CompoundData] = []
        self.selected_compound: Optional[CompoundData] = None

        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        self.logger.addHandler(handler)

    def update_data(self, compounds: List[CompoundData]) -> None:
        """Update dashboard with new compound data.

        Args:
            compounds: List of compounds to display
        """
        try:
            self.logger.info(f"Updating dashboard with {len(compounds)} compounds")
            self.compounds = compounds
            self.list_view.update_data(compounds)
            if self.selected_compound:
                self.detail_view.update_data(self.selected_compound)

            # Clear plot cache when data changes
            self._create_cached_plot.cache_clear()

        except Exception as e:
            self.logger.error(f"Error updating dashboard data: {str(e)}")

    def layout(self) -> html.Div:
        """Get dashboard layout.

        Returns:
            Dash layout
        """
        return html.Div(
            [
                # Header
                html.Div(
                    [
                        html.H1("ChemData Dashboard"),
                        html.P(
                            "Explore and analyze psychopharmacological compounds",
                            className="lead",
                        ),
                    ],
                    className="p-5 bg-light rounded-3",
                ),
                # Main content
                dbc.Container(
                    [
                        # Navigation tabs
                        dcc.Tabs(
                            id="main-tabs",
                            value="list",
                            children=[
                                dcc.Tab(
                                    label="Compound List",
                                    value="list",
                                    children=self.list_view.layout(),
                                ),
                                dcc.Tab(
                                    label="Compound Details",
                                    value="details",
                                    children=self.detail_view.layout(),
                                ),
                                dcc.Tab(
                                    label="ML Predictions",
                                    value="predictions",
                                    children=self._create_predictions_layout(),
                                ),
                                dcc.Tab(
                                    label="Web Data",
                                    value="web_data",
                                    children=self._create_web_data_layout(),
                                ),
                            ],
                        ),
                        # Progress and status
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Progress(
                                            id="progress-bar",
                                            style={"visibility": "hidden"},
                                        ),
                                        html.P(
                                            id="progress-text",
                                            style={"visibility": "hidden"},
                                        ),
                                    ],
                                    width=12,
                                ),
                            ],
                            className="mt-3",
                        ),
                        # Status and info
                        html.Div(
                            [
                                html.P(id="status-message"),
                                html.P(id="info-message", className="text-muted"),
                            ],
                            className="mt-3",
                        ),
                    ],
                    fluid=True,
                    className="py-3",
                ),
            ]
        )

    def _create_predictions_layout(self) -> html.Div:
        """Create ML predictions tab layout.

        Returns:
            Dash layout
        """
        return html.Div(
            [
                # Plot type selection
                html.Div(
                    [
                        html.H4("Select Prediction Types"),
                        dcc.Checklist(
                            id="plot-types",
                            options=[
                                {"label": "Binding Affinity", "value": "binding"},
                                {"label": "Abuse Potential", "value": "abuse"},
                                {"label": "Toxicity", "value": "toxicity"},
                                {"label": "Activity", "value": "activity"},
                                {"label": "Psychoactive", "value": "psychoactive"},
                                {"label": "Nootropic", "value": "nootropic"},
                                {"label": "Interactions", "value": "interactions"},
                                {"label": "Similarity", "value": "similarity"},
                                {"label": "Uncertainty", "value": "uncertainty"},
                            ],
                            value=["binding", "abuse", "toxicity"],
                        ),
                    ],
                    className="mb-4",
                ),
                # Layout type selection
                html.Div(
                    [
                        html.H4("Layout Type"),
                        dcc.RadioItems(
                            id="layout-type",
                            options=[
                                {"label": "Grid", "value": "grid"},
                                {"label": "Tabs", "value": "tabs"},
                                {"label": "Combined", "value": "combined"},
                            ],
                            value="grid",
                        ),
                    ],
                    className="mb-4",
                ),
                # Plot container
                html.Div(id="plot-container"),
            ]
        )

    def _create_web_data_layout(self) -> html.Div:
        """Create web data tab layout.

        Returns:
            Dash layout
        """
        return html.Div(
            [
                # Data source selection
                html.Div(
                    [
                        html.H4("Select Data Sources"),
                        dcc.Checklist(
                            id="data-sources",
                            options=[
                                {"label": "Community Data", "value": "community"},
                                {"label": "Social Media", "value": "social"},
                                {"label": "Swiss Data", "value": "swiss"},
                                {"label": "Experience Reports", "value": "reports"},
                                {"label": "Sentiment Analysis", "value": "sentiment"},
                            ],
                            value=["community", "social"],
                        ),
                    ],
                    className="mb-4",
                ),
                # Visualization type
                html.Div(
                    [
                        html.H4("Visualization Type"),
                        dcc.RadioItems(
                            id="viz-type",
                            options=[
                                {"label": "Timeline", "value": "timeline"},
                                {"label": "Network", "value": "network"},
                                {"label": "Heatmap", "value": "heatmap"},
                                {"label": "Word Cloud", "value": "wordcloud"},
                            ],
                            value="timeline",
                        ),
                    ],
                    className="mb-4",
                ),
                # Visualization container
                html.Div(id="viz-container"),
            ]
        )

    def register_callbacks(self) -> None:
        """Register dashboard callbacks."""
        super().register_callbacks()
        self.list_view.register_callbacks()
        self.detail_view.register_callbacks()

        @self.app.callback(
            Output("plot-container", "children"),
            Input("plot-types", "value"),
            Input("layout-type", "value"),
            prevent_initial_call=True,
        )
        def update_plots(selected_types: List[str], layout_type: str) -> html.Div:
            """Update prediction plots based on selection.

            Args:
                selected_types: List of selected plot types
                layout_type: Type of layout (grid, tabs, or combined)

            Returns:
                Updated plot container
            """
            if not selected_types:
                return html.Div("Please select at least one prediction type")

            if not self.selected_compound:
                return html.Div("Please select a compound to view predictions")

            try:
                return self._create_cached_plot(
                    tuple(selected_types),
                    layout_type,
                    self.selected_compound.name,
                )

            except Exception as e:
                self.logger.error(f"Error updating plots: {str(e)}")
                return html.Div(f"Error creating plots: {str(e)}")

        @self.app.callback(
            Output("viz-container", "children"),
            Input("data-sources", "value"),
            Input("viz-type", "value"),
            prevent_initial_call=True,
        )
        def update_web_viz(selected_sources: List[str], viz_type: str) -> html.Div:
            """Update web data visualization based on selection.

            Args:
                selected_sources: List of selected data sources
                viz_type: Type of visualization

            Returns:
                Updated visualization container
            """
            if not selected_sources:
                return html.Div("Please select at least one data source")

            if not self.selected_compound:
                return html.Div("Please select a compound to view web data")

            try:
                return self._create_cached_web_viz(
                    tuple(selected_sources),
                    viz_type,
                    self.selected_compound.name,
                )

            except Exception as e:
                self.logger.error(f"Error updating web visualization: {str(e)}")
                return html.Div(f"Error creating visualization: {str(e)}")

        @self.app.callback(
            [
                Output("detail-view", "data"),
                Output("status-message", "children"),
                Output("info-message", "children"),
            ],
            Input("compound-table", "selected_rows"),
            State("compound-table", "data"),
            prevent_initial_call=True,
        )
        def update_selected_compound(
            selected_rows: List[int], table_data: List[Dict]
        ) -> tuple:
            """Update selected compound when table selection changes.

            Args:
                selected_rows: List of selected row indices
                table_data: Current table data

            Returns:
                Tuple of (None, status message, info message)
            """
            if not selected_rows or not table_data:
                raise PreventUpdate

            try:
                row_idx = selected_rows[0]
                compound_data = table_data[row_idx]
                compound = next(
                    (c for c in self.compounds if c.name == compound_data["name"]),
                    None,
                )

                if compound:
                    self.selected_compound = compound
                    self.detail_view.update_data(compound)

                    # Clear visualization caches when compound changes
                    self._create_cached_plot.cache_clear()
                    self._create_cached_web_viz.cache_clear()

                    return None, f"Selected compound: {compound.name}", ""

                return None, "Error: Compound not found", "Please try selecting again"

            except Exception as e:
                self.logger.error(f"Error updating selected compound: {str(e)}")
                return None, f"Error: {str(e)}", "Please try again"

        @self.app.callback(
            [
                Output("progress-bar", "value"),
                Output("progress-bar", "style"),
                Output("progress-text", "children"),
                Output("progress-text", "style"),
            ],
            [Input("plot-container", "children"), Input("viz-container", "children")],
            prevent_initial_call=True,
        )
        def update_progress(plot_children, viz_children):
            """Update progress indicators.

            Args:
                plot_children: Plot container children
                viz_children: Visualization container children

            Returns:
                Tuple of progress bar value, style, text, and text style
            """
            ctx = dash.callback_context
            if not ctx.triggered:
                raise PreventUpdate

            triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

            if triggered_id == "plot-container":
                if isinstance(plot_children, str) and "Error" in plot_children:
                    return (
                        100,
                        {"visibility": "visible"},
                        plot_children,
                        {
                            "visibility": "visible",
                            "color": "red",
                        },
                    )
                return 100, {"visibility": "hidden"}, "", {"visibility": "hidden"}

            if triggered_id == "viz-container":
                if isinstance(viz_children, str) and "Error" in viz_children:
                    return (
                        100,
                        {"visibility": "visible"},
                        viz_children,
                        {
                            "visibility": "visible",
                            "color": "red",
                        },
                    )
                return 100, {"visibility": "hidden"}, "", {"visibility": "hidden"}

            raise PreventUpdate

    @lru_cache(maxsize=32)
    def _create_cached_plot(
        self,
        selected_types: tuple,
        layout_type: str,
        compound_name: str,
    ) -> html.Div:
        """Create cached prediction plots.

        Args:
            selected_types: Tuple of selected plot types
            layout_type: Type of layout
            compound_name: Name of selected compound for cache key

        Returns:
            Plot container
        """
        if not self.selected_compound:
            return html.Div("No compound selected")

        predictions = {
            "binding": self.selected_compound.binding_predictions,
            "abuse": self.selected_compound.abuse_predictions,
            "toxicity": self.selected_compound.toxicity_predictions,
            "activity": self.selected_compound.activity_predictions,
            "psychoactive": self.selected_compound.psychoactive_predictions,
            "nootropic": self.selected_compound.nootropic_predictions,
        }

        return self.plot_manager.create_combined_plot(
            predictions=predictions,
            selected_types=list(selected_types),
            layout_type=layout_type,
        )

    @lru_cache(maxsize=32)
    def _create_cached_web_viz(
        self,
        selected_sources: tuple,
        viz_type: str,
        compound_name: str,
    ) -> html.Div:
        """Create cached web data visualization.

        Args:
            selected_sources: Tuple of selected data sources
            viz_type: Type of visualization
            compound_name: Name of selected compound for cache key

        Returns:
            Visualization container
        """
        if not self.selected_compound:
            return html.Div("No compound selected")

        web_data = {
            "community": self.selected_compound.community_data,
            "social": self.selected_compound.social_data,
            "swiss": self.selected_compound.swiss_data,
            "reports": getattr(self.selected_compound, "experience_reports", None),
            "sentiment": getattr(self.selected_compound, "sentiment_data", None),
        }

        return self.plot_manager.create_web_visualization(
            web_data=web_data,
            selected_sources=list(selected_sources),
            viz_type=viz_type,
        )
