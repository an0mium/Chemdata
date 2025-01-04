"""Enhanced dashboard for chemical compound analysis and visualization.

This module provides:
1. Component initialization and management
2. State management and persistence 
3. Event handling and routing
4. Layout management and configuration
5. ML model integration
6. Data enrichment and analysis
7. Interactive visualization
8. Export capabilities
"""

import logging
from typing import Dict, List, Optional, Any
from pathlib import Path
import json
import pandas as pd

from dash import Dash, html, dcc, Input, Output, State, callback
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from .data_processor import DataProcessor
from .visualization import PlotManager
from .ml_integration import MLIntegrator
from ..components.input import InputComponent
from ..components.analysis import AnalysisComponent
from ..components.visualization import VisualizationComponent
from ..components.export import ExportComponent
from ...processors.structure.depiction import StructureDepiction
from ...processors.structure.ml.activity import ActivityPredictor
from ...processors.structure.ml.ensemble import EnsemblePredictor
from ...web_enrichment.data_sources.swiss import SwissClient
from ...web_enrichment.data_sources.chembl import ChEMBLClient
from ...web_enrichment.data_sources.pubchem import PubChemClient


class Dashboard:
    """Enhanced dashboard for chemical compound analysis."""

    # Default configuration
    DEFAULT_CONFIG = {
        "theme": "bootstrap",
        "title": "Chemical Data Analysis",
        "debug": False,
        "cache_dir": "cache",
        "max_compounds": 10000,
        "plot_height": 600,
        "enrichment_types": ["web", "swiss", "chembl", "pubchem"],
        "ml_models": {
            "binding": True,
            "activity": True,
            "toxicity": True,
            "abuse": True,
        },
        "visualization_options": {
            "2d_depiction": True,
            "3d_depiction": True,
            "interactive": True,
        },
    }

    def __init__(
        self,
        config_path: Optional[str] = None,
        data_processor: Optional[DataProcessor] = None,
        plot_manager: Optional[PlotManager] = None,
        ml_integrator: Optional[MLIntegrator] = None,
        structure_depiction: Optional[StructureDepiction] = None,
    ):
        """Initialize enhanced dashboard.

        Args:
            config_path: Path to config file
            data_processor: Data processing component
            plot_manager: Plot management component
            ml_integrator: ML integration component
            structure_depiction: Structure depiction handler
        """
        self.logger = logging.getLogger(__name__)

        # Load configuration
        self.config = self._load_config(config_path)

        # Initialize core components
        self.data_processor = data_processor or DataProcessor()
        self.plot_manager = plot_manager or PlotManager()
        self.ml_integrator = ml_integrator or MLIntegrator()
        self.structure_depiction = structure_depiction or StructureDepiction()

        # Initialize ML models
        self.activity_predictor = ActivityPredictor()
        self.ensemble_predictor = EnsemblePredictor()

        # Initialize data source clients
        self.swiss_client = SwissClient()
        self.chembl_client = ChEMBLClient()
        self.pubchem_client = PubChemClient()

        # Create Dash app
        self.app = Dash(
            __name__,
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            suppress_callback_exceptions=True,
        )
        self.app.title = self.config["title"]

        # Initialize UI components
        self.input_component = InputComponent(self.data_processor)
        self.analysis_component = AnalysisComponent(
            self.data_processor, self.ml_integrator
        )
        self.visualization_component = VisualizationComponent(
            self.plot_manager, self.structure_depiction
        )
        self.export_component = ExportComponent(self.data_processor)

        # Set up layout
        self.app.layout = self._create_layout()

        # Register callbacks
        self._register_callbacks()

    def _load_config(self, config_path: Optional[str]) -> Dict:
        """Load dashboard configuration."""
        config = self.DEFAULT_CONFIG.copy()

        if config_path:
            try:
                with open(config_path) as f:
                    user_config = json.load(f)
                config.update(user_config)
            except Exception as e:
                self.logger.error(f"Error loading config: {str(e)}")

        return config

    def _create_layout(self) -> html.Div:
        """Create enhanced dashboard layout."""
        return html.Div(
            [
                # Navigation
                dbc.NavbarSimple(
                    children=[
                        dbc.NavItem(dbc.NavLink("Input", href="#input")),
                        dbc.NavItem(dbc.NavLink("Analysis", href="#analysis")),
                        dbc.NavItem(dbc.NavLink("ML Predictions", href="#predictions")),
                        dbc.NavItem(
                            dbc.NavLink("Visualization", href="#visualization")
                        ),
                        dbc.NavItem(dbc.NavLink("Export", href="#export")),
                    ],
                    brand=self.config["title"],
                    color="primary",
                    dark=True,
                ),
                # Main content
                dbc.Container(
                    [
                        # Input section
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H2("Data Input", id="input"),
                                        self.input_component.layout,
                                        dbc.Card(
                                            [
                                                dbc.CardHeader("Data Sources"),
                                                dbc.CardBody(
                                                    [
                                                        dbc.Checklist(
                                                            id="data-sources",
                                                            options=[
                                                                {
                                                                    "label": "BindingDB",
                                                                    "value": "bindingdb",
                                                                },
                                                                {
                                                                    "label": "ChEMBL",
                                                                    "value": "chembl",
                                                                },
                                                                {
                                                                    "label": "PubChem",
                                                                    "value": "pubchem",
                                                                },
                                                                {
                                                                    "label": "Swiss",
                                                                    "value": "swiss",
                                                                },
                                                            ],
                                                            value=["bindingdb"],
                                                        )
                                                    ]
                                                ),
                                            ],
                                            className="mt-3",
                                        ),
                                    ],
                                    width=12,
                                )
                            ],
                            className="mb-4",
                        ),
                        # Analysis section
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H2("Analysis", id="analysis"),
                                        self.analysis_component.layout,
                                        dbc.Card(
                                            [
                                                dbc.CardHeader("Analysis Options"),
                                                dbc.CardBody(
                                                    [
                                                        dbc.Checklist(
                                                            id="analysis-options",
                                                            options=[
                                                                {
                                                                    "label": "Structure",
                                                                    "value": "structure",
                                                                },
                                                                {
                                                                    "label": "Activity",
                                                                    "value": "activity",
                                                                },
                                                                {
                                                                    "label": "Binding",
                                                                    "value": "binding",
                                                                },
                                                                {
                                                                    "label": "Properties",
                                                                    "value": "properties",
                                                                },
                                                            ],
                                                            value=[
                                                                "structure",
                                                                "activity",
                                                            ],
                                                        )
                                                    ]
                                                ),
                                            ],
                                            className="mt-3",
                                        ),
                                    ],
                                    width=12,
                                )
                            ],
                            className="mb-4",
                        ),
                        # ML Predictions section
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H2("ML Predictions", id="predictions"),
                                        dbc.Card(
                                            [
                                                dbc.CardHeader("Prediction Models"),
                                                dbc.CardBody(
                                                    [
                                                        dbc.Checklist(
                                                            id="ml-models",
                                                            options=[
                                                                {
                                                                    "label": "Binding Affinity",
                                                                    "value": "binding",
                                                                },
                                                                {
                                                                    "label": "Activity Profile",
                                                                    "value": "activity",
                                                                },
                                                                {
                                                                    "label": "Toxicity Risk",
                                                                    "value": "toxicity",
                                                                },
                                                                {
                                                                    "label": "Abuse Potential",
                                                                    "value": "abuse",
                                                                },
                                                            ],
                                                            value=[
                                                                "binding",
                                                                "toxicity",
                                                            ],
                                                        ),
                                                        html.Div(
                                                            id="prediction-output"
                                                        ),
                                                    ]
                                                ),
                                            ],
                                            className="mt-3",
                                        ),
                                    ],
                                    width=12,
                                )
                            ],
                            className="mb-4",
                        ),
                        # Visualization section
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H2("Visualization", id="visualization"),
                                        self.visualization_component.layout,
                                        dbc.Card(
                                            [
                                                dbc.CardHeader("Plot Options"),
                                                dbc.CardBody(
                                                    [
                                                        dbc.Select(
                                                            id="plot-type",
                                                            options=[
                                                                {
                                                                    "label": "Structure Grid",
                                                                    "value": "grid",
                                                                },
                                                                {
                                                                    "label": "Activity Heatmap",
                                                                    "value": "heatmap",
                                                                },
                                                                {
                                                                    "label": "Property Distribution",
                                                                    "value": "dist",
                                                                },
                                                                {
                                                                    "label": "3D Viewer",
                                                                    "value": "3d",
                                                                },
                                                            ],
                                                            value="grid",
                                                        ),
                                                        html.Div(id="plot-output"),
                                                    ]
                                                ),
                                            ],
                                            className="mt-3",
                                        ),
                                    ],
                                    width=12,
                                )
                            ],
                            className="mb-4",
                        ),
                        # Export section
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H2("Export", id="export"),
                                        self.export_component.layout,
                                    ],
                                    width=12,
                                )
                            ],
                            className="mb-4",
                        ),
                        # Status and storage
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Div(id="status-area"),
                                        dcc.Store(id="data-store"),
                                        dcc.Store(id="prediction-store"),
                                        dcc.Store(id="analysis-store"),
                                        dcc.Store(id="visualization-store"),
                                    ]
                                )
                            ]
                        ),
                    ],
                    fluid=True,
                    className="py-4",
                ),
            ]
        )

    def _register_callbacks(self):
        """Register enhanced dashboard callbacks."""
        self._register_data_callbacks()
        self._register_analysis_callbacks()
        self._register_prediction_callbacks()
        self._register_visualization_callbacks()
        self._register_export_callbacks()

    def _register_data_callbacks(self):
        """Register data handling callbacks."""

        @self.app.callback(
            Output("data-store", "data"),
            [Input("upload-data", "contents"), Input("data-sources", "value")],
            [State("upload-data", "filename")],
        )
        def handle_data_input(contents, sources, filename):
            """Process data input from multiple sources."""
            if not contents and not sources:
                return {}

            try:
                data = {}

                # Process uploaded data
                if contents:
                    uploaded_data = self.data_processor.load_data(contents, filename)
                    data.update(uploaded_data)

                # Gather data from selected sources
                if sources:
                    for source in sources:
                        if source == "bindingdb":
                            binding_data = self.data_processor.get_bindingdb_data()
                            data.update(binding_data)
                        elif source == "chembl":
                            chembl_data = self.chembl_client.get_compound_data()
                            data.update(chembl_data)
                        elif source == "pubchem":
                            pubchem_data = self.pubchem_client.get_compound_data()
                            data.update(pubchem_data)
                        elif source == "swiss":
                            swiss_data = self.swiss_client.get_compound_data()
                            data.update(swiss_data)

                return data

            except Exception as e:
                self.logger.error(f"Error processing data input: {str(e)}")
                return {}

    def _register_analysis_callbacks(self):
        """Register analysis callbacks."""

        @self.app.callback(
            Output("analysis-store", "data"),
            [Input("data-store", "data"), Input("analysis-options", "value")],
        )
        def run_analysis(data, options):
            """Run selected analyses."""
            if not data or not options:
                return {}

            try:
                results = {}

                # Convert data to DataFrame
                df = pd.DataFrame(data)

                # Run selected analyses
                for option in options:
                    if option == "structure":
                        results["structure"] = (
                            self.analysis_component.analyze_structures(df)
                        )
                    elif option == "activity":
                        results["activity"] = (
                            self.analysis_component.analyze_activities(df)
                        )
                    elif option == "binding":
                        results["binding"] = self.analysis_component.analyze_binding(df)
                    elif option == "properties":
                        results["properties"] = (
                            self.analysis_component.analyze_properties(df)
                        )

                return results

            except Exception as e:
                self.logger.error(f"Error running analysis: {str(e)}")
                return {}

    def _register_prediction_callbacks(self):
        """Register ML prediction callbacks."""

        @self.app.callback(
            [
                Output("prediction-store", "data"),
                Output("prediction-output", "children"),
            ],
            [Input("data-store", "data"), Input("ml-models", "value")],
        )
        def run_predictions(data, models):
            """Generate ML predictions."""
            if not data or not models:
                return {}, "No predictions available"

            try:
                predictions = {}
                outputs = []

                # Convert data to DataFrame
                df = pd.DataFrame(data)

                # Run selected models
                for model in models:
                    if model == "binding":
                        binding_pred = self.activity_predictor.predict_binding(df)
                        predictions["binding"] = binding_pred
                        outputs.append(
                            dbc.Card(
                                dbc.CardBody(
                                    [
                                        html.H5("Binding Affinity Predictions"),
                                        html.Pre(json.dumps(binding_pred, indent=2)),
                                    ]
                                )
                            )
                        )
                    elif model == "activity":
                        activity_pred = self.activity_predictor.predict_activity(df)
                        predictions["activity"] = activity_pred
                        outputs.append(
                            dbc.Card(
                                dbc.CardBody(
                                    [
                                        html.H5("Activity Profile Predictions"),
                                        html.Pre(json.dumps(activity_pred, indent=2)),
                                    ]
                                )
                            )
                        )
                    elif model == "toxicity":
                        toxicity_pred = self.ensemble_predictor.predict_toxicity(df)
                        predictions["toxicity"] = toxicity_pred
                        outputs.append(
                            dbc.Card(
                                dbc.CardBody(
                                    [
                                        html.H5("Toxicity Risk Predictions"),
                                        html.Pre(json.dumps(toxicity_pred, indent=2)),
                                    ]
                                )
                            )
                        )
                    elif model == "abuse":
                        abuse_pred = self.ensemble_predictor.predict_abuse_potential(df)
                        predictions["abuse"] = abuse_pred
                        outputs.append(
                            dbc.Card(
                                dbc.CardBody(
                                    [
                                        html.H5("Abuse Potential Predictions"),
                                        html.Pre(json.dumps(abuse_pred, indent=2)),
                                    ]
                                )
                            )
                        )

                return predictions, html.Div(outputs)

            except Exception as e:
                self.logger.error(f"Error generating predictions: {str(e)}")
                return {}, f"Error: {str(e)}"

    def _register_visualization_callbacks(self):
        """Register visualization callbacks."""

        @self.app.callback(
            [Output("visualization-store", "data"), Output("plot-output", "children")],
            [
                Input("data-store", "data"),
                Input("analysis-store", "data"),
                Input("prediction-store", "data"),
                Input("plot-type", "value"),
            ],
        )
        def update_visualization(data, analysis, predictions, plot_type):
            """Update visualization."""
            if not data:
                return {}, "No data to visualize"

            try:
                # Combine all data
                combined_data = {
                    "raw_data": data,
                    "analysis": analysis or {},
                    "predictions": predictions or {},
                }

                # Create requested plot
                if plot_type == "grid":
                    fig = self.visualization_component.create_structure_grid(data)
                elif plot_type == "heatmap":
                    fig = self.visualization_component.create_activity_heatmap(
                        combined_data
                    )
                elif plot_type == "dist":
                    fig = self.visualization_component.create_property_distribution(
                        combined_data
                    )
                elif plot_type == "3d":
                    fig = self.visualization_component.create_3d_viewer(data)
                else:
                    raise ValueError(f"Unknown plot type: {plot_type}")

                return combined_data, dcc.Graph(
                    figure=fig, style={"height": self.config["plot_height"]}
                )

            except Exception as e:
                self.logger.error(f"Error updating visualization: {str(e)}")
                return {}, f"Error: {str(e)}"

    def _register_export_callbacks(self):
        """Register export callbacks."""

        @self.app.callback(
            Output("status-area", "children"),
            [Input("export-button", "n_clicks")],
            [
                State("data-store", "data"),
                State("analysis-store", "data"),
                State("prediction-store", "data"),
                State("export-format", "value"),
                State("export-path", "value"),
            ],
        )
        def export_data(
            n_clicks, data, analysis, predictions, export_format, export_path
        ):
            """Export processed data."""
            if not n_clicks or not data or not export_path:
                return ""

            try:
                # Combine all data
                combined_data = {
                    "compounds": data,
                    "analysis": analysis or {},
                    "predictions": predictions or {},
                }

                # Export data
                success = self.export_component.export_data(
                    combined_data, export_path, export_format
                )

                if success:
                    return html.Div(
                        f"Data exported successfully to {export_path}",
                        className="text-success",
                    )
                else:
                    return html.Div("Error exporting data", className="text-danger")

            except Exception as e:
                self.logger.error(f"Error exporting data: {str(e)}")
                return html.Div(f"Error: {str(e)}", className="text-danger")

    def run(self, debug: Optional[bool] = None, **kwargs) -> None:
        """Run dashboard server.

        Args:
            debug: Debug mode flag
            **kwargs: Additional arguments for Dash server
        """
        debug = debug if debug is not None else self.config["debug"]
        self.app.run_server(debug=debug, **kwargs)
