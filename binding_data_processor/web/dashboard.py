"""Interactive dashboard for compound analysis and visualization.

This module provides:
1. Interactive compound data visualization and analysis
2. Machine learning-powered predictions
3. Real-time data filtering and processing
4. Dynamic chart generation
5. Comparative analysis tools
6. Export capabilities
"""

import logging
from typing import Dict, List, Optional, Union, Tuple
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import plotly.express as px
from rdkit import Chem
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
import torch
import torch.nn as nn

from .components.input import InputOptions
from .components.analysis import AnalysisOptions
from .components.visualization import VisualizationOptions
from .components.export import ExportOptions
from ..processors.structure.depiction import StructureDepiction
from ..processors.structure.ml.advanced_visualization import AdvancedVisualizer
from ..processors.structure.descriptors import DescriptorCalculator
from ..processors.activity_analysis import ActivityAnalyzer
from ..processors.structure.ml.activity import ActivityPredictor
from ..processors.structure.ml.ensemble import EnsemblePredictor
from ..models.compound import CompoundData


class Dashboard:
    """Interactive dashboard for compound analysis and visualization."""

    def __init__(
        self,
        app: dash.Dash,
        structure_depiction: Optional[StructureDepiction] = None,
        visualizer: Optional[AdvancedVisualizer] = None,
        descriptor_calc: Optional[DescriptorCalculator] = None,
        activity_analyzer: Optional[ActivityAnalyzer] = None,
        activity_predictor: Optional[ActivityPredictor] = None,
        ensemble_predictor: Optional[EnsemblePredictor] = None,
    ):
        """Initialize dashboard with ML components.

        Args:
            app: Dash application instance
            structure_depiction: Structure depiction handler
            visualizer: Advanced visualization component
            descriptor_calc: Molecular descriptor calculator
            activity_analyzer: Activity analysis component
            activity_predictor: ML-based activity predictor
            ensemble_predictor: Ensemble prediction model
        """
        self.app = app
        self.logger = logging.getLogger(__name__)

        # Core components
        self.structure_depiction = structure_depiction or StructureDepiction()
        self.visualizer = visualizer or AdvancedVisualizer()
        self.descriptor_calc = descriptor_calc or DescriptorCalculator()
        self.activity_analyzer = activity_analyzer or ActivityAnalyzer()

        # ML components
        self.activity_predictor = activity_predictor or ActivityPredictor()
        self.ensemble_predictor = ensemble_predictor or EnsemblePredictor()

        # UI components
        self.input_options = InputOptions(app, self.structure_depiction)
        self.analysis_options = AnalysisOptions(app, self.structure_depiction)
        self.visualization_options = VisualizationOptions(app, self.structure_depiction)
        self.export_options = ExportOptions(app, self.structure_depiction)

        # Initialize ML models
        self._initialize_ml_models()

        # Setup callbacks
        self._setup_callbacks()

    def _initialize_ml_models(self):
        """Initialize machine learning models."""
        try:
            # Activity prediction models
            self.binding_affinity_model = RandomForestRegressor(
                n_estimators=100, max_depth=10, random_state=42
            )

            self.abuse_potential_model = RandomForestClassifier(
                n_estimators=100, max_depth=10, random_state=42
            )

            # Deep learning model for toxicity prediction
            self.toxicity_model = nn.Sequential(
                nn.Linear(1024, 512),
                nn.ReLU(),
                nn.Dropout(0.3),
                nn.Linear(512, 128),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.Linear(128, 32),
                nn.ReLU(),
                nn.Linear(32, 5),  # 5 toxicity endpoints
            )

            # Load pre-trained models if available
            self.activity_predictor.load_models()
            self.ensemble_predictor.load_models()

        except Exception as e:
            self.logger.error(f"Error initializing ML models: {str(e)}")

    def get_layout(self) -> html.Div:
        """Get dashboard layout."""
        return html.Div(
            [
                # Header
                html.Div(
                    [
                        html.H1("Chemical Data Analysis Dashboard"),
                        html.P("Process, analyze and visualize chemical compound data"),
                    ],
                    className="header",
                    style={
                        "textAlign": "center",
                        "padding": "20px",
                        "backgroundColor": "#f8f9fa",
                    },
                ),
                # Main content
                html.Div(
                    [
                        # Left sidebar - Input and Analysis
                        html.Div(
                            [
                                self.input_options.get_layout(),
                                html.Hr(),
                                self.analysis_options.get_layout(),
                                html.Hr(),
                                html.H3("ML Predictions"),
                                dcc.Checklist(
                                    id="ml-options",
                                    options=[
                                        {
                                            "label": "Binding Affinity",
                                            "value": "binding",
                                        },
                                        {"label": "Abuse Potential", "value": "abuse"},
                                        {"label": "Toxicity Risk", "value": "toxicity"},
                                        {
                                            "label": "Activity Profile",
                                            "value": "activity",
                                        },
                                        {
                                            "label": "Drug Interactions",
                                            "value": "interactions",
                                        },
                                    ],
                                    value=["binding", "toxicity"],
                                ),
                            ],
                            className="sidebar",
                            style={
                                "width": "25%",
                                "float": "left",
                                "padding": "20px",
                                "backgroundColor": "#ffffff",
                                "boxShadow": "2px 0px 5px rgba(0,0,0,0.1)",
                            },
                        ),
                        # Main content area - Visualization
                        html.Div(
                            [
                                self.visualization_options.get_layout(),
                                html.Div(id="ml-predictions-display"),
                                dcc.Graph(id="prediction-plots"),
                            ],
                            className="main-content",
                            style={
                                "width": "50%",
                                "float": "left",
                                "padding": "20px",
                            },
                        ),
                        # Right sidebar - Export and Settings
                        html.Div(
                            [
                                self.export_options.get_layout(),
                                html.Hr(),
                                html.Div(
                                    [
                                        html.H3("Settings"),
                                        dcc.Checklist(
                                            id="settings",
                                            options=[
                                                {
                                                    "label": "Auto-update visualizations",
                                                    "value": "auto_update",
                                                },
                                                {
                                                    "label": "Show advanced options",
                                                    "value": "advanced",
                                                },
                                                {
                                                    "label": "Enable ML predictions",
                                                    "value": "ml_enabled",
                                                },
                                                {
                                                    "label": "Real-time analysis",
                                                    "value": "realtime",
                                                },
                                            ],
                                            value=["auto_update", "ml_enabled"],
                                        ),
                                        html.H4("Model Configuration"),
                                        dcc.Dropdown(
                                            id="model-config",
                                            options=[
                                                {
                                                    "label": "High Accuracy (Slow)",
                                                    "value": "accurate",
                                                },
                                                {
                                                    "label": "Balanced",
                                                    "value": "balanced",
                                                },
                                                {
                                                    "label": "Fast Prediction",
                                                    "value": "fast",
                                                },
                                            ],
                                            value="balanced",
                                        ),
                                    ]
                                ),
                            ],
                            className="sidebar",
                            style={
                                "width": "25%",
                                "float": "right",
                                "padding": "20px",
                                "backgroundColor": "#ffffff",
                                "boxShadow": "-2px 0px 5px rgba(0,0,0,0.1)",
                            },
                        ),
                    ],
                    style={"display": "flex", "minHeight": "calc(100vh - 100px)"},
                ),
                # Hidden data storage
                dcc.Store(id="stored-data"),
                dcc.Store(id="analysis-results"),
                dcc.Store(id="ml-predictions"),
                dcc.Store(id="visualization-state"),
            ]
        )

    def _setup_callbacks(self):
        """Setup dashboard callbacks."""

        @self.app.callback(
            [Output("ml-predictions", "data"), Output("prediction-plots", "figure")],
            [Input("stored-data", "data"), Input("ml-options", "value")],
            [State("model-config", "value")],
        )
        def update_ml_predictions(data, selected_predictions, model_config):
            """Update ML predictions based on data and selected options."""
            if not data or not selected_predictions:
                return None, go.Figure()

            try:
                predictions = {}
                figures = []

                for prediction_type in selected_predictions:
                    if prediction_type == "binding":
                        binding_predictions = (
                            self.activity_predictor.predict_binding_affinity(
                                data, model_config
                            )
                        )
                        predictions["binding"] = binding_predictions
                        figures.append(self._create_binding_plot(binding_predictions))

                    elif prediction_type == "abuse":
                        abuse_predictions = (
                            self.activity_predictor.predict_abuse_potential(
                                data, model_config
                            )
                        )
                        predictions["abuse"] = abuse_predictions
                        figures.append(self._create_abuse_plot(abuse_predictions))

                    elif prediction_type == "toxicity":
                        toxicity_predictions = self.activity_predictor.predict_toxicity(
                            data, model_config
                        )
                        predictions["toxicity"] = toxicity_predictions
                        figures.append(self._create_toxicity_plot(toxicity_predictions))

                # Combine plots
                combined_figure = go.Figure()
                for fig in figures:
                    for trace in fig.data:
                        combined_figure.add_trace(trace)

                combined_figure.update_layout(
                    title="ML Predictions Overview", showlegend=True, height=800
                )

                return predictions, combined_figure

            except Exception as e:
                self.logger.error(f"Error updating ML predictions: {str(e)}")
                return None, go.Figure()

        @self.app.callback(
            [Output("analysis-results", "data"), Output("visualization-state", "data")],
            [Input("stored-data", "data"), Input("ml-predictions", "data")],
            [State("settings", "value")],
        )
        def update_analysis_and_visualization(data, ml_predictions, settings):
            """Update analysis results and visualization state."""
            if not data:
                return None, None

            try:
                auto_update = "auto_update" in (settings or [])
                if not auto_update:
                    return dash.no_update, dash.no_update

                # Combine experimental and predicted data
                combined_data = self._combine_data(data, ml_predictions)

                # Trigger analysis and visualization updates
                analysis_results = self.analysis_options.analyze_data(combined_data)
                visualization_state = self.visualization_options.update_state(
                    combined_data, analysis_results
                )

                return analysis_results, visualization_state

            except Exception as e:
                self.logger.error(f"Error updating analysis/visualization: {str(e)}")
                return None, None

        @self.app.callback(
            Output("stored-data", "data"),
            [Input("upload-data", "contents")],
            [State("upload-data", "filename")],
        )
        def store_uploaded_data(contents, filename):
            """Store and process uploaded data."""
            if not contents:
                return None

            try:
                # Process upload
                data = self.input_options.process_upload(contents, filename)

                # Enrich with web data if available
                if data:
                    enriched_data = self._enrich_compound_data(data)
                    return enriched_data

                return data

            except Exception as e:
                self.logger.error(f"Error processing upload: {str(e)}")
                return None

    def _create_binding_plot(self, predictions: Dict) -> go.Figure:
        """Create binding affinity prediction plot."""
        try:
            fig = go.Figure()

            # Add predicted binding affinities
            receptors = list(predictions.keys())
            affinities = [predictions[r]["affinity"] for r in receptors]
            confidence = [predictions[r]["confidence"] for r in receptors]

            # Create heatmap
            fig.add_trace(
                go.Heatmap(
                    z=[affinities],
                    x=receptors,
                    y=["Binding Affinity"],
                    colorscale="RdBu",
                    showscale=True,
                )
            )

            # Add confidence intervals
            fig.add_trace(
                go.Scatter(
                    x=receptors,
                    y=confidence,
                    mode="markers",
                    name="Prediction Confidence",
                    marker=dict(size=10),
                )
            )

            fig.update_layout(
                title="Predicted Binding Affinities",
                xaxis_title="Receptor",
                yaxis_title="Affinity (pKi)",
                showlegend=True,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating binding plot: {str(e)}")
            return go.Figure()

    def _create_abuse_plot(self, predictions: Dict) -> go.Figure:
        """Create abuse potential prediction plot."""
        try:
            fig = go.Figure()

            # Add radar plot of abuse-related predictions
            categories = list(predictions.keys())
            values = [predictions[c]["score"] for c in categories]
            confidence = [predictions[c]["confidence"] for c in categories]

            fig.add_trace(
                go.Scatterpolar(
                    r=values, theta=categories, fill="toself", name="Abuse Potential"
                )
            )

            # Add confidence intervals
            fig.add_trace(
                go.Scatterpolar(
                    r=confidence, theta=categories, fill="tonext", name="Confidence"
                )
            )

            fig.update_layout(
                title="Predicted Abuse Potential",
                polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
                showlegend=True,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating abuse plot: {str(e)}")
            return go.Figure()

    def _create_toxicity_plot(self, predictions: Dict) -> go.Figure:
        """Create toxicity prediction plot."""
        try:
            fig = go.Figure()

            # Add bar plot of toxicity predictions
            endpoints = list(predictions.keys())
            risks = [predictions[e]["risk"] for e in endpoints]
            uncertainty = [predictions[e]["uncertainty"] for e in endpoints]

            # Add risk bars
            fig.add_trace(
                go.Bar(x=endpoints, y=risks, name="Risk Level", marker_color="red")
            )

            # Add error bars for uncertainty
            fig.add_trace(
                go.Bar(
                    x=endpoints,
                    y=uncertainty,
                    name="Uncertainty",
                    marker_color="gray",
                    opacity=0.5,
                )
            )

            fig.update_layout(
                title="Predicted Toxicity Risks",
                xaxis_title="Toxicity Endpoint",
                yaxis_title="Risk Level",
                barmode="overlay",
                showlegend=True,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating toxicity plot: {str(e)}")
            return go.Figure()

    def _combine_data(self, experimental_data: Dict, ml_predictions: Dict) -> Dict:
        """Combine experimental data with ML predictions."""
        try:
            combined = experimental_data.copy()

            if ml_predictions:
                combined["predictions"] = ml_predictions

                # Calculate consensus values where applicable
                if (
                    "binding" in ml_predictions
                    and "experimental_binding" in experimental_data
                ):
                    combined["consensus_binding"] = self._calculate_consensus(
                        experimental_data["experimental_binding"],
                        ml_predictions["binding"],
                    )

            return combined

        except Exception as e:
            self.logger.error(f"Error combining data: {str(e)}")
            return experimental_data

    def _calculate_consensus(
        self, experimental: Dict, predicted: Dict, weight_experimental: float = 0.7
    ) -> Dict:
        """Calculate consensus between experimental and predicted values."""
        try:
            consensus = {}

            # Find overlapping keys
            common_keys = set(experimental.keys()) & set(predicted.keys())

            for key in common_keys:
                exp_value = experimental[key].get("value", 0)
                pred_value = predicted[key].get("value", 0)

                # Calculate weighted average
                consensus[key] = {
                    "value": weight_experimental * exp_value
                    + (1 - weight_experimental) * pred_value,
                    "experimental_weight": weight_experimental,
                    "experimental_value": exp_value,
                    "predicted_value": pred_value,
                }

            return consensus

        except Exception as e:
            self.logger.error(f"Error calculating consensus: {str(e)}")
            return {}

    def _enrich_compound_data(self, data: Dict) -> Dict:
        """Enrich compound data with additional web-sourced information."""
        try:
            enriched = data.copy()

            # Add web-sourced data if available
            if "smiles" in data:
                web_data = self.activity_predictor.get_web_data(data["smiles"])
                if web_data:
                    enriched["web_data"] = web_data

            return enriched

        except Exception as e:
            self.logger.error(f"Error enriching data: {str(e)}")
            return data

    def run_server(self, debug: bool = False, port: int = 8050):
        """Run the dashboard server.

        Args:
            debug: Enable debug mode
            port: Server port
        """
        self.app.run_server(debug=debug, port=port)
