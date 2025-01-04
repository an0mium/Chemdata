"""Analysis options component for the web interface."""

import logging
from typing import Callable, Dict, List, Optional
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State

from ...processors.activity import ActivityProcessor
from ...processors.structure.descriptors import DescriptorCalculator
from ...processors.structure.similarity import SimilarityCalculator


class AnalysisOptions:
    """Component for analysis options and settings."""

    def __init__(
        self,
        parent_app: dash.Dash,
        activity_processor: Optional[ActivityProcessor] = None,
        descriptor_calc: Optional[DescriptorCalculator] = None,
        similarity_calc: Optional[SimilarityCalculator] = None,
        on_settings_change: Optional[Callable[[Dict], None]] = None,
    ):
        """Initialize analysis options component.

        Args:
            parent_app: Parent Dash application
            activity_processor: Optional activity processor
            descriptor_calc: Optional descriptor calculator
            similarity_calc: Optional similarity calculator
            on_settings_change: Optional callback for settings changes
        """
        self.app = parent_app
        self.activity_processor = activity_processor or ActivityProcessor()
        self.descriptor_calc = descriptor_calc or DescriptorCalculator()
        self.similarity_calc = similarity_calc or SimilarityCalculator()
        self.on_settings_change = on_settings_change
        self.logger = logging.getLogger(__name__)
        self._setup_callbacks()

    def get_layout(self) -> html.Div:
        """Get component layout."""
        return html.Div(
            [
                html.H3("Analysis Options"),
                dcc.Tabs(
                    [
                        dcc.Tab(
                            label="Activity",
                            children=[
                                html.Label("Activity Type"),
                                dcc.Dropdown(
                                    id="activity-type",
                                    options=[
                                        {"label": "Agonist", "value": "agonist"},
                                        {"label": "Antagonist", "value": "antagonist"},
                                        {"label": "Allosteric", "value": "allosteric"},
                                        {"label": "Unknown", "value": "unknown"},
                                    ],
                                    value="unknown",
                                ),
                                html.Label("Activity Threshold (nM)"),
                                dcc.Input(
                                    id="activity-threshold",
                                    type="number",
                                    value=100,
                                    min=0,
                                ),
                            ],
                        ),
                        dcc.Tab(
                            label="Structure",
                            children=[
                                html.Label("Similarity Method"),
                                dcc.Dropdown(
                                    id="similarity-method",
                                    options=[
                                        {"label": "Tanimoto", "value": "tanimoto"},
                                        {"label": "Dice", "value": "dice"},
                                        {"label": "Cosine", "value": "cosine"},
                                    ],
                                    value="tanimoto",
                                ),
                                html.Label("Similarity Threshold"),
                                dcc.Slider(
                                    id="similarity-threshold",
                                    min=0,
                                    max=1,
                                    step=0.05,
                                    value=0.7,
                                    marks={
                                        i / 10: str(i / 10) for i in range(0, 11, 2)
                                    },
                                ),
                            ],
                        ),
                        dcc.Tab(
                            label="Descriptors",
                            children=[
                                html.Label("Descriptor Types"),
                                dcc.Checklist(
                                    id="descriptor-types",
                                    options=[
                                        {
                                            "label": "Topological",
                                            "value": "topological",
                                        },
                                        {
                                            "label": "Constitutional",
                                            "value": "constitutional",
                                        },
                                        {"label": "Electronic", "value": "electronic"},
                                        {"label": "Geometric", "value": "geometric"},
                                    ],
                                    value=["topological", "constitutional"],
                                ),
                                html.Label("Normalization"),
                                dcc.Dropdown(
                                    id="normalization",
                                    options=[
                                        {"label": "None", "value": "none"},
                                        {"label": "Min-Max", "value": "minmax"},
                                        {"label": "Z-Score", "value": "zscore"},
                                    ],
                                    value="none",
                                ),
                            ],
                        ),
                    ]
                ),
                html.Button("Apply", id="apply-settings"),
                html.Div(id="settings-error"),
            ],
            style={"padding": "20px"},
        )

    def _setup_callbacks(self):
        """Setup component callbacks."""

        @self.app.callback(
            Output("settings-error", "children"),
            [Input("apply-settings", "n_clicks")],
            [
                State("activity-type", "value"),
                State("activity-threshold", "value"),
                State("similarity-method", "value"),
                State("similarity-threshold", "value"),
                State("descriptor-types", "value"),
                State("normalization", "value"),
            ],
        )
        def update_settings(
            n_clicks: int,
            activity_type: str,
            activity_threshold: float,
            similarity_method: str,
            similarity_threshold: float,
            descriptor_types: List[str],
            normalization: str,
        ) -> str:
            """Update analysis settings."""
            if not n_clicks:
                return ""

            try:
                settings = {
                    "activity": {
                        "type": activity_type,
                        "threshold": float(activity_threshold),
                    },
                    "structure": {
                        "similarity_method": similarity_method,
                        "similarity_threshold": float(similarity_threshold),
                    },
                    "descriptors": {
                        "types": descriptor_types,
                        "normalization": normalization,
                    },
                }

                # Validate settings
                if activity_threshold < 0:
                    return "Activity threshold must be non-negative"

                if not 0 <= similarity_threshold <= 1:
                    return "Similarity threshold must be between 0 and 1"

                if not descriptor_types:
                    return "At least one descriptor type must be selected"

                # Update processors
                self.activity_processor.set_threshold(activity_threshold)
                self.similarity_calc.set_method(similarity_method)
                self.descriptor_calc.set_types(descriptor_types)

                # Notify parent
                if self.on_settings_change:
                    self.on_settings_change(settings)

                return ""

            except Exception as e:
                self.logger.error(f"Error updating settings: {str(e)}")
                return f"Error: {str(e)}"
