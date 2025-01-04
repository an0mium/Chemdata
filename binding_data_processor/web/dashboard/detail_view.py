"""Compound detail view with comprehensive data display."""

import json
from typing import Dict, List, Optional

import dash
import plotly.graph_objects as go
from dash import html, dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State

from binding_data_processor.models.compound import CompoundData
from binding_data_processor.web.components.visualization import (
    StructureViewer,
    ActivityPlot,
    PredictionPlot,
)
from binding_data_processor.web.dashboard.base import BaseDashboard


class CompoundDetailView(BaseDashboard):
    """Dashboard for viewing detailed compound information."""

    def __init__(self):
        """Initialize the detail view dashboard."""
        super().__init__()
        self.compound: Optional[CompoundData] = None

    def update_data(self, compound: CompoundData) -> None:
        """Update dashboard with new compound data.

        Args:
            compound: Compound to display
        """
        self.compound = compound

    def layout(self) -> html.Div:
        """Get dashboard layout.

        Returns:
            Dash layout
        """
        if not self.compound:
            return html.Div("No compound selected")

        return html.Div(
            [
                # Header
                html.Div(
                    [
                        html.H2(self.compound.name),
                        html.H4(f"Source: {self.compound.source}"),
                    ],
                    className="mb-4",
                ),
                # Structure and basic info
                dbc.Row(
                    [
                        # Structure viewer
                        dbc.Col(
                            [
                                html.H3("Structure"),
                                StructureViewer(
                                    id="structure-viewer",
                                    smiles=self.compound.smiles,
                                ),
                                html.Pre(self.compound.smiles),
                            ],
                            width=6,
                        ),
                        # Basic info
                        dbc.Col(
                            [
                                html.H3("Properties"),
                                self._render_properties(),
                            ],
                            width=6,
                        ),
                    ],
                    className="mb-4",
                ),
                # Binding data
                html.Div(
                    [
                        html.H3("Binding Data"),
                        self._render_binding_data(),
                    ],
                    className="mb-4",
                ),
                # Predictions
                html.Div(
                    [
                        html.H3("Predictions"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H4("Toxicity"),
                                        PredictionPlot(
                                            id="toxicity-plot",
                                            data=self.compound.toxicity_predictions,
                                        ),
                                    ],
                                    width=6,
                                ),
                                dbc.Col(
                                    [
                                        html.H4("Abuse Potential"),
                                        PredictionPlot(
                                            id="abuse-plot",
                                            data=self.compound.abuse_predictions,
                                        ),
                                    ],
                                    width=6,
                                ),
                            ]
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H4("Activity"),
                                        ActivityPlot(
                                            id="activity-plot",
                                            data=self.compound.activity_predictions,
                                        ),
                                    ]
                                )
                            ]
                        ),
                    ],
                    className="mb-4",
                ),
                # Web data
                html.Div(
                    [
                        html.H3("Web Data"),
                        self._render_web_data(),
                    ],
                    className="mb-4",
                ),
                # Export
                html.Div(
                    [
                        dbc.Button(
                            "Export JSON",
                            id="export-json",
                            color="primary",
                        ),
                        dcc.Download(id="download"),
                    ],
                    className="mt-4",
                ),
            ]
        )

    def _render_properties(self) -> html.Div:
        """Render chemical properties section.

        Returns:
            Dash layout
        """
        properties = []
        if hasattr(self.compound, "molecular_weight"):
            properties.append(
                html.P(f"Molecular Weight: {self.compound.molecular_weight:.2f}")
            )
        if hasattr(self.compound, "logp"):
            properties.append(html.P(f"LogP: {self.compound.logp:.2f}"))
        if hasattr(self.compound, "hbd"):
            properties.append(html.P(f"H-Bond Donors: {self.compound.hbd}"))
        if hasattr(self.compound, "hba"):
            properties.append(html.P(f"H-Bond Acceptors: {self.compound.hba}"))
        if hasattr(self.compound, "tpsa"):
            properties.append(html.P(f"TPSA: {self.compound.tpsa:.2f}"))
        if hasattr(self.compound, "rotatable_bonds"):
            properties.append(
                html.P(f"Rotatable Bonds: {self.compound.rotatable_bonds}")
            )
        return html.Div(properties)

    def _render_binding_data(self) -> html.Div:
        """Render binding data section.

        Returns:
            Dash layout
        """
        if not self.compound.binding_data:
            return html.P("No binding data available")

        rows = []
        for data in self.compound.binding_data:
            rows.append(
                html.Tr(
                    [
                        html.Td(data.get("target", "")),
                        html.Td(data.get("activity_type", "")),
                        html.Td(
                            f"{data.get('activity_value', '')} "
                            f"{data.get('activity_unit', '')}"
                        ),
                        html.Td(data.get("reference", "")),
                    ]
                )
            )

        return html.Table(
            [
                html.Thead(
                    html.Tr(
                        [
                            html.Th("Target"),
                            html.Th("Activity Type"),
                            html.Th("Value"),
                            html.Th("Reference"),
                        ]
                    )
                ),
                html.Tbody(rows),
            ],
            className="table",
        )

    def _render_web_data(self) -> html.Div:
        """Render web data section.

        Returns:
            Dash layout
        """
        if not self.compound.web_data:
            return html.P("No web data available")

        sections = []
        for source, data in self.compound.web_data.items():
            sections.append(
                html.Div(
                    [
                        html.H4(source),
                        html.Pre(
                            json.dumps(data, indent=2),
                            style={"maxHeight": "300px", "overflow": "auto"},
                        ),
                    ]
                )
            )
        return html.Div(sections)

    def register_callbacks(self) -> None:
        """Register dashboard callbacks."""

        @self.app.callback(
            Output("download", "data"),
            Input("export-json", "n_clicks"),
            prevent_initial_call=True,
        )
        def export_data(n_clicks):
            """Export compound data as JSON."""
            if not n_clicks or not self.compound:
                return None

            return dict(
                content=json.dumps(self.compound.to_dict(), indent=2),
                filename=f"{self.compound.name}.json",
            )
