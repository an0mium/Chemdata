"""Compound list view with filtering and sorting."""

import pandas as pd
import json
from typing import Dict, List, Optional

import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State

from binding_data_processor.models.compound import CompoundData
from binding_data_processor.web.components.input import FilterInput
from binding_data_processor.web.components.visualization import CompoundTable
from binding_data_processor.web.dashboard.base import BaseDashboard


class CompoundListView(BaseDashboard):
    """Dashboard for viewing compound lists with filtering and sorting."""

    def __init__(self):
        """Initialize the list view dashboard."""
        super().__init__()
        self.compounds: List[CompoundData] = []
        self.df: Optional[pd.DataFrame] = None

    def update_data(self, compounds: List[CompoundData]) -> None:
        """Update dashboard with new compound data.

        Args:
            compounds: List of compounds to display
        """
        self.compounds = compounds
        self.df = pd.DataFrame([c.to_dict() for c in compounds])

    def layout(self) -> html.Div:
        """Get dashboard layout.

        Returns:
            Dash layout
        """
        return html.Div(
            [
                # Filters
                html.Div(
                    [
                        html.H3("Filters"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Label("Target"),
                                        FilterInput(
                                            id="target-filter",
                                            options=self._get_target_options(),
                                        ),
                                    ]
                                ),
                                dbc.Col(
                                    [
                                        html.Label("Activity Type"),
                                        FilterInput(
                                            id="activity-filter",
                                            options=self._get_activity_options(),
                                        ),
                                    ]
                                ),
                                dbc.Col(
                                    [
                                        html.Label("Source"),
                                        FilterInput(
                                            id="source-filter",
                                            options=self._get_source_options(),
                                        ),
                                    ]
                                ),
                            ]
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Label("Toxicity Score Range"),
                                        dcc.RangeSlider(
                                            id="toxicity-range",
                                            min=0,
                                            max=1,
                                            step=0.1,
                                            value=[0, 1],
                                            marks={
                                                i / 10: str(i / 10) for i in range(11)
                                            },
                                        ),
                                    ]
                                ),
                                dbc.Col(
                                    [
                                        html.Label("Abuse Potential Range"),
                                        dcc.RangeSlider(
                                            id="abuse-range",
                                            min=0,
                                            max=1,
                                            step=0.1,
                                            value=[0, 1],
                                            marks={
                                                i / 10: str(i / 10) for i in range(11)
                                            },
                                        ),
                                    ]
                                ),
                            ]
                        ),
                    ],
                    className="mb-4",
                ),
                # Results
                html.Div(
                    [
                        html.H3("Results"),
                        html.Div(id="result-stats"),
                        CompoundTable(id="compound-table"),
                    ]
                ),
                # Export
                html.Div(
                    [
                        html.H3("Export"),
                        dbc.Button(
                            "Export TSV",
                            id="export-tsv",
                            color="primary",
                            className="mr-2",
                        ),
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

    def _get_target_options(self) -> List[Dict[str, str]]:
        """Get target filter options.

        Returns:
            List of {label, value} dicts
        """
        if self.df is None:
            return []
        targets = set()
        for binding_data in self.df["binding_data"]:
            if binding_data:
                for data in binding_data:
                    if "target" in data:
                        targets.add(data["target"])
        return [{"label": t, "value": t} for t in sorted(targets)]

    def _get_activity_options(self) -> List[Dict[str, str]]:
        """Get activity type filter options.

        Returns:
            List of {label, value} dicts
        """
        if self.df is None:
            return []
        activities = set()
        for binding_data in self.df["binding_data"]:
            if binding_data:
                for data in binding_data:
                    if "activity_type" in data:
                        activities.add(data["activity_type"])
        return [{"label": a, "value": a} for a in sorted(activities)]

    def _get_source_options(self) -> List[Dict[str, str]]:
        """Get source filter options.

        Returns:
            List of {label, value} dicts
        """
        if self.df is None:
            return []
        sources = set(self.df["source"].dropna())
        return [{"label": s, "value": s} for s in sorted(sources)]

    def register_callbacks(self) -> None:
        """Register dashboard callbacks."""

        @self.app.callback(
            Output("compound-table", "data"),
            Output("result-stats", "children"),
            Input("target-filter", "value"),
            Input("activity-filter", "value"),
            Input("source-filter", "value"),
            Input("toxicity-range", "value"),
            Input("abuse-range", "value"),
        )
        def update_table(
            target: Optional[str],
            activity: Optional[str],
            source: Optional[str],
            toxicity_range: List[float],
            abuse_range: List[float],
        ) -> tuple:
            if self.df is None:
                return [], "No data loaded"

            # Apply filters
            mask = pd.Series([True] * len(self.df))

            if target:
                target_mask = self.df["binding_data"].apply(
                    lambda x: (
                        any(d.get("target") == target for d in x if d) if x else False
                    )
                )
                mask &= target_mask

            if activity:
                activity_mask = self.df["binding_data"].apply(
                    lambda x: (
                        any(d.get("activity_type") == activity for d in x if d)
                        if x
                        else False
                    )
                )
                mask &= activity_mask

            if source:
                mask &= self.df["source"] == source

            if toxicity_range:
                tox_mask = self.df["toxicity_predictions"].apply(
                    lambda x: (
                        (toxicity_range[0] <= x.get("score", 0) <= toxicity_range[1])
                        if x
                        else False
                    )
                )
                mask &= tox_mask

            if abuse_range:
                abuse_mask = self.df["abuse_predictions"].apply(
                    lambda x: (
                        (abuse_range[0] <= x.get("score", 0) <= abuse_range[1])
                        if x
                        else False
                    )
                )
                mask &= abuse_mask

            filtered_df = self.df[mask]

            # Update stats
            stats = (
                f"Showing {len(filtered_df)} of {len(self.df)} compounds "
                f"({(len(filtered_df) / len(self.df)) * 100:.1f}%)"
            )

            return filtered_df.to_dict("records"), stats

        @self.app.callback(
            Output("download", "data"),
            Input("export-tsv", "n_clicks"),
            Input("export-json", "n_clicks"),
            State("compound-table", "data"),
            prevent_initial_call=True,
        )
        def export_data(tsv_clicks, json_clicks, data):
            """Export filtered data."""
            if not data:
                return None

            ctx = dash.callback_context
            if not ctx.triggered:
                return None

            trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

            if trigger_id == "export-tsv":
                df = pd.DataFrame(data)
                return dcc.send_data_frame(
                    df.to_csv, "compounds.tsv", sep="\t", index=False
                )
            elif trigger_id == "export-json":
                return dict(
                    content=json.dumps(data, indent=2), filename="compounds.json"
                )
