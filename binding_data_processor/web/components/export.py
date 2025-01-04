"""Export options component for the web interface."""

import logging
from typing import Dict, List, Optional
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import pandas as pd
import json
import io
import base64

from ...models.compound import CompoundData
from ...processors.structure.depiction import StructureDepiction


class ExportOptions:
    """Component for data export and download options."""

    def __init__(
        self,
        parent_app: dash.Dash,
        structure_depiction: Optional[StructureDepiction] = None,
    ):
        """Initialize export options component.

        Args:
            parent_app: Parent Dash application
            structure_depiction: Optional structure depiction handler
        """
        self.app = parent_app
        self.structure_depiction = structure_depiction or StructureDepiction()
        self.logger = logging.getLogger(__name__)
        self._setup_callbacks()

    def get_layout(self) -> html.Div:
        """Get component layout."""
        return html.Div(
            [
                html.H3("Export Options"),
                html.Label("Export Format"),
                dcc.Dropdown(
                    id="export-format",
                    options=[
                        {"label": "TSV", "value": "tsv"},
                        {"label": "CSV", "value": "csv"},
                        {"label": "Excel", "value": "excel"},
                        {"label": "JSON", "value": "json"},
                        {"label": "SDF", "value": "sdf"},
                    ],
                    value="tsv",
                ),
                html.Label("Include Columns"),
                dcc.Checklist(
                    id="export-columns",
                    options=[
                        {"label": "Basic Info", "value": "basic"},
                        {"label": "Activity Data", "value": "activity"},
                        {"label": "Structure Data", "value": "structure"},
                        {"label": "Descriptors", "value": "descriptors"},
                        {"label": "References", "value": "references"},
                    ],
                    value=["basic", "activity"],
                ),
                html.Label("Structure Format"),
                dcc.RadioItems(
                    id="structure-format",
                    options=[
                        {"label": "SMILES", "value": "smiles"},
                        {"label": "InChI", "value": "inchi"},
                        {"label": "Both", "value": "both"},
                    ],
                    value="smiles",
                ),
                html.Button(
                    "Export Data",
                    id="export-button",
                    style={"margin-top": "20px"},
                ),
                dcc.Download(id="download-data"),
                html.Div(id="export-error"),
            ],
            style={"padding": "20px"},
        )

    def _setup_callbacks(self):
        """Setup component callbacks."""

        @self.app.callback(
            Output("download-data", "data"),
            [Input("export-button", "n_clicks")],
            [
                State("export-format", "value"),
                State("export-columns", "value"),
                State("structure-format", "value"),
                State("stored-data", "data"),
            ],
        )
        def export_data(
            n_clicks: int,
            export_format: str,
            columns: List[str],
            structure_format: str,
            data: List[Dict],
        ) -> Dict:
            """Export data in selected format."""
            if not n_clicks or not data:
                return None

            try:
                # Convert to DataFrame
                df = pd.DataFrame(data)

                # Filter columns based on selection
                selected_columns = []
                if "basic" in columns:
                    selected_columns.extend(
                        ["name", "cas_number", "chembl_id", "drugbank_id"]
                    )
                if "activity" in columns:
                    selected_columns.extend(
                        [
                            "activity_type",
                            "activity_value",
                            "activity_unit",
                            "target",
                            "assay_type",
                        ]
                    )
                if "structure" in columns:
                    if structure_format == "smiles":
                        selected_columns.append("smiles")
                    elif structure_format == "inchi":
                        selected_columns.append("inchi")
                    else:
                        selected_columns.extend(["smiles", "inchi"])
                if "descriptors" in columns:
                    descriptor_cols = [
                        col for col in df.columns if col.startswith("descriptors.")
                    ]
                    selected_columns.extend(descriptor_cols)
                if "references" in columns:
                    selected_columns.extend(
                        ["doi", "pmid", "patent_number", "reference_urls"]
                    )

                # Filter DataFrame
                df = df[selected_columns]

                # Export in selected format
                if export_format == "tsv":
                    buffer = io.StringIO()
                    df.to_csv(buffer, sep="\t", index=False)
                    content = buffer.getvalue()
                    filename = "compounds.tsv"
                    mimetype = "text/tab-separated-values"

                elif export_format == "csv":
                    buffer = io.StringIO()
                    df.to_csv(buffer, index=False)
                    content = buffer.getvalue()
                    filename = "compounds.csv"
                    mimetype = "text/csv"

                elif export_format == "excel":
                    buffer = io.BytesIO()
                    df.to_excel(buffer, index=False)
                    content = base64.b64encode(buffer.getvalue()).decode()
                    filename = "compounds.xlsx"
                    mimetype = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"

                elif export_format == "json":
                    content = df.to_json(orient="records")
                    filename = "compounds.json"
                    mimetype = "application/json"

                elif export_format == "sdf":
                    # Convert compounds to SDF format
                    buffer = io.StringIO()
                    compounds = [CompoundData(**d) for d in data]
                    self.structure_depiction.write_sdf(compounds, buffer)
                    content = buffer.getvalue()
                    filename = "compounds.sdf"
                    mimetype = "chemical/x-mdl-sdfile"

                return dict(
                    content=content,
                    filename=filename,
                    type=mimetype,
                )

            except Exception as e:
                self.logger.error(f"Error exporting data: {str(e)}")
                return None

        @self.app.callback(
            Output("export-error", "children"),
            [Input("export-button", "n_clicks")],
            [State("stored-data", "data")],
        )
        def show_export_error(n_clicks: int, data: List[Dict]) -> str:
            """Show export error message."""
            if not n_clicks:
                return ""
            if not data:
                return "No data available to export"
            return ""
