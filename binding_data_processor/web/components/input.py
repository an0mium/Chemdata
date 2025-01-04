"""Compound input component for the web interface."""

import logging
from typing import Callable, Dict, Optional, List, Tuple
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import pandas as pd
from rdkit import Chem

from ...models.compound import CompoundData
from ...processors.structure.descriptors import DescriptorCalculator


class CompoundInput:
    """Component for compound structure input."""

    def __init__(
        self,
        parent_app: dash.Dash,
        descriptor_calc: Optional[DescriptorCalculator] = None,
        on_compound_change: Optional[Callable[[CompoundData], None]] = None,
    ):
        """Initialize compound input component.

        Args:
            parent_app: Parent Dash application
            descriptor_calc: Optional descriptor calculator
            on_compound_change: Optional callback for compound changes
        """
        self.app = parent_app
        self.descriptor_calc = descriptor_calc or DescriptorCalculator()
        self.on_compound_change = on_compound_change
        self.logger = logging.getLogger(__name__)
        self._setup_callbacks()

    def get_layout(self) -> html.Div:
        """Get component layout."""
        return html.Div(
            [
                html.H3("Input Compound"),
                dcc.Tabs(
                    [
                        dcc.Tab(
                            label="Structure",
                            children=[
                                dcc.Input(
                                    id="smiles-input",
                                    placeholder="Enter SMILES",
                                    style={"width": "100%"},
                                ),
                                html.Button("Analyze", id="analyze-btn"),
                                html.Div(id="structure-error"),
                            ],
                        ),
                        dcc.Tab(
                            label="File Upload",
                            children=[
                                dcc.Upload(
                                    id="file-upload",
                                    children=html.Div(
                                        [
                                            "Drag and Drop or ",
                                            html.A("Select Files"),
                                        ]
                                    ),
                                    style={
                                        "width": "100%",
                                        "height": "60px",
                                        "lineHeight": "60px",
                                        "borderWidth": "1px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "5px",
                                        "textAlign": "center",
                                    },
                                    multiple=True,
                                ),
                                html.Div(id="upload-error"),
                            ],
                        ),
                    ]
                ),
                html.Div(id="compound-preview"),
            ],
            style={"padding": "20px"},
        )

    def _setup_callbacks(self):
        """Setup component callbacks."""

        @self.app.callback(
            [
                Output("structure-error", "children"),
                Output("compound-preview", "children"),
            ],
            [Input("analyze-btn", "n_clicks")],
            [State("smiles-input", "value")],
        )
        def validate_structure(n_clicks: int, smiles: str) -> Tuple[str, html.Div]:
            """Validate input structure and update preview."""
            if not n_clicks or not smiles:
                return "", html.Div()

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return "Invalid SMILES string", html.Div()

                # Create compound data
                compound = CompoundData(
                    smiles=smiles,
                    mol=mol,
                    descriptors=self.descriptor_calc.calculate_descriptors(mol),
                )

                # Notify parent
                if self.on_compound_change:
                    self.on_compound_change(compound)

                # Show preview
                return "", self._create_preview(compound)

            except Exception as e:
                self.logger.error(f"Error validating structure: {str(e)}")
                return f"Error: {str(e)}", html.Div()

        @self.app.callback(
            [
                Output("upload-error", "children"),
                Output("smiles-input", "value"),
            ],
            [Input("file-upload", "contents")],
            [State("file-upload", "filename")],
        )
        def process_upload(
            contents: List[str], filenames: List[str]
        ) -> Tuple[str, str]:
            """Process uploaded files."""
            if not contents or not filenames:
                return "", ""

            try:
                # Process first file only for now
                content = contents[0]
                filename = filenames[0]

                if filename.endswith(".csv") or filename.endswith(".tsv"):
                    # Parse CSV/TSV
                    df = pd.read_csv(
                        content,
                        sep="\t" if filename.endswith(".tsv") else ",",
                    )
                    if "SMILES" not in df.columns:
                        return "File must contain SMILES column", ""
                    return "", str(df["SMILES"].iloc[0])

                elif filename.endswith(".sdf"):
                    # Parse SDF
                    suppl = Chem.SDMolSupplier(content)
                    if len(suppl) == 0:
                        return "No valid structures found in SDF", ""
                    mol = suppl[0]
                    return "", Chem.MolToSmiles(mol)

                else:
                    return "Unsupported file format", ""

            except Exception as e:
                self.logger.error(f"Error processing upload: {str(e)}")
                return f"Error: {str(e)}", ""

    def _create_preview(self, compound: CompoundData) -> html.Div:
        """Create compound preview element.

        Args:
            compound: Compound data

        Returns:
            Preview div element
        """
        try:
            # Get 2D depiction
            img = Chem.Draw.MolToImage(compound.mol)

            return html.Div(
                [
                    html.Img(
                        src=img,
                        style={
                            "maxWidth": "300px",
                            "maxHeight": "300px",
                        },
                    ),
                    html.Table(
                        [
                            html.Tr(
                                [
                                    html.Td("Molecular Weight"),
                                    html.Td(f"{compound.descriptors['MW']:.2f}"),
                                ]
                            ),
                            html.Tr(
                                [
                                    html.Td("LogP"),
                                    html.Td(f"{compound.descriptors['LogP']:.2f}"),
                                ]
                            ),
                            html.Tr(
                                [
                                    html.Td("TPSA"),
                                    html.Td(f"{compound.descriptors['TPSA']:.2f}"),
                                ]
                            ),
                        ]
                    ),
                ]
            )

        except Exception as e:
            self.logger.error(f"Error creating preview: {str(e)}")
            return html.Div()
