"""Web application serving the ChemData dashboard and API.

This application provides:
1. Interactive dashboard for data exploration
2. REST API for programmatic access
3. ML-powered predictions and analysis
4. Advanced compound search and filtering
5. Structure visualization
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Any
from io import BytesIO, StringIO
from datetime import datetime
import logging
from concurrent.futures import ThreadPoolExecutor

import dash
from dash import html
import dash_bootstrap_components as dbc
from flask import Flask, render_template, jsonify, request, send_file, Response
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, SDWriter
from tqdm import tqdm

from binding_data_processor.models.compound import CompoundData
from binding_data_processor.web.dashboard.main import MainDashboard
from binding_data_processor.pipeline import PipelineManager
from binding_data_processor.processors.structure.ml.predictors.abuse import (
    AbusePotentialPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.toxicity import (
    ToxicityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.activity import (
    ActivityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.affinity import (
    AffinityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.psychoactive import (
    PsychoactivePredictor,
)
from binding_data_processor.processors.structure.ml.predictors.nootropic import (
    NootropicPredictor,
)
from web_enrichment.data_sources.community import CommunityClient
from web_enrichment.data_sources.swiss import SwissClient
from web_enrichment.data_sources.social import SocialDataHarvester
from web_enrichment.http_client import HttpClient


class ChemDataApp:
    """Web application combining dashboard and API functionality."""

    def __init__(
        self,
        data_dir: str = "../data",
        model_dir: str = "../models",
        debug: bool = False,
        n_workers: int = 4,
        batch_size: int = 100,
        checkpoint_interval: int = 1000,
        max_retries: int = 3,
        cache_dir: Optional[str] = None,
        reddit_client_id: Optional[str] = None,
        reddit_client_secret: Optional[str] = None,
        twitter_api_key: Optional[str] = None,
        twitter_api_secret: Optional[str] = None,
        discord_token: Optional[str] = None,
        bluesky_handle: Optional[str] = None,
        bluesky_password: Optional[str] = None,
    ):
        """Initialize web application.

        Args:
            data_dir: Directory containing data files
            model_dir: Directory containing ML models
            debug: Whether to run in debug mode
            n_workers: Number of worker threads
            batch_size: Batch size for processing
            checkpoint_interval: Save checkpoint every N compounds
            max_retries: Maximum number of retries for failed operations
            cache_dir: Optional custom cache directory
            reddit_client_id: Optional Reddit API client ID
            reddit_client_secret: Optional Reddit API client secret
            twitter_api_key: Optional Twitter API key
            twitter_api_secret: Optional Twitter API secret
            discord_token: Optional Discord bot token
            bluesky_handle: Optional Bluesky handle
            bluesky_password: Optional Bluesky password
        """
        # Initialize Flask server
        self.server = Flask(__name__)
        self.server.config.update(
            {
                "DATA_DIR": Path(data_dir).resolve(),
                "MODEL_DIR": Path(model_dir).resolve(),
                "DEBUG": debug,
            }
        )

        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        self.logger.addHandler(handler)

        # Initialize Dash app
        self.app = dash.Dash(
            __name__,
            server=self.server,
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            suppress_callback_exceptions=True,
            url_base_pathname="/dashboard/",
        )
        self.app.title = "ChemData Dashboard"

        # Initialize dashboard
        self.dashboard = MainDashboard()
        self.app.layout = self.dashboard.layout

        # Initialize pipeline manager
        self.pipeline = PipelineManager(
            data_dir=str(self.server.config["DATA_DIR"]),
            model_dir=str(self.server.config["MODEL_DIR"]),
            n_workers=n_workers,
            batch_size=batch_size,
            checkpoint_interval=checkpoint_interval,
            max_retries=max_retries,
            cache_dir=cache_dir,
            reddit_client_id=reddit_client_id,
            reddit_client_secret=reddit_client_secret,
            twitter_api_key=twitter_api_key,
            twitter_api_secret=twitter_api_secret,
            discord_token=discord_token,
            bluesky_handle=bluesky_handle,
            bluesky_password=bluesky_password,
        )

        # Initialize data browser
        self.browser = DataBrowser(
            data_dir=str(self.server.config["DATA_DIR"]),
            model_dir=str(self.server.config["MODEL_DIR"]),
            pipeline=self.pipeline,
            n_workers=n_workers,
            cache_dir=cache_dir,
        )

        # Register routes and callbacks
        self._register_routes()
        self.dashboard.register_callbacks()

        # Load initial data
        self._load_data()

    def _register_routes(self) -> None:
        """Register Flask routes."""

        @self.server.route("/")
        def index():
            """Render main page."""
            return render_template(
                "index.html",
                target_types=self.browser.get_target_types(),
                activity_types=self.browser.get_activity_types(),
                legal_statuses=self.browser.get_legal_statuses(),
                column_groups=self.browser.get_column_groups(),
                stats=self.browser.get_stats(),
            )

        @self.server.route("/api/compounds")
        def get_compounds():
            """API endpoint for compound data."""
            try:
                # Get filter parameters
                search = request.args.get("search", "")
                target = request.args.get("target", "")
                structure = request.args.get("structure", "")
                similarity = float(request.args.get("similarity", 0.7))
                activity = request.args.get("activity", "")
                legal = request.args.get("legal", "")
                high_toxicity = (
                    request.args.get("high_toxicity", "false").lower() == "true"
                )
                high_abuse = request.args.get("high_abuse", "false").lower() == "true"
                columns = request.args.getlist("columns")
                page = int(request.args.get("page", 1))
                per_page = int(request.args.get("per_page", 50))

                # Get compounds with progress tracking
                with tqdm(
                    desc="Getting compounds",
                    unit="compounds",
                    disable=not self.server.config["DEBUG"],
                ) as pbar:
                    data = self.browser.get_compounds(
                        search=search,
                        target_type=target,
                        structure_search=structure,
                        similarity_threshold=similarity,
                        activity_type=activity,
                        legal_status=legal,
                        high_toxicity=high_toxicity,
                        high_abuse=high_abuse,
                        selected_columns=columns,
                        page=page,
                        per_page=per_page,
                        progress_callback=pbar.update,
                    )

                # Format response for DataTables
                response = {
                    "draw": int(request.args.get("draw", 1)),
                    "recordsTotal": data["total"],
                    "recordsFiltered": data["total"],
                    "data": data["compounds"],
                }

                return jsonify(response)

            except Exception as e:
                self.logger.error(f"Error getting compounds: {str(e)}")
                return jsonify(
                    {
                        "draw": int(request.args.get("draw", 1)),
                        "recordsTotal": 0,
                        "recordsFiltered": 0,
                        "data": [],
                        "error": str(e),
                    }
                )

        @self.server.route("/api/structure/<smiles>")
        def get_structure_image(smiles: str):
            """API endpoint for structure images."""
            try:
                img_data = self.browser.structure_to_image(smiles)
                if img_data:
                    return send_file(BytesIO(img_data), mimetype="image/png")
                return "", 404
            except Exception as e:
                self.logger.error(f"Error generating structure image: {str(e)}")
                return "", 500

        @self.server.route("/api/targets")
        def get_targets():
            """API endpoint for target types."""
            try:
                return jsonify(self.browser.get_target_types())
            except Exception as e:
                self.logger.error(f"Error getting target types: {str(e)}")
                return jsonify([])

        @self.server.route("/api/activities")
        def get_activities():
            """API endpoint for activity types."""
            try:
                return jsonify(self.browser.get_activity_types())
            except Exception as e:
                self.logger.error(f"Error getting activity types: {str(e)}")
                return jsonify([])

        @self.server.route("/api/legal-statuses")
        def get_legal_statuses():
            """API endpoint for legal statuses."""
            try:
                return jsonify(self.browser.get_legal_statuses())
            except Exception as e:
                self.logger.error(f"Error getting legal statuses: {str(e)}")
                return jsonify([])

        @self.server.route("/api/compounds/<compound_id>/sdf")
        def export_sdf(compound_id: str):
            """Export compound as SDF file."""
            try:
                compound = self.browser.get_compound_by_id(compound_id)
                if not compound:
                    return "", 404

                sio = StringIO()
                writer = SDWriter(sio)
                mol = Chem.MolFromSmiles(compound.smiles)
                if mol:
                    for key, value in compound.to_dict().items():
                        if value:
                            mol.SetProp(key, str(value))
                    writer.write(mol)
                    writer.close()

                    return Response(
                        sio.getvalue(),
                        mimetype="chemical/x-mdl-sdfile",
                        headers={
                            "Content-Disposition": f"attachment;filename={compound.name}.sdf"
                        },
                    )
                return "", 404

            except Exception as e:
                self.logger.error(f"Error exporting SDF: {str(e)}")
                return "", 500

        @self.server.route("/api/compounds/<compound_id>/csv")
        def export_csv(compound_id: str):
            """Export compound as CSV file."""
            try:
                compound = self.browser.get_compound_by_id(compound_id)
                if not compound:
                    return "", 404

                df = pd.DataFrame([compound.to_dict()])
                csv_data = df.to_csv(index=False)

                return Response(
                    csv_data,
                    mimetype="text/csv",
                    headers={
                        "Content-Disposition": f"attachment;filename={compound.name}.csv"
                    },
                )

            except Exception as e:
                self.logger.error(f"Error exporting CSV: {str(e)}")
                return "", 500

        @self.server.route("/api/compounds/<compound_id>/report")
        def export_report(compound_id: str):
            """Export compound as detailed report."""
            try:
                compound = self.browser.get_compound_by_id(compound_id)
                if not compound:
                    return "", 404

                # Generate HTML report
                report = self.browser.generate_compound_report(compound)

                return Response(
                    report,
                    mimetype="text/html",
                    headers={
                        "Content-Disposition": f"attachment;filename={compound.name}_report.html"
                    },
                )

            except Exception as e:
                self.logger.error(f"Error generating report: {str(e)}")
                return "", 500

        @self.server.route("/api/compounds/export/csv")
        def export_filtered_csv():
            """Export filtered compounds as CSV."""
            try:
                compounds = self._get_filtered_compounds()
                if not compounds:
                    return "", 404

                df = pd.DataFrame([c.to_dict() for c in compounds])
                csv_data = df.to_csv(index=False)

                return Response(
                    csv_data,
                    mimetype="text/csv",
                    headers={
                        "Content-Disposition": "attachment;filename=compounds.csv"
                    },
                )

            except Exception as e:
                self.logger.error(f"Error exporting filtered CSV: {str(e)}")
                return "", 500

        @self.server.route("/api/compounds/export/sdf")
        def export_filtered_sdf():
            """Export filtered compounds as SDF."""
            try:
                compounds = self._get_filtered_compounds()
                if not compounds:
                    return "", 404

                sio = StringIO()
                writer = SDWriter(sio)

                for compound in compounds:
                    mol = Chem.MolFromSmiles(compound.smiles)
                    if mol:
                        for key, value in compound.to_dict().items():
                            if value:
                                mol.SetProp(key, str(value))
                        writer.write(mol)

                writer.close()

                return Response(
                    sio.getvalue(),
                    mimetype="chemical/x-mdl-sdfile",
                    headers={
                        "Content-Disposition": "attachment;filename=compounds.sdf"
                    },
                )

            except Exception as e:
                self.logger.error(f"Error exporting filtered SDF: {str(e)}")
                return "", 500

        @self.server.route("/api/compounds/export/json")
        def export_filtered_json():
            """Export filtered compounds as JSON."""
            try:
                compounds = self._get_filtered_compounds()
                if not compounds:
                    return "", 404

                data = [c.to_dict() for c in compounds]

                return Response(
                    jsonify(data).get_data(),
                    mimetype="application/json",
                    headers={
                        "Content-Disposition": "attachment;filename=compounds.json"
                    },
                )

            except Exception as e:
                self.logger.error(f"Error exporting filtered JSON: {str(e)}")
                return "", 500

    def _get_filtered_compounds(self) -> Optional[List[CompoundData]]:
        """Get compounds based on current filters."""
        try:
            search = request.args.get("search", "")
            target = request.args.get("target", "")
            structure = request.args.get("structure", "")
            similarity = float(request.args.get("similarity", 0.7))
            activity = request.args.get("activity", "")
            legal = request.args.get("legal", "")
            high_toxicity = request.args.get("high_toxicity", "false").lower() == "true"
            high_abuse = request.args.get("high_abuse", "false").lower() == "true"

            data = self.browser.get_compounds(
                search=search,
                target_type=target,
                structure_search=structure,
                similarity_threshold=similarity,
                activity_type=activity,
                legal_status=legal,
                high_toxicity=high_toxicity,
                high_abuse=high_abuse,
                page=1,
                per_page=10000,  # Large number to get all compounds
            )

            return data.get("compounds", [])

        except Exception as e:
            self.logger.error(f"Error getting filtered compounds: {str(e)}")
            return None

    def _load_data(self) -> None:
        """Load and process compound data."""
        try:
            # Process all data sources with progress tracking
            with tqdm(
                desc="Processing data sources",
                unit="sources",
                disable=not self.server.config["DEBUG"],
            ) as pbar:
                stats = self.pipeline.run_pipeline(
                    skip_predictions=False,
                    skip_web_data=False,
                    use_cache=True,
                    progress_callback=pbar.update,
                )

            self.logger.info(
                f"Loaded {stats['total_compounds']} compounds:\n"
                f"- From BindingDB: {stats['bindingdb_compounds']}\n"
                f"- From web sources: {stats['web_compounds']}\n"
                f"- From social media: {stats['social_compounds']}\n"
                f"- With predictions: {stats['with_predictions']}\n"
                f"- With web data: {stats['with_web_data']}"
            )

            # Update dashboard with processed data
            self.dashboard.update_data(self.pipeline.compounds)

        except Exception as e:
            self.logger.error(f"Error loading data: {str(e)}")

    def run(
        self,
        host: str = "0.0.0.0",
        port: int = 8050,
        debug: Optional[bool] = None,
    ) -> None:
        """Run the web application.

        Args:
            host: Host to run on
            port: Port to run on
            debug: Whether to run in debug mode
        """
        if debug is None:
            debug = self.server.config["DEBUG"]

        self.app.run_server(
            host=host,
            port=port,
            debug=debug,
        )


class DataBrowser:
    """Browser for compound data with advanced search capabilities."""

    def __init__(
        self,
        data_dir: str = "../data",
        model_dir: str = "../models",
        pipeline: Optional[PipelineManager] = None,
        n_workers: int = 4,
        cache_dir: Optional[str] = None,
    ):
        """Initialize browser.

        Args:
            data_dir: Directory containing data files
            model_dir: Directory containing ML models
            pipeline: Optional pipeline manager instance
            n_workers: Number of worker threads
            cache_dir: Optional cache directory
        """
        self.data_dir = data_dir
        self.model_dir = model_dir
        self.pipeline = pipeline
        self.n_workers = n_workers
        self.cache_dir = cache_dir

        # Initialize logger
        self.logger = logging.getLogger(__name__)

        # Get latest data file
        self.current_file = self._get_latest_tsv()
        self.df = pd.read_csv(self.current_file, sep="\t")

        # Convert SMILES and InChI to strings
        self.df["smiles"] = self.df["smiles"].fillna("").astype(str)
        self.df["inchi"] = self.df["inchi"].fillna("").astype(str)

        # Initialize predictors
        self.logger.info("Initializing ML predictors...")
        self.predictors = {
            "abuse": AbusePotentialPredictor(
                model_dir=os.path.join(model_dir, "abuse")
            ),
            "toxicity": ToxicityPredictor(
                model_dir=os.path.join(model_dir, "toxicity")
            ),
            "activity": ActivityPredictor(
                model_dir=os.path.join(model_dir, "activity")
            ),
            "affinity": AffinityPredictor(
                model_dir=os.path.join(model_dir, "affinity")
            ),
            "psychoactive": PsychoactivePredictor(
                model_dir=os.path.join(model_dir, "psychoactive")
            ),
            "nootropic": NootropicPredictor(
                model_dir=os.path.join(model_dir, "nootropic")
            ),
        }

        # Initialize thread pool
        self.executor = ThreadPoolExecutor(max_workers=n_workers)

        # Get column info
        self.columns = self.df.columns.tolist()
        self.url_columns = [col for col in self.columns if col.endswith("_url")]

        # Define column groups
        self.column_groups = {
            "identifiers": [
                "name",
                "cas_number",
                "smiles",
                "inchi",
                "molecular_formula",
            ],
            "properties": [
                "molecular_weight",
                "logp",
                "hbd",
                "hba",
                "tpsa",
                "rotatable_bonds",
            ],
            "activity": [
                "target_type",
                "activity_type",
                "activity_value",
                "activity_unit",
            ],
            "pharmacology": ["mechanism", "primary_target", "secondary_targets"],
            "predictions": [
                "toxicity_predictions",
                "abuse_potential",
                "binding_predictions",
                "activity_predictions",
                "psychoactive_predictions",
                "nootropic_predictions",
            ],
            "legal": ["legal_status", "scheduling", "controlled_status"],
            "social": [
                "community_data",
                "social_data",
                "experience_reports",
                "sentiment_data",
            ],
            "references": self.url_columns,
        }

    def _get_latest_tsv(self) -> str:
        """Get path of most recent TSV file."""
        tsv_files = [
            f
            for f in os.listdir(self.data_dir)
            if f.startswith("receptor_compounds_") and f.endswith(".tsv")
        ]
        if not tsv_files:
            # Create an empty TSV file with required columns
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"receptor_compounds_{timestamp}.tsv"
            filepath = os.path.join(self.data_dir, filename)
            pd.DataFrame(columns=self._get_default_columns()).to_csv(
                filepath, sep="\t", index=False
            )
            return filepath

        latest = max(tsv_files)
        return os.path.join(self.data_dir, latest)

    def _get_default_columns(self) -> List[str]:
        """Get default columns for new TSV files."""
        return [
            # Basic identifiers
            "name",
            "smiles",
            "inchi",
            "cas_number",
            "molecular_formula",
            # Properties
            "molecular_weight",
            "logp",
            "hbd",
            "hba",
            "tpsa",
            "rotatable_bonds",
            # Activity data
            "target_type",
            "activity_type",
            "activity_value",
            "activity_unit",
            # Pharmacology
            "mechanism",
            "primary_target",
            "secondary_targets",
            # Legal status
            "legal_status",
            "scheduling",
            "controlled_status",
            # References
            "reference_urls",
            "patent_data",
            "literature_data",
            # ML predictions
            "toxicity_predictions",
            "abuse_potential",
            "binding_predictions",
            "activity_predictions",
            "psychoactive_predictions",
            "nootropic_predictions",
            # Web-enriched data
            "community_data",
            "social_data",
            "experience_reports",
            "sentiment_data",
            "swiss_data",
            "regulatory_data",
        ]

    def get_compounds(
        self,
        search: str = "",
        target_type: str = "",
        structure_search: str = "",
        similarity_threshold: float = 0.7,
        activity_type: str = "",
        legal_status: str = "",
        high_toxicity: bool = False,
        high_abuse: bool = False,
        selected_columns: Optional[List[str]] = None,
        page: int = 1,
        per_page: int = 50,
        progress_callback: Optional[callable] = None,
    ) -> Dict[str, Any]:
        """Get paginated compound data with optional filtering.

        Args:
            search: Text search query
            target_type: Filter by target type
            structure_search: SMILES for structure similarity search
            similarity_threshold: Similarity threshold for structure search
            activity_type: Filter by activity type
            legal_status: Filter by legal status
            high_toxicity: Filter for high toxicity predictions
            high_abuse: Filter for high abuse potential
            selected_columns: Optional list of columns to include
            page: Page number
            per_page: Items per page
            progress_callback: Optional callback for progress updates

        Returns:
            Dictionary with compounds and pagination info
        """
        try:
            df = self.df

            # Apply filters
            if search:
                conditions = []
                for col in [
                    "name",
                    "cas_number",
                    "molecular_formula",
                    "smiles",
                    "inchi",
                ]:
                    if col in df.columns:
                        conditions.append(
                            df[col].str.contains(search, case=False, na=False)
                        )
                if conditions:
                    df = df[pd.concat(conditions, axis=1).any(axis=1)]

            if target_type and "target_type" in df.columns:
                df = df[df["target_type"] == target_type]

            if structure_search and "smiles" in df.columns:
                query_mol = Chem.MolFromSmiles(structure_search)
                if query_mol:
                    query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2)

                    def calc_similarity(smiles):
                        try:
                            mol = Chem.MolFromSmiles(smiles)
                            if mol:
                                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
                                return DataStructs.TanimotoSimilarity(query_fp, fp)
                        except Exception as e:
                            self.logger.debug(f"Error calculating similarity: {e}")
                        return 0

                    df["similarity"] = df["smiles"].apply(calc_similarity)
                    df = df[df["similarity"] >= similarity_threshold]

            if activity_type and "activity_type" in df.columns:
                df = df[df["activity_type"] == activity_type]

            if legal_status and "legal_status" in df.columns:
                df = df[df["legal_status"] == legal_status]

            # Apply ML filters
            if high_toxicity and "toxicity_predictions" in df.columns:
                try:
                    df = df[
                        df["toxicity_predictions"].apply(
                            lambda x: isinstance(x, dict)
                            and any(
                                pred.get("probability", 0) > 0.7 for pred in x.values()
                            )
                        )
                    ]
                except Exception as e:
                    self.logger.error(f"Error applying toxicity filter: {e}")

            if high_abuse and "abuse_potential" in df.columns:
                try:
                    df = df[
                        df["abuse_potential"].apply(
                            lambda x: isinstance(x, dict)
                            and any(
                                pred.get("probability", 0) > 0.7 for pred in x.values()
                            )
                        )
                    ]
                except Exception as e:
                    self.logger.error(f"Error applying abuse filter: {e}")

            # Calculate pagination
            total = len(df)
            start = (page - 1) * per_page
            end = start + per_page

            # Get page of data
            page_df = df.iloc[start:end]

            # Convert to CompoundData objects and then to dicts
            compounds = []
            for _, row in tqdm(
                page_df.iterrows(),
                desc="Converting compounds",
                total=len(page_df),
                disable=not progress_callback,
            ):
                try:
                    compound = CompoundData(
                        name=row["name"],
                        smiles=row["smiles"],
                        inchi=row.get("inchi"),
                        cas_number=row.get("cas_number"),
                    )
                    compound_dict = compound.to_dict()

                    # Add predictions
                    for pred_type in self.predictors.keys():
                        pred_col = f"{pred_type}_predictions"
                        if pred_col in row and pd.notna(row[pred_col]):
                            compound_dict[pred_col] = row[pred_col]

                    # Add web data
                    for data_type in ["community", "social", "swiss"]:
                        data_col = f"{data_type}_data"
                        if data_col in row and pd.notna(row[data_col]):
                            compound_dict[data_col] = row[data_col]

                    # Add URL links
                    compound_dict["urls"] = {
                        col: url
                        for col, url in row.items()
                        if col in self.url_columns and pd.notna(url)
                    }

                    compounds.append(compound_dict)

                    if progress_callback:
                        progress_callback(1)

                except Exception as e:
                    self.logger.error(
                        f"Error converting compound {row.get('name')}: {e}"
                    )
                    continue

            return {
                "compounds": compounds,
                "total": total,
                "page": page,
                "per_page": per_page,
                "pages": (total + per_page - 1) // per_page,
            }

        except Exception as e:
            self.logger.error(f"Error getting compounds: {e}")
            return {
                "compounds": [],
                "total": 0,
                "page": page,
                "per_page": per_page,
                "pages": 0,
                "error": str(e),
            }

    def get_compound_by_id(self, compound_id: str) -> Optional[CompoundData]:
        """Get compound by ID."""
        try:
            row = self.df[self.df["name"] == compound_id].iloc[0]
            compound = CompoundData(
                name=row["name"],
                smiles=row["smiles"],
                inchi=row.get("inchi"),
                cas_number=row.get("cas_number"),
            )

            # Add predictions
            for pred_type in self.predictors.keys():
                pred_col = f"{pred_type}_predictions"
                if pred_col in row and pd.notna(row[pred_col]):
                    setattr(compound, f"{pred_type}_predictions", row[pred_col])

            # Add web data
            for data_type in ["community", "social", "swiss"]:
                data_col = f"{data_type}_data"
                if data_col in row and pd.notna(row[data_col]):
                    setattr(compound, f"{data_type}_data", row[data_col])

            return compound

        except (IndexError, KeyError) as e:
            self.logger.error(f"Error getting compound {compound_id}: {e}")
            return None

    def generate_compound_report(self, compound: CompoundData) -> str:
        """Generate detailed HTML report for compound."""
        from web.report_generator import ReportGenerator

        # Get current versions
        data_version = self._get_data_version()
        model_version = self._get_model_version()

        # Generate report
        generator = ReportGenerator(
            template_dir=os.path.join(os.path.dirname(__file__), "templates")
        )
        return generator.generate_report(
            compound=compound,
            data_version=data_version,
            model_version=model_version,
        )

    def _get_data_version(self) -> str:
        """Get current data version from latest TSV file."""
        try:
            filename = os.path.basename(self.current_file)
            # Extract timestamp from filename (format: receptor_compounds_YYYYMMDD_HHMMSS.tsv)
            timestamp = filename.split("_")[2:4]  # Get YYYYMMDD and HHMMSS parts
            return "_".join(timestamp).replace(".tsv", "")
        except Exception:
            return "Unknown"

    def _get_model_version(self) -> str:
        """Get current ML model versions."""
        versions = []
        try:
            # Get version from model files if available
            for predictor in [
                self.abuse_predictor,
                self.toxicity_predictor,
                self.activity_predictor,
                self.affinity_predictor,
            ]:
                if hasattr(predictor, "version"):
                    versions.append(
                        f"{predictor.__class__.__name__}: {predictor.version}"
                    )

            if versions:
                return "; ".join(versions)
        except Exception:
            pass
        return "Unknown"

    def get_target_types(self) -> List[str]:
        """Get list of unique target types."""
        if "target_type" in self.df.columns:
            return sorted(self.df["target_type"].unique().tolist())
        return []

    def get_activity_types(self) -> List[str]:
        """Get list of unique activity types."""
        if "activity_type" in self.df.columns:
            return sorted(self.df["activity_type"].unique().tolist())
        return []

    def get_legal_statuses(self) -> List[str]:
        """Get list of unique legal statuses."""
        if "legal_status" in self.df.columns:
            return sorted(self.df["legal_status"].unique().tolist())
        return []

    def get_column_groups(self) -> Dict[str, List[str]]:
        """Get column groups."""
        return self.column_groups

    def get_stats(self) -> Dict[str, Any]:
        """Get summary statistics."""
        stats = {
            "total_compounds": len(self.df),
            "with_activity": 0,
            "by_target": {},
            "by_activity": {},
            "by_legal_status": {},
        }

        if "activity_data" in self.df.columns:
            stats["with_activity"] = len(
                self.df[self.df["activity_data"].str.len() > 0]
            )

        if "target_type" in self.df.columns:
            stats["by_target"] = self.df["target_type"].value_counts().to_dict()

        if "activity_type" in self.df.columns:
            stats["by_activity"] = self.df["activity_type"].value_counts().to_dict()

        if "legal_status" in self.df.columns:
            stats["by_legal_status"] = self.df["legal_status"].value_counts().to_dict()

        return stats

    def structure_to_image(self, smiles: str) -> Optional[bytes]:
        """Convert SMILES to PNG image."""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            img_io = BytesIO()
            img.save(img_io, "PNG")
            img_io.seek(0)
            return img_io.getvalue()
        return None


def create_app(
    data_dir: Optional[str] = None,
    model_dir: Optional[str] = None,
    debug: bool = False,
) -> ChemDataApp:
    """Create ChemData web application.

    Args:
        data_dir: Optional directory containing data files
        model_dir: Optional directory containing ML models
        debug: Whether to run in debug mode

    Returns:
        ChemData web application
    """
    if data_dir is None:
        data_dir = os.getenv("CHEMDATA_DATA_DIR", "../data")
    if model_dir is None:
        model_dir = os.getenv("CHEMDATA_MODEL_DIR", "../models")

    return ChemDataApp(
        data_dir=data_dir,
        model_dir=model_dir,
        debug=debug,
    )


def main():
    """Run ChemData web application."""
    import argparse

    parser = argparse.ArgumentParser(description="Run ChemData web application")
    parser.add_argument(
        "--data-dir",
        type=str,
        help="Directory containing data files",
    )
    parser.add_argument(
        "--model-dir",
        type=str,
        help="Directory containing ML models",
    )
    parser.add_argument(
        "--host",
        type=str,
        default="0.0.0.0",
        help="Host to run on",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8050,
        help="Port to run on",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run in debug mode",
    )
    args = parser.parse_args()

    app = create_app(
        data_dir=args.data_dir,
        model_dir=args.model_dir,
        debug=args.debug,
    )
    app.run(
        host=args.host,
        port=args.port,
    )


if __name__ == "__main__":
    main()
