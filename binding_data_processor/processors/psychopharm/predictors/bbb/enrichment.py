"""BBB permeability predictor with web data enrichment.

This module provides the BBBPredictorWebEnriched class that adds:
1. Web data enrichment from scientific databases (ChEMBL, PubChem)
2. Community data enrichment (PsychonautWiki, Erowid, TripSit)
3. Social media monitoring (Reddit, Twitter, Bluesky)
4. LLM-based data extraction and analysis
5. Data validation and standardization
"""

import logging
from typing import Dict, List, Optional, Set, Any
import pandas as pd

from .....models.core import CompoundData
from ..types import PredictionResult
from .integration import BBBPredictorEnhanced

from binding_data_processor.data_sources.chembl import ChEMBLClient
from binding_data_processor.data_sources.pubchem import PubChemClient
from binding_data_processor.web_enrichment.swiss_client_enhanced import SwissClientEnhanced
from binding_data_processor.web_enrichment.community_client_enhanced import CommunityClientEnhanced
from binding_data_processor.web_enrichment.social_client_enhanced import SocialClientEnhanced
from binding_data_processor.web_enrichment.clients.web_search import WebSearchClientEnhanced
from binding_data_processor.web_enrichment.llm_utils import LLMProcessor


class BBBPredictorWebEnriched(BBBPredictorEnhanced):
    """BBB permeability predictor with web data enrichment."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        transporters: Optional[Dict[str, Set[str]]] = None,
        receptor_transporters: Optional[Dict[str, Set[str]]] = None,
        predictors: Optional[Dict[str, Any]] = None,
        web_clients: Optional[Dict[str, Any]] = None,
        from_pt: bool = True,  # Added parameter for loading from PyTorch weights
    ):
        """Initialize web-enriched BBB predictor.

        Args:
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            log_level: Logging level
            transporters: Optional custom transporter definitions
            receptor_transporters: Optional custom receptor transporter definitions
            predictors: Optional custom predictor instances
            web_clients: Optional custom web client instances
            from_pt: Whether to load models from PyTorch weights
        """
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
            transporters=transporters,
            receptor_transporters=receptor_transporters,
            predictors=predictors,
            from_pt=from_pt,  # Pass through to parent class
        )

        # Initialize web clients
        self.web_clients = web_clients or {
            "chembl": ChEMBLClient(cache_dir=cache_dir),
            "pubchem": PubChemClient(cache_dir=cache_dir),
            "swiss": SwissClientEnhanced(cache_dir=cache_dir),
            "community": CommunityClientEnhanced(cache_dir=cache_dir),
            "social": SocialClientEnhanced(cache_dir=cache_dir),
            "web_search": WebSearchClientEnhanced(cache_dir=cache_dir),
        }

        # Initialize LLM processor
        self.llm_processor = LLMProcessor(cache_dir=cache_dir)

        # Add web data columns to prediction history
        self.prediction_history = pd.concat(
            [
                self.prediction_history,
                pd.DataFrame(
                    columns=[
                        "chembl_data",
                        "pubchem_data",
                        "swiss_data",
                        "community_data",
                        "social_data",
                        "web_mentions",
                        "literature_data",
                        "patent_data",
                        "regulatory_data",
                        "clinical_data",
                        "mechanism_data",
                        "pharmacology_data",
                        "safety_data",
                        "interaction_data",
                    ]
                ),
            ],
            axis=1,
        )

        self.logger.info("BBBPredictorWebEnriched initialized successfully")

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate comprehensive BBB permeability predictions with web enrichment.

        Args:
            compound: Compound to predict BBB permeability for

        Returns:
            PredictionResult containing:
            - BBB permeability class
            - Confidence score
            - Supporting data including:
                - Permeability score
                - Transporter interactions
                - Receptor-mediated transport
                - Duration metrics
                - Abuse potential
                - Toxicity risks
                - Receptor binding
                - Psychoactive effects
                - Nootropic activity
                - Web data enrichment
        """
        self.logger.debug(f"Generating predictions for {compound.name}")
        try:
            # Get base predictions
            result = super().predict(compound)

            # Get web data
            web_data = self._get_web_data(compound)

            # Process web data with LLM
            processed_data = self._process_web_data(web_data, compound)

            # Update prediction result
            result = self._update_prediction_with_web_data(result, processed_data, compound)

            # Update prediction history
            self._update_web_prediction_history(processed_data, compound)

            return result

        except Exception as e:
            self.logger.error("Error predicting BBB properties with web enrichment", f"Compound {compound.name}: {str(e)}", exc_info=True)
            return result

    def _get_web_data(self, compound: CompoundData) -> Dict[str, Any]:
        """Get data from web sources."""
        web_data = {}

        # Get ChEMBL data
        web_data["chembl"] = self.web_clients["chembl"].get_compound_data(compound.smiles, compound.name)

        # Get PubChem data
        web_data["pubchem"] = self.web_clients["pubchem"].get_compound_data(compound.smiles, compound.name)

        # Get Swiss data
        web_data["swiss"] = self.web_clients["swiss"].get_compound_data(compound.smiles, compound.name)

        # Get community data
        web_data["community"] = self.web_clients["community"].get_compound_data(compound.smiles, compound.name)

        # Get social data
        web_data["social"] = self.web_clients["social"].get_compound_data(compound.smiles, compound.name)

        # Get web search data
        web_data["web_search"] = self.web_clients["web_search"].get_compound_data(compound.smiles, compound.name)

        return web_data

    def _process_web_data(self, web_data: Dict[str, Any], compound: CompoundData) -> Dict[str, Any]:
        """Process web data using LLM."""
        processed_data = {}

        # Extract BBB-related information
        processed_data["bbb_data"] = self.llm_processor.extract_bbb_data(web_data, compound)

        # Extract mechanism information
        processed_data["mechanism_data"] = self.llm_processor.extract_mechanism_data(web_data, compound)

        # Extract pharmacology information
        processed_data["pharmacology_data"] = self.llm_processor.extract_pharmacology_data(web_data, compound)

        # Extract safety information
        processed_data["safety_data"] = self.llm_processor.extract_safety_data(web_data, compound)

        # Extract interaction information
        processed_data["interaction_data"] = self.llm_processor.extract_interaction_data(web_data, compound)

        return processed_data

    def _update_prediction_with_web_data(
        self,
        result: PredictionResult,
        web_data: Dict[str, Any],
        compound: CompoundData,
    ) -> PredictionResult:
        """Update prediction result with web data."""
        # Add web data to supporting data
        result.supporting_data.update(
            {
                "web_data": {
                    "bbb_data": web_data["bbb_data"],
                    "mechanism_data": web_data["mechanism_data"],
                    "pharmacology_data": web_data["pharmacology_data"],
                    "safety_data": web_data["safety_data"],
                    "interaction_data": web_data["interaction_data"],
                }
            }
        )

        # Adjust confidence based on web data support
        if web_data["bbb_data"].get("supports_prediction", False):
            result.confidence = min(1.0, result.confidence * 1.2)

        return result

    def _update_web_prediction_history(
        self,
        web_data: Dict[str, Any],
        compound: CompoundData,
    ) -> None:
        """Update prediction history with web data."""
        # Find the latest prediction for this compound
        mask = self.prediction_history["compound_name"] == compound.name
        if not mask.any():
            return

        latest_idx = self.prediction_history.index[-1]

        # Update web data columns
        self.prediction_history.at[latest_idx, "chembl_data"] = str(web_data.get("chembl", {}))
        self.prediction_history.at[latest_idx, "pubchem_data"] = str(web_data.get("pubchem", {}))
        self.prediction_history.at[latest_idx, "swiss_data"] = str(web_data.get("swiss", {}))
        self.prediction_history.at[latest_idx, "community_data"] = str(web_data.get("community", {}))
        self.prediction_history.at[latest_idx, "social_data"] = str(web_data.get("social", {}))
        self.prediction_history.at[latest_idx, "web_mentions"] = str(web_data.get("web_search", {}))
        self.prediction_history.at[latest_idx, "mechanism_data"] = str(web_data.get("mechanism_data", {}))
        self.prediction_history.at[latest_idx, "pharmacology_data"] = str(web_data.get("pharmacology_data", {}))
        self.prediction_history.at[latest_idx, "safety_data"] = str(web_data.get("safety_data", {}))
        self.prediction_history.at[latest_idx, "interaction_data"] = str(web_data.get("interaction_data", {}))

    def export_predictions(
        self,
        output_path: str,
        columns: Optional[List[str]] = None,
        include_supporting_data: bool = False,
        include_web_data: bool = False,
    ) -> None:
        """Export prediction history to TSV file with web data options.

        Args:
            output_path: Path to save TSV file
            columns: Optional list of columns to include
            include_supporting_data: Whether to include supporting data
            include_web_data: Whether to include web data
        """
        if columns is None:
            columns = self.prediction_history.columns

        # Create export DataFrame
        export_df = self.prediction_history[columns].copy()

        # Add supporting data if requested
        if include_supporting_data:
            for predictor in self.predictors:
                predictor_history = self.predictors[predictor].prediction_history
                if not predictor_history.empty:
                    # Merge on compound_name and timestamp
                    export_df = pd.merge(
                        export_df,
                        predictor_history,
                        on=["compound_name", "timestamp"],
                        how="left",
                        suffixes=("", f"_{predictor}"),
                    )

        # Add web data if requested
        if include_web_data:
            web_columns = [
                "chembl_data",
                "pubchem_data",
                "swiss_data",
                "community_data",
                "social_data",
                "web_mentions",
                "literature_data",
                "patent_data",
                "regulatory_data",
                "clinical_data",
                "mechanism_data",
                "pharmacology_data",
                "safety_data",
                "interaction_data",
            ]
            export_df = export_df.join(self.prediction_history[web_columns], how="left")

        # Export to TSV
        export_df.to_csv(
            output_path,
            sep="\t",
            index=False,
        )
        self.logger.info(f"Exported predictions to {output_path}")
