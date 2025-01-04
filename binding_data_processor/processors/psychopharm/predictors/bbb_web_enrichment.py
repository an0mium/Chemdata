"""Web data enrichment for BBB permeability prediction.

This module provides functionality to:
1. Harvest additional BBB permeability data from web sources
2. Extract relevant information using LLMs
3. Validate and standardize web data
4. Integrate web data with ML predictions
5. Track amino acid transporter data
"""

import logging
from typing import Dict, List, Optional, Set, Any
import pandas as pd

from ....models.core import CompoundData
from ..base import PredictionResult
from .bbb_enhanced import BBBPredictorEnhanced
from web_enrichment.data_sources.chembl import ChemblClient
from web_enrichment.data_sources.pubchem import PubchemClient
from web_enrichment.data_sources.swiss import SwissClient
from web_enrichment.data_sources.community import CommunityDataClient
from web_enrichment.data_sources.social import SocialDataClient
from web_enrichment.data_sources.web_search import WebSearchClient
from web_enrichment.llm_utils import LLMProcessor


class BBBPredictorWebEnriched(BBBPredictorEnhanced):
    """BBB permeability predictor with web data enrichment."""

    # Keywords for amino acid transporter detection
    AMINO_ACID_KEYWORDS = {
        "lat1": {
            "lat1", "lat-1", "l-type amino acid transporter 1", 
            "slc7a5", "solute carrier family 7 member 5"
        },
        "lat2": {
            "lat2", "lat-2", "l-type amino acid transporter 2",
            "slc7a8", "solute carrier family 7 member 8"
        },
        "asct1": {
            "asct1", "asct-1", "alanine serine cysteine transporter 1",
            "slc1a4", "solute carrier family 1 member 4"
        },
        "asct2": {
            "asct2", "asct-2", "alanine serine cysteine transporter 2",
            "slc1a5", "solute carrier family 1 member 5"
        },
        "b0at1": {
            "b0at1", "b0at-1", "slc6a19", 
            "solute carrier family 6 member 19"
        },
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        transporters: Optional[Dict[str, Set[str]]] = None,
        receptor_transporters: Optional[Dict[str, Set[str]]] = None,
        predictors: Optional[Dict[str, Any]] = None,
        web_clients: Optional[Dict[str, Any]] = None,
    ):
        """Initialize web-enriched BBB predictor."""
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
            transporters=transporters,
            receptor_transporters=receptor_transporters,
            predictors=predictors,
        )

        # Initialize web clients
        self.web_clients = web_clients or {
            "chembl": ChemblClient(cache_dir=cache_dir),
            "pubchem": PubchemClient(cache_dir=cache_dir),
            "swiss": SwissClient(cache_dir=cache_dir),
            "community": CommunityDataClient(cache_dir=cache_dir),
            "social": SocialDataClient(cache_dir=cache_dir),
            "web_search": WebSearchClient(cache_dir=cache_dir),
        }

        # Initialize LLM processor
        self.llm_processor = LLMProcessor(cache_dir=cache_dir)

        # Add web data columns to prediction history
        self.prediction_history = pd.concat([
            self.prediction_history,
            pd.DataFrame(columns=[
                'chembl_data',
                'pubchem_data',
                'swiss_data',
                'community_data',
                'social_data',
                'web_mentions',
                'literature_data',
                'patent_data',
                'regulatory_data',
                'clinical_data',
                'mechanism_data',
                'pharmacology_data',
                'safety_data',
                'interaction_data',
                # Add amino acid transporter columns
                'lat1_literature',
                'lat1_evidence',
                'lat1_confidence_web',
                'amino_acid_transporters',
                'transporter_interactions',
                'substrate_evidence',
            ]),
        ], axis=1)

        self.logger.info("BBBPredictorWebEnriched initialized successfully")

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate BBB permeability predictions with web enrichment."""
        self.logger.debug(f"Generating predictions for {compound.name}")
        try:
            # Get base predictions
            base_result = super().predict(compound)
            
            # Enrich with web data
            web_data = self._get_web_data(compound)
            
            # Process web data with LLM
            processed_data = self._process_web_data(web_data, compound)
            
            # Update prediction result
            result = self._update_prediction_with_web_data(
                base_result, processed_data, compound
            )
            
            # Update prediction history
            self._update_web_prediction_history(processed_data, compound)
            
            return result

        except Exception as e:
            self.logger.error(
                "Error predicting BBB properties with web enrichment",
                f"Compound {compound.name}: {str(e)}",
                exc_info=True
            )
            return base_result

    def _get_web_data(self, compound: CompoundData) -> Dict[str, Any]:
        """Get data from web sources."""
        web_data = {}
        
        # Get ChEMBL data
        web_data["chembl"] = self.web_clients["chembl"].get_compound_data(
            compound.smiles, compound.name
        )
        
        # Get PubChem data
        web_data["pubchem"] = self.web_clients["pubchem"].get_compound_data(
            compound.smiles, compound.name
        )
        
        # Get Swiss data
        web_data["swiss"] = self.web_clients["swiss"].get_compound_data(
            compound.smiles, compound.name
        )
        
        # Get community data
        web_data["community"] = self.web_clients["community"].get_compound_data(
            compound.smiles, compound.name
        )
        
        # Get social data
        web_data["social"] = self.web_clients["social"].get_compound_data(
            compound.smiles, compound.name
        )
        
        # Get web search data with amino acid transporter focus
        web_data["web_search"] = self._get_transporter_web_data(compound)
        
        return web_data

    def _get_transporter_web_data(self, compound: CompoundData) -> Dict[str, Any]:
        """Get web data focused on amino acid transporters."""
        web_data = {}
        
        # Search for each transporter type
        for transporter, keywords in self.AMINO_ACID_KEYWORDS.items():
            # Build search query
            query = f"{compound.name} {' OR '.join(keywords)}"
            
            # Get web search results
            results = self.web_clients["web_search"].get_compound_data(
                compound.smiles,
                compound.name,
                additional_query=query
            )
            
            web_data[transporter] = results
            
        return web_data

    def _process_web_data(
        self, web_data: Dict[str, Any], compound: CompoundData
    ) -> Dict[str, Any]:
        """Process web data using LLM."""
        processed_data = {}
        
        # Extract BBB-related information
        processed_data["bbb_data"] = self.llm_processor.extract_bbb_data(
            web_data, compound
        )
        
        # Extract amino acid transporter information
        processed_data["transporter_data"] = (
            self.llm_processor.extract_transporter_data(web_data, compound)
        )
        
        # Extract mechanism information
        processed_data["mechanism_data"] = self.llm_processor.extract_mechanism_data(
            web_data, compound
        )
        
        # Extract pharmacology information
        processed_data["pharmacology_data"] = (
            self.llm_processor.extract_pharmacology_data(web_data, compound)
        )
        
        # Extract safety information
        processed_data["safety_data"] = self.llm_processor.extract_safety_data(
            web_data, compound
        )
        
        # Extract interaction information
        processed_data["interaction_data"] = (
            self.llm_processor.extract_interaction_data(web_data, compound)
        )
        
        return processed_data

    def _update_prediction_with_web_data(
        self,
        base_result: PredictionResult,
        web_data: Dict[str, Any],
        compound: CompoundData,
    ) -> PredictionResult:
        """Update prediction result with web data."""
        # Add web data to supporting data
        base_result.supporting_data.update({
            "web_data": {
                "bbb_data": web_data["bbb_data"],
                "transporter_data": web_data["transporter_data"],
                "mechanism_data": web_data["mechanism_data"],
                "pharmacology_data": web_data["pharmacology_data"],
                "safety_data": web_data["safety_data"],
                "interaction_data": web_data["interaction_data"],
            }
        })
        
        # Adjust confidence based on web data support
        confidence_boost = 1.0
        
        # Boost if BBB data supports prediction
        if web_data["bbb_data"].get("supports_prediction", False):
            confidence_boost += 0.2
            
        # Boost if transporter data supports prediction
        if web_data["transporter_data"].get("supports_prediction", False):
            confidence_boost += 0.2
            
        # Apply confidence boost
        base_result.confidence = min(
            1.0,
            base_result.confidence * confidence_boost
        )
        
        return base_result

    def _update_web_prediction_history(
        self, web_data: Dict[str, Any], compound: CompoundData
    ) -> None:
        """Update prediction history with web data."""
        # Find the latest prediction for this compound
        mask = self.prediction_history['compound_name'] == compound.name
        if not mask.any():
            return
        
        latest_idx = self.prediction_history[mask].index[-1]
        
        # Update standard web data columns
        self.prediction_history.loc[latest_idx, 'chembl_data'] = str(
            web_data.get("chembl", {})
        )
        self.prediction_history.loc[latest_idx, 'pubchem_data'] = str(
            web_data.get("pubchem", {})
        )
        self.prediction_history.loc[latest_idx, 'swiss_data'] = str(
            web_data.get("swiss", {})
        )
        self.prediction_history.loc[latest_idx, 'community_data'] = str(
            web_data.get("community", {})
        )
        self.prediction_history.loc[latest_idx, 'social_data'] = str(
            web_data.get("social", {})
        )
        self.prediction_history.loc[latest_idx, 'web_mentions'] = str(
            web_data.get("web_search", {})
        )
        
        # Update amino acid transporter columns
        transporter_data = web_data.get("transporter_data", {})
        self.prediction_history.loc[latest_idx, 'lat1_literature'] = str(
            transporter_data.get("lat1_literature", [])
        )
        self.prediction_history.loc[latest_idx, 'lat1_evidence'] = str(
            transporter_data.get("lat1_evidence", {})
        )
        self.prediction_history.loc[latest_idx, 'lat1_confidence_web'] = float(
            transporter_data.get("lat1_confidence", 0.0)
        )
        self.prediction_history.loc[latest_idx, 'amino_acid_transporters'] = str(
            transporter_data.get("transporters", [])
        )
        self.prediction_history.loc[latest_idx, 'transporter_interactions'] = str(
            transporter_data.get("interactions", {})
        )
        self.prediction_history.loc[latest_idx, 'substrate_evidence'] = str(
            transporter_data.get("substrate_evidence", {})
        )
        
        # Update other data columns
        self.prediction_history.loc[latest_idx, 'mechanism_data'] = str(
            web_data.get("mechanism_data", {})
        )
        self.prediction_history.loc[latest_idx, 'pharmacology_data'] = str(
            web_data.get("pharmacology_data", {})
        )
        self.prediction_history.loc[latest_idx, 'safety_data'] = str(
            web_data.get("safety_data", {})
        )
        self.prediction_history.loc[latest_idx, 'interaction_data'] = str(
            web_data.get("interaction_data", {})
        )

    def export_predictions(
        self,
        output_path: str,
        columns: Optional[List[str]] = None,
        include_supporting_data: bool = False,
        include_web_data: bool = False,
    ) -> None:
        """Export prediction history to TSV file with web data options."""
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
                        on=['compound_name', 'timestamp'],
                        how='left',
                        suffixes=('', f'_{predictor}'),
                    )
        
        # Add web data if requested
        if include_web_data:
            web_columns = [
                'chembl_data',
                'pubchem_data',
                'swiss_data',
                'community_data',
                'social_data',
                'web_mentions',
                'literature_data',
                'patent_data',
                'regulatory_data',
                'clinical_data',
                'mechanism_data',
                'pharmacology_data',
                'safety_data',
                'interaction_data',
                'lat1_literature',
                'lat1_evidence', 
                'lat1_confidence_web',
                'amino_acid_transporters',
                'transporter_interactions',
                'substrate_evidence',
            ]
            export_df = export_df.join(
                self.prediction_history[web_columns],
                how='left'
            )
        
        # Export to TSV
        export_df.to_csv(
            output_path,
            sep='\t',
            index=False,
        )
        self.logger.info(f"Exported predictions to {output_path}")
