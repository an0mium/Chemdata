"""Data enrichment for BBB permeability prediction.

This module provides functionality to:
1. Enrich compound data from multiple sources
2. Integrate web data with ML predictions
3. Combine data from different predictors
4. Validate and merge enriched data
5. Generate comprehensive compound profiles
"""

import logging
from typing import Dict, List, Optional, Any, Set
import pandas as pd
from dataclasses import dataclass
from datetime import datetime

from ....models.core import CompoundData
from .web_validation import WebDataValidator
from .data_standardization import DataStandardizer, StandardizedData
from web_enrichment.data_sources.chembl import ChemblClient
from web_enrichment.data_sources.pubchem import PubchemClient
from web_enrichment.data_sources.swiss import SwissClient
from web_enrichment.data_sources.community import CommunityDataClient
from web_enrichment.data_sources.social import SocialDataClient
from web_enrichment.data_sources.web_search import WebSearchClient
from web_enrichment.llm_utils import LLMProcessor


@dataclass
class EnrichedData:
    """Container for enriched compound data."""
    
    # Basic data
    compound: CompoundData
    standardized: StandardizedData
    
    # Web data
    web_data: Dict[str, Any]
    web_validation: Dict[str, Any]
    
    # ML predictions
    bbb_predictions: Dict[str, Any]
    activity_predictions: Dict[str, Any]
    toxicity_predictions: Dict[str, Any]
    abuse_predictions: Dict[str, Any]
    
    # Enriched properties
    properties: Dict[str, Any]
    activities: Dict[str, List[Dict[str, Any]]]
    targets: Dict[str, List[Dict[str, Any]]]
    mechanisms: Dict[str, List[Dict[str, Any]]]
    
    # Safety data
    safety: Dict[str, Any]
    warnings: List[str]
    contraindications: List[str]
    
    # References
    sources: Set[str]
    references: Dict[str, List[str]]
    timestamps: Dict[str, datetime]


class DataEnricher:
    """Enricher for compound data from multiple sources."""

    def __init__(
        self,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        web_clients: Optional[Dict[str, Any]] = None,
        predictors: Optional[Dict[str, Any]] = None,
    ):
        """Initialize data enricher."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.cache_dir = cache_dir

        # Initialize components
        self.standardizer = DataStandardizer(cache_dir=cache_dir)
        self.validator = WebDataValidator(cache_dir=cache_dir)
        self.llm_processor = LLMProcessor(cache_dir=cache_dir)

        # Initialize web clients
        self.web_clients = web_clients or {
            "chembl": ChemblClient(cache_dir=cache_dir),
            "pubchem": PubchemClient(cache_dir=cache_dir),
            "swiss": SwissClient(cache_dir=cache_dir),
            "community": CommunityDataClient(cache_dir=cache_dir),
            "social": SocialDataClient(cache_dir=cache_dir),
            "web_search": WebSearchClient(cache_dir=cache_dir),
        }

        # Initialize predictors
        self.predictors = predictors or {}

    def enrich_compound_data(
        self, compound: CompoundData
    ) -> EnrichedData:
        """Enrich compound data from multiple sources."""
        self.logger.debug(f"Enriching data for {compound.name}")
        
        try:
            # Get web data
            web_data = self._get_web_data(compound)
            
            # Validate web data
            web_validation = self.validator.validate_web_data(
                web_data, compound.name
            )
            
            # Standardize data
            standardized = self.standardizer.standardize_compound_data(
                web_data,
                compound.compound_id,
                compound.name,
                compound.smiles,
            )
            
            # Get ML predictions
            predictions = self._get_predictions(compound)
            
            # Process data with LLM
            llm_data = self._process_with_llm(web_data, predictions)
            
            # Combine all data
            enriched = self._combine_data(
                compound=compound,
                standardized=standardized,
                web_data=web_data,
                web_validation=web_validation.__dict__,
                predictions=predictions,
                llm_data=llm_data,
            )
            
            return enriched

        except Exception as e:
            self.logger.error(
                f"Error enriching data for {compound.name}: {str(e)}",
                exc_info=True
            )
            raise

    def _get_web_data(
        self, compound: CompoundData
    ) -> Dict[str, Dict[str, Any]]:
        """Get data from web sources."""
        web_data = {}
        
        for source, client in self.web_clients.items():
            try:
                data = client.get_compound_data(
                    compound.smiles,
                    compound.name,
                )
                if data:
                    web_data[source] = data
            except Exception as e:
                self.logger.error(
                    f"Error getting {source} data for {compound.name}: {str(e)}"
                )
        
        return web_data

    def _get_predictions(
        self, compound: CompoundData
    ) -> Dict[str, Any]:
        """Get predictions from ML models."""
        predictions = {
            "bbb": {},
            "activity": {},
            "toxicity": {},
            "abuse": {},
        }
        
        # Get predictions from each predictor
        for predictor_type in ["bbb", "activity", "toxicity", "abuse"]:
            predictions[predictor_type] = self._get_predictor_result(
                predictor_type, compound
            )
        
        return predictions

    def _get_predictor_result(
        self, predictor_type: str, compound: CompoundData
    ) -> Dict[str, Any]:
        """Get prediction result from a specific predictor."""
        if predictor_type not in self.predictors:
            return {}
            
        try:
            result = self.predictors[predictor_type].predict(compound)
            
            # Format result based on predictor type
            if predictor_type == "activity":
                return {
                    "predictions": result.value,
                    "confidence": result.confidence,
                    "supporting_data": result.supporting_data,
                }
            elif predictor_type == "abuse":
                return {
                    "potential": result.value,
                    "confidence": result.confidence,
                    "supporting_data": result.supporting_data,
                }
            else:  # bbb and toxicity use same format
                return {
                    "class": result.value,
                    "confidence": result.confidence,
                    "supporting_data": result.supporting_data,
                }
                
        except Exception as e:
            self.logger.error(
                f"Error getting {predictor_type} predictions for "
                f"{compound.name}: {str(e)}"
            )
            return {}

    def _process_with_llm(
        self,
        web_data: Dict[str, Dict[str, Any]],
        predictions: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Process data using LLM."""
        llm_data = {}
        
        try:
            # Extract mechanism descriptions
            llm_data["mechanisms"] = self.llm_processor.extract_mechanism_data(
                web_data
            )
            
            # Extract safety information
            llm_data["safety"] = self.llm_processor.extract_safety_data(
                web_data
            )
            
            # Extract pharmacology information
            llm_data["pharmacology"] = (
                self.llm_processor.extract_pharmacology_data(web_data)
            )
            
            # Extract interaction information
            llm_data["interactions"] = (
                self.llm_processor.extract_interaction_data(web_data)
            )
            
            # Analyze predictions
            llm_data["prediction_analysis"] = (
                self.llm_processor.analyze_predictions(predictions)
            )
            
        except Exception as e:
            self.logger.error(
                f"Error in LLM processing: {str(e)}",
                exc_info=True
            )
        
        return llm_data

    def _combine_data(
        self,
        compound: CompoundData,
        standardized: StandardizedData,
        web_data: Dict[str, Dict[str, Any]],
        web_validation: Dict[str, Any],
        predictions: Dict[str, Any],
        llm_data: Dict[str, Any],
    ) -> EnrichedData:
        """Combine all data sources."""
        # Combine properties
        properties = self._combine_properties(
            standardized, web_data, predictions
        )
        
        # Combine activities
        activities = self._combine_activities(
            standardized, web_data, predictions
        )
        
        # Combine targets
        targets = self._combine_targets(
            standardized, web_data, predictions
        )
        
        # Combine mechanisms
        mechanisms = self._combine_mechanisms(
            standardized, web_data, predictions, llm_data
        )
        
        # Combine safety data
        safety = self._combine_safety_data(
            standardized, web_data, predictions, llm_data
        )
        
        # Get all warnings
        warnings = self._get_all_warnings(
            standardized, web_data, predictions, llm_data
        )
        
        # Get contraindications
        contraindications = self._get_contraindications(
            standardized, web_data, predictions, llm_data
        )
        
        # Combine sources and references
        sources = standardized.sources | set(web_data.keys())
        references = {**standardized.references}
        for source, data in web_data.items():
            if "references" in data:
                references[source] = data["references"]
        
        # Combine timestamps
        timestamps = {**standardized.timestamps}
        for source, data in web_data.items():
            if "timestamp" in data:
                try:
                    timestamps[source] = pd.to_datetime(
                        data["timestamp"]
                    ).to_pydatetime()
                except (ValueError, TypeError):
                    pass
        
        return EnrichedData(
            compound=compound,
            standardized=standardized,
            web_data=web_data,
            web_validation=web_validation,
            bbb_predictions=predictions["bbb"],
            activity_predictions=predictions["activity"],
            toxicity_predictions=predictions["toxicity"],
            abuse_predictions=predictions["abuse"],
            properties=properties,
            activities=activities,
            targets=targets,
            mechanisms=mechanisms,
            safety=safety,
            warnings=warnings,
            contraindications=contraindications,
            sources=sources,
            references=references,
            timestamps=timestamps,
        )

    def _combine_properties(
        self,
        standardized: StandardizedData,
        web_data: Dict[str, Dict[str, Any]],
        predictions: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Combine property data from all sources."""
        properties = {}
        
        # Add standardized properties
        properties["molecular_weight"] = standardized.molecular_weight
        properties["logp"] = standardized.logp
        properties["psa"] = standardized.psa
        properties["hba"] = standardized.hba
        properties["hbd"] = standardized.hbd
        properties["rotatable_bonds"] = standardized.rotatable_bonds
        
        # Add predicted properties
        if "properties" in predictions.get("bbb", {}).get("supporting_data", {}):
            properties.update(
                predictions["bbb"]["supporting_data"]["properties"]
            )
        
        return properties

    def _combine_activities(
        self,
        standardized: StandardizedData,
        web_data: Dict[str, Dict[str, Any]],
        predictions: Dict[str, Any],
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Combine activity data from all sources."""
        activities = {**standardized.activities}
        
        # Add predicted activities
        if "predictions" in predictions.get("activity", {}):
            for target, pred in predictions["activity"]["predictions"].items():
                if target not in activities:
                    activities[target] = []
                activities[target].append({
                    "value": pred["value"],
                    "confidence": pred["confidence"],
                    "source": "ml_prediction",
                })
        
        return activities

    def _combine_targets(
        self,
        standardized: StandardizedData,
        web_data: Dict[str, Dict[str, Any]],
        predictions: Dict[str, Any],
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Combine target data from all sources."""
        targets = {**standardized.targets}
        
        # Add predicted targets
        if "targets" in predictions.get("activity", {}).get("supporting_data", {}):
            pred_targets = predictions["activity"]["supporting_data"]["targets"]
            for target, data in pred_targets.items():
                if target not in targets:
                    targets[target] = []
                targets[target].append({
                    "type": data.get("type", "unknown"),
                    "confidence": data.get("confidence", 1.0),
                    "source": "ml_prediction",
                })
        
        return targets

    def _combine_mechanisms(
        self,
        standardized: StandardizedData,
        web_data: Dict[str, Dict[str, Any]],
        predictions: Dict[str, Any],
        llm_data: Dict[str, Any],
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Combine mechanism data from all sources."""
        mechanisms = {**standardized.mechanisms}
        
        # Add LLM-extracted mechanisms
        if "mechanisms" in llm_data:
            for mech_type, mech_data in llm_data["mechanisms"].items():
                if mech_type not in mechanisms:
                    mechanisms[mech_type] = []
                mechanisms[mech_type].append({
                    "description": mech_data.get("description", ""),
                    "confidence": mech_data.get("confidence", 1.0),
                    "source": "llm_extraction",
                })
        
        return mechanisms

    def _combine_safety_data(
        self,
        standardized: StandardizedData,
        web_data: Dict[str, Dict[str, Any]],
        predictions: Dict[str, Any],
        llm_data: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Combine safety data from all sources."""
        safety = {
            "toxicity": {**standardized.toxicity},
            "side_effects": {**standardized.side_effects},
        }
        
        # Add predicted toxicity
        if "class" in predictions.get("toxicity", {}):
            safety["toxicity"]["predicted_class"] = {
                "value": predictions["toxicity"]["class"],
                "confidence": predictions["toxicity"]["confidence"],
                "source": "ml_prediction",
            }
        
        # Add LLM-extracted safety data
        if "safety" in llm_data:
            if "toxicity" in llm_data["safety"]:
                safety["toxicity"].update(llm_data["safety"]["toxicity"])
            if "side_effects" in llm_data["safety"]:
                safety["side_effects"].update(llm_data["safety"]["side_effects"])
        
        return safety

    def _get_all_warnings(
        self,
        standardized: StandardizedData,
        web_data: Dict[str, Dict[str, Any]],
        predictions: Dict[str, Any],
        llm_data: Dict[str, Any],
    ) -> List[str]:
        """Get all warnings from all sources."""
        warnings = set(standardized.warnings)
        
        # Add predicted warnings
        if "warnings" in predictions.get("toxicity", {}).get("supporting_data", {}):
            warnings.update(
                predictions["toxicity"]["supporting_data"]["warnings"]
            )
        
        # Add LLM-extracted warnings
        if "safety" in llm_data and "warnings" in llm_data["safety"]:
            warnings.update(llm_data["safety"]["warnings"])
        
        return list(warnings)

    def _get_contraindications(
        self,
        standardized: StandardizedData,
        web_data: Dict[str, Dict[str, Any]],
        predictions: Dict[str, Any],
        llm_data: Dict[str, Any],
    ) -> List[str]:
        """Get contraindications from all sources."""
        contraindications = set()
        
        # Add from safety data
        if "contraindications" in standardized.side_effects:
            contraindications.update(
                standardized.side_effects["contraindications"]
            )
        
        # Add from predictions
        if "contraindications" in predictions.get("toxicity", {}).get("supporting_data", {}):
            contraindications.update(
                predictions["toxicity"]["supporting_data"]["contraindications"]
            )
        
        # Add from LLM data
        if "safety" in llm_data and "contraindications" in llm_data["safety"]:
            contraindications.update(llm_data["safety"]["contraindications"])
        
        return list(contraindications)
