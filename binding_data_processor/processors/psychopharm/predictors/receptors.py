"""Receptor binding profile predictor.

This module provides predictors for analyzing receptor binding profiles:
1. Target-specific binding affinity prediction
2. Binding type classification (agonist/antagonist)
3. Receptor selectivity analysis
4. Cross-target interaction prediction
5. Web data enrichment and validation
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor

from .base import PredictorBase, PredictorConfig, WebEnrichedPredictor
from ....models.compound import PsychoactiveCompound
from ....models.predictions import PredictionResult
from ...structure.ml.features import EnhancedFeatureExtractor
from ...structure.ml.fingerprints import FingerprintGenerator
from ....web_enrichment.manager import WebEnrichmentManager


class ReceptorBindingConfig(PredictorConfig):
    """Configuration for receptor binding prediction."""
    
    # Key receptor families to predict
    RECEPTOR_FAMILIES = {
        "serotonin": {
            "5-HT1A", "5-HT1B", "5-HT1D", "5-HT2A", "5-HT2B", "5-HT2C",
            "5-HT3", "5-HT4", "5-HT5A", "5-HT6", "5-HT7"
        },
        "dopamine": {
            "D1", "D2", "D3", "D4", "D5"
        },
        "norepinephrine": {
            "α1A", "α1B", "α1D", "α2A", "α2B", "α2C",
            "β1", "β2", "β3"
        },
        "glutamate": {
            "NMDA", "AMPA", "Kainate", "mGluR1", "mGluR2", "mGluR3",
            "mGluR4", "mGluR5", "mGluR6", "mGluR7", "mGluR8"
        },
        "gaba": {
            "GABA-A", "GABA-B", "GABA-C"
        },
        "opioid": {
            "μ", "κ", "δ", "NOP"
        },
        "cannabinoid": {
            "CB1", "CB2"
        },
        "histamine": {
            "H1", "H2", "H3", "H4"
        },
        "sigma": {
            "σ1", "σ2"
        },
    }

    # Activity types to predict
    ACTIVITY_TYPES = [
        "full_agonist",
        "partial_agonist", 
        "antagonist",
        "inverse_agonist",
        "allosteric_modulator",
        "unknown"
    ]

    def __init__(
        self,
        model_dir: Optional[str] = None,
        receptor_families: Optional[Dict[str, Set[str]]] = None,
        activity_types: Optional[List[str]] = None,
        **kwargs
    ):
        """Initialize config.
        
        Args:
            model_dir: Optional directory containing trained models
            receptor_families: Optional custom receptor families
            activity_types: Optional custom activity types
            **kwargs: Additional config parameters
        """
        super().__init__(**kwargs)
        self.model_dir = Path(model_dir) if model_dir else None
        self.receptor_families = receptor_families or self.RECEPTOR_FAMILIES
        self.activity_types = activity_types or self.ACTIVITY_TYPES
        
        # Flatten receptor list
        self.target_receptors = []
        for family in self.receptor_families.values():
            self.target_receptors.extend(family)


class ReceptorBindingPredictor(WebEnrichedPredictor):
    """Predicts receptor binding profiles for compounds."""

    def __init__(
        self,
        config: Optional[ReceptorBindingConfig] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize predictor.
        
        Args:
            config: Optional predictor configuration
            log_level: Logging level
        """
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        
        # Add file handler if model_dir provided
        if config and config.model_dir:
            log_path = config.model_dir / "receptor_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

        super().__init__(config or ReceptorBindingConfig())
        
        # Feature extractors
        self.fingerprint_gen = FingerprintGenerator()
        self.feature_gen = EnhancedFeatureExtractor()
        
        # Web enrichment
        self.web_manager = WebEnrichmentManager()

        # Initialize prediction history
        self.prediction_history = pd.DataFrame(
            columns=[
                'compound_name',
                'receptor',
                'affinity_value',
                'affinity_confidence',
                'activity_type',
                'activity_confidence',
                'timestamp',
            ]
        )

        self.logger.info("ReceptorBindingPredictor initialized successfully")

    def _initialize(self) -> None:
        """Initialize predictor components."""
        self.logger.debug("Initializing predictor components")
        
        # Initialize affinity prediction models
        self.affinity_models = {}
        for receptor in self.config.target_receptors:
            self.affinity_models[receptor] = self._load_affinity_model(receptor)

        # Initialize binding type classification models
        self.type_models = {}
        for receptor in self.config.target_receptors:
            self.type_models[receptor] = self._load_type_model(receptor)

    def _load_affinity_model(self, receptor: str) -> RandomForestRegressor:
        """Load affinity prediction model for receptor."""
        if self.config.model_dir:
            # Load trained model from disk
            model_path = self.config.model_dir / f"{receptor}_affinity.pkl"
            try:
                self.logger.debug(f"Loading affinity model for {receptor}")
                return np.load(model_path, allow_pickle=True)
            except Exception as e:
                self.logger.warning(
                    f"Failed to load affinity model for {receptor}: {str(e)}"
                )
                
        # Create new model if loading fails
        self.logger.debug(f"Creating new affinity model for {receptor}")
        return RandomForestRegressor(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1,
            verbose=0
        )

    def _load_type_model(self, receptor: str) -> GradientBoostingRegressor:
        """Load binding type classification model for receptor."""
        if self.config.model_dir:
            # Load trained model from disk
            model_path = self.config.model_dir / f"{receptor}_type.pkl"
            try:
                self.logger.debug(f"Loading type model for {receptor}")
                return np.load(model_path, allow_pickle=True)
            except Exception as e:
                self.logger.warning(
                    f"Failed to load type model for {receptor}: {str(e)}"
                )
                
        # Create new model if loading fails
        self.logger.debug(f"Creating new type model for {receptor}")
        return GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=0
        )

    def _extract_features(self, compound: PsychoactiveCompound) -> np.ndarray:
        """Extract features for binding prediction."""
        try:
            # Generate fingerprints
            fp = self.fingerprint_gen.generate(compound.smiles)
            
            # Generate additional features
            features = self.feature_gen.generate(compound.smiles)
            
            # Combine features
            return np.concatenate([fp, features])
            
        except Exception as e:
            self.logger.error(
                f"Feature extraction error for {compound.name}: {str(e)}"
            )
            return np.array([])

    def _predict_raw(self, features: np.ndarray) -> Tuple[Dict, float]:
        """Generate raw binding predictions."""
        predictions = {}
        confidences = []
        
        # Predict for each receptor
        for receptor in self.config.target_receptors:
            try:
                # Predict binding affinity
                affinity_model = self.affinity_models[receptor]
                affinity = float(affinity_model.predict([features])[0])
                affinity_conf = 1.0 - np.std([
                    tree.predict([features])[0]
                    for tree in affinity_model.estimators_
                ])
                
                # Predict binding type
                type_model = self.type_models[receptor]
                type_probs = type_model.predict_proba([features])[0]
                binding_type = self.config.activity_types[np.argmax(type_probs)]
                type_conf = float(np.max(type_probs))
                
                # Store predictions
                predictions[receptor] = {
                    "affinity": affinity,
                    "type": binding_type,
                }
                confidences.extend([affinity_conf, type_conf])
                
            except Exception as e:
                self.logger.error(
                    f"Prediction error for {receptor}: {str(e)}"
                )
                continue
            
        # Calculate overall confidence
        confidence = float(np.mean(confidences)) if confidences else 0.0
        
        return predictions, confidence

    def _process_prediction(
        self,
        prediction: Dict,
        confidence: float,
        compound: PsychoactiveCompound
    ) -> PredictionResult:
        """Process raw predictions into result."""
        try:
            # Format predictions
            results = {}
            for receptor, data in prediction.items():
                results[receptor] = {
                    "affinity": float(data["affinity"]),
                    "affinity_unit": "Ki (nM)",
                    "type": str(data["type"]),
                    "confidence": float(confidence),
                }
                
                # Update prediction history
                self.prediction_history = pd.concat([
                    self.prediction_history,
                    pd.DataFrame([{
                        'compound_name': compound.name,
                        'receptor': receptor,
                        'affinity_value': data["affinity"],
                        'affinity_confidence': confidence,
                        'activity_type': data["type"],
                        'activity_confidence': confidence,
                        'timestamp': pd.Timestamp.now(),
                    }])
                ], ignore_index=True)
                
            # Add metadata
            metadata = {
                "target_receptors": self.config.target_receptors,
                "activity_types": self.config.activity_types,
                "model_version": self.config.model_version,
            }
                
            return PredictionResult(
                value=results,
                confidence=confidence,
                metadata=metadata
            )
            
        except Exception as e:
            self.logger.error(
                f"Error processing predictions for {compound.name}: {str(e)}"
            )
            return PredictionResult(
                value={},
                confidence=0.0,
                metadata={"error": str(e)}
            )

    def _get_web_data(self, compound: PsychoactiveCompound) -> Optional[Dict]:
        """Get binding data from web sources."""
        if not self.config.use_web_data:
            return None
            
        binding_data = {}
        
        try:
            # Get ChEMBL data
            chembl_data = self.web_manager.get_chembl_data(
                compound.smiles,
                targets=self.config.target_receptors
            )
            if chembl_data:
                binding_data["chembl"] = chembl_data
                
            # Get PubChem data
            pubchem_data = self.web_manager.get_pubchem_data(
                compound.smiles,
                targets=self.config.target_receptors
            )
            if pubchem_data:
                binding_data["pubchem"] = pubchem_data
                
            # Get patent data
            patent_data = self.web_manager.get_patent_data(
                compound.smiles,
                targets=self.config.target_receptors
            )
            if patent_data:
                binding_data["patents"] = patent_data
                
            return binding_data if binding_data else None
            
        except Exception as e:
            self.logger.error(
                f"Error getting web data for {compound.name}: {str(e)}"
            )
            return None

    def _combine_predictions(
        self,
        ml_result: PredictionResult,
        web_data: Dict,
        compound: PsychoactiveCompound
    ) -> PredictionResult:
        """Combine ML and web-based predictions."""
        try:
            combined_results = {}
            
            # Process each receptor
            for receptor in self.config.target_receptors:
                ml_data = ml_result.value.get(receptor, {})
                
                # Collect web data
                web_values = []
                web_confs = []
                
                for source, data in web_data.items():
                    if receptor in data:
                        web_values.append(data[receptor].get("affinity"))
                        web_confs.append(data[receptor].get("confidence", 0.5))
                        
                # Combine predictions if web data exists
                if web_values:
                    # Weight ML and web predictions
                    ml_weight = ml_result.confidence
                    web_weight = np.mean(web_confs)
                    
                    # Calculate weighted average
                    combined_affinity = (
                        ml_data["affinity"] * ml_weight +
                        np.mean(web_values) * web_weight
                    ) / (ml_weight + web_weight)
                    
                    combined_conf = (ml_weight + web_weight) / 2
                    
                    combined_results[receptor] = {
                        "affinity": float(combined_affinity),
                        "affinity_unit": "Ki (nM)",
                        "type": ml_data["type"],  # Keep ML binding type
                        "confidence": float(combined_conf),
                    }
                else:
                    # Use ML prediction if no web data
                    combined_results[receptor] = ml_data
                    
            return PredictionResult(
                value=combined_results,
                confidence=float(np.mean([
                    r["confidence"] for r in combined_results.values()
                ])),
                metadata=ml_result.metadata
            )
            
        except Exception as e:
            self.logger.error(
                f"Error combining predictions for {compound.name}: {str(e)}"
            )
            return ml_result

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save models to disk with versioning."""
        save_dir = Path(save_dir) if save_dir else self.config.model_dir
        if not save_dir:
            self.logger.warning("No save directory specified")
            return

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving models to {save_dir}")

        # Save models with versioning
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        version_dir = save_dir / f"version_{timestamp}"
        version_dir.mkdir(exist_ok=True)

        try:
            # Save affinity models
            for receptor, model in self.affinity_models.items():
                # Save current version
                model_path = save_dir / f"{receptor}_affinity.pkl"
                np.save(model_path, model)
                
                # Save versioned copy
                version_path = version_dir / f"{receptor}_affinity.pkl"
                np.save(version_path, model)
                
                self.logger.debug(f"Saved affinity model for {receptor}")

            # Save type models
            for receptor, model in self.type_models.items():
                # Save current version
                model_path = save_dir / f"{receptor}_type.pkl"
                np.save(model_path, model)
                
                # Save versioned copy
                version_path = version_dir / f"{receptor}_type.pkl"
                np.save(version_path, model)
                
                self.logger.debug(f"Saved type model for {receptor}")

            # Save prediction history
            history_path = save_dir / "prediction_history.csv"
            self.prediction_history.to_csv(history_path, index=False)
            self.logger.debug("Saved prediction history")

        except Exception as e:
            self.logger.error(f"Error saving models: {str(e)}")

    def get_prediction_statistics(self) -> pd.DataFrame:
        """Get statistics about predictions made so far."""
        stats = pd.DataFrame()
        
        # Affinity value distribution
        stats['affinity_mean'] = (
            self.prediction_history.groupby('receptor')['affinity_value'].mean()
        )
        stats['affinity_std'] = (
            self.prediction_history.groupby('receptor')['affinity_value'].std()
        )
        
        # Activity type distribution
        stats['activity_dist'] = (
            self.prediction_history.groupby('receptor')['activity_type']
            .value_counts(normalize=True)
        )
        
        # Average confidence scores
        stats['affinity_confidence'] = (
            self.prediction_history.groupby('receptor')['affinity_confidence'].mean()
        )
        stats['activity_confidence'] = (
            self.prediction_history.groupby('receptor')['activity_confidence'].mean()
        )
        
        return stats
