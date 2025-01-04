"""Enhanced BBB permeability predictor.

This module provides the BBBPredictorEnhanced class that extends BBBPredictor with:
1. Transporter predictions (P-gp, BCRP, etc.)
2. Receptor-mediated transport predictions
3. Duration and pharmacokinetic predictions
4. Integration with other predictors (abuse, toxicity)
5. Enhanced amino acid transport detection
6. Web data enrichment
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Any
import numpy as np  # Required for array operations in training methods
import pandas as pd

from ....models.core import CompoundData
from ....models.psychopharm import BBBPermeability
from ..base import PredictionResult
from .predictors import BBBPredictor
from .abuse import AbusePotentialPredictor
from .toxicity import ToxicityPredictor
from .receptors import ReceptorPredictor
from .psychoactive import PsychoactivePredictor
from .nootropic import NootropicPredictor


class BBBPredictorEnhanced(BBBPredictor):
    """Enhanced BBB permeability predictor with additional capabilities."""

    # BBB-related transporters
    TRANSPORTERS = {
        "efflux": {
            "p_glycoprotein", "bcrp", "mrp1", "mrp2", "mrp4",
        },
        "uptake": {
            "lat1", "mct1", "glut1", "oatp1a2", "oat3",
        },
        "specialized": {
            "organic_cation_transporter", "amino_acid_transporter",
            "peptide_transporter", "fatty_acid_transporter",
        },
    }

    # Receptor-mediated transport mechanisms
    RECEPTOR_TRANSPORTERS = {
        "transferrin_receptor": {"iron_transport", "antibody_transport"},
        "insulin_receptor": {"insulin_transport", "protein_transport"},
        "leptin_receptor": {"peptide_transport"},
        "ldl_receptor": {"lipid_transport"},
    }

    # Duration-based features
    DURATION_FEATURES = {
        "rapid": {"onset": "<15min", "half_life": "<2h"},
        "intermediate": {"onset": "15-60min", "half_life": "2-12h"},
        "slow": {"onset": ">60min", "half_life": ">12h"},
    }

    # Additional amino acid patterns
    AMINO_ACID_PATTERNS = {
        **BBBPredictor.LAT1_PATTERNS,
        "carboxyl_group": r"C\(=O\)O",  # Carboxyl group
        "amine_group": r"CN|NC",  # Primary/secondary amine
        "peptide_bond": r"NC\(=O\)",  # Peptide bond
    }

    # Additional model filenames
    MODEL_FILENAMES = {
        **BBBPredictor.MODEL_FILENAMES,
        "pgp_classifier": "pgp_classifier.pkl",
        "transporter_classifier": "transporter_classifier.pkl",
        "ensemble_classifier": "ensemble_classifier.pkl",
        "receptor_classifier": "receptor_classifier.pkl",
        "duration_classifier": "duration_classifier.pkl",
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        transporters: Optional[Dict[str, Set[str]]] = None,
        receptor_transporters: Optional[Dict[str, Set[str]]] = None,
        predictors: Optional[Dict[str, Any]] = None,
    ):
        """Initialize enhanced BBB predictor.
        
        Args:
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            log_level: Logging level
            transporters: Optional custom transporter definitions
            receptor_transporters: Optional custom receptor transporter definitions
            predictors: Optional custom predictor instances
        """
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
            transporters=transporters,
            receptor_transporters=receptor_transporters,
        )

        # Initialize other predictors
        self.predictors = predictors or {
            "abuse": AbusePotentialPredictor(
                model_dir=str(Path(model_dir) / "abuse") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
            "toxicity": ToxicityPredictor(
                model_dir=str(Path(model_dir) / "toxicity") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
            "receptors": ReceptorPredictor(
                model_dir=str(Path(model_dir) / "receptors") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
            "psychoactive": PsychoactivePredictor(
                model_dir=str(Path(model_dir) / "psychoactive") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
            "nootropic": NootropicPredictor(
                model_dir=str(Path(model_dir) / "nootropic") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
        }

        # Configure model weights for ensemble
        self.model_weights = {
            "rf_classifier": 0.25,
            "gb_classifier": 0.15,
            "rf_regressor": 0.15,
            "amino_acid_features": 0.15,  # Add weight for amino acid features
            "ensemble_classifier": 0.10,  # Reduced from 0.25
            "receptor_classifier": 0.1,
            "duration_classifier": 0.1,
        }

        # Initialize enhanced prediction history
        self.prediction_history = pd.DataFrame(
            columns=[
                'compound_name',
                'permeability_class',
                'confidence',
                'permeability_score',
                'p_gp_substrate',
                'transporter_category',
                'transporter',
                'is_substrate',
                'transporter_confidence',
                'lat1_substrate',  # Added for amino acid transport
                'lat1_confidence',  # Added for amino acid transport
                'amino_acid_features',  # Added for amino acid patterns
                'receptor_mediated',
                'receptor_type',
                'duration_class',
                'half_life',
                'cross_tolerance',
                'abuse_potential',
                'toxicity_class',
                'psychoactive_class',
                'nootropic_mechanisms',
                'timestamp',
            ]
        )

        self.logger.info("BBBPredictorEnhanced initialized successfully")

    def _extract_compound_features(self, compound: CompoundData) -> Dict[str, np.ndarray]:
        """Extract compound features including amino acid patterns."""
        # Get base features
        features = super()._extract_compound_features(compound)
        
        # Add amino acid pattern features
        amino_features = []
        for pattern_name, pattern in self.AMINO_ACID_PATTERNS.items():
            match = bool(re.search(pattern, compound.smiles))
            amino_features.append(float(match))
        
        features["amino_acid_features"] = np.array(amino_features)
        
        return features

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate comprehensive BBB permeability predictions."""
        self.logger.debug(f"Generating predictions for {compound.name}")
        try:
            # Extract features and get base predictions
            features = self._extract_compound_features(compound)
            predictions = self._get_all_predictions(features, compound)
            
            # Get predictions from other predictors
            predictions.update(self._get_integrated_predictions(compound))
            
            # Format and return result
            result = self._format_prediction_result(predictions, compound)
            self._update_prediction_history(predictions, compound)
            
            self.logger.info(
                f"Generated predictions for {compound.name}: "
                f"{result.value.value} (confidence: {result.confidence:.3f})"
            )
            
            return result

        except Exception as e:
            self.logger.error(
                "Error predicting BBB properties",
                f"Compound {compound.name}: {str(e)}",
                exc_info=True
            )
            return PredictionResult(
                value=BBBPermeability.UNKNOWN,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _get_integrated_predictions(
        self, compound: CompoundData
    ) -> Dict[str, Any]:
        """Get predictions from integrated predictors."""
        predictions = {}
        
        # Get abuse potential predictions
        abuse_result = self.predictors["abuse"].predict(compound)
        predictions["abuse"] = {
            "potential": abuse_result.value,
            "confidence": abuse_result.confidence,
            "supporting_data": abuse_result.supporting_data,
        }
        
        # Get toxicity predictions
        toxicity_result = self.predictors["toxicity"].predict(compound)
        predictions["toxicity"] = {
            "class": toxicity_result.value,
            "confidence": toxicity_result.confidence,
            "supporting_data": toxicity_result.supporting_data,
        }
        
        # Get receptor predictions
        receptor_result = self.predictors["receptors"].predict(compound)
        predictions["receptors"] = {
            "binding": receptor_result.value,
            "confidence": receptor_result.confidence,
            "supporting_data": receptor_result.supporting_data,
        }
        
        # Get psychoactive predictions
        psychoactive_result = self.predictors["psychoactive"].predict(compound)
        predictions["psychoactive"] = {
            "class": psychoactive_result.value,
            "confidence": psychoactive_result.confidence,
            "supporting_data": psychoactive_result.supporting_data,
        }
        
        # Get nootropic predictions
        nootropic_result = self.predictors["nootropic"].predict(compound)
        predictions["nootropic"] = {
            "mechanisms": nootropic_result.value,
            "confidence": nootropic_result.confidence,
            "supporting_data": nootropic_result.supporting_data,
        }
        
        return predictions

    def _format_prediction_result(
        self, predictions: Dict[str, Any], compound: CompoundData
    ) -> PredictionResult:
        """Format predictions into PredictionResult."""
        # Get base BBB predictions
        base_result = super()._format_prediction_result(predictions, compound)
        
        # Add integrated predictions
        base_result.supporting_data.update({
            "abuse_potential": predictions["abuse"],
            "toxicity": predictions["toxicity"],
            "receptor_binding": predictions["receptors"],
            "psychoactive": predictions["psychoactive"],
            "nootropic": predictions["nootropic"],
        })
        
        return base_result

    def _update_prediction_history(
        self, predictions: Dict[str, Any], compound: CompoundData
    ) -> None:
        """Update prediction history with new predictions."""
        # Get LAT1 substrate prediction
        is_lat1, lat1_conf = self._is_lat1_substrate(compound.smiles)
        
        # Get amino acid pattern matches
        amino_patterns = [
            k for k, v in self.AMINO_ACID_PATTERNS.items() 
            if re.search(v, compound.smiles)
        ]
        
        # Create history entry
        self.prediction_history = pd.concat([
            self.prediction_history,
            pd.DataFrame([{
                'compound_name': compound.name,
                'permeability_class': predictions["permeability"]["class"],
                'confidence': predictions["permeability"]["class_confidence"],
                'permeability_score': predictions["permeability"]["score"],
                'p_gp_substrate': predictions["transporters"]["p_gp"]["is_substrate"],
                'transporter_category': None,
                'transporter': None,
                'is_substrate': None,
                'transporter_confidence': None,
                'lat1_substrate': is_lat1,
                'lat1_confidence': lat1_conf,
                'amino_acid_features': ','.join(amino_patterns),
                'receptor_mediated': bool(predictions["receptors"]),
                'receptor_type': next(iter(predictions["receptors"]), None),
                'duration_class': predictions["duration"]["class"],
                'half_life': predictions["duration"]["half_life"],
                'cross_tolerance': ','.join(predictions["cross_tolerance"]),
                'abuse_potential': predictions["abuse"]["potential"].value,
                'toxicity_class': predictions["toxicity"]["class"].value,
                'psychoactive_class': predictions["psychoactive"]["class"].value,
                'nootropic_mechanisms': ','.join(
                    [m.value for m in predictions["nootropic"]["mechanisms"]]
                ),
                'timestamp': pd.Timestamp.now(),
            }])
        ], ignore_index=True)

    def export_predictions(
        self,
        output_path: str,
        columns: Optional[List[str]] = None,
        include_supporting_data: bool = False,
    ) -> None:
        """Export prediction history to TSV file with enhanced options."""
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
        
        # Export to TSV
        export_df.to_csv(
            output_path,
            sep='\t',
            index=False,
        )
        self.logger.info(f"Exported predictions to {output_path}")

    def retrain(
        self,
        compounds: List[CompoundData],
        labels: List[BBBPermeability],
        scores: Optional[List[float]] = None,
        pgp_labels: Optional[List[bool]] = None,
        transporter_data: Optional[Dict[str, List[bool]]] = None,
        receptor_data: Optional[Dict[str, List[bool]]] = None,
        duration_labels: Optional[List[str]] = None,
        **kwargs,
    ) -> Dict[str, float]:
        """Retrain models with new data.
        
        Args:
            compounds: List of compounds to train on
            labels: BBB permeability class labels
            scores: Optional permeability scores for regression
            pgp_labels: Optional P-gp substrate labels
            transporter_data: Optional transporter substrate labels
            receptor_data: Optional receptor transport labels
            duration_labels: Optional duration class labels
            **kwargs: Additional training parameters
            
        Returns:
            Dictionary of training metrics
        """
        # Retrain base models
        metrics = super().retrain(
            compounds=compounds,
            labels=labels,
            scores=scores,
            **kwargs,
        )
        
        # Extract features once for all models
        X = self._extract_training_features(compounds)
        
        # Train transporter models
        if transporter_data:
            metrics.update(self._train_transporter_models(X, transporter_data))
        
        # Train receptor models
        if receptor_data:
            metrics.update(self._train_receptor_models(X, receptor_data))
        
        # Train duration model
        if duration_labels:
            metrics.update(self._train_duration_model(X, duration_labels))
        
        # Train P-gp model
        if pgp_labels:
            metrics.update(self._train_pgp_model(X, pgp_labels))
        
        # Save updated models
        self.save_models()
        
        self.logger.info("Model retraining completed successfully")
        return metrics

    def _train_transporter_models(
        self, X: np.ndarray, transporter_data: Dict[str, List[bool]]
    ) -> Dict[str, float]:
        """Train transporter models."""
        metrics = {}
        self.logger.debug("Training transporter models")
        
        for transporter, labels in transporter_data.items():
            self.models["transporter_classifier"].fit(X, labels)
            score = self.models["transporter_classifier"].score(X, labels)
            metrics[f"transporter_{transporter}_score"] = score
            self.logger.debug(
                f"Transporter classifier for {transporter} "
                f"training score: {score:.3f}"
            )
        
        return metrics

    def _train_receptor_models(
        self, X: np.ndarray, receptor_data: Dict[str, List[bool]]
    ) -> Dict[str, float]:
        """Train receptor models."""
        metrics = {}
        self.logger.debug("Training receptor models")
        
        for receptor, labels in receptor_data.items():
            self.models["receptor_classifier"].fit(X, labels)
            score = self.models["receptor_classifier"].score(X, labels)
            metrics[f"receptor_{receptor}_score"] = score
            self.logger.debug(
                f"Receptor classifier for {receptor} "
                f"training score: {score:.3f}"
            )
        
        return metrics

    def _train_duration_model(
        self, X: np.ndarray, duration_labels: List[str]
    ) -> Dict[str, float]:
        """Train duration model."""
        metrics = {}
        self.logger.debug("Training duration model")
        
        self.models["duration_classifier"].fit(X, duration_labels)
        score = self.models["duration_classifier"].score(X, duration_labels)
        metrics["duration_classifier_score"] = score
        self.logger.debug(f"Duration classifier training score: {score:.3f}")
        
        return metrics

    def _train_pgp_model(
        self, X: np.ndarray, pgp_labels: List[bool]
    ) -> Dict[str, float]:
        """Train P-gp substrate model."""
        metrics = {}
        self.logger.debug("Training P-gp substrate model")
        
        self.models["pgp_classifier"].fit(X, pgp_labels)
        score = self.models["pgp_classifier"].score(X, pgp_labels)
        metrics["pgp_classifier_score"] = score
        self.logger.debug(f"P-gp classifier training score: {score:.3f}")
        
        return metrics
