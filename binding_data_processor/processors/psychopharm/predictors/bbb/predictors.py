"""BBB permeability predictors with transporter analysis.

This module extends the base BBB permeability prediction with:
1. Transporter substrate prediction (P-gp, BCRP, etc.)
2. Receptor-mediated transport prediction
3. Duration and pharmacokinetic predictions
4. Ensemble prediction methods
"""

import logging
from typing import Any, Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
import re

from .....models.core import CompoundData
from .....models.psychopharm import BBBPermeability
from .base import BBBPredictorBase, PredictionResult


class BBBPredictor(BBBPredictorBase):
    """BBB permeability predictor with transporter analysis."""

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

    # Amino acid patterns for LAT1 substrates
    LAT1_PATTERNS = {
        "alpha_amino_acid": r"NC\([R]\)C\(=O\)O",  # Basic α-amino acid pattern
        "aromatic_side_chain": r"c1[cn][cnh][cnh]c1",  # Aromatic rings (phenyl, indole)
        "hydroxyl_groups": r"[OH]",  # Hydroxyl groups common in LAT1 substrates
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

    # Additional model filenames
    MODEL_FILENAMES = {
        **BBBPredictorBase.MODEL_FILENAMES,
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
    ):
        """Initialize BBB predictor.
        
        Args:
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            log_level: Logging level
            transporters: Optional custom transporter definitions
            receptor_transporters: Optional custom receptor transporter definitions
        """
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
        )

        # Use custom transporters if provided
        self.transporters = transporters or self.TRANSPORTERS
        self.receptor_transporters = receptor_transporters or self.RECEPTOR_TRANSPORTERS

        # Configure model weights for ensemble
        self.model_weights.update({
            "ensemble_classifier": 0.25,
            "receptor_classifier": 0.1,
            "duration_classifier": 0.1,
        })

        # Add transporter columns to prediction history
        self.prediction_history = pd.concat([
            self.prediction_history,
            pd.DataFrame(columns=[
                'p_gp_substrate',
                'transporter_category',
                'transporter',
                'is_substrate',
                'transporter_confidence',
                'receptor_mediated',
                'receptor_type',
                'duration_class',
                'half_life',
                'cross_tolerance',
            ]),
        ], axis=1)

        self.logger.info("BBBPredictor initialized successfully")

    def _is_lat1_substrate(self, smiles: str) -> Tuple[bool, float]:
        """Check if compound matches LAT1 substrate patterns.
        
        Args:
            smiles: SMILES string of compound
            
        Returns:
            Tuple of (is_substrate, confidence)
        """
        matches = []
        confidences = []

        # Check each pattern
        for pattern_name, pattern in self.LAT1_PATTERNS.items():
            match = bool(re.search(pattern, smiles))
            matches.append(match)
            
            # Weight patterns differently
            if pattern_name == "alpha_amino_acid":
                confidences.append(0.6 if match else 0)
            elif pattern_name == "aromatic_side_chain":
                confidences.append(0.3 if match else 0)
            elif pattern_name == "hydroxyl_groups":
                confidences.append(0.1 if match else 0)

        # Compound is LAT1 substrate if it has α-amino acid pattern
        # plus at least one other feature
        is_substrate = matches[0] and any(matches[1:])
        
        # Confidence is sum of matched pattern confidences
        confidence = sum(confidences)

        return is_substrate, min(confidence, 1.0)

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate comprehensive BBB permeability predictions."""
        self.logger.debug(f"Generating predictions for {compound.name}")
        try:
            # Get base predictions
            features = self._extract_compound_features(compound)
            base_predictions = self._get_all_predictions(features)
            
            # Get transporter predictions
            transporter_predictions = self._predict_transporters(features, compound)
            
            # Get receptor predictions
            receptor_predictions = self._predict_receptors(features)
            
            # Get duration predictions
            duration_predictions = self._predict_duration(features)
            
            # Combine all predictions
            result = self._format_enhanced_prediction_result(
                base_predictions=base_predictions,
                transporter_predictions=transporter_predictions,
                receptor_predictions=receptor_predictions,
                duration_predictions=duration_predictions,
                compound=compound,
            )
            
            # Update history
            self._update_enhanced_prediction_history(
                base_predictions=base_predictions,
                transporter_predictions=transporter_predictions,
                receptor_predictions=receptor_predictions,
                duration_predictions=duration_predictions,
                compound=compound,
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

    def _predict_transporters(
        self, 
        features: Dict[str, np.ndarray],
        compound: CompoundData,
    ) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Predict transporter interactions."""
        predictions = {}
        
        # Check LAT1 substrate status first
        is_lat1, lat1_conf = self._is_lat1_substrate(compound.smiles)
        if is_lat1:
            predictions["uptake"] = {
                "lat1": {
                    "is_substrate": True,
                    "confidence": lat1_conf,
                    "mechanism": "amino_acid_transport",
                }
            }
        
        # Predict for each transporter category
        for category, transporters in self.transporters.items():
            if category not in predictions:
                predictions[category] = {}
                
            for transporter in transporters:
                # Skip LAT1 if already predicted
                if transporter == "lat1" and is_lat1:
                    continue
                    
                # Predict substrate status
                is_substrate, conf = self._predict_transporter_substrate(
                    self.models["transporter_classifier"],
                    features,
                    transporter,
                )
                
                if conf > 0.5:  # Confidence threshold
                    predictions[category][transporter] = {
                        "is_substrate": is_substrate,
                        "confidence": conf,
                    }
                    
        return predictions

    def _predict_transporter_substrate(
        self,
        model,
        features: Dict[str, np.ndarray],
        transporter: str,
    ) -> Tuple[bool, float]:
        """Predict if compound is a substrate for a transporter."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"transporter_{transporter}"].transform(X)

        # Get prediction probabilities
        probs = model.predict_proba(X_scaled)[0]
        is_substrate = bool(np.argmax(probs))
        confidence = float(probs[int(is_substrate)])

        return is_substrate, confidence

    def _predict_receptors(
        self, features: Dict[str, np.ndarray]
    ) -> Dict[str, Dict[str, Any]]:
        """Predict receptor-mediated transport."""
        predictions = {}
        
        for receptor, mechanisms in self.receptor_transporters.items():
            is_substrate, conf = self._predict_receptor_transport(
                self.models["receptor_classifier"],
                features,
                receptor,
            )
            
            if conf > 0.5:  # Confidence threshold
                predictions[receptor] = {
                    "is_substrate": is_substrate,
                    "confidence": conf,
                    "mechanisms": mechanisms,
                }
                
        return predictions

    def _predict_receptor_transport(
        self,
        model,
        features: Dict[str, np.ndarray],
        receptor: str,
    ) -> Tuple[bool, float]:
        """Predict if compound uses receptor-mediated transport."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"receptor_{receptor}"].transform(X)

        # Get prediction probabilities
        probs = model.predict_proba(X_scaled)[0]
        is_substrate = bool(np.argmax(probs))
        confidence = float(probs[int(is_substrate)])

        return is_substrate, confidence

    def _predict_duration(
        self, features: Dict[str, np.ndarray]
    ) -> Dict[str, str]:
        """Predict duration class and half-life."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers["duration"].transform(X)

        # Get prediction probabilities
        probs = self.models["duration_classifier"].predict_proba(X_scaled)[0]
        pred_idx = np.argmax(probs)
        
        # Map to duration class
        duration_class = list(self.DURATION_FEATURES.keys())[pred_idx]
        half_life = self.DURATION_FEATURES[duration_class]["half_life"]

        return {
            "class": duration_class,
            "half_life": half_life,
        }

    def _format_enhanced_prediction_result(
        self,
        base_predictions: Dict[str, Tuple[str, float]],
        transporter_predictions: Dict[str, Dict[str, Dict[str, float]]],
        receptor_predictions: Dict[str, Dict[str, Any]],
        duration_predictions: Dict[str, str],
        compound: CompoundData,
    ) -> PredictionResult:
        """Format all predictions into enhanced PredictionResult."""
        # Get base result
        result = super()._format_prediction_result(base_predictions, compound)
        
        # Add enhanced predictions
        result.supporting_data.update({
            "transporters": transporter_predictions,
            "receptor_transport": receptor_predictions,
            "duration_metrics": duration_predictions,
        })
        
        return result

    def _update_enhanced_prediction_history(
        self,
        base_predictions: Dict[str, Tuple[str, float]],
        transporter_predictions: Dict[str, Dict[str, Dict[str, float]]],
        receptor_predictions: Dict[str, Dict[str, Any]],
        duration_predictions: Dict[str, str],
        compound: CompoundData,
    ) -> None:
        """Update prediction history with enhanced predictions."""
        # Update base predictions
        super()._update_prediction_history(base_predictions, compound)
        
        # Get latest prediction index
        latest_idx = self.prediction_history.index[-1]
        
        # Update transporter predictions
        for category, transporters in transporter_predictions.items():
            for transporter, data in transporters.items():
                self.prediction_history.at[latest_idx, 'transporter_category'] = category
                self.prediction_history.at[latest_idx, 'transporter'] = transporter
                self.prediction_history.at[latest_idx, 'is_substrate'] = (
                    data["is_substrate"]
                )
                self.prediction_history.at[latest_idx, 'transporter_confidence'] = (
                    data["confidence"]
                )
        
        # Update receptor predictions
        if receptor_predictions:
            self.prediction_history.at[latest_idx, 'receptor_mediated'] = True
            self.prediction_history.at[latest_idx, 'receptor_type'] = (
                next(iter(receptor_predictions))
            )
        else:
            self.prediction_history.at[latest_idx, 'receptor_mediated'] = False
        
        # Update duration predictions
        self.prediction_history.at[latest_idx, 'duration_class'] = (
            duration_predictions["class"]
        )
        self.prediction_history.at[latest_idx, 'half_life'] = (
            duration_predictions["half_life"]
        )

    def retrain(
        self,
        compounds: List[CompoundData],
        labels: List[BBBPermeability],
        scores: Optional[List[float]] = None,
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
            transporter_data: Optional transporter substrate labels
            receptor_data: Optional receptor transport labels
            duration_labels: Optional duration class labels
            **kwargs: Additional training parameters
            
        Returns:
            Dictionary of training metrics
        """
        # Train base models
        metrics = super().retrain(
            compounds=compounds,
            labels=labels,
            scores=scores,
            **kwargs,
        )
        
        # Extract features once for all models
        X = []
        for compound in compounds:
            features = []
            for feature_type in self.feature_types:
                feat = self._extract_features(compound, feature_type)
                features.append(feat)
            X.append(np.hstack(features))
        X = np.vstack(X)
        
        # Train transporter models
        if transporter_data:
            self.logger.debug("Training transporter models")
            for transporter, labels in transporter_data.items():
                self.models["transporter_classifier"].fit(X, labels)
                score = self.models["transporter_classifier"].score(X, labels)
                metrics[f"transporter_{transporter}_score"] = score
                self.logger.debug(
                    f"Transporter classifier for {transporter} "
                    f"training score: {score:.3f}"
                )
        
        # Train receptor models
        if receptor_data:
            self.logger.debug("Training receptor models")
            for receptor, labels in receptor_data.items():
                self.models["receptor_classifier"].fit(X, labels)
                score = self.models["receptor_classifier"].score(X, labels)
                metrics[f"receptor_{receptor}_score"] = score
                self.logger.debug(
                    f"Receptor classifier for {receptor} "
                    f"training score: {score:.3f}"
                )
        
        # Train duration model
        if duration_labels:
            self.logger.debug("Training duration model")
            self.models["duration_classifier"].fit(X, duration_labels)
            score = self.models["duration_classifier"].score(X, duration_labels)
            metrics["duration_classifier_score"] = score
            self.logger.debug(f"Duration classifier training score: {score:.3f}")
        
        # Save updated models
        self.save_models()
        
        self.logger.info("Model retraining completed successfully")
        return metrics
